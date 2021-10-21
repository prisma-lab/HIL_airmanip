#include "ros/ros.h"
#include "nav_msgs/Odometry.h"
#include "boost/thread.hpp"
#include "mav_msgs/Actuators.h"
#include <tf/tf.h>
#include <tf_conversions/tf_eigen.h>
#include <tf/transform_listener.h>
#include <eigen3/Eigen/Dense>
#include <csignal>

#include <nav_msgs/Path.h>
#include <std_msgs/Bool.h>
#include <std_msgs/Float64MultiArray.h>
#include <std_msgs/Float32MultiArray.h>
#include <geometry_msgs/WrenchStamped.h>

#include "../include/aerialmanip_control/planner_spline.h"
#include "../include/aerialmanip_control/LowPassFilter.hpp"
#include <aerialmanip_control/drone_waypointsAction.h>
#include <actionlib/server/simple_action_server.h>

#define GEOMETRIC

double g = 9.8;

void rotateYaw(const geometry_msgs::PoseStamped& init, geometry_msgs::PoseStamped& final, double incyaw) {
  tf::Matrix3x3 initR;
  tf::Quaternion initq(init.pose.orientation.x,init.pose.orientation.y,init.pose.orientation.z,init.pose.orientation.w);
  initR.setRotation(initq);

  double roll,pitch,yaw;
  initR.getRPY(roll,pitch,yaw);
  yaw+=incyaw;
  initR.setRPY(roll,pitch,yaw);
  initR.getRotation(initq);

  final.pose.orientation.x = initq.x();
  final.pose.orientation.y = initq.y();
  final.pose.orientation.z = initq.z();
  final.pose.orientation.w = initq.w();
}

class DERIV {
	public:
		DERIV(double freq=200,double gain=5) {_f=freq;_dt=1.0/_f;_integral=Eigen::Vector3d::Zero();_gain=gain;};
		void update(Eigen::Vector3d x) {_xd=_gain*(x-_integral); _integral+=_gain*_dt*(x-_integral);};

		Eigen::Vector3d _xd;
	private:
		double _f,_dt;
		Eigen::Vector3d _integral;
		double _gain;

};


class QUAD_CTRL {
    public:
        QUAD_CTRL(double freq);
        ~QUAD_CTRL() {_shutdown=true; ROS_INFO("Shutting down...");};
        void odom_cb( nav_msgs::OdometryConstPtr );
        void run();
        void ctrl_loop();
        void setTraj(std::vector<geometry_msgs::PoseStamped> p, std::vector<geometry_msgs::TwistStamped> v, std::vector<geometry_msgs::AccelStamped> a);
        void updateFilter();
        bool isTrajReady() {return _traj_ready;};
        bool isTrajCompleted() {return _traj_completed;};
        bool isLanded() {return _landed;};
        void land() {_isLanding=true;};
        void getPose(double * Pose);
        void getPose(geometry_msgs::PoseStamped &pose);
        void takeoff();
    private:
        void updateError();
        void updatePassiveError();
        bool correctW(Vector4d & w);
        void actionCB(const aerialmanip_control::drone_waypointsGoalConstPtr &goal);
        void fbPos_loop();
        ros::NodeHandle _nh;
        ros::Subscriber _odom_sub;
        ros::Publisher _cmd_vel_pub, _fbackPos_pub;
        ros::Publisher _error_pub, _est_wrench_pub;
        ros::Publisher _forces_pub,_path_pub;
        double _freq;
        Vector3d _P;
        Vector3d _P_dot;
        Vector3d _Eta;
        Vector3d _Eta_dot;
        Vector3d _etaDot_des, _etadd_des;
        Vector3d _wbb;
        Matrix3d _Rb, _Rb_des;
        Matrix3d _RNed;
        Vector3d _Ep;
        Vector3d _Ev;
        Vector3d _Er, _Ew;
        Matrix3d _Q, _Qdot;
        Matrix4d _G;
        double _c_T;
        double _c_a;
        double _l;
        double _m;
        Matrix3d _I_b;
        std::vector<geometry_msgs::PoseStamped> _poses;
        std::vector<geometry_msgs::TwistStamped> _velocities;
        std::vector<geometry_msgs::AccelStamped> _accelerations;
        bool _traj_ready;
        bool _traj_completed;
        bool _isLanding;
        bool _landed;
        bool _odomOk;
        Vector3d _P_des, _Pd_des, _Pdd_des;
        Vector3d _wbb_des, _wbbd_des;
        int _iter;
        Matrix3d _Kp, _Kv, _Kr, _Kw;
        double _uT;
        Vector3d _tau_b;
        bool _isSettingTraj;
        tf::Quaternion _q_des;
        mav_msgs::Actuators _comm;
        bool _shutdown;
        VectorXd _Fe, _Fe_integral, _Fe_integral_out;
        bool _filterReady;
        DERIV _first,_second;
        actionlib::SimpleActionServer<aerialmanip_control::drone_waypointsAction> _actionServer;
		    aerialmanip_control::drone_waypointsFeedback _actionFeedback;
  		  aerialmanip_control::drone_waypointsResult _actionResult;
        bool _fbackPos_ready;
};

void QUAD_CTRL::takeoff() {

  CARTESIAN_PLANNER	cplanner(_freq);

  std::vector<geometry_msgs::PoseStamped> waypoints;
  geometry_msgs::PoseStamped final;
  getPose(final);
  waypoints.push_back(final);
  final.pose.position.z=0.6;
  rotateYaw(final,final,-(M_PI/4 + M_PI/2) );
  waypoints.push_back(final);

  std::vector<double> times;
  times.push_back(0);
  times.push_back(10);

	cplanner.set_waypoints(waypoints,times);
	cplanner.compute();

  setTraj(cplanner._x, cplanner._xd, cplanner._xdd);
  while(!isTrajCompleted()) usleep(1000);
  ROS_WARN("TAKEOFF COMPLETED!");
  ROS_WARN("UPDATING WRENCH FILTER...");
  ros::Duration d = ros::Duration(10, 0);
  d.sleep();
  ROS_WARN("FILTER READY!");
  _filterReady=true;
  d = ros::Duration(5, 0);
  d.sleep();
  ROS_WARN("SENDING POSITION!");
  _fbackPos_ready = true;
}

bool QUAD_CTRL::correctW(Vector4d & w2) {

  for (int i=0; i<4; i++) {
    if(w2(i)<0)
      return false;
  }
  return true;

}

void QUAD_CTRL::getPose(double * Pose) {
  Vector3d P;
  Matrix3d Rb;
  tf::Matrix3x3 RbTf;
  double roll, pitch, yaw;

  while(!_odomOk) usleep(10);

  P = _RNed.transpose()*_P;

  Rb = _RNed.transpose()*_Rb*_RNed;
  matrixEigenToTF(Rb, RbTf);
  RbTf.getRPY(roll, pitch, yaw);
  Pose[0] = P(0); Pose[1] = P(1); Pose[2] = P(2); Pose[3] = yaw;
}

void QUAD_CTRL::getPose(geometry_msgs::PoseStamped &pose) {
  Vector3d P;
  Matrix3d Rb;
  tf::Matrix3x3 RbTf;
  tf::Quaternion q;

  while(!_odomOk) usleep(10);

  P = _RNed.transpose()*_P;

  Rb = _RNed.transpose()*_Rb*_RNed;
  matrixEigenToTF(Rb, RbTf);
  RbTf.getRotation(q);

  pose.pose.position.x = P(0);
  pose.pose.position.y = P(1);
  pose.pose.position.z = P(2);
  pose.pose.orientation.x = q.x();
  pose.pose.orientation.y = q.y();
  pose.pose.orientation.z = q.z();
  pose.pose.orientation.w = q.w();

}

void QUAD_CTRL::updateFilter() {
  double zita_=1.0, omega_lin=20, omega_tor=20, omega_z=omega_tor, omega_yaw=10;
  MatrixXd zita(6,6), omega(6,6);
  MatrixXd K1(6,6),K2(6,6);

  zita = MatrixXd::Identity(6,6);
  omega = MatrixXd::Identity(6,6);
  zita.diagonal() << zita_,zita_,zita_,zita_,zita_,zita_;
  omega.diagonal() << omega_lin,omega_lin,omega_z,omega_tor,omega_tor,omega_yaw;

  K1 = 2*zita*omega;
  K2 = omega*omega*K1.inverse();

  Vector3d e3(0,0,1);

  MatrixXd M_xi(6,6);
  M_xi << _m*MatrixXd::Identity(3,3) , MatrixXd::Zero(3,3),
      MatrixXd::Zero(3,3) , _I_b;

  VectorXd alpha(6);
  alpha.head(3)=_Rb.transpose()*_P_dot;
  //alpha.head(3)=_P_dot;
  alpha.tail(3)=_wbb;

  VectorXd internal(6);
  internal.head(3) = -_uT*e3 + _m*g*_Rb.transpose()*e3;
  //internal.head(3) = -_uT*_Rb*e3 + _m*g*e3;
  internal.tail(3) = _tau_b - Skew(_wbb)*_I_b*_wbb;

  _Fe_integral += ( internal + _Fe )*(1.0/_freq);
  _Fe_integral_out += ( -_Fe + K2*( M_xi*alpha - _Fe_integral ) )*(1.0/_freq);
  _Fe = K1*_Fe_integral_out;

  //Vector3d _Fe_b = _Rb.transpose()*_Fe.head(3);

  geometry_msgs::WrenchStamped est_w;
  est_w.header.stamp=ros::Time::now();
  est_w.wrench.force.x = _Fe(0);
  est_w.wrench.force.y = _Fe(1);
  est_w.wrench.force.z = _Fe(2);
  est_w.wrench.torque.x = _Fe(3);
  est_w.wrench.torque.y = _Fe(4);
  est_w.wrench.torque.z = _Fe(5);
  _est_wrench_pub.publish(est_w);

}

void QUAD_CTRL::updatePassiveError() {
  Vector3d mu_d;
  Vector3d mu_s;
  Vector3d xb_des;
  Vector3d e3(0,0,1);
  Matrix3d Rb_des;
  Matrix3d Q_des;
  Matrix3d Qd_des;
  Vector3d wbb_des;
  Vector3d etaDot_des;
  Vector3d eta_des;
  Vector3d wbbd_des;
  Vector3d etadd_des;

  double ni = 1;

  MatrixXd M(3,3);
  M << (_Q.transpose()*_I_b*_Q);

  MatrixXd C(3,3);
  C << (_Q.transpose()*_I_b*_Qdot + _Q.transpose()*Skew(_Q*_Eta_dot)*_I_b*_Q);

  _P_des[0] = _poses[_iter].pose.position.x;
  _P_des[1] = _poses[_iter].pose.position.y;
  _P_des[2] = _poses[_iter].pose.position.z;
  _P_des = _RNed * _P_des;
  _Pd_des[0] = _velocities[_iter].twist.linear.x;
  _Pd_des[1] = _velocities[_iter].twist.linear.y;
  _Pd_des[2] = _velocities[_iter].twist.linear.z;
  _Pd_des = _RNed*_Pd_des;
  _Pdd_des[0] = _accelerations[_iter].accel.linear.x;
  _Pdd_des[1] = _accelerations[_iter].accel.linear.y;
  _Pdd_des[2] = _accelerations[_iter].accel.linear.z;
  _Pdd_des = _RNed*_Pdd_des;

  _Ep = _P - _P_des;
  _Ev = _P_dot - _Pd_des;

  mu_d = _Pdd_des - (1.0/_m) * ( _Kp*_Ep + _Kv*_Ev );

  mu_s = mu_d;// - (1.0/_m)*_Fe.head(3);
  if(_filterReady) {
    //mu_s -= (1.0/_m)*_Fe.head(3);
  }
  _uT = _m*sqrt( mu_s(0)*mu_s(0) + mu_s(1)*mu_s(1) + (mu_s(2)-g)*(mu_s(2)-g));

  tf::Quaternion q_des(_poses[_iter].pose.orientation.x, _poses[_iter].pose.orientation.y, _poses[_iter].pose.orientation.z,  _poses[_iter].pose.orientation.w);
  _q_des = q_des;
  tf::Matrix3x3 Rb_des_tf(_q_des);

  tf::matrixTFToEigen(Rb_des_tf,Rb_des);
  Rb_des = _RNed*Rb_des*(_RNed.transpose()); //Rb NED transform
  tf::matrixEigenToTF(Rb_des,Rb_des_tf);

  Rb_des_tf.getRPY(eta_des(0),eta_des(1),eta_des(2));
  eta_des(0) = asin(_m*( mu_s(1)*cos(eta_des(2)) - mu_s(0)*sin(eta_des(2)) ) / _uT );
  eta_des(1) = atan2( mu_s(0)*cos(eta_des(2)) + mu_s(1)*sin(eta_des(2)) , mu_s(2)-g );
  if(eta_des(1)<0) eta_des(1)+=M_PI;
  else eta_des(1)-=M_PI;
  
  //Derivate
  _first.update(eta_des);
  etaDot_des = _first._xd;
  _second.update(etaDot_des);
  etadd_des = _second._xd;

  Vector3d e_eta, edot_eta, v_eta, etadot_r, etadd_r;
  e_eta = _Eta - eta_des;
  edot_eta = _Eta_dot - etaDot_des;
  v_eta = edot_eta + ni*e_eta;
  etadot_r = etaDot_des - ni*e_eta;
  etadd_r = etadd_des - ni*edot_eta;

  Vector3d controltau = M*etadd_des + C*etadot_r - _Kw*v_eta - _Kr*e_eta;
  if(_filterReady) {
    //controltau -= _Fe.tail(3);
  }
  _tau_b = (_Q.transpose()).inverse() * controltau;
  

    if(_isLanding) {
    _uT=0;
    _tau_b=Vector3d::Zero();
  }


  if (_iter < (_poses.size() - 1)) _iter++;
  else if (!_traj_completed) {
    _traj_completed = true;
    _traj_ready = false;
  }

  _odomOk =true;

  std_msgs::Float64MultiArray errors_data;
  std_msgs::Float64MultiArray forces_data;
  forces_data.data.resize(4);
  errors_data.data.resize(12);
  for (int i=0; i<3; i++) {
    errors_data.data[i] = _Ep(i);
    errors_data.data[i+3] = _Ev(i);
    errors_data.data[i+6] = e_eta(i);
    errors_data.data[i+9] = edot_eta(i);
  }

  forces_data.data[0] = _uT;
  for (int i=0; i<3; i++)
    forces_data.data[1+i] = _tau_b(i);

  _forces_pub.publish(forces_data);
  _error_pub.publish(errors_data);


}

void QUAD_CTRL::updateError() {
  Vector3d zb_des;
  Vector3d yb_des;
  Vector3d xb_des;
  Vector3d e3(0,0,1);
  Matrix3d Rb_des;
  Matrix3d Q_des;
  Matrix3d Qd_des;
  Vector3d wbb_des;
  Vector3d etaDot_des;
  Vector3d eta_des;
  Vector3d wbbd_des;
  Vector3d etadd_des;

  _P_des[0] = _poses[_iter].pose.position.x;
  _P_des[1] = _poses[_iter].pose.position.y;
  _P_des[2] = _poses[_iter].pose.position.z;
  _P_des = _RNed * _P_des;
  _Pd_des[0] = _velocities[_iter].twist.linear.x;
  _Pd_des[1] = _velocities[_iter].twist.linear.y;
  _Pd_des[2] = _velocities[_iter].twist.linear.z;
  _Pd_des = _RNed*_Pd_des;
  _Pdd_des[0] = _accelerations[_iter].accel.linear.x;
  _Pdd_des[1] = _accelerations[_iter].accel.linear.y;
  _Pdd_des[2] = _accelerations[_iter].accel.linear.z;
  _Pdd_des = _RNed*_Pdd_des;

  _Ep = _P - _P_des;
  _Ev = _P_dot - _Pd_des;

  zb_des = _Kp*_Ep + _Kv*_Ev + _m*g*e3 - _m*_Pdd_des;
  if(_filterReady) {
    Vector3d worldForce = _Rb*_Fe.head(3);
    zb_des(0) += worldForce(0);
    zb_des(1) += worldForce(1);
    zb_des(2) += worldForce(2);
  }

  //Projection on a cone
  double cone_a =10*M_PI/180;
  double cone_tan = tan(cone_a);
  double angle_a = acos(e3.dot(zb_des)/(e3.norm()*zb_des.norm()));
  if( (angle_a>cone_a) || isnan(angle_a)) {
    double fnorm = sqrt(zb_des(0)*zb_des(0) + zb_des(1)*zb_des(1));
    Vector3d zb_proj(zb_des(0),zb_des(1),fnorm/cone_tan);
    zb_proj = zb_proj*(zb_des(2)/zb_proj(2));
    ROS_WARN("Correggo");
    cout<<"proj: "<<zb_proj.transpose()<<endl;
    cout<<"real: "<<zb_des.transpose()<<endl;
    zb_des = zb_proj;
  }
  
  

  _uT = zb_des.transpose() * _Rb * e3;

  if(_isLanding) {
    _Ep(2)=_Ev(2)=0;
    zb_des = _Kp*_Ep + _Kv*_Ev + _m*g*e3 - _m*_Pdd_des;
    _uT = zb_des.transpose() * _Rb * e3;
    _uT = 0.95*_uT;
  }
  zb_des = zb_des/zb_des.norm();

  tf::Quaternion q_des(_poses[_iter].pose.orientation.x, _poses[_iter].pose.orientation.y, _poses[_iter].pose.orientation.z,  _poses[_iter].pose.orientation.w);
  _q_des = q_des;
  wbb_des << _velocities[_iter].twist.angular.x,_velocities[_iter].twist.angular.y,_velocities[_iter].twist.angular.z;
  wbbd_des << _accelerations[_iter].accel.angular.x,_accelerations[_iter].accel.angular.y,_accelerations[_iter].accel.angular.z;

  tf::Matrix3x3 Rb_des_tf(_q_des);
  //Rb_des_tf.getRPY(eta_des(0), eta_des(1), eta_des(2)); //phi theta psi
  tf::matrixTFToEigen(Rb_des_tf,Rb_des);

  Rb_des = _RNed*Rb_des*(_RNed.transpose()); //Rb NED transform
  
    xb_des = Rb_des.col(0);
    yb_des = zb_des.cross(xb_des);
    yb_des = yb_des / yb_des.norm();
    xb_des = yb_des.cross(zb_des);
    Rb_des << xb_des(0), yb_des(0), zb_des(0),
              xb_des(1), yb_des(1), zb_des(1),
              xb_des(2), yb_des(2), zb_des(2);
  

  _Rb_des = Rb_des;

  Vector3d appo(wbb_des(1),wbb_des(0),-wbb_des(2));
  wbb_des = appo;
  _wbb_des = wbb_des;

  Vector3d appo1(wbbd_des(1),wbbd_des(0),-wbbd_des(2));
  wbbd_des = appo1;
  _wbbd_des = wbbd_des;

  _Er = 0.5*Vee(_Rb_des.transpose()*_Rb - _Rb.transpose()*_Rb_des);
  _Ew = _wbb - _Rb.transpose()*_Rb_des*_wbb_des;

  _tau_b = -_Kr*_Er - _Kw*_Ew + Skew(_wbb)*_I_b*_wbb - _I_b*( Skew(_wbb)*_Rb.transpose()*_Rb_des*_wbb_des - _Rb.transpose()*_Rb_des*_wbbd_des );
  if(_filterReady) {
    Vector3d estTorq = _Fe.tail(3);
    _tau_b -= estTorq;
  }

  double xy_sat = 3, yaw_sat = 1.5;
  if(fabs(_tau_b(2))>yaw_sat) {
    if(_tau_b(2)>0) _tau_b(2) = yaw_sat;
    else _tau_b(2) = -yaw_sat;
  }

  if(_isLanding) {
    _uT=0;
    _tau_b=Vector3d::Zero();
  }
  if (_iter < (_poses.size() - 1)) _iter++;
  else if (!_traj_completed) {
    _traj_completed = true;
    _traj_ready = false;
  }

  _odomOk =true;

}

void QUAD_CTRL::setTraj(std::vector<geometry_msgs::PoseStamped> p, std::vector<geometry_msgs::TwistStamped> v, std::vector<geometry_msgs::AccelStamped> a) {

  nav_msgs::Path path;
  path.header.frame_id = "world";
  path.header.stamp=ros::Time::now();
  path.poses = p;
  _path_pub.publish(path);
  cout<<"traj size: "<<p.size()<<endl;
  for (int i=0; i<p.size(); i++) {
    _poses.push_back(p[i]);
    _velocities.push_back(v[i]);
    _accelerations.push_back(a[i]);
  }

  _traj_ready = true;
  _traj_completed = false;
  //_iter = 0;
}

QUAD_CTRL::QUAD_CTRL(double freq) : _actionServer(_nh, "droneActionServer", boost::bind(&QUAD_CTRL::actionCB, this, _1), false) {
    _odom_sub = _nh.subscribe("/hummingbird/ground_truth/odometry", 0, &QUAD_CTRL::odom_cb, this, ros::TransportHints().tcpNoDelay());

    _cmd_vel_pub = _nh.advertise< mav_msgs::Actuators>("/hummingbird/command/motor_speed", 0);
    _error_pub = _nh.advertise<std_msgs::Float64MultiArray>("/controller/errors", 0);
    _forces_pub = _nh.advertise<std_msgs::Float64MultiArray>("/controller/forces", 0);
    _path_pub = _nh.advertise<nav_msgs::Path>("/controller/plannedPath", 0);
    _est_wrench_pub = _nh.advertise<geometry_msgs::WrenchStamped>("/controller/estimatedWrench", 0);
    _fbackPos_pub = _nh.advertise<std_msgs::Float64MultiArray>("/controller/posFeedback", 0);

    _traj_ready = false;
    _traj_completed = false;
    _isLanding = false;
    _landed = false;
    _odomOk = false;
    _iter = 0;
    _shutdown = false;
    _filterReady = false;
    _fbackPos_ready = false;

    _freq = freq;

    _P.resize(3);
    _Eta.resize(3);

    _l = 0.17; //meters
    _c_T = 11.54802e-05;
    _c_a = -0.016*_c_T;//-0.000001/7.311;
    cout<<_c_a<<endl;
    //_c_a = -0.016;
    _m = 11.68 + 0.05;//0.68; //Kg
 
    //Inertia matrix
    _I_b << 0.777, 0, 0,
           0, 0.777, 0,
           0, 0, 1.112;

    double mRobot = 0.97;//Kg

    Matrix3d IRobot;
    IRobot<< 0.34, 0,0,
             0, 0.37, 0,
             0, 0, 0.2;

    //_I_b+=_I_b;
    _m+=mRobot;

    _G(0,0) = _c_T;    _G(0,1) = _c_T;    _G(0,2) = _c_T; _G(0,3) = _c_T;
    _G(1,0) = 0;       _G(1,1) = _l*_c_T; _G(1,2) = 0;    _G(1,3) = -_l*_c_T;
    _G(2,0) = -_l*_c_T; _G(2,1) = 0;       _G(2,2) = _l*_c_T; _G(2,3) = 0;
    _G(3,0) = -_c_a;    _G(3,1) = _c_a;    _G(3,2) = -_c_a; _G(3,3) = _c_a;

    //RPY CONTROLLER
    #ifdef RPY
      ROS_WARN("Starting RPY controller!");
      _Kp = 0.2 *Vector3d(20,20,50).asDiagonal();
      _Kv = _Kp/6;
      _Kr = 0.4 *Vector3d(40,40,20).asDiagonal();
      _Kw = _Kr/12;
    #else
    //GEOMETRIC CONTROLLER
      ROS_WARN("Starting GEOMETRIC controller!");
      _Kp = Vector3d(20,20,50).asDiagonal();
      _Kv = _Kp/1.2;
      _Kv(2) = _Kp(2)/6;
      _Kr = 4.0*Vector3d(30,30,20).asDiagonal();
      _Kw = _Kr/4;
      //_Kw(2) = _Kr(2)/2;
    #endif

    _Fe.resize(6);
    _Fe = Eigen::VectorXd::Zero(6);
    _Fe_integral.resize(6);
    _Fe_integral = Eigen::VectorXd::Zero(6);
    _Fe_integral_out.resize(6);
    _Fe_integral_out = Eigen::VectorXd::Zero(6);
    _Qdot.resize(3,3);
    _uT = 0;
    _tau_b = Vector3d::Zero();

    _actionServer.start();

}

void QUAD_CTRL::odom_cb( nav_msgs::OdometryConstPtr odom ) {
    //static double secs;
    //cout << ros::Time::now().toSec() - secs <<endl;
    //secs = ros::Time::now().toSec();
    //fflush(NULL);

    tf::Matrix3x3 RNed;
    RNed.setEulerYPR(M_PI/2,0,M_PI);
    //RNed = RNed.transpose();
    tf::Vector3 p;
    tf::Vector3 pDot;
    tf::Vector3 wbb;
    tf::Vector3 wb;
    tf::Vector3 etaDot;

    tf::Quaternion q(odom->pose.pose.orientation.x, odom->pose.pose.orientation.y, odom->pose.pose.orientation.z,  odom->pose.pose.orientation.w);
    tf::Matrix3x3 Rb(q);
    tf::Matrix3x3 RbNed = RNed*Rb*(RNed.transpose());

    RbNed.getRPY(_Eta(0), _Eta(1), _Eta(2)); //phi theta psi
    tf::matrixTFToEigen (RbNed, _Rb);
    tf::matrixTFToEigen (RNed, _RNed);

    p[0] = odom->pose.pose.position.x;
    p[1] = odom->pose.pose.position.y;
    p[2] = odom->pose.pose.position.z;

    p = RNed*p;
    _P(0) = p[0];
    _P(1) = p[1];
    _P(2) = p[2];

    _Q(0,0) = 1; _Q(0,1) = 0;             _Q(0,2) = -sin(_Eta(1));
    _Q(1,0) = 0; _Q(1,1) = cos(_Eta(0));  _Q(1,2) = cos(_Eta(1))*sin(_Eta(0));
    _Q(2,0) = 0; _Q(2,1) = -sin(_Eta(0)); _Q(2,2) = cos(_Eta(1))*cos(_Eta(0));

    pDot[0] = odom->twist.twist.linear.x;
    pDot[1] = odom->twist.twist.linear.y;
    pDot[2] = odom->twist.twist.linear.z;
    //cout<<odom->twist.twist.linear.x;
    pDot = RNed*Rb*pDot;
    //cout<<pDot[0];
    //pDot = RNed*Rb*pDot*RNed.transpose();

    _P_dot(0) = pDot[0];
    _P_dot(1) = pDot[1];
    _P_dot(2) = pDot[2];

    wbb[0] = odom->twist.twist.angular.y;
    wbb[1] = odom->twist.twist.angular.x;
    wbb[2] = -odom->twist.twist.angular.z;

    _wbb(0) = wbb[0];
    _wbb(1) = wbb[1];
    _wbb(2) = wbb[2];

    _Eta_dot = _Q.inverse() * _wbb;

    _Qdot(0,0) = 0; _Qdot(0,1) = 0;                         _Qdot(0,2) = -cos(_Eta(1))*_Eta_dot(1);
    _Qdot(1,0) = 0; _Qdot(1,1) = -sin(_Eta(0))*_Eta_dot(0); _Qdot(1,2) = -sin(_Eta(1))*sin(_Eta(0))*_Eta_dot(1) + cos(_Eta(1))*cos(_Eta(0))*_Eta_dot(0);
    _Qdot(2,0) = 0; _Qdot(2,1) = -cos(_Eta(0))*_Eta_dot(0); _Qdot(2,2) = -sin(_Eta(1))*cos(_Eta(0))*_Eta_dot(1) + -cos(_Eta(1))*sin(_Eta(0))*_Eta_dot(0);

    _odomOk = true;
    //ROS_INFO("x: %f  y: %f z: %f",_P_dot(0),_P_dot(1),_P_dot(2));
    //ROS_INFO("x: %f  y: %f z: %f",_P(0),_P(1),_P(2));
    //ROS_INFO("phi: %f  tetha: %f psi: %f",_Eta(0)*180.0/M_PI,_Eta(1)*180.0/M_PI,_Eta(2)*180.0/M_PI);
}

void QUAD_CTRL::fbPos_loop() {
  ros::Rate r(100);

  while(!_fbackPos_ready) usleep(10);
  geometry_msgs::PoseStamped initPose;
  getPose(initPose);

  while(ros::ok() && !_isLanding) {
    geometry_msgs::PoseStamped actual;
    getPose(actual);
    std_msgs::Float64MultiArray msg;
    msg.data.resize(3);
    msg.data[0] = actual.pose.position.x - initPose.pose.position.x;
    msg.data[1] = actual.pose.position.y - initPose.pose.position.y;
    msg.data[2] = actual.pose.position.z - initPose.pose.position.z;
    _fbackPos_pub.publish(msg);

    r.sleep();
  }
}

void QUAD_CTRL::ctrl_loop() {

  ros::Rate r(_freq);
  _comm.angular_velocities.resize(4);
  _comm.angular_velocities[0] = 0;
  _comm.angular_velocities[1] = 0;
  _comm.angular_velocities[2] = 0;
  _comm.angular_velocities[3] = 0;
  Vector4d w2, controlInput;


  /* Red: 0 poi gli 1-2-3 in senso antiorario
     0 frontale asse x - 1 frontale asse y
     NED: 1 frontale asse x - 0 frontale asse y
          w1 = 1, w4 = 2, w3 = 3, w2 = 0*/

  ROS_INFO("Controllo attivo");

  while(!_odomOk) {
    usleep(1e3);
  }
  while(!isTrajReady() || !_odomOk) {
    _cmd_vel_pub.publish (_comm);
    r.sleep();
  }

  while( ros::ok() && !_shutdown ) {

    #ifdef RPY
      updatePassiveError();
    #else
      updateError();
    #endif

    controlInput(0) = _uT      ;
    controlInput(1) = _tau_b(0);
    controlInput(2) = _tau_b(1);
    controlInput(3) = _tau_b(2);
    w2 = _G.inverse() * controlInput;
    _comm.header.stamp = ros::Time::now();

    if (correctW(w2)) {
      _comm.angular_velocities[0] = sqrt(w2(3));
      _comm.angular_velocities[1] = sqrt(w2(2));
      _comm.angular_velocities[2] = sqrt(w2(1));
      _comm.angular_velocities[3] = sqrt(w2(0));
    }
    else {
      ROS_WARN("w problem");
      w2 << _comm.angular_velocities[3]*_comm.angular_velocities[3], _comm.angular_velocities[2]*_comm.angular_velocities[2], _comm.angular_velocities[1]*_comm.angular_velocities[1], _comm.angular_velocities[0]*_comm.angular_velocities[0];
      controlInput = _G * w2;
      _uT = controlInput(0);
      _tau_b(0) = controlInput(1);
      _tau_b(1) = controlInput(2);
      _tau_b(2) = controlInput(3);
    }
    
    if(_filterReady && !_isLanding)
      updateFilter();
    
    if(_isLanding && _P(2)>-0.065) _landed=true;

    if(_landed) {
      _comm.angular_velocities[0] = 0;
      _comm.angular_velocities[1] = 0;
      _comm.angular_velocities[2] = 0;
      _comm.angular_velocities[3] = 0;
    }

   // cout<<"gain:"<<((_m*9.8)/(w2.sum())*1e6) << endl;
    _cmd_vel_pub.publish (_comm);

    std_msgs::Float64MultiArray errors_data;
    std_msgs::Float64MultiArray forces_data;
    forces_data.data.resize(4);
    errors_data.data.resize(12);
    for (int i=0; i<3; i++) {
      errors_data.data[i] = _Ep(i);
      errors_data.data[i+3] = _Ev(i);
      errors_data.data[i+6] = _Er(i);
      errors_data.data[i+9] = _Ew(i);
    }

    forces_data.data[0] = _uT;
    for (int i=0; i<3; i++)
      forces_data.data[1+i] = _tau_b(i);

    _forces_pub.publish(forces_data);
    _error_pub.publish(errors_data);
    r.sleep();
  }

}

void QUAD_CTRL::run() {
    boost::thread ctrl_loop_t ( &QUAD_CTRL::ctrl_loop, this);
    boost::thread fbPos_loop_t ( &QUAD_CTRL::fbPos_loop, this);
}

void QUAD_CTRL::actionCB(const aerialmanip_control::drone_waypointsGoalConstPtr &goal) {
	bool result;
	
  std::vector<geometry_msgs::PoseStamped> waypoints = goal->waypoints.poses;
	std::vector<double> times = goal->times;
	//Eigen::VectorXd xdi,xdf,xddi,xddf;
	//twist2Vector(goal->initVel,xdi);
	//twist2Vector(goal->finalVel,xdf);
	//accel2Vector(goal->initAcc,xddi);
	//accel2Vector(goal->finalAcc,xddf);
	//result = newTrajectory(waypoints,times,xdi,xdf,xddi,xddf);

  CARTESIAN_PLANNER	cplanner(_freq);
  cplanner.set_waypoints(waypoints,times);
	cplanner.compute();

  setTraj(cplanner._x, cplanner._xd, cplanner._xdd);
  while(!isTrajCompleted() && ros::ok()) usleep(1000);
  result = true;
	if(result) {
      _actionResult.ok = true;
      ROS_INFO("ACTION: Succeeded");
      // set the action state to succeeded
      _actionServer.setSucceeded(_actionResult);
  } else {
		_actionResult.ok = false;
		ROS_INFO("ACTION: Aborted. Probably already following another trajectory.");
		// set the action state to succeeded
		_actionServer.setAborted(_actionResult);
	}
}

bool shut=false;
QUAD_CTRL* cpoint;

void signalHandler( int signum ) {
    ROS_WARN("LANDING!");
    cpoint->land();
    sleep(1);
    shut = true;
    ROS_INFO("PATH COMPLETED!");
    exit(0);
}

int main( int argc, char** argv) {

  ros::init(argc, argv, "aerialmanip_controller" );
  ros::AsyncSpinner spinner(1); // Use 1 thread
	spinner.start();

  QUAD_CTRL c(200);
  cpoint = &c;

  c.run();
  signal(SIGINT, signalHandler);

  CARTESIAN_PLANNER	cplanner(200);

  c.takeoff();
  sleep(5);

  while(ros::ok()) sleep(1);

  // Uncomment to use trajectory planning
  /*
  std::vector<geometry_msgs::PoseStamped> waypoints;
  geometry_msgs::PoseStamped init,final;
  c.getPose(init);
  final=init;
  final.pose.position.x+=0.0;
  final.pose.position.y-=3.0;
  final.pose.position.z+=6.0;
  rotateYaw(init,final,M_PI/2);
  waypoints.push_back(init);
  waypoints.push_back(final);
  std::vector<double> times;
  times.push_back(0);
  times.push_back(5);
  cplanner.set_waypoints(waypoints,times);
	cplanner.compute();
  c.setTraj(cplanner._x, cplanner._xd, cplanner._xdd);
  while(!c.isTrajCompleted() && ros::ok()) usleep(1000);
  ROS_WARN("WAYPOINT REACHED!");
  */

  return 0;

}
