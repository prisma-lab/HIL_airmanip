#include "ros/ros.h"
#include "boost/thread.hpp"
#include "sensor_msgs/JointState.h"
#include "geometry_msgs/PoseStamped.h"
#include "geometry_msgs/TwistStamped.h"
#include "geometry_msgs/AccelStamped.h"
#include "gazebo_msgs/ContactsState.h"
#include <std_msgs/Float64.h>
#include <std_msgs/Float64MultiArray.h>
#include <geometry_msgs/WrenchStamped.h>
#include <geometry_msgs/PointStamped.h>

#include <kdl_parser/kdl_parser.hpp>
#include <kdl/chainfksolvervel_recursive.hpp>
#include <kdl/chainfksolverpos_recursive.hpp>
#include <kdl/chainiksolvervel_pinv.hpp>
#include <kdl/chainfksolverpos_recursive.hpp>
#include <kdl/chainiksolverpos_nr.hpp>
#include <kdl/chainjnttojacsolver.hpp>
#include <kdl/chainjnttojacdotsolver.hpp>
#include <kdl/chaindynparam.hpp>

#include "../include/hardware_controller/planner.h"
#include <hardware_controller/waypointsAction.h>
#include <actionlib/server/simple_action_server.h>

#include "../include/hardware_controller/LowPassFilter.hpp"

using namespace std;

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

class ETankGen {
	public:
		ETankGen(double Einit, double Emin, double Emax, double dt, int inputSize);
		void update(const std::vector<double> inputs, const std::vector<double> dissInputs);
		double getEt() {return _Et;};
		std::vector<double> _alpha;
	private:
		double _Et, _Emin, _Emax, _xt;
		double _beta;
		double _dt;
};

ETankGen::ETankGen(double Einit, double Emin, double Emax, double dt, int inputSize) {
	_Et=Einit;
	_Emin=Emin;
	_Emax=Emax;
	_xt=sqrt(2*_Et);
	_dt=dt;
	for (int i=0; i<inputSize; i++) {
		_alpha.push_back(1);
	}
}

void ETankGen::update(const std::vector<double> inputs, const std::vector<double> dissInputs) {
	if(_Et<=_Emax) _beta=1;
	else _beta=0;

	double f_energy = 0.5*( 1 - cos(M_PI*(_Et-_Emin)/(_Emax-_Emin)) );
	double eta=0.8;
	double Diss=0;
	double wtot = 0;

	for (int i=0; i<dissInputs.size(); i++)
		Diss+=dissInputs[i];

	for (int i=0; i<_alpha.size(); i++) {
		double w = inputs[i];
		double g_input;
		if(w<=0) g_input=0;
		else	g_input=1;

		_alpha[i] = f_energy*g_input + (1-g_input);
		double gamma;
		if( (_Et>=_Emin) && w>=0 ) gamma=_alpha[i];
		else gamma=0;

		//if(i==1) cout<<"Power: "<<gamma*w<<endl;

		wtot+=gamma*w;
		//if(w!=0) cout<<"Power: "<<w<<"  Diss: "<<Diss<<endl;
	}

	double ut = -(1.0/_xt)*wtot;
	double xt_dot = (_beta*eta/_xt)*Diss + ut;
	_xt += xt_dot*_dt;
	if (_xt>sqrt(_Emax*2)) _xt=sqrt(_Emax*2);
	_Et = 0.5*_xt*_xt;
}

enum diverterState {IMPACT, DETACHED, HOOKED, NORMAL};

class KUKA_INVDYN {
	public:
		KUKA_INVDYN(double sampleTime);
		void run();
		bool init_robot_model();
		void get_dirkin();

		void joint_states_cb( sensor_msgs::JointState );
		void real_interaction_wrench_cb(const geometry_msgs::WrenchStampedConstPtr&);
		void drone_posfb_cb(const std_msgs::Float64MultiArrayConstPtr& message);
		void ctrl_loop();
		void compute_compliantFrame(const geometry_msgs::PoseStamped& p_des, const geometry_msgs::TwistStamped& v_des, const geometry_msgs::AccelStamped& a_des);
		bool newTrajectory(const std::vector<geometry_msgs::PoseStamped> waypoints, const std::vector<double> times);
		bool newTrajectory(const std::vector<geometry_msgs::PoseStamped> waypoints, const std::vector<double> times, const Eigen::VectorXd xdi, const Eigen::VectorXd xdf, const Eigen::VectorXd xddi, const Eigen::VectorXd xddf);
		bool getDesPose(geometry_msgs::PoseStamped& p_des);
		const diverterState getState() {return _state;};
		void actionCB(const hardware_controller::waypointsGoalConstPtr &goal);
		void setDone(bool done) {_mainDone=done;};
	private:
		void updatePose();
		void updateState();
		ros::NodeHandle _nh;
		KDL::Tree iiwa_tree;

		KDL::ChainFkSolverPos_recursive *_fksolver; //Forward position solver
		KDL::ChainFkSolverVel_recursive *_fk_solver_pos_vel; //Forward position and velocity solver
		KDL::ChainIkSolverVel_pinv *_ik_solver_vel;   	//Inverse velocity solver
		KDL::ChainIkSolverPos_NR *_ik_solver_pos;
		KDL::ChainJntToJacSolver *_J_solver;

		KDL::Chain _k_chain;

		ros::Subscriber _js_sub;
		ros::Publisher _js_pub;
		ros::Subscriber _real_wrench_sub, _dronePosFb_sub;
		ros::Publisher _cartpose_pub, _cartvel_pub, _desPose_pub, _extWrench_pub;
		ros::Publisher _plannedpose_pub,_plannedtwist_pub,_plannedacc_pub,_plannedwrench_pub;
		ros::Publisher _tankEnergy_pub, _kpvalue_pub;
		KDL::JntArray *_initial_q;
		KDL::JntArray *_q_in;
		KDL::JntArray *_q_out;
		KDL::JntArray *_q_in_old;
		KDL::JntArray *_dq_in;
		ros::Publisher _cmd_pub[7];
		bool _first_js;
		bool _first_fk;
		bool _sync, _first_wrench;
		KDL::FrameVel _dirkin_out;
		KDL::Frame _p_out;
		KDL::Twist _v_out;
		KDL::ChainDynParam *_dyn_param;
		geometry_msgs::PoseStamped _pose;
		geometry_msgs::TwistStamped _vel;
		Eigen::VectorXd x_t;
		Eigen::VectorXd xDot_t;
		Eigen::VectorXd xDotDot;
		Eigen::MatrixXd _J;
		Eigen::VectorXd _extWrench, _wrenchBias;
		int _wrenchCount;
		Eigen::VectorXd z_t,zDot_t,zDotDot_t;
		geometry_msgs::PoseStamped _complPose;
		geometry_msgs::TwistStamped _complVel;
		geometry_msgs::AccelStamped _complAcc;
		geometry_msgs::PoseStamped _desPose;
		geometry_msgs::TwistStamped _desVel;
		geometry_msgs::AccelStamped _desAcc;
		bool _trajEnd;
		bool _newPosReady;
		geometry_msgs::PoseStamped _nextdesPose;
		geometry_msgs::TwistStamped _nextdesVel;
		geometry_msgs::AccelStamped _nextdesAcc;
		Eigen::MatrixXd _Mt;
		Eigen::MatrixXd _Kdt;
		Eigen::MatrixXd _Kpt;
		Eigen::VectorXd xf,xf_dot,xf_dotdot;
		double _sTime,_freq;
		actionlib::SimpleActionServer<hardware_controller::waypointsAction> _kukaActionServer;
		hardware_controller::waypointsFeedback _actionFeedback;
  		hardware_controller::waypointsResult _actionResult;
		LowPassFilter lpf[6];
		double _contTime;
		diverterState _state;
		bool _firstCompliant, _mainDone, _dronePos_ready;
		Eigen::Vector3d _dronePos;
		NotchFilter notchF[6];
};

bool KUKA_INVDYN::init_robot_model() {

    if (!kdl_parser::treeFromFile("/home/eugenio/kuka_ws/src/hardware_controller/urdf/iiwa7.urdf", iiwa_tree)){
    	ROS_ERROR("Failed to construct kdl tree");
      	return false;
    }

	std::string base_link = "iiwa_link_0";
	std::string tip_link  = "iiwa_link_sensor_kuka";
	if ( !iiwa_tree.getChain(base_link, tip_link, _k_chain) ) return false;

	_fksolver = new KDL::ChainFkSolverPos_recursive( _k_chain );
	_fk_solver_pos_vel = new KDL::ChainFkSolverVel_recursive( _k_chain );
	_ik_solver_vel = new KDL::ChainIkSolverVel_pinv( _k_chain );
	_ik_solver_pos = new KDL::ChainIkSolverPos_NR( _k_chain, *_fksolver, *_ik_solver_vel, 500, 1e-6 );
	_J_solver = new KDL::ChainJntToJacSolver( _k_chain );

	_q_in = new KDL::JntArray( _k_chain.getNrOfJoints() );
	_q_out = new KDL::JntArray( _k_chain.getNrOfJoints() );
	_q_in_old = new KDL::JntArray( _k_chain.getNrOfJoints() );
	_dq_in = new KDL::JntArray( _k_chain.getNrOfJoints() );
	_initial_q = new KDL::JntArray( _k_chain.getNrOfJoints() );
	_dyn_param = new KDL::ChainDynParam(_k_chain,KDL::Vector(0,0,-9.81));

	return true;
}

KUKA_INVDYN::KUKA_INVDYN(double sampleTime) :
    _kukaActionServer(_nh, "kukaActionServer", boost::bind(&KUKA_INVDYN::actionCB, this, _1), false) {

	_sTime=sampleTime;
	_freq = 1.0/_sTime;

	if (!init_robot_model()) exit(1);
	ROS_INFO("Robot tree correctly loaded from parameter server!");

	cout << "Joints and segments: " << iiwa_tree.getNrOfJoints() << " - " << iiwa_tree.getNrOfSegments() << endl;

	_js_sub = _nh.subscribe("/iiwa/joint_states", 0, &KUKA_INVDYN::joint_states_cb, this);
	_js_pub = _nh.advertise<std_msgs::Float64MultiArray>("/iiwa/jointsCommand", 0);
	_real_wrench_sub = _nh.subscribe("/netft_data", 0, &KUKA_INVDYN::real_interaction_wrench_cb, this);
	_dronePosFb_sub = _nh.subscribe("/controller/posFeedback", 0, &KUKA_INVDYN::drone_posfb_cb, this);

	_cartpose_pub = _nh.advertise<geometry_msgs::PoseStamped>("/iiwa/eef_pose", 0);
	_cartvel_pub = _nh.advertise<geometry_msgs::TwistStamped>("/iiwa/eef_twist", 0);
	_plannedpose_pub = _nh.advertise<geometry_msgs::PoseStamped>("/iiwa/planned_pose", 0);
	_plannedtwist_pub = _nh.advertise<geometry_msgs::TwistStamped>("/iiwa/planned_twist", 0);
	_plannedacc_pub = _nh.advertise<geometry_msgs::AccelStamped>("/iiwa/planned_acc", 0);
	_plannedwrench_pub = _nh.advertise<std_msgs::Float64>("/iiwa/planned_wrench", 0);
	_desPose_pub = _nh.advertise<geometry_msgs::PoseStamped>("/iiwa/eef_des_pose", 0);
	_extWrench_pub = _nh.advertise<geometry_msgs::WrenchStamped>("/iiwa/eef_ext_wrench", 0);
	_tankEnergy_pub = _nh.advertise<std_msgs::Float64>("/iiwa/tank_energy", 0);
	_kpvalue_pub = _nh.advertise<std_msgs::Float64MultiArray>("/iiwa/admit_gains", 0);

	x_t.resize(6);
	xDot_t.resize(6);
	xDotDot.resize(6);
	_extWrench.resize(6);
	_extWrench = Eigen::VectorXd::Zero(6);
	_J.resize(6,_k_chain.getNrOfJoints());
	_J = MatrixXd::Zero(6,_k_chain.getNrOfJoints());

	z_t.resize(6); z_t=Eigen::VectorXd::Zero(6);
	zDot_t.resize(6); zDot_t=Eigen::VectorXd::Zero(6);
	zDotDot_t.resize(6); zDotDot_t=Eigen::VectorXd::Zero(6);

	xf.resize(7); xf=Eigen::VectorXd::Zero(7);
	xf_dot.resize(6); xf_dot=Eigen::VectorXd::Zero(6);
	xf_dotdot.resize(6); xf_dotdot=Eigen::VectorXd::Zero(6);

	_Mt =  1*Eigen::MatrixXd::Identity(6,6);
	_Kdt = 10*Eigen::MatrixXd::Identity(6,6);
	_Kpt = 5*Eigen::MatrixXd::Identity(6,6);

	_wrenchCount = 0;
	_wrenchBias = Eigen::VectorXd::Zero(6);

	_dronePos = Eigen::Vector3d::Zero();

	_state = NORMAL;

	_first_js = false;
	_first_fk = false;
	_trajEnd = true;
	_newPosReady = false;
	_first_wrench = false;
	_firstCompliant = false;
	_mainDone = false;

	_contTime=0;

	_kukaActionServer.start();
}

void KUKA_INVDYN::drone_posfb_cb(const std_msgs::Float64MultiArrayConstPtr& message) {
	_dronePos(0) = message->data[0];
	_dronePos(1) = message->data[1];
	_dronePos(2) = message->data[2];

	_dronePos_ready = true;
}

bool KUKA_INVDYN::getDesPose(geometry_msgs::PoseStamped& p_des) {
	if(!_first_fk) return false;

	p_des = _desPose;
	return true;
}

void KUKA_INVDYN::real_interaction_wrench_cb(const geometry_msgs::WrenchStampedConstPtr& message) {
	
	Eigen::VectorXd localWrench, outWrench;
	int nSamples = 500;
	localWrench.resize(6);
	outWrench.resize(6);

	localWrench(0)=message->wrench.force.x;
	localWrench(1)=message->wrench.force.y;
	localWrench(2)=message->wrench.force.z;
	localWrench(3)=message->wrench.torque.x;
	localWrench(4)=message->wrench.torque.y;
	localWrench(5)=message->wrench.torque.z;

	if(_wrenchCount<nSamples) {
		_wrenchBias += localWrench;
		_wrenchCount++;
		return;
	}
	else if (_wrenchCount==nSamples) {
		_wrenchBias/=_wrenchCount;
		_wrenchCount++;
		ROS_WARN("Force biased!");
		cout<<_wrenchBias.transpose()<<endl;
	}
	
	for(int i=0; i<6; i++)
		localWrench(i) -= _wrenchBias(i);

	outWrench = localWrench;

	for (int i=0; i<6; i++) {
		localWrench(i) = lpf[i].update(localWrench(i),0.002,15.0/(2.0*M_PI));
		localWrench(i) = notchF[i].update(localWrench(i));
	}

	tf::Quaternion qe(_pose.pose.orientation.x,_pose.pose.orientation.y,_pose.pose.orientation.z,_pose.pose.orientation.w);
    Vector3d pe(_pose.pose.position.x,_pose.pose.position.y,_pose.pose.position.z);
    tf::Matrix3x3 Re_tf;
    Eigen::Matrix3d Re;
    Re_tf.setRotation(qe);
    tf::matrixTFToEigen(Re_tf,Re);
    MatrixXd staticTransf(6,6);
    staticTransf << Re, MatrixXd::Zero(3,3),
                    Skew(pe)*Re, Re;
    localWrench = staticTransf*localWrench;
    outWrench = staticTransf*outWrench;
	outWrench = localWrench;
	geometry_msgs::WrenchStamped wrenchstamp;
	wrenchstamp.header.stamp = ros::Time::now();
	wrenchstamp.wrench.force.x = outWrench(0);
	wrenchstamp.wrench.force.y = outWrench(1);
	wrenchstamp.wrench.force.z = outWrench(2);
	wrenchstamp.wrench.torque.x = outWrench(3);
	wrenchstamp.wrench.torque.y = outWrench(4);
	wrenchstamp.wrench.torque.z = outWrench(5);

	_extWrench = localWrench;
	
	if(outWrench.norm()<1000)
		_extWrench_pub.publish(wrenchstamp);

	_first_wrench=true;
}

void KUKA_INVDYN::joint_states_cb( sensor_msgs::JointState js ) {
	//cout<<"Joint states"<<endl;
	_q_in_old->data=_q_in->data;

	for(int i=0; i<7; i++ ) {
		_q_in->data[i] = js.position[i];
		_dq_in->data[i] = js.velocity[i];
		if( !_first_js ) {
			_initial_q->data[i] = js.position[i];
			_q_out->data[i] = js.position[i];
		}
	}

	get_dirkin();

	_first_js = true;
	_sync = true;
}

void KUKA_INVDYN::updateState() {
	double uptresh = 0.5, imptresh = 3, timetresh = 1.0;
	bool up_cond = (_extWrench.norm()>uptresh) && (_extWrench.norm()<imptresh);
	bool down_cond = ( _extWrench(2)>(-2) ) && ( _extWrench(2)<(-1) );

	switch(_state) {
        case NORMAL:
          if( up_cond ) {
			  _contTime += 1.0/_freq;
			  ROS_WARN("Entrato");
			  if(_contTime>(timetresh)) {
				_state = HOOKED;
				ROS_WARN("STATE: HOOKED!");
				_contTime = 0;
			  }
		  } else
		  {
			  _contTime = 0;
		  }
          break;
        case HOOKED:
          if( down_cond ) {
			  _contTime += 1.0/_freq;
			  ROS_WARN("Uscendo");
			  if(_contTime>(timetresh*4)) {
				_state = DETACHED;
				ROS_WARN("STATE: DETACHED!");
				_contTime = 0;
			  }
          } else {
			  _contTime = 0;
		  }
          break;
        case DETACHED:
			if( _trajEnd ) {
			  _contTime += 1.0/_freq;
			  if(_contTime>(timetresh*4)) {
				_state = IMPACT;
				ROS_WARN("STATE: IMPACT!");
				_contTime = 0;
			  }
          	} else {
				_contTime = 0;
		  	}
        	break;
		case IMPACT:
			if( _mainDone ) {
			  _contTime += 1.0/_freq;
			  if(_contTime>(timetresh*4)) {
				_state = NORMAL;
				ROS_WARN("STATE: NORMAL!");
				_contTime = 0;
			  }
          	} else {
				_contTime = 0;
		  	}
        	break;
      }
}

void KUKA_INVDYN::ctrl_loop() {

	std_msgs::Float64 cmd[7];
	std_msgs::Float64MultiArray jcmd;
	jcmd.data.resize(7);
	KDL::JntArray q_out_new(_k_chain.getNrOfJoints());
	
	KDL::JntArray qd_out(_k_chain.getNrOfJoints());

	ros::Rate r(_freq);

	ETankGen stiffnessTank(0.15,0.01,0.15,_sTime,1);

	Eigen::MatrixXd finalKp = 1*_Kpt;
	Eigen::MatrixXd initialKp = _Kpt;
	Eigen::MatrixXd KpDot = Eigen::MatrixXd::Zero(6,6);
	Eigen::MatrixXd finalKd = 3*_Kdt; //4
	Eigen::MatrixXd initialKd = _Kdt;
	Eigen::MatrixXd KdDot = Eigen::MatrixXd::Zero(6,6);
	Eigen::MatrixXd finalM = 3*_Mt; //4
	Eigen::MatrixXd initialM = _Mt;
	Eigen::MatrixXd MDot = Eigen::MatrixXd::Zero(6,6);
	double finalT = 0.5; //transition in 0.5 seconds

	bool emergencyShut = false;
	LowPassFilter* dronepos_filter[3];
	for (int i=0; i<3; i++)
		dronepos_filter[i] = new LowPassFilter(15.0/(2.0*M_PI),(1.0/_freq));

	while( !_first_js ) usleep(0.1);
	while( !_first_wrench ) usleep(0.1);

	while( ros::ok() && (!emergencyShut)) {

    while( !_sync ) usleep(0.1);

		updatePose();

		if( (_state == HOOKED) || (_state == IMPACT) ) {
			if(_state == HOOKED) {
				finalKp = 20*Eigen::MatrixXd::Identity(6,6); 
				finalKd = 10*Eigen::MatrixXd::Identity(6,6); //2
				finalM =  2*Eigen::MatrixXd::Identity(6,6); //3
			}
			else if(_state == IMPACT) {
				finalKp = 400*Eigen::MatrixXd::Identity(6,6); 
				finalKd = 800*Eigen::MatrixXd::Identity(6,6);
				finalM =  30*Eigen::MatrixXd::Identity(6,6);
				finalKp(1,1) = 20;
				finalKd(1,1) = 30;
				finalM(1,1) = 1;
			}

			for(int i=0; i<6; i++) {
				if(_Kpt(i,i)<finalKp(i,i)) {
					KpDot(i,i) = ((finalKp(i,i)-initialKp(i,i))/finalT);
				}
				else
					KpDot(i,i) = 0;

				if(_Kdt(i,i)<finalKd(i,i)) {
					KdDot(i,i) = ((finalKd(i,i)-initialKd(i,i))/finalT);
				}
				else
					KdDot(i,i) = 0;

				if(_Mt(i,i)<finalM(i,i)) {
					MDot(i,i) = ((finalM(i,i)-initialM(i,i))/finalT);
				}
				else
					MDot(i,i) = 0;
			}
		}
		else if( _state == NORMAL || (_state == DETACHED)) {
			if(_state == DETACHED) {
				initialKp = 20*Eigen::MatrixXd::Identity(6,6);
				initialKd = 10*Eigen::MatrixXd::Identity(6,6);
				initialM =  2*Eigen::MatrixXd::Identity(6,6);
			}
			for(int i=0; i<6; i++) {
				if(_Kpt(i,i)>initialKp(i,i)) {
					KpDot(i,i) = -((finalKp(i,i)-initialKp(i,i))/finalT);
				}
				else
					KpDot(i,i) = 0;

				if(_Kdt(i,i)>initialKd(i,i)) {
					KdDot(i,i) = -((finalKd(i,i)-initialKd(i,i))/finalT);
				}
				else
					KdDot(i,i) = 0;

				if(_Mt(i,i)>initialM(i,i)) {
					MDot(i,i) = -((finalM(i,i)-initialM(i,i))/finalT);
				}
				else
					MDot(i,i) = 0;
			}
		}

		std::vector<double> tankInputs;
		std::vector<double> tankDiss;

		tankDiss.push_back(zDot_t.dot(_Kdt*zDot_t));
		if (stiffnessTank._alpha[0] != 1)
			cout<<stiffnessTank._alpha[0] << " " << _Kpt(0,0)<<endl;
		KpDot *= stiffnessTank._alpha[0];
		MDot *= stiffnessTank._alpha[0];
		tankInputs.push_back(0.5*z_t.dot(KpDot*z_t) + 0.5*zDot_t.dot(MDot*zDot_t));
		stiffnessTank.update(tankInputs,tankDiss);
	
		_Kpt += _sTime * KpDot;	
		_Mt += _sTime * MDot;
		_Kdt += _sTime * KdDot;
		if(_Kpt.norm()>finalKp.norm())
			_Kpt = finalKp;
		else if(_Kpt.norm()<initialKp.norm())
			_Kpt = initialKp;

		if(_Kdt.norm()>finalKd.norm())
			_Kdt = finalKd;
		else if(_Kdt.norm()<initialKd.norm())
			_Kdt = initialKd;

		if(_Mt.norm()>finalM.norm())
			_Mt = finalM;
		else if(_Mt.norm()<initialM.norm())
			_Mt = initialM;

		std_msgs::Float64MultiArray msg;
		msg.data.resize(3);
		msg.data[0] = _Kpt(0,0);
		msg.data[1] = _Kdt(0,0);
		msg.data[2] = _Mt(0,0);
		_kpvalue_pub.publish(msg);

		std_msgs::Float64 msgtank;
		msgtank.data = stiffnessTank.getEt();
		_tankEnergy_pub.publish(msgtank);

		_desPose_pub.publish(_desPose);
		updateState();
		compute_compliantFrame(_desPose,_desVel,_desAcc);

		_complPose.header.stamp = ros::Time::now();
		_complVel.header.stamp = _complPose.header.stamp;
		_complAcc.header.stamp = _complPose.header.stamp;
		_plannedpose_pub.publish(_complPose);
		_plannedtwist_pub.publish(_complVel);
		_plannedacc_pub.publish(_complAcc);

		KDL::Frame F_dest;
		tf::Quaternion qdes(_complPose.pose.orientation.x,_complPose.pose.orientation.y,_complPose.pose.orientation.z,_complPose.pose.orientation.w);
		tf::Matrix3x3 R(qdes);
		F_dest.M.data[0] = R[0][0];
		F_dest.M.data[1] = R[0][1];
		F_dest.M.data[2] = R[0][2];
		F_dest.M.data[3] = R[1][0];
		F_dest.M.data[4] = R[1][1];
		F_dest.M.data[5] = R[1][2];
		F_dest.M.data[6] = R[2][0];
		F_dest.M.data[7] = R[2][1];
		F_dest.M.data[8] = R[2][2];

		if(!_dronePos_ready) {
			F_dest.p.data[0] = _complPose.pose.position.x;
			F_dest.p.data[1] = _complPose.pose.position.y;
			F_dest.p.data[2] = _complPose.pose.position.z;
		} else {
			Vector3d diff = _dronePos;
			Vector3d actualPos(_pose.pose.position.x,_pose.pose.position.y,_pose.pose.position.z);
			Vector3d actualVel;
			double tresh = 0.05;//bounding box
			for(int i=0; i<3; i++) {
				if(diff(i)>tresh) diff(i) = tresh;
				else if(diff(i)<(-tresh)) diff(i) = (-tresh);
			}

			actualVel(0) =_pose.pose.position.x-_complPose.pose.position.x-diff(0);
			actualVel(1) =_pose.pose.position.y-_complPose.pose.position.y-diff(1);
			actualVel(2) =_pose.pose.position.z-_complPose.pose.position.z-diff(2);
			actualVel = _freq*actualVel;

			F_dest.p.data[0] = _complPose.pose.position.x + (dronepos_filter[0]->update(diff(0)) );
			F_dest.p.data[1] = _complPose.pose.position.y + (dronepos_filter[1]->update(diff(1)) );
			F_dest.p.data[2] = _complPose.pose.position.z + (dronepos_filter[2]->update(diff(2)) );
			
		}

		if(F_dest.p.data[1]>0.75) F_dest.p.data[1]=0.75; //workspace saturation

		if( _ik_solver_pos->CartToJnt(*_q_out, F_dest, q_out_new) != KDL::SolverI::E_NOERROR )
			cout << "failing in ik!" << endl;
		else {
			_q_out->data = q_out_new.data;
		}
		
		for(int i=0; i<7; i++) jcmd.data[i]=_q_out->data[i];
		_js_pub.publish(jcmd);

		_sync = false;

		r.sleep();
	}

}

void KUKA_INVDYN::get_dirkin() {
	KDL::JntArrayVel q_qdot(*_q_in,*_dq_in);
	_fk_solver_pos_vel->JntToCart(q_qdot, _dirkin_out);
	_p_out = _dirkin_out.GetFrame();
	_v_out = _dirkin_out.GetTwist();
	_pose.pose.position.x = _p_out.p.x();
	_pose.pose.position.y = _p_out.p.y();
	_pose.pose.position.z = _p_out.p.z();

	double qx, qy, qz, qw;
	_p_out.M.GetQuaternion( qx, qy, qz, qw);
	_pose.pose.orientation.w = qw;
	_pose.pose.orientation.x = qx;
	_pose.pose.orientation.y = qy;
	_pose.pose.orientation.z = qz;

	if(!_first_fk) _desPose = _pose;

	KDL::Jacobian Jac(_k_chain.getNrOfJoints());
	if( _J_solver->JntToJac(*_q_in, Jac) != KDL::ChainJntToJacSolver::E_NOERROR )
		cout << "failing in Jacobian computation!" << endl;

	_J = Jac.data;

	Eigen::VectorXd vel(6);
	vel = _J*(_dq_in->data);

	_vel.twist.linear.x = vel(0);
	_vel.twist.linear.y = vel(1);
	_vel.twist.linear.z = vel(2);
	_vel.twist.angular.x = vel(3);
	_vel.twist.angular.y = vel(4);
	_vel.twist.angular.z = vel(5);

	_pose.header.stamp = ros::Time::now();
	_vel.header.stamp = _pose.header.stamp;
	_cartpose_pub.publish( _pose );
	_cartvel_pub.publish( _vel );
	_first_fk = true;
}

void KUKA_INVDYN::compute_compliantFrame(const geometry_msgs::PoseStamped& p_des, const geometry_msgs::TwistStamped& v_des, const geometry_msgs::AccelStamped& a_des) {

	if(!(_extWrench.norm()<1000000))
		_extWrench = Eigen::VectorXd::Zero(6);
	zDotDot_t = _Mt.inverse() * ( _extWrench - _Kdt*zDot_t - _Kpt*z_t);
	zDotDot_t.tail(3) = Eigen::VectorXd::Zero(3);

	if(!_first_wrench) {
		zDotDot_t = Eigen::VectorXd::Zero(6);
	}
	zDot_t += zDotDot_t*_sTime;
	z_t += zDot_t*_sTime;

	_complAcc.accel.linear.x = a_des.accel.linear.x + zDotDot_t(0);
	_complAcc.accel.linear.y = a_des.accel.linear.y + zDotDot_t(1);
	_complAcc.accel.linear.z = a_des.accel.linear.z + zDotDot_t(2);
	_complAcc.accel.angular.x = a_des.accel.angular.x + zDotDot_t(3);
	_complAcc.accel.angular.y = a_des.accel.angular.y + zDotDot_t(4);
	_complAcc.accel.angular.z = a_des.accel.angular.z + zDotDot_t(5);

	_complVel.twist.linear.x = v_des.twist.linear.x + zDot_t(0);
	_complVel.twist.linear.y = v_des.twist.linear.y + zDot_t(1);
	_complVel.twist.linear.z = v_des.twist.linear.z + zDot_t(2);
	_complVel.twist.angular.x = v_des.twist.angular.x + zDot_t(3);
	_complVel.twist.angular.y = v_des.twist.angular.y + zDot_t(4);
	_complVel.twist.angular.z = v_des.twist.angular.z + zDot_t(5);

	_complPose.pose.position.x = p_des.pose.position.x + z_t(0);
	_complPose.pose.position.y = p_des.pose.position.y + z_t(1);
	_complPose.pose.position.z = p_des.pose.position.z + z_t(2);

	tf::Quaternion qe(p_des.pose.orientation.x,p_des.pose.orientation.y,p_des.pose.orientation.z,p_des.pose.orientation.w);
	tf::Quaternion qd(p_des.pose.orientation.x,p_des.pose.orientation.y,p_des.pose.orientation.z,p_des.pose.orientation.w);
	tf::Matrix3x3 Re_tf, Rd_tf;
	Eigen::Matrix3d Re,Rd;
	Re_tf.setRotation(qe);
	tf::matrixTFToEigen(Re_tf,Re);
	Eigen::Vector3d eps;
	eps << z_t(3),z_t(4),z_t(5);
	eps = Re.transpose()*eps;
	double eta = sqrt(1-eps(0)*eps(0)-eps(1)*eps(1)-eps(2)*eps(2));
	if(eta>1) eta=1;
	else if (eta<-1) eta=-1;
	double theta = 2*acos(eta);
	if(theta!=0) { //qd actually different from qe
		Eigen::Vector3d axis = (1.0/sin(theta*0.5))*eps;
		tf::Vector3 axis_tf;
		tf::vectorEigenToTF(axis,axis_tf);
		tf::Quaternion qerr(axis_tf,theta);
		tf::Matrix3x3 Rerr_tf(qerr);
		Eigen::Matrix3d Rerr;
		tf::matrixTFToEigen(Rerr_tf,Rerr);
		Rd = Re*Rerr;
		tf::matrixEigenToTF(Rd,Rd_tf);
		Rd_tf.getRotation(qd);
	}

	_complPose.pose.orientation.x = qd.x();
	_complPose.pose.orientation.y = qd.y();
	_complPose.pose.orientation.z = qd.z();
	_complPose.pose.orientation.w = qd.w();

	VectorXd vel;
	vel.resize(6);
	vel << _complVel.twist.linear.x ,
		_complVel.twist.linear.y ,
		_complVel.twist.linear.z ,
		_complVel.twist.angular.x,
		_complVel.twist.angular.y,
		_complVel.twist.angular.z;

	_firstCompliant = true;
}

void KUKA_INVDYN::updatePose() {

	if(!_newPosReady)
		return;

	_desPose=_nextdesPose;
	_desVel = _nextdesVel;
	_desAcc = _nextdesAcc;

	_newPosReady = false;
}

bool KUKA_INVDYN::newTrajectory(const std::vector<geometry_msgs::PoseStamped> waypoints, const std::vector<double> times, const Eigen::VectorXd xdi, const Eigen::VectorXd xdf, const Eigen::VectorXd xddi, const Eigen::VectorXd xddf) {
	if(!_trajEnd) return false;

	_trajEnd=false;
	CARTESIAN_PLANNER	cplanner(_freq);
	cplanner.set_waypoints(waypoints,times,xdi,xdf,xddi,xddf);
	cplanner.compute();

	int trajsize = cplanner._x.size();
	int trajpoint = 0;
	double status = 0;

	while(cplanner.isReady() && ros::ok()) {
		while(_newPosReady && ros::ok()) usleep(1);
		cplanner.getNext(_nextdesPose,_nextdesVel,_nextdesAcc);
		_newPosReady=true;
		trajpoint++;
		status = 100.0*((double)(trajpoint))/trajsize;
		if(_kukaActionServer.isActive()) {
			_actionFeedback.completePerc=status;
			_kukaActionServer.publishFeedback(_actionFeedback);
			if (_kukaActionServer.isPreemptRequested())
      			{
      			  ROS_INFO("ACTION Preempted");
      			  // set the action state to preempted
      			  _kukaActionServer.setPreempted();
							_trajEnd=true;
      			  return false;
      			}
		}
	}

	_trajEnd=true;
	return true;
}

bool KUKA_INVDYN::newTrajectory(const std::vector<geometry_msgs::PoseStamped> waypoints, const std::vector<double> times) {
	Eigen::VectorXd dummy(6);
	dummy = Eigen::VectorXd::Zero(6);
	newTrajectory(waypoints,times,dummy,dummy,dummy,dummy);
}

void KUKA_INVDYN::run() {
	boost::thread ctrl_loop_t( &KUKA_INVDYN::ctrl_loop, this);
	//ros::spin();
}

void KUKA_INVDYN::actionCB(const hardware_controller::waypointsGoalConstPtr &goal) {
	bool result;
	std::vector<geometry_msgs::PoseStamped> waypoints = goal->waypoints.poses;
	std::vector<double> times = goal->times;
	Eigen::VectorXd xdi,xdf,xddi,xddf;
	twist2Vector(goal->initVel,xdi);
	twist2Vector(goal->finalVel,xdf);
	accel2Vector(goal->initAcc,xddi);
	accel2Vector(goal->finalAcc,xddf);
	result = newTrajectory(waypoints,times,xdi,xdf,xddi,xddf);

	if(result) {
      _actionResult.ok = true;
      ROS_INFO("ACTION: Succeeded");
      // set the action state to succeeded
      _kukaActionServer.setSucceeded(_actionResult);
  } else {
		_actionResult.ok = false;
		ROS_INFO("ACTION: Aborted. Probably already following another trajectory.");
		// set the action state to succeeded
		_kukaActionServer.setAborted(_actionResult);
	}
}


int main(int argc, char** argv) {
	ros::init(argc, argv, "iiwa_kdl");

	ros::AsyncSpinner spinner(1); // Use 1 thread
	spinner.start();

	KUKA_INVDYN iiwa(0.01);
	iiwa.run();
	ros::Rate r(50);
	diverterState state,oldState;
	state = iiwa.getState();
	oldState = state;

	while(ros::ok()) { 
		state = iiwa.getState();
		if((state == DETACHED) && (oldState==HOOKED)) {
			ROS_WARN("Transition from HOOKED to DETACHED.");
			sleep(2);
			exit(0);
			std::vector<geometry_msgs::PoseStamped> waypoints;
			geometry_msgs::PoseStamped p;
			iiwa.getDesPose(p);
			waypoints.push_back(p); //Initial
			p.pose.position.x = 0.5;
			p.pose.position.y = 0.0;
			p.pose.position.z = 0.4;
  			rotateYaw(p,p,M_PI/2);
  			waypoints.push_back(p); //medio
			p.pose.position.x = -0.041;
			p.pose.position.y = 0.655;//0.65
			p.pose.position.z = 0.50-0.30;//0.50-0.30
  			tf::Quaternion qinit(0.617,0.784,0.041,-0.038);
  			qinit.normalize();
			p.pose.orientation.z = qinit.z();
			p.pose.orientation.w = qinit.w();
			p.pose.orientation.x = qinit.x();
			p.pose.orientation.y = qinit.y();
  			waypoints.push_back(p); //finale
			std::vector<double> times;
  			times.push_back(0);
  			times.push_back(12);//15
			times.push_back(22);//30
			iiwa.newTrajectory(waypoints,times);
		}
		else if((state == IMPACT) && (oldState==DETACHED)) {
			ROS_WARN("Transition from DETACHED to IMPACT.");
			iiwa.setDone(false);
			//char c;
			//cin>>c;
			sleep(4);
			if(!ros::ok()) exit(0);
			std::vector<geometry_msgs::PoseStamped> waypoints;
			geometry_msgs::PoseStamped p;
			iiwa.getDesPose(p);
			waypoints.push_back(p); //Initial
			p.pose.position.z += 0.30;//0.30
			waypoints.push_back(p); //Initial
			std::vector<double> times;
			times.push_back(0);
  			times.push_back(1);
			
			iiwa.newTrajectory(waypoints,times);
			sleep(1);
			exit(0);
			//iiwa.setDone(true);
		}

		oldState=state;
		r.sleep();
	}

	ros::waitForShutdown();

	return 0;
}
