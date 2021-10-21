#include <ros/ros.h>
#include <gazebo/gazebo.hh>
#include <gazebo/physics/physics.hh>
#include <gazebo/common/Plugin.hh>
#include "geometry_msgs/WrenchStamped.h"
#include "geometry_msgs/PoseStamped.h"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <tf/tf.h>
#include <tf_conversions/tf_eigen.h>

using namespace Eigen;
using namespace std;

namespace gazebo
{
  class ConnectionPlugin : public WorldPlugin
  {
	private: ros::NodeHandle* _nh;
  private:
     ros::Publisher _force_pub;
     ros::Subscriber _force_sub;
     ros::Subscriber _pose_sub;
     ros::Subscriber _dronepose_sub;
     Matrix3d _Rsensor, _Rb;
     event::ConnectionPtr updateConnection;
     physics::WorldPtr _world;
     physics::ModelPtr _model;
    gazebo::physics::LinkPtr _link;
    bool _firstcb, _firstcb_drone;
  public:
   void pose_cb( geometry_msgs::PoseStampedConstPtr );
   void dronepose_cb( geometry_msgs::PoseConstPtr );
   void force_cb( geometry_msgs::WrenchStampedConstPtr );
   void compute_cubic();
   void ctrl_loop();
   void OnUpdate();
	 void Load(physics::WorldPtr world, sdf::ElementPtr _sdf) {
		_nh = new ros::NodeHandle();
    _world = world;
    ROS_WARN("Connection force plugin started!");
    _pose_sub = _nh->subscribe("/arm/eef_pose", 0, &ConnectionPlugin::pose_cb, this);
    _dronepose_sub = _nh->subscribe("/hummingbird/ground_truth/pose", 0, &ConnectionPlugin::dronepose_cb, this);
    _force_sub = _nh->subscribe("/iiwa/eef_ext_wrench", 0, &ConnectionPlugin::force_cb, this);
    _force_pub = _nh->advertise<geometry_msgs::WrenchStamped>("/arm/meas_ext_wrench", 0);

    _firstcb=_firstcb_drone=false;
    ROS_WARN("Connection plugin started!");
	 }

  };

 void ConnectionPlugin::force_cb( geometry_msgs::WrenchStampedConstPtr message) {
  Eigen::VectorXd localWrench;
	localWrench.resize(6);

	localWrench(0)=message->wrench.force.x;
	localWrench(1)=message->wrench.force.y;
	localWrench(2)=message->wrench.force.z;
	localWrench(3)=message->wrench.torque.x;
	localWrench(4)=message->wrench.torque.y;
	localWrench(5)=message->wrench.torque.z;

  if(_firstcb && _firstcb_drone) {
    Eigen::Vector3d ef = _Rsensor.transpose() * _Rb.transpose() * localWrench.head(3); //local
    Eigen::Vector3d ef_global = localWrench.head(3); //global
    Eigen::Vector3d et_global = localWrench.tail(3); //global
    ignition::math::Vector3d link_force(ef_global(0),ef_global(1),ef_global(2));
    ignition::math::Vector3d link_torque(et_global(0),et_global(1),et_global(2));
    _link->AddForce(link_force); //for global force
    _link->AddTorque(link_torque); //for global force
    geometry_msgs::WrenchStamped out;
    out.wrench.force.x = ef(0);
    out.wrench.force.y = ef(1);
    out.wrench.force.z = ef(2);
    out.wrench.torque.x = 0;
		out.wrench.torque.y = 0;
		out.wrench.torque.z = 0;
    out.header.stamp = ros::Time::now();
    _force_pub.publish(out);
  }
 }

  void ConnectionPlugin::dronepose_cb( geometry_msgs::PoseConstPtr pos ) {
    tf::Quaternion qsensor(pos->orientation.x,pos->orientation.y,pos->orientation.z,pos->orientation.w);
    tf::Matrix3x3 Rsensor(qsensor);
    tf::matrixTFToEigen(Rsensor,_Rb);

    if(!_firstcb) {
      // Get model
      _model = _world->ModelByName("hummingbird");
      if(_model == NULL) {
        ROS_WARN("Model not found!");
        return;
      }
      _link = _model->GetChildLink("Roll_1");
      if(_link == NULL) {
        ROS_WARN("Link not found!");
        return;
      }

      ROS_WARN("Link obtained!");
    }
    _firstcb=true;
  }

  void ConnectionPlugin::pose_cb( geometry_msgs::PoseStampedConstPtr pos ) {
      tf::Quaternion qsensor(pos->pose.orientation.x,pos->pose.orientation.y,pos->pose.orientation.z,pos->pose.orientation.w);
      tf::Matrix3x3 Rsensor(qsensor);
      tf::matrixTFToEigen(Rsensor,_Rsensor);

      _firstcb_drone = true;
  }

  // Register this plugin with the simulator
  GZ_REGISTER_WORLD_PLUGIN(ConnectionPlugin)
}