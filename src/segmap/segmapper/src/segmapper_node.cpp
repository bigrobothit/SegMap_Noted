#include <thread>

#include <ros/ros.h>

#include "segmapper/segmapper.hpp"

/// \brief 程序入口
int main(int argc, char **argv)
{
  ros::init(argc, argv, "SegMapper");
  ros::NodeHandle node_handle("~");

  // 一个好的程序架构就是所有的启动程序统一到一个类里面
  // 然后在main函数中通过开启线程来启动整个系统
  // 并通过线程间用指针交互数据 eg:ORB-SLAM2
  SegMapper mapper(node_handle);

  std::thread publish_map_thread(&SegMapper::publishMapThread, &mapper);
  std::thread publish_tf_thread(&SegMapper::publishTfThread, &mapper);
  std::thread segmatch_thread(&SegMapper::segMatchThread, &mapper);

  try
  {
    ros::spin();
  }
  catch (const std::exception &e)
  {
    ROS_ERROR_STREAM("Exception: " << e.what());
    return 1;
  }
  catch (...)
  {
    ROS_ERROR_STREAM("Unknown Exception");
    return 1;
  }

  publish_map_thread.join();
  publish_tf_thread.join();
  segmatch_thread.join();

  return 0;
}
