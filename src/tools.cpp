#include "tools.h"
#include <iostream>
#include <assert.h>
#include <math.h>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
   VectorXd rmse(4);
   rmse << 0, 0, 0, 0;
   if (estimations.size() == 0) {
      return rmse;
   }
   if (estimations.size() != ground_truth.size()) {
      return rmse;
   }
   unsigned int vsize = estimations.size();
   for (unsigned int i = 0; i < vsize; ++i) {
      VectorXd diff = estimations[i] - ground_truth[i];
      diff  = diff.array() * diff.array();
      rmse += diff; 
   }
   rmse /= vsize;
   rmse = rmse.cwiseSqrt();
   return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
   float px = x_state[0];
   float py = x_state[1];
   float vx = x_state[2];
   float vy = x_state[3];

   MatrixXd jocobian(3,4);
   jocobian << 0,0,0,0,
               0,0,0,0,
               0,0,0,0;
   float squred_sum = px * px + py * py;
   if (squred_sum < 0.00001) {
      std::cout << "Divided by zero!" << std::endl;
      return jocobian; 
   }
   jocobian(0,0) = px / sqrt(squred_sum);
   jocobian(0,1) = py / sqrt(squred_sum);
   jocobian(1,0) = -py / squred_sum;
   jocobian(1,1) = px / squred_sum;
   jocobian(2,0) = py * (vx * py - vy * px) / pow(squred_sum, 3.0/2.0);
   jocobian(2,1) = px * (vy * px - vx * py) / pow(squred_sum, 3.0/2.0);
   jocobian(2,2) = jocobian(0, 0);
   jocobian(2,3) = jocobian(0, 1);
   return jocobian;
}

VectorXd Tools::MapPolar(const VectorXd& x_state) {
   float px = x_state[0];
   float py = x_state[1];
   float vx = x_state[2];
   float vy = x_state[3];
   float squred_sum = px * px + py * py;
   VectorXd polar_coor(3);
   polar_coor << 0, 0, 0;
   if (squred_sum < 0.00001) {
      std::cout << "Divided by zero!" << std::endl;
      return polar_coor; 
   }
   polar_coor(0) = sqrt(squred_sum);
   auto angle = atan2(py, px);
   if (angle < -M_PI) {
      angle += 2 * M_PI;
   }
   if (angle > M_PI) {
      angle -= 2 * M_PI;
   }
   polar_coor(1) = angle;
   polar_coor(2) = (px*vx + py*vy) / polar_coor(0);
   return polar_coor;
}

VectorXd Tools::MapCartesian(const VectorXd& x_polar) {
   float rho = x_polar[0];
   float phi = x_polar[1];
   while (phi < -M_PI) {
      phi += 2 * M_PI;
   }
   while (phi > M_PI) {
      phi -= 2 * M_PI;
   }
   float dot_rho = x_polar[2];
   VectorXd cartesian_coor(4);
   cartesian_coor << 0, 0, 0, 0;
   cartesian_coor[0] = rho * sin(phi);
   cartesian_coor[1] = rho * cos(phi);
   cartesian_coor[2] = dot_rho * sin(phi);
   cartesian_coor[3] = dot_rho * cos(phi);
   return cartesian_coor;
}