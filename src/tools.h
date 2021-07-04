#ifndef TOOLS_H_
#define TOOLS_H_

#include <vector>
#include "Eigen/Dense"

class Tools {
 public:
  /**
   * Constructor.
   */
  Tools();

  /**
   * Destructor.
   */
  virtual ~Tools();

  /**
   * A helper method to calculate RMSE.
   */
  static Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, 
                                const std::vector<Eigen::VectorXd> &ground_truth);

  /**
   * A helper method to calculate Jacobians.
   */
  static Eigen::MatrixXd CalculateJacobian(const Eigen::VectorXd& x_state);

  /**
   * A helper method to map Cartesian to polar coordinates .
   */

  static Eigen::VectorXd MapPolar(const Eigen::VectorXd& x_state);

  static Eigen::VectorXd MapCartesian(const Eigen::VectorXd& x_polar);
};

#endif  // TOOLS_H_
