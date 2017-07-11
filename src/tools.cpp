#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

    Eigen::VectorXd rmse(4);
    rmse << 0.0, 0.0, 0.0, 0.0;
    
    // check the validity of the following inputs:
    //  * the estimation vector size should not be zero
    //  * the estimation vector size should equal ground truth vector size
    if(estimations.size() != ground_truth.size()
       || estimations.size() == 0){
        cout << "Invalid estimation or ground_truth data" << endl;
        return rmse;
    }
    
    //cout << "Tools step: " << estimations.size() << endl;

    //accumulate squared residuals
    for(unsigned int i=0; i < estimations.size(); ++i){
        
        VectorXd residual = estimations[i] - ground_truth[i];
        
        //coefficient-wise multiplication
        residual = residual.array()*residual.array();
        //cout << "rmse: " << rmse << endl;
        //cout << "residual: " << residual << endl;
        
        rmse += residual;
        
    }
    
    //calculate the mean
    rmse = rmse/estimations.size();
    
    //calculate the squared root
    rmse = rmse.array().sqrt();
    //if (estimations.size()==9)
    //{
       /* cout << rmse << endl;
        double x = 0.0;
        cout << "i, px, py, gtx, gty" << endl;
        for(unsigned int i=0; i < estimations.size(); ++i)
        {
            cout << i << "," << estimations[i][0] << "," << estimations[i][1] << "," << ground_truth[i][0] << "," << ground_truth[i][1] << endl;
        }
        x = 1.0; */
    //}

    return rmse;
}

        
// simplified, non-vectorized computation
VectorXd Tools::CalculateRMSE2(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    
    //
    const int numSamples = estimations.size();
    double data[4] = {0.0, 0.0, 0.0, 0.0};
    // compute the rmse of each component,
    for (int i=0; i<numSamples; ++i)
    {
        for (int j=0; j<4; ++j)
            data[j] += pow((estimations[i][j] - ground_truth[i][j]),2.0);
    }
    
    for (int j=0; j<4; ++j)
    {
        data[j] /= numSamples;
        data[j] = sqrt(data[j]);
    }
    
    Eigen::VectorXd rmse(4);
    rmse(0) = data[0];
    rmse(1) = data[1];
    rmse(2) = data[2];
    rmse(3) = data[3];
    
    //cout << rmse << endl;
    if (rmse[0]>0.15)
    {
        double x = 0.0;
        cout << "i, px, py, gtx, gty" << endl;
        for(unsigned int i=0; i < estimations.size(); ++i)
        {
            cout << i << "," << estimations[i][0] << "," << estimations[i][1] << "," << ground_truth[i][0] << "," << ground_truth[i][1] << endl;
        }
        x = 1.0;
    }
    return rmse;
}


MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
    MatrixXd Hj(3,4);
    //recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);
    
    //pre-compute a set of terms to avoid repeated calculation
    float c1 = px*px+py*py;
    float c2 = sqrt(c1);
    float c3 = (c1*c2);
    
    //check division by zero
    if(fabs(c1) < 0.0001){
        cout << "CalculateJacobian () - Error - Division by Zero" << endl;
        return Hj;
    }
    
    //compute the Jacobian matrix
    Hj << (px/c2), (py/c2), 0, 0,
		  -(py/c1), (px/c1), 0, 0,
		  py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;
    
    return Hj;
}
