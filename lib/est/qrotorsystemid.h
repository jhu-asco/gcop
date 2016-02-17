// This file is part of libgcop, a library for Geometric Control, Optimization, and Planning (GCOP)
//
// Copyright (C) 2004-2014 Marin Kobilarov <marin(at)jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef QROTORSYSTEMID_H
#define QROTORSYSTEMID_H
/**
 * @brief The QRotorSystemID class
 * This class does a MLE estimation of quadrotor system parameters. In particular this uses
 * the Quad rotor model and identify parameters based on measured position and rpy data.
 */

#include "qrotoridmodel.h"
#include "so3.h"
#include <Eigen/Dense>

#include "ceres/ceres.h" //Nonlinear Least Squares Optimization Package

namespace gcop {

using namespace std;
using namespace Eigen;

struct QRotorSystemIDMeasurement {
    Vector3d position;///< Position
    Vector3d rpy;///< Measured rpy
    Vector4d control;///< Input rpyt command
    double t;///< Measured time
};

class QRotorSystemID
{
public:
    QRotorSystemID();

    void EstimateParameters(const vector<QRotorSystemIDMeasurement> &inputs, QRotorIDState &initial_state);
public:
    MatrixXd offsets_matrix;///< Offsets matrix from Optimization
    Vector7d qrotor_gains;///< Parameter vector of kt, kp, kd
    Matrix<double,7,7> qrotor_gains_residualgain;///< Inverse Square root of covariance of qrotor_gains
    Vector6d offsets_prior;///< Prior on offsets
    Matrix<double,6,6> offsets_prior_residualgain;///< Inverse square root of covariance of offsets
    Vector15d stdev_initial_state_prior;///< stdeviation on initial state prior
    double stdev_position;///< Stdeviation on position measurement
    double stdev_rpy;///< Stdeviation on rpy measurement
    double offsets_timeperiod;///< Time periof for offsets
    static QRotorIDModel sys_;///< Reference system for which parameters are optimized over

protected:
    ceres::Solver::Options options;
    ceres::Solver::Summary summary;
    /**
      * Standard stereo camera residual error
      */
    struct MeasurementError {
      MeasurementError(const vector<QRotorSystemIDMeasurement> &inputs, const QRotorIDState &initial_state, const QRotorSystemID &parent_class)
        : inputs_(inputs), parent_(parent_class)
          ,so3(SO3::Instance()), initial_state_prior_(initial_state) {

          number_offsets = std::ceil((inputs_.back().t -inputs_.front().t)/parent_.offsets_timeperiod);//Offsets change every 0.1 secs
          //DEBUG:
          //cout<<"Number of Offsets: "<<number_offsets<<endl;
      }

      inline void setInitialState(double const* p, QRotorIDState &initial_state) const
      {
          initial_state.p = Vector3d(p[0], p[1], p[2]);
          initial_state.v = Vector3d(p[3],p[4],p[5]);
          Vector3d rpy(p[6],p[7],p[8]);
          so3.q2g(initial_state.R, rpy);
          initial_state.w = Vector3d(p[9],p[10],p[11]);
          initial_state.u = Vector3d(p[12],p[13],p[14]);
      }

      /**
       * @brief operator ()
       * @param p_ Input Parameters (kt, kp, kd, offsets)
       * @param res Output Residual based on measurements and parameters
       * @return  true
       */
      bool operator()(double const* const* p_, double* res) const
      {
        //Sample a trajectory using the given parameters from the system
        int number_of_controls = inputs_.size();
        const double *first_block = p_[0];// Initial State
        const double *second_block = p_[1];//kt, kp, etc
        const double *third_block = p_[2];//Offsets
        QRotorIDState state, temp_state;// Iterating state for trajectory
        setInitialState(first_block,state);// From parameters
        //DEBUG:
        //cout<<"Initial State"<<endl;
        //cout<<state.p.transpose()<<" "<<state.v.transpose()<<" "<<state.w.transpose()<<" "<<state.u.transpose()<<endl;
        //Find difference between the given state and initial state prior:
        Vector15d prior_res;
        sys_.X.Lift(prior_res,initial_state_prior_,state);
        prior_res = prior_res.cwiseQuotient(parent_.stdev_initial_state_prior);
        Matrix<double,13,1> sys_params;
        for(int i = 0; i < 7; i++)
            sys_params[i] = second_block[i];//kt, kp, kd etc
        double offsets_prevtime = inputs_.at(0).t;
        const double *offsets_ref = third_block;
        double *offsets_dest_ref = sys_params.data()+7;
        memcpy(offsets_dest_ref,offsets_ref,6*sizeof(double));
        //int offset_count = 0;//Offset being used
        //cout<<"Offset count: "<<offset_count<<endl;//DEBUG

        double * res_ref = res;//Residual Reference to Copy to
        Vector3d current_rpy;
        for(int i = 0; i < number_of_controls-1; i++)
        {
            //DEBUG
            //cout<<"Index step: "<<i<<endl;
            const QRotorSystemIDMeasurement &current_meas = inputs_.at(i);
            const QRotorSystemIDMeasurement &next_meas = inputs_.at(i+1);
            //DEBUG
            //cout<<"Current time, prev_offset_time: "<<current_meas.t<<" "<<offsets_prevtime<<endl;
            if((current_meas.t - offsets_prevtime) >= parent_.offsets_timeperiod)
            {
                offsets_ref = offsets_ref + 6;
                memcpy(offsets_dest_ref,offsets_ref,6*sizeof(double));
                offsets_prevtime = current_meas.t;
                //offset_count++;
                //cout<<"Offset count: "<<offset_count<<endl;//DEBUG
            }
            double h = next_meas.t - current_meas.t;
            sys_.Step(temp_state,current_meas.t,state,current_meas.control,h,&sys_params,0,0,0);
            state = temp_state;
            so3.g2q(current_rpy, state.R);
            //Find residual for each measurement:
            for(int j = 0; j < 3; j++)
            {
                res_ref[j] = (1/parent_.stdev_position)*(state.p[j] - next_meas.position[j]);
                double error_angle = current_rpy[j] - next_meas.rpy[j];
                error_angle = (error_angle > M_PI)?(error_angle - 2*M_PI):(error_angle < -M_PI)?(error_angle +2*M_PI):error_angle;
                res_ref[j+3] = (1/parent_.stdev_rpy)*error_angle;
            }
            res_ref = res_ref+6;
        }
        //Residual for initial state
        memcpy(res_ref,prior_res.data(),15*sizeof(double));
        res_ref = res_ref + 15;

        //Residual for gains
        Map<const Vector7d> qrotor_gains(second_block);
        Map<Vector7d> residual_gains(res_ref);
        residual_gains = parent_.qrotor_gains_residualgain*(qrotor_gains - parent_.qrotor_gains);
        res_ref = res_ref+7;

        //Offsets:
        offsets_ref = third_block;
        for(int i = 0; i < number_offsets; i++)
        {
            Map<const Vector6d> current_offsets(offsets_ref);
            Map<Vector6d> offset_res(res_ref);
            offset_res = parent_.offsets_prior_residualgain*(current_offsets - parent_.offsets_prior);
            //if(i != number_offsets)
            //{
            res_ref = res_ref + 6;//Size of offset
            offsets_ref = offsets_ref+ 6;
            //}
        }

        return true;
      }

      static ceres::CostFunction* Create(const vector<QRotorSystemIDMeasurement> &inputs, const QRotorIDState &initial_state, const QRotorSystemID &parent) {

        ceres::DynamicNumericDiffCostFunction<MeasurementError,ceres::CENTRAL> *costFunctor = new ceres::DynamicNumericDiffCostFunction<MeasurementError, ceres::CENTRAL>(
                  new MeasurementError(inputs, initial_state, parent));

        int number_offsets = std::ceil((inputs.back().t -inputs.front().t)/parent.offsets_timeperiod);//Offsets change every 0.1 secs
        int number_of_parameters = 22+6*number_offsets;
        costFunctor->AddParameterBlock(15);//15 for initial state
        costFunctor->AddParameterBlock(7);//7 for kt, kp, etc;
        costFunctor->AddParameterBlock(6*number_offsets);
        int number_of_residuals = 6*(inputs.size()-1) + number_of_parameters;
        costFunctor->SetNumResiduals(number_of_residuals);
        return costFunctor;
      }
      const vector<QRotorSystemIDMeasurement> &inputs_;
      const QRotorIDState &initial_state_prior_;
      SO3 &so3;
      int number_offsets;///<  Number of offsets
      const QRotorSystemID &parent_;
    };
};
}

#endif // QROTORSYSTEMID_H
