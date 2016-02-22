#include "qrotorsystemid.h"
#define ONEDEG M_PI/180

using namespace gcop;
using namespace std;
using namespace Eigen;

QRotorIDModel QRotorSystemID::sys_;//Static member initialization

QRotorSystemID::QRotorSystemID():qrotor_gains(7), offsets_timeperiod(0.1), CSVFormat(StreamPrecision, DontAlignCols, ", ", "\n")
{
    //Parameters for DJI
    //qrotor_gains<<sys_.kt, sys_.kp, sys_.kd;//kt, kp, kd
    qrotor_gains<<0.15,4,4,4, 10,10,10;
    //Vector7d residualgain_vector;
    Vector7d residualgain_stdev;
    residualgain_stdev<<0.2, 0.1,0.1,0.1, 0.1,0.1,0.1;
    //residualgain_vector<<0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01;
    qrotor_gains_residualgain = residualgain_stdev.cwiseInverse().asDiagonal();
    offsets_prior<<-0.4, -0.22, 0, 0, 0, -0.01;
    //offsets_prior<<sys_.a0, sys_.tau0;
    Vector6d offsets_residualgain_vector;
    //offsets_residualgain_vector<<0.1855, 0.2691, 0.2069, 0.1440, 0.0710, 0.0871;
    offsets_residualgain_vector<< 2, 2, 4, 10,10,10;
    //offsets_residualgain_vector<< 2, 2, 2, 10,10,10;
    offsets_prior_residualgain = offsets_residualgain_vector.asDiagonal();
    //stdev_initial_state_prior<<0.05,0.05,0.05, 0.05,0.05,0.05, ONEDEG,ONEDEG,ONEDEG, ONEDEG,ONEDEG,ONEDEG, ONEDEG,ONEDEG,ONEDEG;
    stdev_initial_state_prior<<0.05,0.05,0.05, 0.1,0.1,0.1, 2*ONEDEG,2*ONEDEG,2*ONEDEG, 10*ONEDEG,10*ONEDEG,10*ONEDEG, ONEDEG,ONEDEG,ONEDEG;
    stdev_position = 0.2;
    stdev_rpy = 4*ONEDEG;
    options.linear_solver_type = ceres::DENSE_SCHUR;
    options.minimizer_progress_to_stdout = true;
}

void QRotorSystemID::EstimateParameters(const vector<QRotorSystemIDMeasurement> &inputs, QRotorIDState &initial_state, Matrix7d *stdev_gains, Vector6d *offsets, Matrix6d *stdev_offsets){
  ceres::Problem problem;
  ceres::CostFunction *measurement_error_inst = MeasurementError::Create(inputs,initial_state,*this);
  //Create Prior for params:
  int number_offsets = std::ceil((inputs.back().t -inputs.front().t)/offsets_timeperiod);//Offsets change every 0.1 secs
  //DEBUG:
  cout<<"Number of inputs: "<<inputs.size()<<endl;
  cout<<"Number of Offsets: "<<number_offsets<<endl;

  //Initial State guess = initial state_prior
  double initial_state_prior[15];
  memcpy(initial_state_prior, initial_state.p.data(),3*sizeof(double));
  memcpy(initial_state_prior+3, initial_state.v.data(),3*sizeof(double));
  Vector3d rpy;
  SO3::Instance().g2q(rpy,initial_state.R);
  memcpy(initial_state_prior+6, rpy.data(),3*sizeof(double));
  memcpy(initial_state_prior+9, initial_state.w.data(),3*sizeof(double));
  memcpy(initial_state_prior+12, initial_state.u.data(),3*sizeof(double));
  //DEBUG:
  /*cout<<"Initial State: "<<endl;
  for(int i = 0; i < 15; i++)
      cout<<initial_state_prior[i]<<endl;
      */
  //Guess for gains
  double *qrotor_gains_prior_data = qrotor_gains.data();
  //Guess for Offsets:
  offsets_matrix.resize(6, number_offsets);
  offsets_matrix.colwise() = offsets_prior;
  /*cout<<"Offsets: "<<endl<<offsets_matrix<<endl;
  double *offsets_prior_guess = (double*)malloc(6*number_offsets*sizeof(double));
  double *offsets_prior_ref = offsets_prior_guess;
  for(int j = 0; j < number_offsets; j++)
  {
      memcpy(offsets_prior_ref, offsets_prior.data(),6*sizeof(double));
      offsets_prior_ref  = offsets_prior_ref+6;
  }
  */
  vector<double *> parameter_blocks;
  parameter_blocks.push_back(initial_state_prior);
  parameter_blocks.push_back(qrotor_gains_prior_data);
  parameter_blocks.push_back(offsets_matrix.data());

  problem.AddResidualBlock(measurement_error_inst,NULL,parameter_blocks);

  ceres::Solve(options,&problem,&summary);

  std::cout << summary.FullReport() << "\n";
  std::cout<<"Important Params: "<<std::endl;
  std::cout<<qrotor_gains.transpose().format(CSVFormat)<<std::endl;
  //cout<<"Offsets: "<<endl<<offsets_matrix<<endl;
  //Mean of Offsets:
  //Map<MatrixXd> offsets_matrix(offsets_prior_guess,number_offsets,6);
  if(offsets!= NULL)
  {
    Vector6d mean_offsets = offsets_matrix.rowwise().mean().transpose();
    std::cout<<"Mean offsets: "<<mean_offsets.transpose().format(CSVFormat)<<std::endl;
    *offsets  = mean_offsets;
    if(stdev_offsets != NULL)
    {
      MatrixXd centered = offsets_matrix.colwise() - mean_offsets;
      Matrix6d cov = (centered * centered.adjoint()) / double(offsets_matrix.cols() - 1);
      //cout<<"Covariance Offsets: "<<endl;
      //cout<<cov<<endl;
      //Compute Stdev
      SelfAdjointEigenSolver<Matrix6d> es(cov);
      *stdev_offsets = es.operatorSqrt();
      cout<<"Stdev Offsets"<<endl<<(*stdev_offsets).format(CSVFormat)<<endl;
    }
  }
    //Evaluate Covariance of computed parameters TODO
    vector<pair<const double *,const double *> >covariance_blocks;
    covariance_blocks.push_back(make_pair(parameter_blocks[1], parameter_blocks[1]));//Corresponding to Gains
    ceres::Covariance covariance(cov_options);
    bool cov_res = covariance.Compute(covariance_blocks,&problem);
    if(!cov_res)
    {
      std::cout<<"Covariance Not Found"<<std::endl;
    }
    else
    {
        Matrix<double,7,7,RowMajor> cov_gains;
        covariance.GetCovarianceBlock(parameter_blocks[1],parameter_blocks[1],cov_gains.data());
        //cout<<"Covariance of Gains: "<<endl<<cov_gains<<endl;
        SelfAdjointEigenSolver<Matrix7d> es(cov_gains);
        if(stdev_gains != NULL)
        {
          (*stdev_gains) = es.operatorSqrt();
          cout<<"Stdev_Gains: "<<endl<<(*stdev_gains).format(CSVFormat)<<endl;
          qrotor_gains_residualgain = (*stdev_gains).inverse();//Inverse of stdev is the gains
        }
        else
        {
            qrotor_gains_residualgain = es.operatorInverseSqrt();
        }
        cout<<"Qrotor Residual Gains"<<endl<<qrotor_gains_residualgain<<endl;
    }
}
