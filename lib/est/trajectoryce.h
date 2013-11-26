#ifndef DGC_TRAJECTORYCE_H
#define DGC_TRAJECTORYCE_H

#include "ce.h"
#include "trajectory.h"
#include <vector>

namespace dgc {
  
  class TrajectoryCe : public Ce {
  public:    
    
    //    TrajectoryCe(const std::vector<Trajectory*> &trajs, 
    //                 int k,
    //                 const double *inj = 0);
    
    
    TrajectoryCe(const System &sys,                 
                 int sn,
                 int k,
                 int ni = 0, const int *is = 0,
                 const double *inj = 0);
    
    virtual ~TrajectoryCe();
    
    void TrajToVec(double *z, const Trajectory &traj);
    void VecToTraj(Trajectory &traj, const double *z);
    
    void AddSample(const Trajectory &traj, double c);
    
    double Draw(Trajectory &traj);
    
    void States(std::vector<State>& states, std::vector<double>& ps, int N, int sn = 0);

    State si;
    State sf;
    int sn;

    int ni;    ///< number of indices into the state vector
    int *is;   ///< indices into the state vector

    bool pad;    ///< use padding in the L^2 metric (false by default)

      // bool ext;    ///< extended state space
  };
}


#endif
