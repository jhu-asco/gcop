#ifndef GCOP_FILTER_H
#define GCOP_FILTER_H

#include <cstring>
#include "model.h"

namespace gcop {

  class Filter {
  public:
    /**
     * Basic Filter
     * @param model system model
     */
    Filter(Model &model);
        
    virtual ~Filter();
        
    
    /**
     * Prediction step. The base implementation does not change the state.
     * @param u control (optional - no control by default)
     * @param cov whether to update the covariance as well (true by default)
     * @return predicted state
     */
    virtual bool Predict(const Eigen::VectorXd *u = 0,
                         bool cov = true) = 0;
    
    /**
     * Correction step. The base implementation does not change the state.
     * @param z measurement
     * @param cov whether to update the covariance as well (true by default)
     * @return corrected state
     */
    virtual bool Update(const Eigen::VectorXd &z, 
                        bool cov = true) = 0;
        
    /**
     * Set measurement Chi Gate
     * @param chi squared chi bound
     */
    void SetChiGate(double chiGate) { this->chiGate = chiGate; };
    
    /**
     * Get the preset measurement Chi Gate value
     * @return chi squared chi bound
     */
    double GetChiGate() const {return chiGate; }
    
    /**
     * Returns the internally computed chi-value during the call to Correct
     * @return the chi-value of the most recent measurement
     */
    double GetChi() const {return chi; }
        
    Model &model;           ///< model
    
    Eigen::VectorXd x;            ///< state mean
    Eigen::MatrixXd P;            ///< state covariance
    
    double chi;             ///< computed chi-estimate (computed only if chiGate > 0)
    double chiGate;         ///< chi-gate (default: zero or negative to ignore)

  protected:

    /**
     * Init the filter
     */  
    virtual void Init();
  };
}


#endif
