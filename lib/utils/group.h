#ifndef GCOP_GROUP_H
#define GCOP_GROUP_H

#include <Eigen/Dense>
#include "manifold.h"

namespace gcop {
  
  using namespace Eigen;

  /**
   * General definition of a matrix Lie group and its algebra.
   * Defines basic operations such as identity, inverse, 
   * multiplication, maps and flows from algebra to the group using
   * exponentiation, adjoint transform and operators. This should
   * serve as a base class for implementing specific groups.
   *
   * A Lie group is also a homogeneous manifold, so it subclasses Manifold
   * and provides a default implementation of the operations 
   * Manifold::Lift and Manifold::Retract
   *
   * Subclasses should provide all methods (inv, cay, cayinv are provided for convenience
   * but they can be implemented more efficiently for specific groups)
   *
   * @param n manifold dimension, i.e. algebra alements are vectors \f$ v\in R^n \f$
   * @param m matrix dimension, group elements are m-by-m matrices \f$ g \in G \subset GL(m) \f$
   *
   * Author: Marin Kobilarov -- Copyright (C) 2007
   * 
   */
  template <int n, int m>class Group : Manifold<Matrix<double, m, m>, n > {
    
    typedef Matrix<double, n, 1> Vectornd;
    typedef Matrix<double, m, m> Matrixmd;
    typedef Matrix<double, n, n> Matrixnd;
    
  public:
    /**
     * Compute the inverse of a given group element 
     * @param gi inverse \f$ g^{-1} \in G \f$
     * @param g group element \f$ g \in G \f$
     */   
    virtual void inv(Matrixmd &gi, const Matrixmd &g) const;

    /**
     * Compute the inverse of a given group element 
     * @return inverse \f$ g^{-1} \in G \f$
     * @param g group element \f$ g \in G \f$
     */   
    virtual Matrixmd inv(const Matrixmd &g) const;


    /**
     * Compute the difference b/n two elements \f$ \Delta g=g_a^{-1}g_b \f$.
     * @param dg resulting group element \f$ \Delta g \in G \f$
     * @param ga first group element \f$ g_a \in G \f$
     * @param gb second group element \f$ g_b \in G \f$
     */
    virtual void diff(Matrixmd &dg, const Matrixmd &ga, const Matrixmd &gb) const;


    /**
     * Compute the difference b/n two elements \f$ \Delta g=g_a^{-1}g_b \f$.
     * @return resulting group element \f$ \Delta g \in G \f$
     * @param ga first group element \f$ g_a \in G \f$
     * @param gb second group element \f$ g_b \in G \f$
     */
    virtual Matrixmd diff(const Matrixmd &ga, const Matrixmd &gb) const;


    /**
     * The operator \f$\hat v\f$ such that $\exp(v) \approx I + \hat v$
     * 
     * @param g resulting map
     * @param v given algebra element
     */    
    virtual void hat(Matrixmd &g, const Vectornd &v) const = 0;


    /**
     * The operator \f$\hat v\f$ such that $\exp(v) \approx I + \hat v$
     * 
     * @return resulting map
     * @param v given algebra element
     */    
    virtual Matrixmd hat(const Vectornd &v) const;


    /**
     * Inverse of hat operator
     * @param v resulting algebra element 
     * @param g infinetisamal group element
     */    
    virtual void hatinv(Vectornd &v, const Matrixmd &g) const = 0;


    /**
     * Inverse of hat operator
     * @return resulting algebra element 
     * @param g infinetisamal group element
     */    
    virtual Vectornd hatinv(const Matrixmd &g) const;


    /**
     * Tangent action
     * @param m resulting action as a matrix
     * @param g infinetisamal group element
     */    
    virtual void Tg(Matrixnd &M, const Matrixmd &g) const = 0;
    

    /**
     * Tangent action
     * @return resulting tangent action as a matrix
     * @param g infinetisamal group element
     */    
    virtual Matrixnd Tg(const Matrixmd &g) const;

    /**
     * The operator \f$\operatorname{Ad}_g:\mathfrak{g}\rightarrow\mathfrak{g}\f$. 
     * One can think of this operator as a change-of-frame transformation.
     * @param m resulting transformation
     * @param g group element corresponding to frame change
     */
    virtual void Ad(Matrixnd &M, const Matrixmd &g) const = 0;

    /**
     * The operator \f$\operatorname{Ad}_g:\mathfrak{g}\rightarrow\mathfrak{g}\f$. 
     * One can think of this operator as a change-of-frame transformation.
     * @return resulting transformation
     * @param g group element corresponding to frame change
     */
    virtual Matrixnd Ad(const Matrixmd &g) const;
    

    /**
     * The operator \f$\operatorname{ad}_\xi:\mathfrak{g}\rightarrow\mathfrak{g}\f$. 
     * One can think of this operator as the infinitesimal of
     * the change-of-frame transformation. It can also be written using brackets [*,*].
     * @param m resulting map
     * @param a argument
     */    
    virtual void ad(Matrixnd &M, const Vectornd &v) const = 0;


    /**
     * The operator \f$\operatorname{ad}_\xi:\mathfrak{g}\rightarrow\mathfrak{g}\f$. 
     * One can think of this operator as the infinitesimal of
     * the change-of-frame transformation. It can also be written using brackets [*,*].
     * @return resulting map
     * @param a argument
     */    
    virtual Matrixnd ad(const Vectornd &v) const;
    

    /**
     * Inverse of ad operator
     * @param v resulting algebra element 
     * @param m ad operator matrix
     */    
    virtual void adinv(Vectornd &v, const Matrixnd &M) const = 0;


    /**
     * Inverse of ad operator
     * @return resulting algebra element 
     * @param m ad operator matrix
     */    
    virtual Vectornd adinv(const Matrixnd &M) const;


    /**
     * Exponentiate lie algebra element a \f$ g=\exp(a) \f$.
     * Subclasses should implement this method
     * @param g resulting group element \f$ g \in G \f$
     * @param v Lie algebra element \f$ v\in\mathfrak{g} \f$
     */
    virtual void exp(Matrixmd &g, const Vectornd &v) const = 0;


    /**
     * Exponentiate lie algebra element a \f$ g=\exp(a) \f$.
     * Subclasses should implement this method
     * @return resulting group element \f$ g \in G \f$
     * @param v Lie algebra element \f$ v\in\mathfrak{g} \f$
     */
    virtual Matrixmd exp(const Vectornd &v) const;


    /**
     * Logarithm of a group element \f$ a = \log(g) \f$.
     * Subclasses should implement this method
     * @param v Lie algebra element \f$ v\in\mathfrak{g} \f$
     * @param g group element \f$ g \in G \f$
     */
    virtual void log(Vectornd &v, const Matrixmd &g) const = 0;


    /**
     * Logarithm of a group element \f$ a = \log(g) \f$.
     * Subclasses should implement this method
     * @return v Lie algebra element \f$ v\in\mathfrak{g} \f$
     * @param g group element \f$ g \in G \f$
     */
    virtual Vectornd log(const Matrixmd &g) const;


    /**
     * Cayley map of a lie algebra element a \f$ g=\text{cay}(a) \f$.
     * Subclasses should implement this method
     * @param g resulting group element \f$ g \in G \f$
     * @param v Lie algebra element \f$ v\in\mathfrak{g} \f$
     */
    virtual void cay(Matrixmd &g, const Vectornd &v) const;


    /**
     * Cayley map of a lie algebra element a \f$ g=\text{cay}(a) \f$.
     * Subclasses should implement this method
     * @return resulting group element \f$ g \in G \f$
     * @param v Lie algebra element \f$ v\in\mathfrak{g} \f$
     */
    virtual Matrixmd cay(const Vectornd &v) const;
 

    /**
     * Inverse of cayley\f$ v=\text{cay}^{-1}(g) \f$.
     * Subclasses should implement this method
     * @param v Lie algebra element \f$ v\in\mathfrak{g} \f$
     * @param g group element \f$ g \in G \f$
     */
    virtual void cayinv(Vectornd &v, const Matrixmd &g) const;

    /**
     * Inverse of cayley\f$ v=\text{cay}^{-1}(g) \f$.
     * Subclasses should implement this method
     * @return v Lie algebra element \f$ v\in\mathfrak{g} \f$
     * @param g group element \f$ g \in G \f$
     */
    virtual Vectornd cayinv(const Matrixmd &g) const;


    /**
     * Right-trivialized tangent of cayley
     * @param m resulting linear map
     * @param v Lie algebra element \f$ v\in\mathfrak{g}\f$
     */
    virtual void dcay(Matrixnd &M, const Vectornd &v) const = 0;


    /**
     * Right-trivialized tangent of cayley
     * @return resulting linear map
     * @param v Lie algebra element \f$ v\in\mathfrak{g}\f$
     */
    virtual Matrixnd dcay(const Vectornd &v) const;


    /**
     * Right-trivialized tangent inverse of cayley
     * @param m resulting linear map
     * @param v Lie algebra element \f$ v\in\mathfrak{g}\f$
     */
    virtual void dcayinv(Matrixnd &M, const Vectornd &v) const = 0;


    /**
     * Right-trivialized tangent inverse of cayley
     * @return resulting linear map
     * @param v Lie algebra element \f$ v\in\mathfrak{g}\f$
     */
    virtual Matrixnd dcayinv(const Vectornd &v) const;


    virtual void Lift(Vectornd &v,
                      const Matrixmd &ga,
                      const Matrixmd &gb);
    
    
    virtual void Retract(Matrixmd &gb,
                         const Matrixmd &ga,
                         const Vectornd &v);


    virtual void dtau(Matrixnd &M, const Vectornd &v);

    virtual void Adtau(Matrixnd &M, const Vectornd &v);


    static const Matrixmd Id;
    static const Vectornd e;    
  };


  template <int n, int m> 
    Matrix<double, m, m> const Group<n, m>::Id = Matrix<double, m, m>::Identity();
  
  template <int n, int m> 
    Matrix<double, n, 1> const Group<n, m>::e = Matrix<double, n, 1>::Zero();
  
  
  template <int n, int m> 
    void Group<n, m>::inv(Matrix<double, m, m> &gi, const Matrix<double, m, m> &g) const
    {
      // default uses regular matrix inverse
      // methods typically override this
      gi = g.inverse();
    }

  template <int n, int m> 
    Matrix<double, m, m> Group<n, m>::inv(const Matrix<double, m, m> &g) const
    {
      // default uses regular matrix inverse
      // methods typically override this
      Matrix<double, m, m> gi;
      inv(gi, g);
      return gi;
    }
  
  
  template <int n, int m> 
    void Group<n, m>::diff(Matrix<double, m, m> &dg, const Matrix<double, m, m> &ga, const Matrix<double, m, m> &gb) const
    {
      Matrix<double, m, m> gai;
      inv(gai, ga);
      dg = gai*gb;
    }

  template <int n, int m> 
    Matrix<double, m, m> Group<n, m>::diff(const Matrix<double, m, m> &ga, const Matrix<double, m, m> &gb) const
    {
      Matrix<double, m, m> dg;
      diff(dg, ga, gb);
      return dg;
    }
  
  template <int n, int m> 
    Matrix<double, m, m> Group<n, m>::hat(const Vectornd &v) const
    {
      Matrix<double, m, m> g;
      hat(g, v);
      return g;
    }

  template <int n, int m> 
    Matrix<double, n, 1> Group<n, m>::hatinv(const Matrixmd &g) const
    {
      Matrix<double, n, 1> v;
      hatinv(v, g);
      return v;
    }


  template <int n, int m> 
    Matrix<double, n, n> Group<n, m>::Tg(const Matrixmd &g) const
    {
      Matrixnd M;
      Tg(M, g);
      return M;
    }
  

  /*  
  template <int n, int m> 
    void Group<T, n>::exp(Matrix<double, m, m> &g, const MatrixNd &v) const
    {
      // by default only provide a first-order approximation
      // this function is typically overriden
      g = Id + hat(v);
    }
  */

  template <int n, int m> 
    Matrix<double, m, m> Group<n, m>::exp(const Matrix<double, n, 1> &v) const
    {
      Matrix<double, m, m> g;
      exp(g, v);
      return g;
    }


  /*
  template <int n, int m> 
    void Group<T, n>::log(MatrixNd &v, const Matrix<double, m, m> &g) const
    {
      // by default only provide a first-order approximation
      // this function is typically overriden
      v = hatinv(g - Id);
    }
  */

  template <int n, int m> 
    Matrix<double, n, 1> Group<n, m>::log(const Matrix<double, m, m> &g) const
    {
      Matrix<double, n, 1> v;
      log(v, g);
      return v;
    }

  
  template <int n, int m> 
    void Group<n, m>::cay(Matrix<double, m, m> &g, const Matrix<double, n, 1> &v) const
    {
      Matrix<double, m, m> vhh;
      hat(vhh, v/2);
      g = (Id - vhh).inverse()*(Id + vhh);
    }


  template <int n, int m> 
    Matrix<double, m, m> Group<n, m>::cay(const Matrix<double, n, 1> &v) const
    {
      Matrix<double, m, m> vhh;
      hat(vhh, v/2);
      return (Id - vhh).inverse()*(Id + vhh);
    }

  
  template <int n, int m> 
    void Group<n, m>::cayinv(Matrix<double, n, 1> &v, const Matrix<double, m, m> &g) const
    {
      hatinv(v, (Id + g).inverse()*(Id - g));
      v = -2*v;
    }

  template <int n, int m> 
    Matrix<double, n, 1> Group<n, m>::cayinv(const Matrix<double, m, m> &g) const
    {
      return -2*hatinv((Id + g).inverse()*(Id - g));
    }


  template <int n, int m> 
    Matrix<double, n, n> Group<n, m>::ad(const Matrix<double, n, 1> &v) const
    {
      Matrix<double, n, n> M;
      ad(M, v);
      return M;
    }


  template <int n, int m> 
    Matrix<double, n, 1> Group<n, m>::adinv(const Matrix<double, n, n> &M) const
    {
      Matrix<double, n, 1> v;
      adinv(v, M);
      return v;
    }


  template <int n, int m> 
    Matrix<double, n, n> Group<n, m>::Ad(const Matrix<double, m, m> &g) const
    {
      Matrix<double, n, n> M;
      Ad(M, g);
      return M;
    }


  template <int n, int m> 
    Matrix<double, n, n> Group<n, m>::dcay(const Matrix<double, n, 1> &v) const
    {
      Matrix<double, n, n> M;
      dcay(M, v);
      return M;
    }

  template <int n, int m> 
    Matrix<double, n, n> Group<n, m>::dcayinv(const Matrix<double, n, 1> &v) const
    {
      Matrix<double, n, n> M;
      dcayinv(M, v);
      return M;
    }

  template <int n, int m> 
    void Group<n, m>::Lift(Matrix<double, n, 1> &v,
                           const Matrix<double, m, m> &ga,                           
                           const Matrix<double, m, m> &gb) {
    
    Matrix<double, m, m> gai;
    inv(gai, ga);
    log(v, gai*gb);
  }


  template <int n, int m> 
    void Group<n, m>::Retract(Matrix<double, m, m> &gb,
                              const Matrix<double, m, m> &ga,
                              const Matrix<double, n, 1> &v) {
    
    Matrix<double, m, m> g;
    cay(g, v);
    gb = ga*g;
  }

  template <int n, int m> 
    void Group<n, m>::dtau(Matrixnd &M, const Vectornd &v) {
    dcay(M, v);
  }
  
  template <int n, int m> 
    void Group<n, m>::Adtau(Matrixnd &M, const Vectornd &v) {
    Matrix<double, m, m> g;
    cay(g, v);
    Ad(M, g);
  }



}

#endif

    
