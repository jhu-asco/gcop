#ifndef GCOP_QUAT_H
#define GCOP_QUAT_H

#include <stdio.h>
#include <iostream> 

namespace gcop {

/**
 * The class defines a quaternion and its basic set of operations.  
 * This is a minimalistic implementation that contains only 5 doubles: 
 * 4 numbers describing the quaternion coordinates: \f$ q \in S^3 \subset \mathbb{R}^4 \f$ and a number tol
 * that defines the numerical tolerance for quaternion operations. The implemented
 * functionality is suitable for efficient quaternion algebra; for various 
 * actions on vectors in \f$ \mathbb{R}^3\f$; for transforming b/n 
 * quaternion, matrix, and exponential representations
 * 
 * Author: Marin Kobilarov -- (C) 2003
 */
class Quat {
public:
  /**
   * Initialize identity quaternion
   */
  Quat();


  /**
   * Copy contructor
   * @param q quaternion 
   */
  Quat(const Quat &q);


  /**
   * Initialize using a double array
   * @param q quaternion vector
   */
  Quat(const double q[4]);

  
  /**
   * Initialize using seperate quaternion components \f$ q=[w,x,y,z] \f$
   * @param w w
   * @param x x
   * @param y y
   * @param z z 
   */
  Quat(double w, double x, double y, double z);


  virtual ~Quat();


  /**
   * Rotate a vector v using around axis u by angle w and store resulting
   * vector in c
   * @param c resulting vector
   * @param v initial vector
   * @param u axis of rotation
   * @param w angle of rotation
   * @return a pointer to c
   */
  static double* Rotate(double c[3], 
                        const double v[3], 
                        const double u[3], 
                        double w);

  double* Rotate2(double c[3], const double a[3]) const;
  double* Rotate2(double c[3]) const;


  /**
   * Rotate a vector v and store in c
   * @param c resulting vector
   * @param v initial vector
   * @return a pointer to c
   */  
  double* Rotate(double c[3], const double v[3]) const;


  /**
   * Rotate a vector v and store the resulting vector back in v
   * @param v vector to rotate
   * @return pointer to v
   */
  double* Rotate(double v[3]) const;


  /**
   * Transform this quaternion by exponential coordinates e. 
   * Note: norm(e)<pi so this is suitable for small rotations e
   * @param e exponential coordinates
   */
  void Transform(const double e[3]);


  /**
   * Transform this quaternion by axis u and angle v
   * @param u axis of rotation
   * @param v angle of rotation
   */
  void Transform(const double u[3], double v);


  /**
   * Invert the quaternion
   */
  void Invert();


  /**
   * Create an inverse of this quaternion and store it in q
   * @param q inverse of this quaternion
   */
  void Invert(Quat &q) const;


  /**
   * Set this quaternion to the identity
   */
  void Identity();


  /**
   * Normalize \f$ q = \frac{q}{\|q\|} \f$
   */
  void Normalize();


  /**
   * Multiply with another quaternion
   * @param q other quaterion
   * @return resulting quaternion
   */
  Quat& operator*=(const Quat &q);  


  /**
   * Create from exponential coordinates
   * @param e exponential coordinates
   */
  void FromExp(const double e[3]);
  

  /**
   * To exponential coordinates
   * @param e exponential coordinates
   */
  void ToExp(double e[3]) const;


  /**
   * Create from angle-axis representation
   * @param u axis (unit vector)
   * @param v angle
   */
  void FromAxis(const double u[3], const double v);


  /**
   * To angle-axis (or exponential/canonical) representation
   * @param u resulting axis (unit vector)
   * @param v resulting angle (optional)   
   */
  void ToAxis(double u[3], double &v) const;


  /**
   * Extract this quaternion from an SE(3) matrix
   * @param m 4x4 matrix
   */
  void FromSE3(const double m[16]);


  /**
   * Create a 4x4 SE(3) matrix from this quaternion
   * @param m matrix
   */
  void ToSE3(double m[16]) const;


  /**
   * From roll-pitch-yaw angles
   * @param rpy roll-pitch-yaw angles
   */
  void FromRpy(const double rpy[3]);


  /**
   * To roll-pitch-yaw angles
   * @param rpy roll-pitch-yaw angles
   */
  void ToRpy(double rpy[3]) const;

  
  /**
   * Get a 4-dim array with the quaternion components
   * @param q the array to be filled 
   * @return pointer to q
   */
  double* Q(double q[4]) const;


  /**
   * Check if it's identity up to defined tolerance tol
   */
  bool IsIdentity() const;


  /**
   * Set the operations tolerance (default is 1e-10)
   * @param tol tolerance
   */
  void SetTol(double tol);

  double qw, qx, qy, qz;  ///< quaternion components

 protected:
  double tol;             ///< operations tolerance (default is 1e-10)

 private:
  friend const Quat operator*(const Quat& q0, const Quat& q1);
  friend std::ostream& operator<<(std::ostream &os, const Quat &q);
  friend std::istream& operator>>(std::istream &is, Quat &q);
};
 
/**
 * Multiply two quaternions and return the result \f$ q=q_0*q_1 \f$
 * @param q0 first quat
 * @param q1 second quat
 * @return their product
 */
 const Quat operator*(const Quat& q0, const Quat& q1);
 
 /**
  * Print a quaternion
  * @param os output stream
  * @param q quaternion
  * @return the output stream
  */
 std::ostream& operator<<(std::ostream &os, const gcop::Quat &q);
 
 /**
  * Load a quaternion
  * @param is input stream
  * @param q quaternion
  * @return the input stream
  */
 std::istream& operator>>(std::istream &is, gcop::Quat &q);
 
}

#endif 
