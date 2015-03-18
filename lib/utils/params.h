#ifndef GCOP_PARAMS_H
#define GCOP_PARAMS_H

#include <map>
#include <cstring>
#include <iostream>
#include <Eigen/Dense>

/*
 *  Provides storage and retrieval for arbitrary 
 *  number of parameters of different types (currently only primitives and strings)
 *  TODO: include matrix parameter
 *  Parameters can be loaded/saved from/to a file, or can be inserted at run-time 
 * 
 *
 *  Author: Marin Kobilarov (mkobilar@@robotics.usc.edu)
 */

namespace gcop {

  using namespace Eigen;

#define GCOP_PARAMS_MBL 256

class Params
{
 public:

  typedef Matrix<double, 6, 6> Matrix6d;
  typedef Matrix<double, 6, 1> Vector6d;

  Params();
  Params(FILE *file);
  Params(const char *fileName);
  Params(std::iostream &io);

  virtual ~Params();

  void Load(const char *fileName);
  void Load(FILE *file);
  void Load(std::iostream &io);

  void Save(const char *fileName) const;
  void Save(FILE* file) const;
  void Save(std::iostream &io) const;

  bool Exists(const char *name) const { return valueMap.find(std::string(name)) != valueMap.end(); }

  void SetInt(const char *name, int v);
  bool GetInt(const char *name, int &v) const;
  
  void SetFloat(const char *name, float v);
  bool GetFloat(const char *name, float &v) const;

  void SetDouble(const char *name, double v);
  bool GetDouble(const char *name, double &v) const;

  void SetVectorXd(const char *name, const VectorXd &v);
  bool GetVectorXd(const char *name, VectorXd &v) const;

  void SetVector2d(const char *name, const Vector2d &v);
  bool GetVector2d(const char *name, Vector2d &v) const;

  void SetVector3d(const char *name, const Vector3d &v);
  bool GetVector3d(const char *name, Vector3d &v) const;

  void SetVector4d(const char *name, const Vector4d &v);
  bool GetVector4d(const char *name, Vector4d &v) const;

  void SetVector6d(const char *name, const Vector6d &v);
  bool GetVector6d(const char *name, Vector6d &v) const;

  void SetMatrix3d(const char *name, const Matrix3d &m);
  bool GetMatrix3d(const char *name, Matrix3d &m) const;

  void SetMatrix6d(const char *name, const Matrix6d &m);
  bool GetMatrix6d(const char *name, Matrix6d &m) const;

  void SetDoubleVec(const char *name, const std::vector<double> &v);
  bool GetDoubleVec(const char *name, std::vector<double> &v) const;

  void SetFloatArray(const char *name, int n, const float *v);
  bool GetFloatArray(const char *name, int n, float *v) const;

  void SetDoubleArray(const char *name, int n, const double *v);
  bool GetDoubleArray(const char *name, int n, double *v) const;
  
  void SetString(const char *name, const std::string &v);
  bool GetString(const char *name, std::string &v) const;
  
  void SetBool(const char *name, bool v);
  bool GetBool(const char *name, bool &v) const;

  void SetChar(const char *name, char v);
  char GetChar(const char *name, char &v) const;

  void Print(FILE *file = stdout) const;
  void Print(std::iostream &io) const;

 protected:


  struct RemoveDelimiter
  {
    bool operator()(char c)
    {
      return (c =='\r' || c =='\t' || c == ' ' || c == '\n');
    }
  };
  

  void Parse(char *line);

  std::map<std::string, std::string> valueMap;
  char buf[GCOP_PARAMS_MBL];
};

}

#endif
