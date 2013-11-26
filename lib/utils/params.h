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

  void SetDoubleVec(const char *name, const std::vector<double> &v);
  bool GetDoubleVec(const char *name, std::vector<double> &v) const;

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

  void Parse(char *line);

  std::map<std::string, std::string> valueMap;
  char buf[GCOP_PARAMS_MBL];
};

}

#endif
