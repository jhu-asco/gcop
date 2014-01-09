#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include "params.h"


using namespace gcop;
using namespace std;

static void replace(std::string& str, const std::string& a, const std::string& b)
{
  size_t pos = 0;
  while((pos = str.find(a, pos)) != std::string::npos) {
    str.replace(pos, a.length(), b);
    pos += b.length();
  }
}



static void Trim(string& str)
{
  string::size_type pos = str.find_last_not_of(' ');
  if(pos != string::npos) {
    str.erase(pos + 1);
    pos = str.find_first_not_of(' ');
    if(pos != string::npos) str.erase(0, pos);
  }
  else str.erase(str.begin(), str.end());
}


static void Tokenize(const string& str,
                     vector<string>& tokens,
                     const string& delimiters = " ")
{
  // Skip delimiters at beginning.
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  string::size_type pos     = str.find_first_of(delimiters, lastPos);
  
  while (string::npos != pos || string::npos != lastPos)
    {
      // Found a token, add it to the vector.
      tokens.push_back(str.substr(lastPos, pos - lastPos));
      // Skip delimiters.  Note the "not_of"
      lastPos = str.find_first_not_of(delimiters, pos);
      // Find next "non-delimiter"
      pos = str.find_first_of(delimiters, lastPos);
    }
}




Params::Params()
{
}

Params::Params(const char *fileName)
{
  Load(fileName);
}

Params::Params(FILE *file)
{
  Load(file);
}

Params::Params(iostream &io)
{
  Load(io);
}


Params::~Params()
{
  valueMap.clear();
}

void Params::Load(const char *fileName)
{
  FILE *file = fopen(fileName, "r");
  if (!file) {
    cout << "[W] Params::Load: failed to load file " << fileName << endl;
    return;
  }
  Load(file);
  fclose(file);
}

struct RemoveDelimiter
{
  bool operator()(char c)
  {
    return (c =='\r' || c =='\t' || c == ' ' || c == '\n');
  }
};

void Params::Parse(char *line)
{
  line[strcspn(line, "\n\r\t")] = 0;    
  if (line[0] == '#' || strlen(line) <= 2)
    return;
  char *name = strtok(line, "=");
  assert(name);
  char* value = strtok(NULL, "=");
  if (!value) {
    cerr << "Error:\tParams::Load:\tnull value for param <" << name << "> !" << endl;
    return;
  }
  
  //value = Trim(string(value)).c_str();
  //  name = Trim(string(name)).c_str();
  
  if (value[0] == '\"' || value[0] == ' ')
    value++;
  if (value[strlen(value) - 1] == '\"' || value[strlen(value) - 1] == ' ')
    value[strlen(value)-1] = 0;

  string svalue = string(value);
  svalue.erase( std::remove_if( svalue.begin(), svalue.end(), RemoveDelimiter()), svalue.end()); 

  string sname = string(name);
  sname.erase( std::remove_if( sname.begin(), sname.end(), RemoveDelimiter()), sname.end()); 

  cout << svalue << endl;
  valueMap[sname] = svalue;
}

void Params::Load(FILE *file)
{
  char line[256];
  while(fgets(line, 256, file)) {
    Parse(line);
  }
}

void Params::Load(iostream &io)
{
  char line[256];
  while(io.getline(line, 256)) {
    Parse(line);        
  }
}


void Params::Save(const char *fileName) const
{
  FILE *file = fopen(fileName, "w+");
  if (!file) {
    cout << "[W] Params::Load: failed to save to file " << fileName << endl;
    return;
  }
  Save(file);
  fclose(file);
}

void Params::Save(FILE *file) const
{
  Print(file);
}

void Params::Save(iostream &io) const
{
  Print(io);
}

void Params::SetInt(const char *name, int value)
{
  memset(buf, 0, GCOP_PARAMS_MBL);
  sprintf(buf, "%d", value);
  valueMap[name] = string(buf);
}

void Params::SetFloat(const char *name, float value)
{
  memset(buf, 0, GCOP_PARAMS_MBL);
  sprintf(buf, "%f", value);
  valueMap[name] = string(buf);
}

void Params::SetDouble(const char *name, double value)
{
  memset(buf, 0, GCOP_PARAMS_MBL);
  sprintf(buf, "%lf", value);
  valueMap[name] = string(buf);
}


void Params::SetVector2d(const char *name, const Vector2d &v)
{
  vector<double>::const_iterator it;
  stringstream s;
  unsigned int i = 0;
  for (int i = 0; i < v.size(); ++i) {
    s << v[i];
    if (i < v.size()-1)
      s << ",";
  }
  valueMap[name] = s.str();
}

bool Params::GetVector2d(const char *name, Vector2d &v) const
{
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    // cerr << "Error:\tParams::Get:\tparameter <" << name << "> not found!" << endl;
    return false;
  }

  vector<string> tokens;
  Tokenize(i->second, tokens, ", ");
  vector<string>::iterator it;
  int vi = 0;
  for (it = tokens.begin(); it != tokens.end(); ++it) {
    string str = *it;
    replace(str, string("pi"), string("3.141592"));  
    v[vi] = atof(str.c_str());
    ++vi;
  }
  return true;
}

void Params::SetVector3d(const char *name, const Vector3d &v)
{
  vector<double>::const_iterator it;
  stringstream s;
  unsigned int i = 0;
  for (int i = 0; i < v.size(); ++i) {
    s << v[i];
    if (i < v.size()-1)
      s << ",";
  }
  valueMap[name] = s.str();
}

bool Params::GetVector3d(const char *name, Vector3d &v) const
{
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    // cerr << "Error:\tParams::Get:\tparameter <" << name << "> not found!" << endl;
    return false;
  }

  vector<string> tokens;
  Tokenize(i->second, tokens, ", ");
  vector<string>::iterator it;
  int vi = 0;
  for (it = tokens.begin(); it != tokens.end(); ++it) {
    string str = *it;
    replace(str, string("pi"), string("3.141592"));  
    v[vi] = atof(str.c_str());
    ++vi;
  }
  return true;
}

void Params::SetVector4d(const char *name, const Vector4d &v)
{
  vector<double>::const_iterator it;
  stringstream s;
  unsigned int i = 0;
  for (int i = 0; i < v.size(); ++i) {
    s << v[i];
    if (i < v.size()-1)
      s << ",";
  }
  valueMap[name] = s.str();
}

bool Params::GetVector4d(const char *name, Vector4d &v) const
{
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    // cerr << "Error:\tParams::Get:\tparameter <" << name << "> not found!" << endl;
    return false;
  }

  vector<string> tokens;
  Tokenize(i->second, tokens, ", ");
  vector<string>::iterator it;
  int vi = 0;
  for (it = tokens.begin(); it != tokens.end(); ++it) {
    string str = *it;
    replace(str, string("pi"), string("3.141592"));  
    v[vi] = atof(str.c_str());
    ++vi;
  }
  return true;
}


void Params::SetVectorXd(const char *name, const VectorXd &v)
{
  vector<double>::const_iterator it;
  stringstream s;
  unsigned int i = 0;
  for (int i = 0; i < v.size(); ++i) {
    s << v[i];
    if (i < v.size()-1)
      s << ",";
  }
  valueMap[name] = s.str();
}

bool Params::GetVectorXd(const char *name, VectorXd &v) const
{
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    // cerr << "Error:\tParams::Get:\tparameter <" << name << "> not found!" << endl;
    return false;
  }

  vector<string> tokens;
  Tokenize(i->second, tokens, ", ");
  vector<string>::iterator it;
  int vi = 0;
  for (it = tokens.begin(); it != tokens.end(); ++it) {
    if (vi >= v.size()) {
      cout << "[E] Params::GetVectorXd: mismatched size in element " << name << " with expected v.size=" << v.size() << endl;
      return false;
    }
    string str = *it;
    replace(str, string("pi"), string("3.141592"));  
    v[vi] = atof(str.c_str());
    ++vi;
  }
  return true;
}



void Params::SetDoubleVec(const char *name, const vector<double> &v)
{
  vector<double>::const_iterator it;
  stringstream s;
  unsigned int i = 0;
  for (it = v.begin(); it != v.end(); ++it) {
    s << *it;
    ++i;
    if (i < v.size())
      s << ",";
  }
  valueMap[name] = s.str();
}

void Params::SetDoubleArray(const char *name, int n, const double *v)
{
  stringstream s;
  for (int i = 0; i < n; ++i) {
    s << v[i];
    if (i < n - 1)
      s << ",";
  }
  valueMap[name] = s.str();
}


void Params::SetBool(const char *name, bool value)
{
  memset(buf, 0, GCOP_PARAMS_MBL);
  sprintf(buf, "%ud", value);
  valueMap[name] = string(buf);
}

void Params::SetChar(const char *name, char value)
{
  memset(buf, 0, GCOP_PARAMS_MBL);
  sprintf(buf, "%c", value);
  valueMap[name] = string(buf);
}


void Params::SetString(const char *name, const string &v)
{
  valueMap[name] = v;
}

bool Params::GetInt(const char *name, int &v) const
{
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    //    cerr << "Error:\tParams::Get:\tparameter <" << name << "> not found!" << endl;
    return false;
  }
  v = atoi(i->second.c_str());
  return true;
}


bool Params::GetFloat(const char *name, float &v) const
{
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    //    cerr << "Error:\tParams::Get:\tparameter <" << name << "> not found!" << endl;
    return false;
  }
  string str = i->second;
  replace(str, string("pi"), string("3.141592"));
  v = atof(str.c_str());
  return true;
}


bool Params::GetDouble(const char *name, double &v) const
{
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    // cerr << "Error:\tParams::Get:\tparameter <" << name << "> not found!" << endl;
    return false;
  }
  string str = i->second;
  replace(str, string("pi"), string("3.141592"));  
  v = atof(str.c_str());
  return true;
}



bool Params::GetDoubleVec(const char *name, vector<double> &v) const
{
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    // cerr << "Error:\tParams::Get:\tparameter <" << name << "> not found!" << endl;
    return false;
  }

  vector<string> tokens;
  Tokenize(i->second, tokens, ",");
  vector<string>::iterator it;
  for (it = tokens.begin(); it != tokens.end(); ++it) {
    string str = *it;
    replace(str, string("pi"), string("3.141592"));  
    v.push_back(atof(str.c_str()));
  }
  return true;
}


bool Params::GetDoubleArray(const char *name, int n, double *v) const
{
  std::map<string, string>::const_iterator mi = valueMap.find(name);
  if (mi == valueMap.end()) {
    // cerr << "Error:\tParams::Get:\tparameter <" << name << "> not found!" << endl;
    return false;
  }

  vector<string> tokens;
  Tokenize(mi->second, tokens, ",");
  vector<string>::iterator it;
  int i = 0;
  for (it = tokens.begin(); it != tokens.end(), i < n; ++it, ++i) {    
    string str = *it;
    replace(str, string("pi"), string("3.141592"));
    v[i] = atof(str.c_str());
  }
  return true;
}


bool Params::GetBool(const char *name, bool &v) const
{
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    //    cerr << "Error:\tParams::Get:\tparameter <" << name << "> not found!" << endl;
    return false;
  }
  v = atoi(i->second.c_str());
  return true;
}

char Params::GetChar(const char *name, char &v) const
{
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    //    cerr << "Error:\tParams::Get:\tparameter <" << name << "> not found!" << endl;
    return false;
  }
  v = i->second[0];
  return true;
}


bool Params::GetString(const char *name, string &v) const
{
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    //    cerr << "Error:\tParams::Get:\tparameter <" << name << "> not found!" << endl;
    return false;
  }
  v = i->second;
  return true;
}


void Params::Print(FILE *file) const
{
  fprintf(file, "# Params: count=%d\n\n", (int)valueMap.size());
  std::map<string, string>::const_iterator i;
  for (i = valueMap.begin(); i != valueMap.end(); ++i)
    fprintf(file, "%s=%s\n", i->first.c_str(), i->second.c_str());
}

void Params::Print(iostream &io) const
{
  io << "# Params: count=" << valueMap.size() << "\n" << endl;
  std::map<string, string>::const_iterator i;
  for (i = valueMap.begin(); i != valueMap.end(); ++i)
    io << i->first << " " << i->second << endl;
}
