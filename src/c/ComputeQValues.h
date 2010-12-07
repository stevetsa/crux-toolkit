#ifndef ComputeQValues_H
#define ComputeQValues_H

#include "CruxApplication.h"
#include "DelimitedFileReader.h"

#include <string>

class ComputeQValues: public CruxApplication {

 public:

  ComputeQValues();
  ~ComputeQValues();
  virtual int main(int argc, char** argv);
  virtual std::string getName();
  virtual std::string getDescription();
  virtual std::string getFileString();
};


#endif
