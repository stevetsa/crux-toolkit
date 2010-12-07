#ifndef VERSION_H
#define VERSION_H

#include "CruxApplication.h"
#include "DelimitedFileReader.h"

#include <string>

class Version: public CruxApplication {

 public:

  Version();
  ~Version();
  virtual int main(int argc, char** argv);
  virtual std::string getName();
  virtual std::string getDescription();

};


#endif
