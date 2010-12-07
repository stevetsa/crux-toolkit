#ifndef EXTRACTROWS_H
#define EXTRACTROWS_H

#include "CruxApplication.h"

#include <string>

class ExtractRows: public CruxApplication {

 public:

  ExtractRows();
  ~ExtractRows();
  virtual int main(int argc, char** argv);
  virtual std::string getName();
  virtual std::string getDescription();

};


#endif
