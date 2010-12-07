#ifndef SEARCHFORXLINKS_H
#define SEARCHFORXLINKS_H

#include "CruxApplication.h"
#include "DelimitedFileReader.h"

#include <string>

class SearchForXLinks: public CruxApplication {

 public:

  SearchForXLinks();
  ~SearchForXLinks();
  virtual int main(int argc, char** argv);
  virtual std::string getName();
  virtual std::string getDescription();

};


#endif
