#ifndef EXTRACTCOLUMNS_H
#define EXTRACTCOLUMNS_H

#include "CruxApplication.h"
#include "DelimitedFileReader.h"

#include <string>

class ExtractColumns: public CruxApplication {
 protected:
  void printAvailableColumns(DelimitedFileReader& infile);

 public:

  ExtractColumns();
  ~ExtractColumns();
  virtual int main(int argc, char** argv);
  virtual std::string getName();
  virtual std::string getDescription();

};


#endif
