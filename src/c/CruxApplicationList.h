#ifndef CRUXAPPLICATIONLIST_H
#define CRUXAPPLICATIONLIST_H

#include "CruxApplication.h"

#include <vector>
#include <string>

class CruxApplicationList : std::vector<CruxApplication*> {

 protected:
  std::string list_name_;

 public:

  CruxApplicationList(const char* list_name);
  ~CruxApplicationList();
  void add(CruxApplication* application);
  CruxApplication* find(const std::string& appname);
  CruxApplication* find(const char* appname);
  void usage();

  int main(int argc, char** argv);
};


#endif
