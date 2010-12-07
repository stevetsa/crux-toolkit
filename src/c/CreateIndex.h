#ifndef CREATEINDEX_H
#define CREATEINDEX_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <unistd.h>
#include <signal.h>
#include "objects.h"
#include "carp.h"
#include "peptide.h"
#include "peptide_src.h"
#include "protein.h"
#include "database.h"
#include "index.h"
#include "protein_index.h"
#include "parameter.h"

#include "CruxApplication.h"




#include <string>

class CreateIndex: public CruxApplication {

 public:

  CreateIndex();
  ~CreateIndex();
  virtual int main(int argc, char** argv);
  virtual std::string getName();
  virtual std::string getDescription();
  

};


#endif
