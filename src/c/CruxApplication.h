#ifndef CRUXAPPLICATION_H
#define CRUXAPPLICATION_H

#include <string>

class CruxApplication{
 protected:
  /**
   * \brief Perform the set-up steps common to all crux commands:
   * initialize parameters, parse command line, set verbosity, open
   * output directory, write params file. 
   *
   * Uses the given command name, arguments and options for parsing the
   * command line.
   */
  void initializeRun(
    const char** argument_list, ///< list of required arguments
    int num_arguments,          ///< number of elements in arguments_list
    const char** option_list,   ///< list of optional flags
    int num_options,            ///< number of elements in options_list
    int argc,                   ///< number of tokens on cmd line
    char** argv                 ///< array of command line tokens
  );

 public:
  virtual int main(int argc, char** argv)=0;
  virtual std::string getName()=0;
  virtual std::string getFileString();
  virtual std::string getDescription()=0;


  virtual ~CruxApplication();

};



#endif
