#include "ExtractRows.h"

#include <iostream>

#include "DelimitedFileReader.h"

#include "carp.h"
#include "parameter.h"

using namespace std;

ExtractRows::ExtractRows() {

}

ExtractRows::~ExtractRows() {
}

/*
void ExtractRows::printAvailableRows(DelimitedFileReader& infile) {
}
*/

int ExtractRows::main(int argc, char** argv) {

   /* Define optional command line arguments */
  const char* option_list[] = {
    "verbosity"
  };
  int num_options = sizeof(option_list) / sizeof(char*);

  /* Define required command line arguments */
  const char* argument_list[] = {"tsv file", "column name", "column value"};
  int num_arguments = sizeof(argument_list) / sizeof(char*);

    // Verbosity level for set-up/command line reading 
  set_verbosity_level(CARP_WARNING);

  // Initialize parameter.c and set default values
  initialize_parameters();

  // Define optional and required arguments 
  select_cmd_line_options(option_list, num_options);
  select_cmd_line_arguments(argument_list, num_arguments);

  // Parse the command line, including optional params file
  // Includes syntax, type, and bounds checking, dies on error 
  const char* cmd_name = this->getName().c_str();
  char* full_cmd = cat_string("crux-util ", cmd_name);
  parse_cmd_line_into_params_hash(argc, argv, cmd_name);
  free(full_cmd);

  const char* delimited_filename = 
    get_string_parameter_pointer("tsv file");

  string column_name = 
    string(get_string_parameter_pointer("column name"));
  string column_value = 
    string(get_string_parameter_pointer("column value"));

  DelimitedFileReader delimited_file(delimited_filename, true);
  int column_idx = delimited_file.findColumn(column_name);

  if (column_idx == -1) {
    carp(CARP_FATAL,"column not found:%s\n\n:%s",
      column_name.c_str(),
      delimited_file.getAvailableColumnsString().c_str());
  }

  cout << delimited_file.getHeaderString() << endl;

  while (delimited_file.hasNext()) {
    if (delimited_file.getString(column_idx) == column_value) {
      cout << delimited_file.getString() << endl;
    }
    delimited_file.next();
  }

  return 0;

}

string ExtractRows::getName() {
  return "extract-rows";
}

string ExtractRows::getDescription() {

  return "prints out rows that match a particular column value";

}
