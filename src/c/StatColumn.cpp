/**
 * \file StatColumn.cpp 
 * \brief Give a tab delimited file and a comma-separated list of column names
 * print out a tab delimied file with only those columns
 *****************************************************************************/
#include "StatColumn.h"

#include "DelimitedFile.h"

using namespace std;


/**
 * \returns a blank ExtractRows object
 */
StatColumn::StatColumn() {

}

/**
 * Destructor
 */
StatColumn::~StatColumn() {
}

/**
 * main method for StatColumn
 */
int StatColumn::main(int argc, char** argv) {

   /* Define optional command line arguments */
  const char* option_list[] = {
    "delimiter",
    "header",
    "verbosity"
  };
  int num_options = sizeof(option_list) / sizeof(char*);

  /* Define required command line arguments */
  const char* argument_list[] = {"tsv file", "column name"};
  int num_arguments = sizeof(argument_list) / sizeof(char*);

  /* Initialize the application */
  initialize(argument_list, num_arguments,
    option_list, num_options, argc, argv);

  const char* delimited_filename = get_string_parameter_pointer("tsv file");

  string column_name_string = 
    string(get_string_parameter_pointer("column name"));

  char delimiter = get_delimiter_parameter("delimiter");

  DelimitedFileReader delimited_file(delimited_filename, true, delimiter);
  
  int col_idx = delimited_file.findColumn(column_name_string);

  if (col_idx == -1) {
    carp(CARP_ERROR,"column not found:%s\n\n%s", 
      column_name_string.c_str(),
      delimited_file.getAvailableColumnsString().c_str());
    return(-1);
  }

  vector<FLOAT_T> data;

  FLOAT_T sum = 0;

  while (delimited_file.hasNext()) {

    FLOAT_T current = delimited_file.getFloat(col_idx);
    data.push_back(current);
    sum += current;
    delimited_file.next();

  }
  
  sort(data.begin(), data.end(), less<FLOAT_T>());

  FLOAT_T min = data.front();
  FLOAT_T max = data.back();
  
  FLOAT_T average = sum / (FLOAT_T)data.size();
  
  FLOAT_T median = 0.0;

  if (data.size() >= 2) {
    int half = data.size() / 2;
    if (data.size() % 2 == 0) {
      median = (data[half] + data[half-1]) / 2.0;
    } else {
      median = data[half];
    }
  }

  //print out the header
  cout <<"N\tMin\tMax\tSum\tAverage\tMedian"<<endl;
  
  cout << data.size() << "\t";
  cout << min << "\t";
  cout << max << "\t";
  cout << sum << "\t";
  cout << average << "\t";
  cout << median << endl;

  return 0;
}

/**
 * \returns the command name for StatColumn
 */
string StatColumn::getName() {
  return "stat-column";
}

/**
 * \returns the description for StatColumn
 */
string StatColumn::getDescription() {

  return "Collects statistics from a column of data in "
         "a delimited file";
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
