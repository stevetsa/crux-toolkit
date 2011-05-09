/**
 * \file SortColumn.cpp 
 * \brief Given a delimited file and a column-name, sort the 
 * file.
 *****************************************************************************/
#include "SortColumn.h"

#include "DelimitedFile.h"

#include "fdstream.hpp"

#include <errno.h>

using namespace std;


/**
 * \returns a blank SortColumn object
 */
SortColumn::SortColumn() {

}

/**
 * Destructor
 */
SortColumn::~SortColumn() {
}

/**
 * main method for SortColumn
 */
int SortColumn::main(int argc, char** argv) {

   /* Define optional command line arguments */
  const char* option_list[] = {
    "delimiter",
    "header",
    "column-type",
    "ascending",
    "verbosity"
  };
  int num_options = sizeof(option_list) / sizeof(char*);

  /* Define required command line arguments */
  const char* argument_list[] = {"tsv file", "column name"};
  int num_arguments = sizeof(argument_list) / sizeof(char*);

  /* Initialize the application */
  initialize(argument_list, num_arguments,
    option_list, num_options, argc, argv);

  /* Get parameters */
  delimited_filename_ = 
    string(get_string_parameter_pointer("tsv file"));

  column_name_string_ = 
    string(get_string_parameter_pointer("column name"));

  column_type_ = get_column_type_parameter("column-type");

  ascending_ = get_boolean_parameter("ascending");
  delimiter_ = get_delimiter_parameter("delimiter");
  header_ = get_boolean_parameter("header");


  DelimitedFileReader delimited_file(delimited_filename_, true, delimiter_);
  
  int col_sort_idx = delimited_file.findColumn(column_name_string_);

  if (col_sort_idx == -1) {
    carp(CARP_ERROR,"column not found:%s\n\n%s", 
      column_name_string_.c_str(),
      delimited_file.getAvailableColumnsString().c_str());
    return(-1);
  }

  col_sort_idx_ = (unsigned int)col_sort_idx;

  vector<int> temp_file_descriptors;
  vector<boost::fdostream*> temp_file_streams;
  vector<string> temp_filenames;

  DelimitedFile* temp_delimited = new DelimitedFile();
  temp_delimited->setDelimiter(delimiter_);
  for (unsigned int col_idx=0;col_idx < delimited_file.numCols();col_idx++) {
    temp_delimited->addColumn(delimited_file.getColumnName(col_idx));
  }

  int max_rows = 200000;

  int current_count = 0;

  while(delimited_file.hasNext()) {

    unsigned int new_row = temp_delimited->addRow();
    for (unsigned int col_idx = 0;col_idx < delimited_file.numCols();col_idx++) {
      temp_delimited->setString(col_idx, new_row, delimited_file.getString(col_idx));
    }
  

    current_count++;

    if (current_count >= max_rows) {
      cerr <<"Sorting "<<max_rows<<" rows"<<endl;
      sortDelimited(temp_delimited);

      char ctemp_filename[50] = "SortColumn_XXXXXX";
    
      int fd=mkstemp(ctemp_filename);
      if (fd == -1) {
        cerr <<"Error creating temp file!"<<endl;
        cerr <<string(strerror(errno))<<endl;
        return(-1);
      }
      
      string temp_filename(ctemp_filename);
  
      boost::fdostream* out = new boost::fdostream(fd);
  
      (*out) << (*temp_delimited);
      temp_file_streams.push_back(out);
      temp_file_descriptors.push_back(fd);
      temp_filenames.push_back(temp_filename);
      delete(temp_delimited);
      temp_delimited = new DelimitedFile();
      temp_delimited->setDelimiter(delimiter_);
      for (unsigned int col_idx = 0;col_idx < delimited_file.numCols();col_idx++) {
        temp_delimited->setString(col_idx, new_row, delimited_file.getString(col_idx));
      }
      current_count = 0;
  
    }

    delimited_file.next();
  }

  if (temp_file_descriptors.size() == 0) {

    //no temporary files used, print out the sorted output.
    sortDelimited(temp_delimited);
    cout << (*temp_delimited);
    delete temp_delimited;
  } else {

    //sort the current and save to a file. 
    sortDelimited(temp_delimited);
    char ctemp_filename[50] = "SortColumn_XXXXXX";

    int fd = mkstemp(ctemp_filename);
    if (fd == -1) {
      cerr <<"Error creating temp file!"<<endl;
      cerr <<string(strerror(errno))<<endl;
      exit(-1);
    }
    string temp_filename(ctemp_filename);

    boost::fdostream* out = new boost::fdostream(fd);
    (*out) << (*temp_delimited);
    temp_file_streams.push_back(out);
    temp_file_descriptors.push_back(fd);
    temp_filenames.push_back(temp_filename);
    delete(temp_delimited);

    //flush everything.
    for (unsigned int i=0;i<temp_file_descriptors.size();i++) {
      temp_file_streams[i]->flush();
    }

    //merge the temporary files together
    mergeDelimitedFiles(temp_filenames);
  }
  
  //clean everything up
  for (unsigned int idx=0;idx<temp_file_descriptors.size();idx++) {
    delete temp_file_streams[idx];
    close(temp_file_descriptors[idx]);
    remove(temp_filenames[idx].c_str());
  }

  return 0;

}

/**
 * \returns the command name for SortColumn
 */
string SortColumn::getName() {
  return "sort-by-column";
}

/**
 * \returns the description for SortColumn
 */
string SortColumn::getDescription() {

  return "Sorts a delimited file by a column ";
}

/**
 * sorts the delimited file based upon the column type parameter
 */
void SortColumn::sortDelimited(
  DelimitedFile* delimited ///< delimited file to sort.
  ) {

  switch(column_type_) {

    case COLTYPE_STRING:
      delimited->sortByStringColumn(column_name_string_, ascending_);
      break;
    case COLTYPE_INT:
      delimited->sortByIntegerColumn(column_name_string_, ascending_);
      break;
    case COLTYPE_REAL:
      delimited->sortByFloatColumn(column_name_string_, ascending_);
      break;
    case NUMBER_COLTYPES:
    case COLTYPE_INVALID:
      carp(CARP_FATAL, "Unknow column type");
  }
}

/**
 * merges a list of sorted delimited files and prints 
 * out the resulting sorted file.
 */
void SortColumn::mergeDelimitedFiles(
  vector<string>& temp_filenames
  ) {

  vector<DelimitedFileReader*> delimited_files;
  for (unsigned int idx=0;idx < temp_filenames.size();idx++) {
    delimited_files.push_back(new DelimitedFileReader(temp_filenames[idx], delimiter_));
  }
  
  if (header_) {
    cout << delimited_files[0]->getHeaderString() << endl;
  }

  int best_idx = -1;

  do {
    best_idx = -1;
    for (unsigned int idx = 0;idx < delimited_files.size();idx++) {
      if (delimited_files[idx]->hasNext()) {
        if (best_idx == -1) {
          best_idx = idx;
        } else {
          int compare_result = compare(delimited_files[idx], delimited_files[best_idx]);
          if ((ascending_) && (compare_result == -1)) {
            best_idx = idx;
          } else if ((!ascending_) && (compare_result == 1)) {
            best_idx = idx;
          }
        }
      }
  
      if (best_idx != -1) {
        cout << delimited_files[best_idx]->getString()<<endl;
        delimited_files[best_idx]->next();
      }
    }
  } while (best_idx != -1);

  for (unsigned int idx=0;idx<temp_filenames.size();idx++) {
    delete delimited_files[idx];
  }
}

/**
 * /returns the result of comparing two file(s) current row and column: 
 * 1 : file1 > file2
 * -1 : file1 < file2
 * 0 : file1 = file2
 */
int SortColumn::compare(
  DelimitedFileReader* file1, 
  DelimitedFileReader* file2
  ) {

  bool isless = false;
  bool isequal = false;
  
  switch (column_type_) {
    case COLTYPE_STRING:
      isless = file1->getString(col_sort_idx_) < file2->getString(col_sort_idx_);
      if (!isless) {
        isequal = file1->getString(col_sort_idx_) == file2->getString(col_sort_idx_);
      }
    break;
    case COLTYPE_REAL:
      isless = file1->getDouble(col_sort_idx_) < file2->getDouble(col_sort_idx_);
      if (!isless) {
        isequal = file1->getDouble(col_sort_idx_) == file2->getDouble(col_sort_idx_);
      }
    break;
    case COLTYPE_INT:
      isless = file1->getInteger(col_sort_idx_) < file2->getInteger(col_sort_idx_);
      if (!isless) {
        isequal = file1->getInteger(col_sort_idx_) == file2->getInteger(col_sort_idx_);
      }
    break;
    case NUMBER_COLTYPES:
    case COLTYPE_INVALID:
      carp(CARP_FATAL,"Column type invalid!");

  }

  if (isless) {
    return -1;
  } else if (isequal) {
    return 0;
  } else {
    return 1;
  }
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
