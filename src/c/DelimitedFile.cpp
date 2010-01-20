#include "DelimitedFile.h"

#include <fstream>
#include <sstream>

extern "C" {

#include "carp.h"

}


using namespace std;

/**
 * tokenize a string by delimiter
 */
void DelimitedFile::tokenize(
  const string& str,
  vector<string>& tokens,
  char delimiter
  ) {

  tokens.clear();
  string::size_type lastPos = 0;
  string::size_type pos = str.find(delimiter, lastPos);

  while (string::npos != pos || string::npos != lastPos) {
    //found a token, add to the vector.
    string token = str.substr(lastPos, pos - lastPos);
    tokens.push_back(token);
    lastPos = pos+1;
    if (lastPos >= str.size() || pos >= str.size()) { 
      break;
    }
    pos = str.find(delimiter,lastPos);
  }
}

/**
 * convert string to data type
 */
template<typename TValue>  
bool DelimitedFile::from_string(
  TValue& value,
  const string& s
  ) {

  istringstream iss(s);
  return !(iss >> dec >> value).fail();
} 

/**
 * \returns a DelimitedFile object
 */  
DelimitedFile::DelimitedFile() {
  reset();
}

/**
 * \returns a DelimitedFile object and loads the tab-delimited
 * data specified by file_name.
 */  
DelimitedFile::DelimitedFile(
  const char *file_name, ///< the path of the file to read 
  bool hasHeader ///< indicate whether header exists
  ){

  loadData(file_name, hasHeader);
}

/** 
 * \returns a DelimitedFile object and loads the tab-delimited
 * data specified by file_name.
 */
DelimitedFile::DelimitedFile(
    const string& file_name, ///< the path of the file  to read
    bool hasHeader ///< indicates whether header exists
  ){

  loadData(file_name, hasHeader);
}

/**
 * Destructor
 */
DelimitedFile::~DelimitedFile() {

  //TODO : Do we need to do anything here?
}

/**
 *\returns the number of rows, assuming a square matrix
 */
unsigned int DelimitedFile::numRows() {

  if (data_.size() == 0) {
    carp(CARP_DEBUG, "DelimitedFile::numRows(): 0x0 matrix");
    return 0;
  }

  return data_[0].size();
}

/**
 *\returns the number of rows for a column
 */
unsigned int DelimitedFile::numRows(
  unsigned int col_idx ///<the column index
  ) {

  if (col_idx >= numCols()) {
    return 0;
  }

  return(data_[col_idx].size());
}

/**
 *\returns the number of columns
 */
unsigned int DelimitedFile::numCols() {

  return data_.size();
}

/**
 * clears the current data and column names,
 * parses the header if it exists,
 * reads the file one line at a time while
 * populating the data matrix with the 
 * strings separated by tabs.
 */
void DelimitedFile::loadData(
  const char *file_name, ///< the file path
  bool hasHeader ///< header indicator
  ) {

  data_.clear();
  column_names_.clear();

  fstream file(file_name, ios::in);

  if (!file.is_open()) {
    carp(CARP_ERROR, "Opening %s or reading failed", file_name);
    return;
  }

  std::string line;
  bool hasLine;

  vector<string>tokens;

  if (hasHeader) {
    hasLine = getline(file, line) != NULL;
    if (hasLine) {
      tokenize(line, tokens, '\t');
      for (vector<string>::iterator iter = tokens.begin();
        iter != tokens.end();
        ++iter) {
        addColumn(*iter);
      }
    }
    else {
      carp(CARP_WARNING,"No data/headers found!");
      return;
    }
  }

  hasLine = getline(file, line) != NULL;
  while (hasLine) {
    tokenize(line, tokens, '\t');
    for (unsigned int idx = 0; idx < tokens.size(); idx++) {
      if (numCols() <= idx) {
        addColumn();
      }
      data_[idx].push_back(tokens[idx]);
    }
    hasLine = getline(file, line) != NULL;
  }
  
  file.close();

  //reset the iterator.
  reset();
}

/**
 * loads a tab delimited file
 */
void DelimitedFile::loadData(
  const std::string& file, ///< the file path
  bool hasHeader ///< header indicator
  ) {

  loadData(file.c_str(), hasHeader);
}

/**
 * saves a tab delimited file
 */ 
void DelimitedFile::saveData(
  const std::string& file ///< the file path
  ) {

  saveData(file.c_str());
}

/**
 * saves a tab delimited file
 */ 
void DelimitedFile::saveData(
  const char* file ///< the file path
  ) {
  ofstream fout(file);

  //find the maximum number of rows.
  unsigned int maxRow = 0;
  for (unsigned int col_idx=0; col_idx < numCols(); col_idx++) {
    maxRow = max(maxRow, numRows(col_idx));
  }

  //print out the header if it exists.
  if (column_names_.size() != 0) {
    fout << column_names_[0];
    for (unsigned int col_idx=1; col_idx<column_names_.size(); col_idx++) {
      fout << "\t" << column_names_[col_idx];
    }
    fout << endl;
  }

  //print out all rows, using \t when
  //the row goes past the current column
  //size.
  for (unsigned int row_idx=0; row_idx<maxRow; row_idx++) {
    if (row_idx < numRows(0)) {
      fout << getString(0, row_idx);
    } else {
      fout << "\t";
    }
    for (unsigned int col_idx=1;col_idx<numCols();col_idx++) {
      fout <<"\t";
      if (row_idx < numRows(col_idx))
	fout << getString(col_idx, row_idx);
    }
    fout << endl;
  }

  fout.close();
}

/**
 * adds a column to the delimited file
 *\returns the column index.
 */
unsigned int DelimitedFile::addColumn(
  std::string& column_name ///< the column name
  ) {

  std::vector<std::string> new_col;
  data_.push_back(new_col);
  column_names_.push_back(column_name);
  return data_.size()-1;
}

/**
 * adds a column to the delimited file
 *\returns the new column index.
 */
unsigned int DelimitedFile::addColumn(
  const char* column_name ///< the column name
  ) {

  string string_name(column_name);
  return addColumn(string_name);
}

/**
 * adds a column to the delimited file
 *\returns the new column index.
 */
unsigned int DelimitedFile::addColumn() {
  std::vector<std::string> new_col;
  data_.push_back(new_col);
  return numCols() - 1;
}

/**
 * finds the index of a column
 *\returns the column index, -1 if not found.
 */ 
int DelimitedFile::findColumn(
  std::string& column_name ///< the column name
  ) {

  for (unsigned int col_idx=0;col_idx<column_names_.size();col_idx++) {
    if (column_names_[col_idx] == column_name) {
      return col_idx;
    }
  }
  return -1;
}

/**
 * finds the index of a column
 *\returns the column index, -1 if not found.
 */ 
int DelimitedFile::findColumn(
  const char* column_name ///< the column name
) {
  string sname = string(column_name);
  return findColumn(sname);
}

/**
 *\returns the string vector corresponding to the column
 */
std::vector<std::string>& DelimitedFile::getColumn(
  std::string column ///< the column name 
  ) {

  int col_idx = findColumn(column);

  if (col_idx != -1) {
    return data_[col_idx];
  }

  else {
    carp(CARP_ERROR,"column %s not found, returning column 0", column.c_str());
    return data_[0];
  }
}

/**
 *\returns the string vector corresponding to the column
 */
std::vector<std::string>& DelimitedFile::getColumn(
  unsigned int col_idx ///< the column index
  ) {
  
  return data_.at(col_idx);
}

/**
 *\returns the name of the column
 */
std::string& DelimitedFile::getColumnName(
  unsigned int col_idx ///< the column index
  ) {
  return column_names_.at(col_idx);
}

/**
 *\returns the string value of the cell
 */
std::string& DelimitedFile::getString(
  unsigned int col_idx, ///< the column index
  unsigned int row_idx  ///< the row index
  ) {
  return getColumn(col_idx)[row_idx];
}

/**
 * sets the string value of the cell
 */
void DelimitedFile::setString(
  unsigned int col_idx, ///< the column index
  unsigned int row_idx, ///< the row index
  std::string& value ///< the new value
  ) {

  //ensure there are enough columns
  while (col_idx >= numCols()) {
    addColumn();
  }
  std::vector<std::string>& col = getColumn(col_idx);

  //ensure there are enough rows
  while (row_idx >= col.size()) {
    col.push_back("");
  }

  col[row_idx] = value;
}

/**
 * sets the string value of the cell
 */
void DelimitedFile::setString(
  unsigned int col_idx, ///< the column index
  unsigned int row_idx, ///< the row index
  char* value ///< the new value
  ) {
  string svalue(value);
  setString(col_idx, row_idx, svalue);
}

/**
 *\returns the data type of the cell
 */
template<typename TValue>
TValue DelimitedFile::getValue(
  unsigned int col_idx, ///< the column index 
  unsigned int row_idx  ///< the row index
  ) {
  string& string_ans = getString(col_idx, row_idx);
  TValue type_ans;
  from_string<TValue>(type_ans, string_ans);
  return type_ans;
}

/**
 * gets a double type from cell, checks for infinity. 
 */
double DelimitedFile::getDouble(
  unsigned int col_idx, ///< the column index 
  unsigned int row_idx ///< the row index
  ) {

  string& string_ans = getString(col_idx,row_idx);
  if (string_ans == "Inf") {

    return numeric_limits<double>::infinity();
  } else if (string_ans == "-Inf") {

    return -numeric_limits<double>::infinity();
  }
  else {

    return getValue<double>(col_idx, row_idx);
  }
}

/** 
 * gets a double type from cell, checks for infinity.
 */
double DelimitedFile::getDouble(
  const char* column_name, ///<the column name
  unsigned int row_idx ///<the row index
) {

  int col_idx = findColumn(column_name);
  if (col_idx == -1) {
    carp(CARP_FATAL, "Cannot find column %s", column_name);
  }
  return getDouble(col_idx, row_idx);
}

/**
 * gets a double value from cell, checks for infinity
 * uses the current_row_ as the row index
 */
double DelimitedFile::getDouble(
  const char* column_name ///<the column name
) {

  if (current_row_ >= numRows()) {
    carp(CARP_FATAL, "Iterated past maximum number of rows!");
  }
  return getDouble(column_name, current_row_);
}

/**
 * gets an integer type from cell. 
 */
int DelimitedFile::getInteger(
  unsigned int col_idx, ///< the column index 
  unsigned int row_idx ///< the row index
  ) {
  //TODO : check the string for a valid integer.
  return getValue<int>(col_idx, row_idx);
}

/**
 * get an integer type from cell, checks for infintiy.
 */
int DelimitedFile::getInteger(
  const char* column_name, ///< the column name
  unsigned int row_idx ///<the row index
) {

  int col_idx = findColumn(column_name);
  if (col_idx == -1) {
    carp(CARP_FATAL, "Cannot find column %s", column_name);
  }

  return getInteger(col_idx, row_idx);
}


/**
 * get an integer type from cell, checks for infinity.
 * uses the current_row_ as the row index.
 */
int DelimitedFile::getInteger(
    const char* column_name ///< the column name
  ) {

  if (current_row_ >= numRows()) {
    carp(CARP_FATAL, "Iterated past maximum number of rows!");
  }

  return getInteger(column_name, current_row_);
}


/*Iterator functions.*/
/**
 * resets the current_row_ index to 0.
 */
void DelimitedFile::reset() {

  current_row_ = 0;
}

/**
 * increments the current_row_, 
 */
void DelimitedFile::next() {
  if (current_row_ < numRows())
    current_row_++;
}


/**
 * \returns whether there are more rows to 
 * iterate through
 */
BOOLEAN_T DelimitedFile::hasNext() {
  return current_row_ < numRows();
}
