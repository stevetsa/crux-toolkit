#include "DelimitedFileReader.h"

#include <fstream>

#include "carp.h"

using namespace std;

/**
 * \returns a DelimitedFileReader object
 */  
DelimitedFileReader::DelimitedFileReader() {

  file_ptr = NULL;
}

/**
 * \returns a DelimitedFileReader object and loads the tab-delimited
 * data specified by file_name.
 */  
DelimitedFileReader::DelimitedFileReader(
  const char *file_name, ///< the path of the file to read 
  bool hasHeader ///< indicate whether header exists
  ){

  loadData(file_name, hasHeader);
}

/** 
 * \returns a DelimitedFileReader object and loads the tab-delimited
 * data specified by file_name.
 */
DelimitedFileReader::DelimitedFileReader(
    const string& file_name, ///< the path of the file  to read
    bool hasHeader ///< indicates whether header exists
  ){

  loadData(file_name, hasHeader);
}

/**
 * Destructor
 */
DelimitedFileReader::~DelimitedFileReader() {
  if (file_ptr != null) {

    file_ptr -> close();
    delete file_ptr;
  }
}

/**
 *\returns the number of columns
 */
unsigned int DelimitedFileReader::numCols() {

  return column_names_.size();
}

/**
 * clears the current data and column names,
 * parses the header if it exists,
 * reads the file one line at a time while
 * populating the data matrix with the 
 * strings separated by tabs.
 */
void DelimitedFileReader::loadData(
  const char *file_name, ///< the file path
  bool hasHeader ///< header indicator
  ) {

  clear();
  file_ptr -> new fstream(file_name, ios::in);

  if (!file_ptr -> is_open()) {
    carp(CARP_ERROR, "Opening %s or reading failed", file_name);
    return;
  }

  string line;
  bool hasLine;

  vector<string>tokens;

  if (hasHeader) {
    hasLine = getline(*file_ptr, line) != NULL;
    if (hasLine) {
      tokenize(line, tokens, '\t');
      for (vector<string>::iterator iter = tokens.begin();
        iter != tokens.end();
        ++iter) {
        column_names_.push_back(*iter);
      }
    }
    else {
      carp(CARP_WARNING,"No data/headers found!");
      return;
    }
  }

  hasLine = getline(*file_ptr, line) != NULL;
  while (hasLine) {
    tokenize(line, tokens, '\t');
    int row_idx = addRow();
    for (unsigned int col_idx = 0; col_idx < tokens.size(); col_idx++) {
      if (numCols() <= col_idx) {
        addColumn();
      }
      setString(col_idx, row_idx, tokens[col_idx]);
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
void DelimitedFileReader::loadData(
  const string& file, ///< the file path
  bool hasHeader ///< header indicator
  ) {

  loadData(file.c_str(), hasHeader);
}

/**
 * saves a tab delimited file
 */ 
void DelimitedFileReader::saveData(
  const string& file ///< the file path
  ) {

  saveData(file.c_str());
}

/**
 * saves a tab delimited file
 */ 
void DelimitedFileReader::saveData(
  const char* file ///< the file path
  ) {
  
  ofstream fout(file);
  fout << *this;
  fout.close();
}

/**
 * adds a column to the delimited file
 *\returns the column index.
 */
unsigned int DelimitedFileReader::addColumn(
  string& column_name ///< the column name
  ) {

  vector<string> new_col;
  data_.push_back(new_col);
  column_names_.push_back(column_name);
  return data_.size()-1;
}

/**
 * adds a column to the delimited file
 *\returns the new column index.
 */
unsigned int DelimitedFileReader::addColumn(
  const char* column_name ///< the column name
  ) {

  string string_name(column_name);
  return addColumn(string_name);
}

/**
 * adds a column to the delimited file
 *\returns the new column index.
 */
unsigned int DelimitedFileReader::addColumn() {
  vector<string> new_col;
  data_.push_back(new_col);
  return numCols() - 1;
}

/**
 * adds a vector of columns to the delimited file
 */
void DelimitedFileReader::addColumns(
  vector<string>& column_names
  ) {
  cout <<"Number of columns:"<<column_names.size()<<endl;
  for (unsigned int col_idx = 0;col_idx < column_names.size(); col_idx++) {
    cout <<"Adding :"<<column_names[col_idx]<<endl;
    //addColumn(column_names[col_idx]);
  }
}



/**
 * finds the index of a column
 *\returns the column index, -1 if not found.
 */ 
int DelimitedFileReader::findColumn(
  const string& column_name ///< the column name
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
int DelimitedFileReader::findColumn(
  const char* column_name ///< the column name
) {
  string sname = string(column_name);
  return findColumn(sname);
}

/**
 *\returns the string vector corresponding to the column
 */
vector<string>& DelimitedFileReader::getColumn(
  string column ///< the column name 
  ) {

  int col_idx = findColumn(column);

  if (col_idx != -1) {
    return data_[col_idx];
  } else {
    carp(CARP_ERROR,"column %s not found, returning column 0", column.c_str());
    return data_[0];
  }
}

/**
 *\returns the string vector corresponding to the column
 */
vector<string>& DelimitedFileReader::getColumn(
  unsigned int col_idx ///< the column index
  ) {
  
  return data_.at(col_idx);
}

/**
 *\returns the name of the column
 */
string& DelimitedFileReader::getColumnName(
  unsigned int col_idx ///< the column index
  ) {
  return column_names_.at(col_idx);
}

/**
 *\returns the column_names
 */
vector<string>& DelimitedFileReader::getColumnNames() {

  return column_names_;
}


/**
 * adds a row to the delimited file
 *\returns the new row index
 */
unsigned int DelimitedFileReader::addRow() {
  if (numCols() == 0) {
    carp(CARP_FATAL,"Must have at least one column before calling add row!");
  }

  unsigned int row_idx = numRows();
  for (unsigned int col_idx = 0;col_idx < numCols(); col_idx++) {
      setString(col_idx, row_idx, "");
  }
  return row_idx;
}

/**
 *\returns the string value of the cell
 */
string& DelimitedFileReader::getString(
  unsigned int col_idx, ///< the column index
  unsigned int row_idx  ///< the row index
  ) {
  return getColumn(col_idx).at(row_idx);
}

/** 
 * gets a string value of the cell.
 */
string& DelimitedFileReader::getString(
  const char* column_name, ///<the column name
  unsigned int row_idx ///< the row index
  ) {
  int col_idx = findColumn(column_name);
  if (col_idx == -1) {
    carp(CARP_FATAL, "Cannot find column %s", column_name);
  }
  return getColumn(col_idx)[row_idx];
}

/**
 * gets a string value of the cell
 * uses the current_row_ as the row index
 */
string& DelimitedFileReader::getString(
  const char* column_name ///<the column name
  ) {

  if (current_row_ >= numRows()) {
    carp(CARP_FATAL, "Iterated past maximum number of rows!");
  }
  return getString(column_name, current_row_);
}

/**
 * sets the string value of the cell
 */
void DelimitedFileReader::setString(
  unsigned int col_idx, ///< the column index
  unsigned int row_idx, ///< the row index
  string& value ///< the new value
  ) {

  //ensure there are enough columns
  while (col_idx >= numCols()) {
    addColumn();
  }
  vector<string>& col = getColumn(col_idx);

  //ensure there are enough rows
  while (row_idx >= col.size()) {
    col.push_back("");
  }

  col[row_idx] = value;
}

/**
 * sets the string value of the cell
 */
void DelimitedFileReader::setString(
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
TValue DelimitedFileReader::getValue(
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
FLOAT_T DelimitedFileReader::getFloat(
    unsigned int col_idx, ///< the column index
    unsigned int row_idx ///< the row index
) {
  
  string& string_ans = getString(col_idx,row_idx);
  if (string_ans == "Inf") {

    return numeric_limits<FLOAT_T>::infinity();
  } else if (string_ans == "-Inf") {

    return -numeric_limits<FLOAT_T>::infinity();
  }
  else {

    return getValue<FLOAT_T>(col_idx, row_idx);
  }
}

/** 
 * gets a double type from cell, checks for infinity.
 */
FLOAT_T DelimitedFileReader::getFloat(
    const char* column_name, ///<the column name
    unsigned int row_idx ///< the row index
) {
  
  int col_idx = findColumn(column_name);
  if (col_idx == -1) {
    carp(CARP_FATAL, "Cannot find column %s", column_name);
  }
  return getFloat(col_idx, row_idx);
}

/**
 * gets a double value from cell, checks for infinity
 * uses the current_row_ as the row index
 */
FLOAT_T DelimitedFileReader::getFloat(
  const char* column_name ///<the column name
) {
  
  if (current_row_ >= numRows()) {
    carp(CARP_FATAL, "Iterated past maximum number of rows!");
  }
  return getFloat(column_name, current_row_);
}


FLOAT_T DelimitedFileReader::getFloat(
    const string& column_name, ///<the column name
    unsigned int row_idx ///<the row index
) {

  return getFloat(column_name.c_str(), row_idx);
}

/**
 * gets a double type from cell, checks for infinity. 
 */
double DelimitedFileReader::getDouble(
  unsigned int col_idx, ///< the column index 
  unsigned int row_idx ///< the row index
  ) {

  string& string_ans = getString(col_idx,row_idx);
  if (string_ans == "") {
    return 0.0;
  } else if (string_ans == "Inf") {

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
double DelimitedFileReader::getDouble(
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
double DelimitedFileReader::getDouble(
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
int DelimitedFileReader::getInteger(
  unsigned int col_idx, ///< the column index 
  unsigned int row_idx ///< the row index
  ) {
  //TODO : check the string for a valid integer.
  return getValue<int>(col_idx, row_idx);
}

/**
 * get an integer type from cell, checks for infintiy.
 */
int DelimitedFileReader::getInteger(
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
int DelimitedFileReader::getInteger(
    const char* column_name ///< the column name
  ) {

  if (current_row_ >= numRows()) {
    carp(CARP_FATAL, "Iterated past maximum number of rows!");
  }

  return getInteger(column_name, current_row_);
}

/**
 * gets an vector of strings from cell where the
 * string in the cell has delimiters that are
 * different than the column delimiter. The
 * default delimiter is a comma
 * uses the current_row_ as the row index.
 * clears the integer vector before 
 * populating it.
 */
void DelimitedFileReader::getStringVectorFromCell(
  const char* column_name, ///< the column name
  std::vector<std::string>& string_vector, ///<the vector of integers
  char delimiter ///<the delimiter to use
  ) {

  string& string_ans = getString(column_name);

  //get the list of strings separated by delimiter
  string_vector.clear();
  tokenize(string_ans, string_vector, delimiter);
}

/**
 * gets an vector of integers from cell where the
 * string in the cell are integers which are separated
 * by a delimiter which is differnt than the column
 * delimiter.  The default delimiter is a comma
 * uses the current_row_ as the row index.
 * clears the integer vector before 
 * populating it.
 */
void DelimitedFileReader::getIntegerVectorFromCell(
    const char* column_name, ///< the column name
    vector<int>& int_vector, ///<the vector of integers
    char delimiter ///<the delimiter to use
  ) {
  
  //get the list of strings separated by delimiter
  vector<string> string_vector_ans;

  getStringVectorFromCell(column_name, string_vector_ans, delimiter);

  //convert each string into an integer.
  int_vector.clear();

  for (vector<string>::iterator string_iter = string_vector_ans.begin();
    string_iter != string_vector_ans.end();
    ++string_iter) {

    int int_ans;
    from_string<int>(int_ans, *string_iter);
    int_vector.push_back(int_ans);
  }
}

/**
 * gets an vector of doubles from cell where the
 * string in the cell are integers which are separated
 * by a delimiter which is differnt than the column
 * delimiter.  The default delimiter is a comma
 * uses the current_row_ as the row index.
 * clears the double vector before 
 * populating it.
 */
void DelimitedFileReader::getDoubleVectorFromCell(
  const char* column_name, ///< the column name
  std::vector<double>& double_vector, ///<the vector of integers
  char delimiter ///<the delimiter to use
) {
  
  //get the list of strings separated by delimiter
  vector<string> string_vector_ans;

  getStringVectorFromCell(column_name, string_vector_ans, delimiter);

  //convert each string into an integer.
  double_vector.clear();

  for (vector<string>::iterator string_iter = string_vector_ans.begin();
    string_iter != string_vector_ans.end();
    ++string_iter) {

    double double_ans;
    from_string<double>(double_ans, *string_iter);
    double_vector.push_back(double_ans);
  }
}



template <typename T>
void DelimitedFileReader::reorderRows(multimap<T, unsigned int>& sort_indices, BOOLEAN_T ascending) {
  //cout <<"reorderRows: start"<<endl;
  vector<vector<string> > newData;

  for (unsigned int col_idx = 0;col_idx < numCols();col_idx++) {
    vector<string> current_col;
    unsigned int row_idx;
    if (ascending) {
      typename multimap<T, unsigned int>::iterator sort_iter;
      for (sort_iter = sort_indices.begin();
	   sort_iter != sort_indices.end();
	   ++sort_iter) {
	row_idx = sort_iter -> second;
	//cout <<"Getting data_["<<col_idx<<"]["<<row_idx<<"]"<<endl;
	string current_cell = getString(col_idx, row_idx);
	//cout <<"Ans:"<<current_cell<<endl;
	current_col.push_back(current_cell);
      }
    } else {
      typename multimap<T, unsigned int>::reverse_iterator sort_iter;
      for (sort_iter = sort_indices.rbegin();
	   sort_iter != sort_indices.rend();
	   ++sort_iter) {
	row_idx = sort_iter -> second;
	string current_cell = data_[col_idx][row_idx];
	current_col.push_back(current_cell);
      }
    }
    //cout <<"Adding new column"<<endl;
    newData.push_back(current_col);
  }
  //cout <<"reorderRows(): swapping"<<endl;
  data_.swap(newData);
  //newData.clear();
  //cout <<"reorderRows(): Done."<<endl;
}

void DelimitedFileReader::sortByFloatColumn(
  const string& column_name,
  BOOLEAN_T ascending) {

  multimap<FLOAT_T, unsigned int> sort_indices;
  int sort_col_idx = findColumn(column_name); 
  
  if (sort_col_idx == -1) {
    carp(CARP_FATAL,"column %s doesn't exist",column_name.c_str());
  }

  for (unsigned int row_idx=0;row_idx<numRows();row_idx++) {
    //cout <<" Inserting "<<getFloat(sort_col_idx, row_idx)<<endl;
    sort_indices.insert(pair<FLOAT_T, unsigned int>(getFloat(sort_col_idx, row_idx), row_idx));
  }

  reorderRows(sort_indices, ascending);
}

void DelimitedFileReader::sortByIntegerColumn(
  const string& column_name,
  BOOLEAN_T ascending
  ) {
  
  multimap<int, unsigned int> sort_indices;
  int sort_col_idx = findColumn(column_name); 
  
  if (sort_col_idx == -1) {
    carp(CARP_FATAL,"column %s doesn't exist",column_name.c_str());
  }

  for (unsigned int row_idx=0;row_idx<numRows();row_idx++) {
    //cout <<" Inserting "<<getInteger(sort_col_idx, row_idx)<<endl;
    sort_indices.insert(pair<int, unsigned int>(getInteger(sort_col_idx, row_idx), row_idx));
  }

  reorderRows(sort_indices, ascending);
}

void DelimitedFileReader::sortByStringColumn(
  const string& column_name,
  BOOLEAN_T ascending) {
  
  multimap<string, unsigned int> sort_indices;

  int sort_col_idx = findColumn(column_name); 
  
  if (sort_col_idx == -1) {
    carp(CARP_FATAL,"column %s doesn't exist",column_name.c_str());
  }

  for (unsigned int row_idx=0;row_idx<numRows();row_idx++) {
    cout <<" Inserting "<<getInteger(sort_col_idx, row_idx)<<endl;
    sort_indices.insert(pair<string, unsigned int>(getString(sort_col_idx, row_idx), row_idx));
  }

  reorderRows(sort_indices, ascending);

}


void DelimitedFileReader::copyToRow(
  DelimitedFileReader& dest,
  int src_row_idx,
  int dest_row_idx
  ) {

  for (unsigned int src_col_idx=0;src_col_idx < numCols();src_col_idx++) {

    int dest_col_idx = dest.findColumn(getColumnName(src_col_idx));
    if (dest_col_idx != -1) {
      dest.setString(dest_col_idx, dest_row_idx, getString(src_col_idx, src_row_idx));
    }
  }
}



/*Iterator functions.*/
/**
 * resets the current_row_ index to 0.
 */
void DelimitedFileReader::reset() {

  current_row_ = 0;
}

/**
 * increments the current_row_, 
 */
void DelimitedFileReader::next() {
  if (current_row_ < numRows())
    current_row_++;
}


/**
 * \returns whether there are more rows to 
 * iterate through
 */
BOOLEAN_T DelimitedFileReader::hasNext() {
  return current_row_ < numRows();
}


/**
 *Allows object to be printed to a stream
 */
std::ostream &operator<< (std::ostream& os, DelimitedFileReader& delimited_file) {

  //find the maximum number of rows.
  unsigned int maxRow = 0;
  for (unsigned int col_idx=0; col_idx < delimited_file.numCols(); col_idx++) {
    maxRow = max(maxRow, delimited_file.numRows(col_idx));
  }

  //print out the header if it exists.
  if (delimited_file.column_names_.size() != 0) {
    os << delimited_file.column_names_[0];
    for (unsigned int col_idx=1; col_idx < delimited_file.column_names_.size(); col_idx++) {
      os << "\t" << delimited_file.column_names_[col_idx];
    }
    os << endl;
  }

  //print out all rows, using \t when
  //the row goes past the current column
  //size.
  for (unsigned int row_idx=0; row_idx<maxRow; row_idx++) {
    if (row_idx < delimited_file.numRows(0)) {
      os << delimited_file.getString((unsigned int)0, row_idx);
    } else {
      os << "\t";
    }
    for (unsigned int col_idx=1;col_idx<delimited_file.numCols();col_idx++) {
      os <<"\t";
      if (row_idx < delimited_file.numRows(col_idx))
	os << delimited_file.getString(col_idx, row_idx);
    }
    os << endl;
  }

  return os;
}
