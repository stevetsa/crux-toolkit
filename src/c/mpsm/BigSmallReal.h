/**
 * @file BigSmallReal.h
 *
 *  Created on: Nov 23, 2008
 *      Author: mcilwain
 *  BigSmallReal represents a number as a sum of exponents.
 *  this allows for really small or really big number to be
 *  stored and manipulated.  This is especially useful when
 *  the final result can be actually stored as a real or 
 *  double value, but the intermediate calculations are way
 *  to big.
 *
 * $Id: BigSmallReal.cpp,v
 * $Log: BigSmallReal.h,v $
 * Revision 1.1  2009/02/23 21:56:59  mcilwain
 * *** empty log message ***
 *
 */
#include <map>

class BigSmallReal {
 private:
  /*We are storing the number as list of exponents
   *in this case the exponent is an integer value (key),
   *which is multiplied by the real part (data).
   *we can use a map for this.*/
  std::map<int, double, std::greater<int> > data;

 public:
  
  BigSmallReal();
  BigSmallReal(double a);
  
  void zero();
  
  /*
   *add functions, add a number to the current value which is in
   * either normal or log format
   */
  void add(double a);
  void add(BigSmallReal& a);

  void addLog(double la);

  /*
   *subtract functions, subtract a number from the current value
   */
  void sub(double a);
  void subLog(double la);


  void mult(double a);
  void mult(BigSmallReal& a);
  void multLog(double la);

  /*
  *getLog: returns the log (base 10) of the number.
  */
  double getLog();

  double get();
  
  void print();

  void refactor();
  
  ~BigSmallReal() {;}






};
