#include "BigSmallReal.h"
#include <cmath>
#include <iostream>

using namespace std;

BigSmallReal::BigSmallReal() {
  ;
}

BigSmallReal::BigSmallReal(double a) {
  add(a);
}

void BigSmallReal::zero() {
  data.clear();
}


void BigSmallReal::add(BigSmallReal& a) {

  map<int, double>::iterator iter1;
  map<int, double>::iterator iter2;

  for (iter1 = a.data.begin();
    iter1 != a.data.end();
    ++iter1) {

    iter2 = data.find(iter1 -> first);
    if (iter2 == data.end()) {

      data[iter1 -> first] = iter1 -> second;
    } else {

      iter2 -> second += iter1 -> second;
    }
  }
}

void BigSmallReal::add(double a) {
  if (a != 0) {
    addLog(log10(a));
  }
}

void BigSmallReal::mult(double a) {
  map<int,double>:: iterator iter;

  for (iter = data.begin();
    iter != data.end();
    ++iter) {

    iter -> second *= a;
  }

}

void BigSmallReal::multLog(double la) {
  int exp_part = (int)floor(la);
  double frac_part = pow(10.0, (la - exp_part));

  map<int, double, greater<int> > temp = this -> data;

  zero();

  for (map<int, double>::iterator iter = temp.begin();
    iter != temp.end();
    ++iter) {

    data[iter -> first + exp_part] = iter -> second;
  }

  mult(frac_part);
}

void BigSmallReal::mult(BigSmallReal& a) {

  //cout <<"mult "<<get()<<" by "<<a.get()<<endl;

  map<int, double>::iterator iter1;
  map<int, double>::iterator iter2;

  map<int, double, greater<int> > temp = this -> data;
  zero();

  for (iter1 = a.data.begin(); iter1 != a.data.end(); ++iter1) {
    for (iter2 = temp.begin(); iter2 != temp.end(); ++iter2) {

      int new_exp = iter1 -> first + iter2 -> first;
      double new_frac = iter1 -> second * iter2 -> second;
      /*
      cout << "exp1:"
        << iter1 -> first
        << " exp2:" 
        << iter2->first
        << " frac1:"
        << iter1->second
        << " frac2:"
        << iter2->second;

      cout <<" new_exp:"<<new_exp<<" new_frac:"<<new_frac<<endl;
      */
      data[new_exp] += new_frac;

    }
  }
}



void BigSmallReal::addLog(double la) {
  int exp_part = (int)floor(la);
  double frac_part = pow(10.0, (la - exp_part));

  //cout <<exp(la) <<" = "<< (exp(exp_part) * frac_part) <<endl;
  
  map<int,double>::iterator pos = data.find(exp_part);
  if (pos == data.end()) {
    //not found, create it.
    data.insert(make_pair<int, double>(exp_part, frac_part));
  }
  else {
    //found, add fracs.
    pos -> second += frac_part;
  }
}

void BigSmallReal::print() {
  map<int,double>::iterator iter;
  cout<<"0";
  for (iter = data.begin();iter != data.end();++iter)
    cout <<" + " << iter -> second<<"*10^"<< iter -> first;
  cout<<endl;
}

double BigSmallReal::getLog() {
  refactor();
  //use the equation i deived.
  //1, find the largest exponent, should be the first number.
  int int_part = data.begin() -> first;

  //then divide the smaller exponents, which should be less significant as we go..

  double accumulator = 0;
  int count = 0;
  map<int,double>::iterator iter;
  for (iter=data.begin(); iter != data.end() ;++iter) {
    accumulator += pow(10.0, iter -> first - int_part)*iter -> second;
   
  }
  
  double ans = (double) int_part + log10(accumulator);

  return ans;
}

double BigSmallReal::get() {

  double lans = getLog();
  return pow(10.0, lans);

}


void BigSmallReal::refactor() {
  //start with the smallest exponent and work our way up.
  map<int, double>::reverse_iterator iter;
  for (iter=data.rbegin(); iter != data.rend(); ++iter) {
    int next_exp = iter -> first + 1;
    if (iter -> second > 10) {
      double incr = floor(iter -> second / 10.0);
      iter -> second -= incr * 10.0;
      //find the next exponent.
      //if it exists, just add,
      //if not create it, and restart.
      map<int, double>::iterator iter2 = data.find(next_exp);
      if (iter2 == data.end()) {
	//create it 
	data.insert(make_pair(next_exp, (double)incr));
	iter = data.rbegin();
      }
      else
	{
	  iter2 -> second += incr;
	}
    }
  }

}

/*
int main(void) {

  BigSmallReal a;

  a.add(50);
  a.add(50);
  a.add(50);
  a.print();
  a.refactor();
  a.print();

}
*/
