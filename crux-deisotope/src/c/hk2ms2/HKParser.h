#ifndef _HKPARSER_H
#define _HKPARSER_H

#include <vector>
#include <iostream>
#include <fstream>
#include <cstring>
#include <stdlib.h>

using namespace std;

typedef struct sPep{
  int charge;
  float intensity;
  double monoMass;
	double basePeak;
  double xCorr;
	char mods[32];
} sPep;

typedef struct sScan {
  vector<sPep> *vPep;
	int scanNum;
	char file[256];
	float rTime;

	//Constructors & Destructor
	sScan(){vPep = new vector<sPep>;}
	sScan(const sScan& s){
		vPep = new vector<sPep>;
		for(unsigned int i=0;i<s.vPep->size();i++) vPep->push_back(s.vPep->at(i));
		scanNum = s.scanNum;
		strcpy(file,s.file);
		rTime=s.rTime;
	}
	~sScan(){delete vPep;}

  //Copy operator
	sScan& operator=(const sScan& s){
		if(this!=&s){
			delete vPep;
			vPep = new vector<sPep>;
			for(unsigned int i=0;i<s.vPep->size();i++) vPep->push_back(s.vPep->at(i));
			scanNum = s.scanNum;
			strcpy(file,s.file);
			rTime=s.rTime;
		}
		return *this;
	}

  //Clear
	void clear(){
		delete vPep;
		vPep = new vector<sPep>;
	}

  static int compareIntRev(const void *p1, const void *p2){
    const sPep d1 = *(sPep *)p1;
    const sPep d2 = *(sPep *)p2;
    if(d1.intensity>d2.intensity) return -1;
    else if(d1.intensity<d2.intensity) return 1;
    else return 0;
  }

  void sortIntRev(){
    if(vPep->size()==0) return;
    qsort(&vPep->at(0),vPep->size(),sizeof(sPep),compareIntRev);
  }

} sScan;

class HKParser {
public:

  //Constructors and destructors
  HKParser();
  ~HKParser();

  //Operator overrides
  sScan& operator[ ](const unsigned int& i);

  //Automation
  bool loadHK(char* in);

  //Data manipulation functions
  unsigned int size();

protected:
private:

  //Data Members
  vector<sScan> hkData;

};

#endif

