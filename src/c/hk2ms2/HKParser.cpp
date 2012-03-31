#include "HKParser.h"

//-------------------------------------
//   Constructors and Destructors
//-------------------------------------
HKParser::HKParser(){}
HKParser::~HKParser(){}

sScan& HKParser::operator[ ](const unsigned int& i){
  return hkData[i];
}


bool HKParser::loadHK(char* in){
	FILE *hkr;
	sScan scan;
	sPep pep;
	double td;
	int ti;
	char tag;
	bool firstScan;

	int pepCount=0;

  //Read in the Hardklor results
  firstScan=true;
	hkr = fopen(in,"rt");
	if(hkr==NULL) {
		cout << "Problem reading file." << endl;
		return false;
	}

  hkData.clear();
	while(!feof(hkr)){
    
		tag=fgetc(hkr);

    if(tag=='S') {
			if(firstScan) firstScan=false;
			else hkData.push_back(scan);	
			scan.clear();
      fscanf(hkr,"\t%d\t%f%s\t%lf\t%d\t%lf\n",&scan.scanNum,&scan.rTime,scan.file,&td,&ti,&td);
		} else {
			pepCount++;
			fscanf(hkr,"\t%lf\t%d\t%f\t%lf\t%lf-%lf\t%lf\t%s\t%lf\n", &pep.monoMass,&pep.charge,&pep.intensity,&pep.basePeak,&td,&td,&td,pep.mods,&pep.xCorr);
			scan.vPep->push_back(pep);
		}
	}
  hkData.push_back(scan);
	fclose(hkr);

  cout << pepCount << " peptide signals from " << hkData.size() << " scans." << endl;
	return true;
}

unsigned int HKParser::size(){
	return hkData.size();
}

