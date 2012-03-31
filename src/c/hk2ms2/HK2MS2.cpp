/* 
	HK2MS2 - A tool for converting spectra de-isotoped with Hardklor to MS2 files.
	Author: Michael Hoopmann, Institute for Systems Biology

	Version 1.0: Mar 16, 2011.
	Version 1.11: Mar 6, 2012.

	Copyright 2011 Institute for Systems Biology (ISB). All rights reserved.

	Redistribution and use in source and binary forms, with or without modification, are
	permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of
      conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list
      of conditions and the following disclaimer in the documentation and/or other materials
      provided with the distribution.

	THIS SOFTWARE IS PROVIDED BY ISB ``AS IS'' AND ANY EXPRESS OR IMPLIED
	WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
	FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL ISB OR
	CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
	CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
	SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
	ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
	NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
	ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

	The views and conclusions contained in the software and documentation are those of the
	authors and should not be interpreted as representing official policies, either expressed
	or implied, of ISB.

*/

#include "HKParser.h"
#include "MSReader.h"

#ifdef CRUX
#include "../../external/MSToolkit/include/Spectrum.h"
//using namespace MSToolkit;
#include "CruxHK2MS2Application.h"
#endif

#include <cmath>

#define GAUSSCONST 5.5451774444795623

double calcFWHM(double mz,double res,int iType);
void centroid(MSToolkit::Spectrum& s, MSToolkit::Spectrum& out, int res);


/* Notes from MH:
Right now centroid is hard-coded to assume Orbitrap. The code is there, though, if you
need to switch to FT instrumentation.
*/

#ifdef CRUX
int CruxHK2MS2Application::hk2ms2Main(int argc, char* argv[]){
#else
int main(int argc, char* argv[]) {
#endif	
	if(argc==1){
		cout << "HK2MS2 version 1.1" << endl;
		cout << "Copyright 2011-2012, Michael Hoopmann, Institute for Systems Biology\n" << endl;
		cout << "USAGE: HK2MS2 <Hardklor File> <mzXML/RAW File> <output file base name (no extension)> [options]" << endl;
		cout << "\nOutput is monoisotopic mass by default." << endl;
		cout << "\nOptions:" << endl;
		cout << "\t-h\tOutput m+H" << endl;
		cout << "\t-z\tOutput m/z" << endl;
		return 1;
	}

	int i,j,k;
	int resolution;
	int iOutputType=0;
	bool bCentroid=false;
	char baseName[256];
	char str[256];
  char outName[256];

	strcpy(str,argv[2]);
	i=strrchr(str,'.')-str;
	strncpy(baseName,str,i);
	baseName[i]='\0';
  strcpy(outName,argv[3]);
	//cout << baseName << "\t" << i << endl;

	//check the options
	for(i=4;i<argc;i++){
		if(strcmp(argv[i],"-c")==0) {
			bCentroid=true;
			i++;
			resolution=atoi(argv[1]);
		}
		if(strcmp(argv[i],"-h")==0) iOutputType=1;
		if(strcmp(argv[i],"-z")==0) iOutputType=2;
	}

	//Load the Hardklor results
	HKParser h;
	h.loadHK(argv[1]);

	//The next block has three major functions
	// 1. Iterate through the Hardklor results
	// 2. Load an MS/MS scan and deisotope using the Hardklor results
	// 3. Output the MS/MS scan as a DTA file
	double mz,ppm;
	double d;
	double precursorMz;
	float intensity;
	int z;
	bool bMatch;
	MSToolkit::MSReader r;
	MSToolkit::MSReader w;
	MSToolkit::MSObject o;
	MSToolkit::Spectrum s;
	MSToolkit::Spectrum c;
	MSToolkit::Spectrum deIso;

	//Open the RAW and read the first scan. This
	//commits parts of the RAW file to memory for faster traversal.
	//r.setFilter(MS2);
	r.readFile(argv[2],s);
        //cout << s.getScanNumber()<<endl;
  char outFile1[256];

  sprintf(outFile1,"%s.ms2",outName);
	w.writeFile(outFile1,MSToolkit::ms2,o);	
	
        i = 0;
        while ((s.getScanNumber() != 0) && (i < h.size())) {
          //cerr << s.getScanNumber() << " " << i << " " << h[i].scanNum << " " <<h.size()<<endl;
          while (i < h.size() && (h[i].scanNum < s.getScanNumber())) {
            //cerr <<"incr i"<<endl;
            i++;
          }
          while (s.getScanNumber() != 0 && (s.getScanNumber() < h[i].scanNum)) {
            //cerr <<"incr s"<<endl;
            r.readFile(NULL, s);
          }
          
          if (s.getScanNumber() == h[i].scanNum) {
          //s.setFileType(MS2);
            //cout << s.getScanNumber() << "\t" << h[i].scanNum << endl;
            if(s.sizeZ()<1) {
              cerr << s.getScanNumber() << "(" << i << ") has no high-res precursor information. skipping" << endl;
	      
	    } else if (h[i].vPep->size() <= 0) {
              cerr << s.getScanNumber() << "(" << i << ") has no hardklor peaks. skipping"<<endl;
            } else {
              s.clearPeaks();
	      for(j=0;j<(int)h[i].vPep->size();j++){
			switch(iOutputType){
				case 1:
					s.add(h[i].vPep->at(j).monoMass+1.00727646677,h[i].vPep->at(j).intensity);
					break;
				case 2:
					z=h[i].vPep->at(j).charge;
					mz=(h[i].vPep->at(j).monoMass+1.00727646677*z)/z;
					s.add(mz,h[i].vPep->at(j).intensity);
					break;
				default:
					s.add(h[i].vPep->at(j).monoMass,h[i].vPep->at(j).intensity);
					break;
			}
		}			

		//I think I do these steps so that the precursor mass is correctly calculated when convertin
		//from .ms2 to .mzXML with msconvert
		//s.setMZ((s.atZ(0).mz+(s.atZ(0).z-1)*1.00727649)/s.atZ(0).z);

		//Put these peaks in order
		if(s.size()>0) s.sortMZ();

		o.add(s);
		if(o.size()>1000){
			w.appendFile(outFile1,o);
			o.clear();
		}
              }
            i++;
            r.readFile(NULL, s);
            
            }
          //cerr << "loop"<<endl;
	}


	w.appendFile(outFile1,o);
	o.clear();
        
	return 0;
}

//Calculates the resolution (FWHM) of a peak
double calcFWHM(double mz,double res,int iType){
	double deltaM;
	switch(iType){
	case 0: //Orbitrap
		deltaM = mz * sqrt(mz) / (20*res);  //sqare root of 400
		break;
	case 1: //TOF
		deltaM = mz / res;
		break;
	case 2: //QIT
		deltaM = res;
		break;
	case 3: //FTICR
	default:
		deltaM = mz * mz / (400*res);
		break;
	}
	return deltaM;
}

//First derivative method, returns base peak intensity of the set
void centroid(MSToolkit::Spectrum& s, MSToolkit::Spectrum& out, int res){
  int i,j;
  float maxIntensity;
  int bestPeak;
  bool bLastPos;

	int nextBest;
	double FWHM;
	MSToolkit::Peak_T centroid;

	out.clear();

  bLastPos=false;
	for(i=0;i<s.size()-1;i++){

    if(s[i].intensity<s[i+1].intensity) {
      bLastPos=true;
      continue;
    } else {
      if(bLastPos){
				bLastPos=false;
	
        //Possible ways to improve this:
				//1. check FWHM - arg! what a computational disaster.
				//2. account for noise - another disaster.

				//find max and add peak
				maxIntensity=0;
				for(j=i;j<i+1;j++){
				  if (s[j].intensity>maxIntensity){
				    maxIntensity=s[j].intensity;
				    bestPeak = j;
				  }
				}

				//Best estimate of Gaussian centroid
				//Get 2nd highest point of peak
				if(bestPeak==s.size()){
					nextBest=bestPeak-1;
				} else if(s[bestPeak-1].intensity > s[bestPeak+1].intensity){
					nextBest=bestPeak-1;
				} else {
					nextBest=bestPeak+1;
				}

				//Get FWHM
				FWHM = calcFWHM(s[bestPeak].mz,(double)res,0);

				//Calc centroid MZ (in three lines for easy reading)
				centroid.mz = pow(FWHM,2)*log(s[bestPeak].intensity/s[nextBest].intensity);
				centroid.mz /= GAUSSCONST*(s[bestPeak].mz-s[nextBest].mz);
				centroid.mz += (s[bestPeak].mz+s[nextBest].mz)/2;

				//Calc centroid intensity
				centroid.intensity=(float)(s[bestPeak].intensity/exp(-pow((s[bestPeak].mz-centroid.mz)/FWHM,2)*GAUSSCONST));

				//some peaks are funny shaped and have bad gaussian fit.
				//if error is more than 10%, keep existing intensity
				if( fabs((s[bestPeak].intensity - centroid.intensity) / centroid.intensity * 100) > 10 ||
            //not a good check for infinity
            centroid.intensity>999999999999.9 ||
            centroid.intensity < 0 ) {
					centroid.intensity=s[bestPeak].intensity;
				}

				//Hack until I put in mass ranges
				if(centroid.mz<0 || centroid.mz>2000) {
					//do nothing if invalid mz
				} else {
					out.add(centroid);
				}
			
      }

    }
  }

}

