#include <cstdlib>
#include <string>
#include <iostream>

#include "ChargeIndex.h"
#include "DelimitedFile.h"
#include "BigSmallReal.h"

using namespace std;

void usage() {
  cout<<"usage(): CountCandidates <file> <spec count>"<<endl;
}

BigSmallReal calculateCount(map<ChargeIndex, BigSmallReal>& spsm_count, ChargeIndex& new_charge_index);
BigSmallReal getTotalCounts(map<ChargeIndex, BigSmallReal>& spsm_count, int max_peptides);

BigSmallReal choose(BigSmallReal count, int n);


int main(int argc, char** argv) {


  if (argc != 3) {
    usage();
    exit(-1);
  }
  //cout <<"Reading "<<argv[1]<<endl;
  DelimitedFile infile(argv[1], true);

  long spec_count = 0;
  DelimitedFile::from_string<long>(spec_count, argv[2]);

  int max_peptides = 3;

  int last_spec = infile.getInteger("scan");
  map<ChargeIndex, BigSmallReal> spsm_count;

  //cout <<"scan\tmatches/spectrum"<<endl;

  BigSmallReal total_total_count;

  do {
    int current_spec = infile.getInteger("scan");
    if (current_spec == last_spec) {
      int charge = infile.getInteger("charge");
      int matches = infile.getInteger("matches/spectrum");
      BigSmallReal bmr_matches(matches);
      ChargeIndex charge_index;
      charge_index.add(charge);
      spsm_count.insert(make_pair(charge_index, bmr_matches));
    } else {
      //cout <<"Getting total counts"<<endl;
      BigSmallReal total_count = getTotalCounts(spsm_count, max_peptides);
      total_total_count.add(total_count);
      //cout << last_spec << "\t" << total_count.getLog() << "\t" << total_total_count.getLog()<<endl;
      spsm_count.clear();
      last_spec = current_spec;
      int charge = infile.getInteger("charge");
      int matches = infile.getInteger("matches/spectrum");
      BigSmallReal bmr_matches(matches);
      ChargeIndex charge_index;
      charge_index.add(charge);
      spsm_count.insert(make_pair(charge_index, bmr_matches));
    }
    
    //cout <<"Reading next entry"<<endl;
    infile.next();
  } while (infile.hasNext());
  BigSmallReal total_count = getTotalCounts(spsm_count, max_peptides);
  total_total_count.add(total_count);
  //cout << last_spec << "\t" << total_count << endl;

  //cout <<"total spec:"<<spec_count<<endl;

  cout << (total_total_count.getLog() - log10(spec_count))<<endl;

}


BigSmallReal getTotalCounts(map<ChargeIndex, BigSmallReal>& spsm_count, int max_peptides) {

      //calculate the total candidates for this combination and print it out.
      map<ChargeIndex, BigSmallReal> mpsm_count = spsm_count;
      for (int npeptides=2; npeptides <= max_peptides; npeptides++) {
        map<ChargeIndex, BigSmallReal>::iterator map_iter;
        map<ChargeIndex, BigSmallReal>::iterator map_iter2;
        
        for (map_iter2 = spsm_count.begin();
          map_iter2 != spsm_count.end();
          ++map_iter2) {

            ChargeIndex spsm_charge_index = map_iter2 -> first;
            map<ChargeIndex, BigSmallReal> new_mpsms;
            for (map_iter = mpsm_count.begin();
              map_iter != mpsm_count.end();
              ++map_iter) {
              
                ChargeIndex charge_index = map_iter -> first;
              if ((charge_index.size() == (unsigned int)(npeptides - 1))) {
              ChargeIndex new_charge_index = charge_index;
              new_charge_index.add(spsm_charge_index[0]);
              BigSmallReal count = calculateCount(spsm_count, new_charge_index);
              //cout <<"Count for "<<new_charge_index<<" "<<count<<endl;
              mpsm_count.insert(make_pair(new_charge_index, count));
            }
          }
        }

      }

      BigSmallReal total_count;
      for (map<ChargeIndex, BigSmallReal>::iterator map_iter = 
        mpsm_count.begin();
        map_iter != mpsm_count.end();
        ++map_iter) {
        ChargeIndex index = map_iter -> first;
        //cout<<"Adding "<<index<<" "<<map_iter -> second<<endl;
        total_count.add(map_iter -> second);
        //cout<<"count is now:"<<total_count.getLog()<<endl;
      }
  return total_count;
}

BigSmallReal calculateCount(map<ChargeIndex, BigSmallReal>& spsm_count, ChargeIndex& new_charge_index) {

  BigSmallReal ans;
  ans.add(1);

  map<ChargeIndex, BigSmallReal>::iterator map_iter;

  for (map_iter = spsm_count.begin();
    map_iter != spsm_count.end();
    ++map_iter) {

    int num_charge = new_charge_index.numCharge(map_iter -> first[0]);

    if (num_charge != 0) {
      BigSmallReal c = choose(map_iter->second, num_charge);
      //cout<<"Choose:"<<c.get()<<endl;
      ans.mult(c);
      //cout <<"ans is now:"<<ans.get()<<endl;
    }
  }

  return ans;
}


BigSmallReal choose(BigSmallReal count, int n) {

  if (n == 1) return count;

  double dcount = count.get();

  BigSmallReal ans = count;

  for (long idx=2;idx<=n;idx++) {
    ans.mult((dcount - (double)idx + 1.0) / (double)idx);
  }

  //cout <<"Choose "<<count.get()<<","<<n<<"="<<ans.get()<<endl;

  return ans;
}
