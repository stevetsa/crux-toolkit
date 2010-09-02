#ifndef MPSM_CHARGEMAP_H
#define MPSM_CHARGEMAP_H

#include <map>
#include <vector>

#include "MPSM_MatchCollection.h"

class MPSM_ChargeMap : public std::map<ChargeIndex, std::vector<MPSM_MatchCollection> > {

  public:
    void insert(ChargeIndex& charge_index);
    void insert(std::vector<MPSM_MatchCollection>& match_collections);
    void insert(MPSM_ChargeMap& mpsm_chargemap);
    void insert(MPSM_MatchCollection& match_collection, int match_collection_idx);
    void insert(MPSM_Match& new_match, int match_collection_idx);

    void sortMatches(SCORER_TYPE_T sort_type);

    void clearMap();

    void calcDeltaCN();
    void calcZScores();

    virtual ~MPSM_ChargeMap();


};









#endif
