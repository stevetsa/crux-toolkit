#ifndef MPSM_CHARGEMAP_H
#define MPSM_CHARGEMAP_H

#include <map>
#include <vector>

#include "MPSM_MatchCollection.h"

class MPSM_ZStateMap : public std::map<ZStateIndex, std::vector<MPSM_MatchCollection> > {

  public:
    void insert(ZStateIndex& charge_index);
    void insert(std::vector<MPSM_MatchCollection>& match_collections);
    void insert(MPSM_ZStateMap& mpsm_chargemap);
    void insert(MPSM_MatchCollection& match_collection, int match_collection_idx);
    void insert(MPSM_Match& new_match, int match_collection_idx);

    void sortMatches(SCORER_TYPE_T sort_type);

    void clearMap();

    void calcDeltaCN();
    void calcZScores();
    void calcRanks();
    bool visited(MPSM_Match& match, int match_collection_idx);


    virtual ~MPSM_ZStateMap();


};









#endif
