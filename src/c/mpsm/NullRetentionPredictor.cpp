#include "NullRetentionPredictor.h"

using namespace std;


NullRetentionPredictor::NullRetentionPredictor() {
}

NullRetentionPredictor::~NullRetentionPredictor() {
}

FLOAT_T NullRetentionPredictor::predictRTime(Match* match) {

  return 0.0;
}

FLOAT_T NullRetentionPredictor::predictRTimeS(const char* sequence) {
  return 0.0;
}

FLOAT_T NullRetentionPredictor::predictRTimeS(const char* sequence, int N) {
  return 0.0;
}
