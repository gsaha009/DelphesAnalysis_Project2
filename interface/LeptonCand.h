#ifndef __LeptonCand__hh
#define __LeptonCand__hh

#include "TLorentzVector.h"

struct LeptonCand {

  float PT {-999};
  float Eta {-999};
  float Phi {-999};
  int Charge {-999};
  TLorentzVector P4;
  int Flavour {-111}; // 1:muon, 2:electron, -1:others

};
#endif
