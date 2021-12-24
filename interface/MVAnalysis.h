#ifndef __MVAnalysis__h
#define __MVAnalysis__h

#include <fstream>
#include <string>

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

typedef struct  
{
  float pt_tauh1;
  float met;
  float mt_Wlepmet;
  float pt_bjet1;
  float dphi_Wleptauh;
  float dphi_Xleptauh;
  float pt_Xlep;
  float pt_Wlep;
  float vectorSumpt_XlwpWlep;
  float dr_XlepWlep;
  float dphi_tauhbjet;
  float dr_min_Xlepjets;
  float dr_min_Wlepjets;
  float deta_Wlepbjet;
  float dphi_Xlepmet;
  float dphi_Wlepmet;
  float dr_min_jets;
  float dphi_bjetljet;
  float effectiveMass;
  float HT_Jets;
  
} InputVariables;

class MVAnalysis {
    
public:

  MVAnalysis(const std::string& mva_algo, const std::string& xmlfile);
  virtual ~MVAnalysis() {}

  double evaluate(const std::string& tag, const InputVariables& varList);

  InputVariables varList_;
  //  std::unique_ptr<TMVA::Reader> reader_;
  TMVA::Reader* reader_;
};
#endif
