#ifndef __MVAnalysis__h
#define __MVAnalysis__h

#include <fstream>
#include <string>

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

typedef struct  
{
  /*
  float dphi_met_tauh;
  float dphi_wlep_leadbj;
  float dphi_xlep_wlep;
  float dr_xlep_tauh;
  float effectivemass;
  float ht_jets;
  float mt_wlep_met;
  float mt_x;
  float met;
  float scalarsumpt;
  float dphi_leadbj_leadlj;
  float m_coll_xlep;
  float m_coll_wlep;
  */
  float pt_tauh;
  float met;
  float pt_leadbj;
  float pt_xlep;
  float pt_wlep;
  float vectorsumpt_xlep_wlep;
  float dr_xlep_wlep;
  float dphi_wlep_tauh;
  float dphi_wlep_leadbj;
  float ht_jets;
  float dr_min_xlep_jets;
  float dr_min_wlep_jets;
  float dphi_xlep_tauh;
  float effectivemass;
  float dr_min_jets;
  float mt_wlep_met;
  float dphi_tauh_leadbj;
  float dphi_met_wlep;
  float dphi_leadbj_leadlj;
  
} InputVariablesDL;

typedef struct  
{
  float pt_tauh;
  float met;
  
} InputVariablesSL;

class MVAnalysis {
    
public:

  MVAnalysis(const std::string& mva_algo, const std::string& xmlfile, bool SL, bool DL);
  virtual ~MVAnalysis() {}

  double evaluate(const std::string& tag, const InputVariablesDL& varList);
  double evaluate(const std::string& tag, const InputVariablesSL& varList);

  InputVariablesDL varListDL_;
  InputVariablesSL varListSL_;
  //  std::unique_ptr<TMVA::Reader> reader_;
  TMVA::Reader* readerDL_;
  TMVA::Reader* readerSL_;
};
#endif
