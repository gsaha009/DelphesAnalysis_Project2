#ifndef __MVASkim__h
#define __MVASkim__h

#include <fstream>
#include <string>

class TTree;
class TFile;

typedef struct  
{
  int event;
  
  int nLeptons;
  int nJets;
  int nbJets;
  int nlJets;
  int nTauh;

  float px_tauh1;
  float py_tauh1;
  float pz_tauh1;
  float energy_tauh1;
  float pt_tauh1;
  float eta_tauh1;
  float phi_tauh1;

  float px_met;
  float py_met;
  float pz_met;
  float energy_met;
  float met;
  float phi_met;

  float px_bjet1;
  float py_bjet1;
  float pz_bjet1;
  float energy_bjet1;
  float pt_bjet1;
  float eta_bjet1;
  float phi_bjet1;

  float px_ljet1;
  float py_ljet1;
  float pz_ljet1;
  float energy_ljet1;
  float pt_ljet1;
  float eta_ljet1;
  float phi_ljet1;

  float px_Xlep;
  float py_Xlep;
  float pz_Xlep;
  float energy_Xlep;
  float pt_Xlep;
  float eta_Xlep;
  float phi_Xlep;

  float px_Wlep;
  float py_Wlep;
  float pz_Wlep;
  float energy_Wlep;
  float pt_Wlep;
  float eta_Wlep;
  float phi_Wlep;

  float scalarSumpt_XlwpWlep;
  float vectorSumpt_XlwpWlep;
  float dr_XlepWlep;
  float dphi_XlepWlep;
  float deta_XlepWlep;

  float dr_Xleptauh;
  float dphi_Xleptauh;
  float deta_Xleptauh;
  float dr_Wleptauh;
  float dphi_Wleptauh;
  float deta_Wleptauh;

  float dr_tauhjet;
  float dphi_tauhjet;
  float deta_tauhjet;
  float dr_tauhljet;
  float dphi_tauhljet;
  float deta_tauhljet;
  float dr_tauhbjet;
  float dphi_tauhbjet;
  float deta_tauhbjet;

  float dr_min_Xlepjets;   
  float dphi_min_Xlepjets;
  float dr_max_Xlepjets;   
  float dphi_max_Xlepjets;
  float dr_min_Wlepjets;   
  float dphi_min_Wlepjets;
  float dr_max_Wlepjets;   
  float dphi_max_Wlepjets;
  
  float dr_Xlepbjet;
  float dr_Wlepbjet;
  float dphi_Xlepbjet;
  float dphi_Wlepbjet;
  float deta_Xlepbjet;
  float deta_Wlepbjet;

  float dr_Xlepljet;
  float dr_Wlepljet;
  float dphi_Xlepljet;
  float dphi_Wlepljet;
  float deta_Xlepljet;
  float deta_Wlepljet;

  float dphi_Xlepmet;
  float dphi_Wlepmet;
  float mt_Xlepmet;
  float mt_Wlepmet;
  
  float dr_min_jets;
  float dr_max_jets;
  float dr_bjetljet;
  float dphi_bjetljet;
  float deta_bjetljet;

  float scalarSumPtVis;
  float scalarSumPt;
  float effectiveMassVis;
  float effectiveMass;
  float M_coll_IC_Xlep;
  float M_coll_IC_Wlep;
  float HT_Jets;
  
} TreeVariables;

class MVASkim {
    
public:

  MVASkim(const std::string& filename);
  virtual ~MVASkim();

  void fill(const TreeVariables& varList);
  void close();

  TFile* _mvaFile;
  TTree* _tree;

  TreeVariables _varList;
};
#endif
