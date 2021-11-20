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

  float pt_tauh1;
  float eta_tauh1;
  float phi_tauh1;
  float met;
  float pt_bjet1;
  float eta_bjet1;
  float phi_bjet1;
  float pt_ljet1;
  float eta_ljet1;
  float phi_ljet1;
  float pt_lep1;
  float pt_lep2;
  float eta_lep1;
  float eta_lep2;
  float phi_lep1;
  float phi_lep2;

  float dr_lep1lep2;
  float dphi_lep1lep2;
  float deta_lep1lep2;

  float dr_lep1tauh;
  float dphi_lep1tauh;
  float deta_lep1tauh;
  float dr_lep2tauh;
  float dphi_lep2tauh;
  float deta_lep2tauh;

  float dr_tauhjet;
  float dphi_tauhjet;
  float deta_tauhjet;
  float dr_tauhljet;
  float dphi_tauhljet;
  float deta_tauhljet;
  float dr_tauhbjet;
  float dphi_tauhbjet;
  float deta_tauhbjet;

  float dr_lep1bjet;
  float dr_lep2bjet;
  float dphi_lep1bjet;
  float dphi_lep2bjet;
  float deta_lep1bjet;
  float deta_lep2bjet;

  float dr_lep1ljet;
  float dr_lep2ljet;
  float dphi_lep1ljet;
  float dphi_lep2ljet;
  float deta_lep1ljet;
  float deta_lep2ljet;

  float dphi_lep1met;
  float dphi_lep2met;
  float mt_lep1met;
  float mt_lep2met;
  float mt_muon1met;
  
  float dr_min_jets;
  float dr_max_jets;
  float dr_bjetljet;
  float dphi_bjetljet;
  float deta_bjetljet;
  
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
