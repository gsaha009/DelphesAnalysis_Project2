#ifndef __MVASkim__h
#define __MVASkim__h

#include <fstream>
#include <string>

class TTree;
class TFile;

typedef struct  
{
  int event;
  float event_wt;
  
  int nleptons;
  int njets;
  int nbjets;
  int nljets;
  int ntauh;

  float px_jet1;
  float py_jet1;
  float pz_jet1;
  float energy_jet1;
  float mass_jet1;
  float pt_jet1;
  float eta_jet1;

  float px_jet2;
  float py_jet2;
  float pz_jet2;
  float energy_jet2;
  float mass_jet2;
  float pt_jet2;
  float eta_jet2;

  float px_jet3;
  float py_jet3;
  float pz_jet3;
  float energy_jet3;
  float mass_jet3;
  float pt_jet3;
  float eta_jet3;

  float px_jet4;
  float py_jet4;
  float pz_jet4;
  float energy_jet4;
  float mass_jet4;
  float pt_jet4;
  float eta_jet4;

  float px_jet5;
  float py_jet5;
  float pz_jet5;
  float energy_jet5;
  float mass_jet5;
  float pt_jet5;
  float eta_jet5;

  float px_jet6;
  float py_jet6;
  float pz_jet6;
  float energy_jet6;
  float mass_jet6;
  float pt_jet6;
  float eta_jet6;

  float px_tauh;
  float py_tauh;
  float pz_tauh;
  float energy_tauh;
  float pt_tauh;
  float eta_tauh;
  float phi_tauh;

  float px_met;
  float py_met;
  float pz_met;
  float energy_met;
  float met;
  float phi_met;

  float px_leadbj;
  float py_leadbj;
  float pz_leadbj;
  float energy_leadbj;
  float pt_leadbj;
  float eta_leadbj;
  float phi_leadbj;

  float px_leadlj;
  float py_leadlj;
  float pz_leadlj;
  float energy_leadlj;
  float pt_leadlj;
  float eta_leadlj;
  float phi_leadlj;

  float px_xlep;
  float py_xlep;
  float pz_xlep;
  float energy_xlep;
  float pt_xlep;
  float eta_xlep;
  float phi_xlep;

  float px_wlep;
  float py_wlep;
  float pz_wlep;
  float energy_wlep;
  float pt_wlep;
  float eta_wlep;
  float phi_wlep;

  float met_along_tauh;
  float frac_met_along_tauh;

  float dr_xlep_tauh;
  float dphi_xlep_tauh;
  float deta_xlep_tauh;

  float dr_xlep_leadlj;
  float dphi_xlep_leadlj;
  float deta_xlep_leadlj;

  float dr_xlep_leadj;
  float dphi_xlep_leadj;
  float deta_xlep_leadj;

  float dr_tauh_leadj;
  float dphi_tauh_leadj;
  float deta_tauh_leadj;

  float dr_tauh_leadlj;
  float dphi_tauh_leadlj;
  float deta_tauh_leadlj;

  float dphi_met_tauh;
  float dphi_met_leadlj;
  float dphi_met_xlep;

  float dr_min_xlep_jets;
  float dphi_min_xlep_jets;
  float dr_max_xlep_jets;
  float dphi_max_xlep_jets;

  float mt_xlep_met;
  float m_coll_xlep;
  float effectiveMass_t1;
  float invMass_x_test;
  float mt_x;

  float dphi_met_leadbj;
  float dphi_met_wlep;

  float dr_wlep_leadbj;
  float dphi_wlep_leadbj;
  float deta_wlep_leadbj;

  float mt_wlep_met;
  float m_coll_wlep;
  float effectiveMass_t2;

  float dr_tauh_leadbj;
  float dphi_tauh_leadbj;
  float deta_tauh_leadbj;

  float dr_xlep_leadbj;
  float dphi_xlep_leadbj;
  float deta_xlep_leadbj;

  float dr_xlep_wlep;
  float dphi_xlep_wlep;
  float deta_xlep_wlep;
  float scalarsumpt_xlep_wlep;
  float vectorsumpt_xlep_wlep;

  float dr_wlep_leadlj;
  float dphi_wlep_leadlj;
  float deta_wlep_leadlj;

  float dr_wlep_tauh;
  float dphi_wlep_tauh;
  float deta_wlep_tauh;

  float dr_min_wlep_jets;   
  float dphi_min_wlep_jets;
  float dr_max_wlep_jets;   
  float dphi_max_wlep_jets;

  float dr_min_jets;
  float dr_max_jets;
  float dr_leadbj_leadlj;
  float dphi_leadbj_leadlj;
  float deta_leadbj_leadlj;

  float scalarsumpt_vis;
  float scalarsumpt;
  float effectivemass_vis;
  float effectivemass;
  float ht_jets;  
  float smin;
  float costhetaS_xlep_tau;
  float costhetaS_xlep_tauh;

} TreeVariablesDL;


typedef struct  
{
  int event;
  float event_wt;

  int nleptons;
  int njets;
  int nbjets;
  int nljets;
  int ntauh;

  float px_jet1;
  float py_jet1;
  float pz_jet1;
  float energy_jet1;
  float mass_jet1;
  float pt_jet1;
  float eta_jet1;

  float px_jet2;
  float py_jet2;
  float pz_jet2;
  float energy_jet2;
  float mass_jet2;
  float pt_jet2;
  float eta_jet2;

  float px_jet3;
  float py_jet3;
  float pz_jet3;
  float energy_jet3;
  float mass_jet3;
  float pt_jet3;
  float eta_jet3;

  float px_jet4;
  float py_jet4;
  float pz_jet4;
  float energy_jet4;
  float mass_jet4;
  float pt_jet4;
  float eta_jet4;

  float px_jet5;
  float py_jet5;
  float pz_jet5;
  float energy_jet5;
  float mass_jet5;
  float pt_jet5;
  float eta_jet5;

  float px_jet6;
  float py_jet6;
  float pz_jet6;
  float energy_jet6;
  float mass_jet6;
  float pt_jet6;
  float eta_jet6;

  float px_tauh;
  float py_tauh;
  float pz_tauh;
  float energy_tauh;
  float pt_tauh;
  float eta_tauh;
  float phi_tauh;

  float px_met;
  float py_met;
  float pz_met;
  float energy_met;
  float pt_met;
  float phi_met;

  float px_leadbj;
  float py_leadbj;
  float pz_leadbj;
  float energy_leadbj;
  float pt_leadbj;
  float eta_leadbj;
  float phi_leadbj;

  float px_toplj;
  float py_toplj;
  float pz_toplj;
  float energy_toplj;
  float pt_toplj;
  float eta_toplj;
  float phi_toplj;

  float px_wjet1;
  float py_wjet1;
  float pz_wjet1;
  float energy_wjet1;
  float pt_wjet1;
  float eta_wjet1;
  float phi_wjet1;

  float px_wjet2;
  float py_wjet2;
  float pz_wjet2;
  float energy_wjet2;
  float pt_wjet2;
  float eta_wjet2;
  float phi_wjet2;

  float px_xlep;
  float py_xlep;
  float pz_xlep;
  float energy_xlep;
  float pt_xlep;
  float eta_xlep;
  float phi_xlep;

  float px_tau;
  float py_tau;
  float pz_tau;
  float energy_tau;
  float pt_tau;
  float eta_tau;
  float phi_tau;

  float px_nu;
  float py_nu;
  float pz_nu;
  float energy_nu;
  float pt_nu;
  float eta_nu;
  float phi_nu;

  float px_w;
  float py_w;
  float pz_w;
  float energy_w;
  float pt_w;
  float eta_w;
  float phi_w;

  float px_chi;
  float py_chi;
  float pz_chi;
  float energy_chi;
  float pt_chi;
  float eta_chi;
  float phi_chi;

  float px_tSM;
  float py_tSM;
  float pz_tSM;
  float energy_tSM;
  float pt_tSM;
  float eta_tSM;
  float phi_tSM;

  float px_tBSM;
  float py_tBSM;
  float pz_tBSM;
  float energy_tBSM;
  float pt_tBSM;
  float eta_tBSM;
  float phi_tBSM;

  // from one leg
  float invm_wj1_wj2;
  float dr_wj1_wj2;
  float dphi_wj1_wj2;
  float deta_wj1_wj2;
  float dr_w_leadbj;
  float dphi_w_leadbj;
  float dr_wj1_leadbj;
  float dphi_wj1_leadbj;
  float dr_wj2_leadbj;
  float dphi_wj2_leadbj;
  float invm_w_leadbj;
  float dphi_met_leadbj;
  float dphi_met_wj1;
  float dphi_met_wj2;
  float dphi_met_w;

  // from other leg
  float m_coll_x;
  float invm_x_toplj;
  float dr_xlep_toplj;
  float dphi_xlep_toplj;
  float dr_tauh_toplj;
  float dphi_tauh_toplj;
  float dr_x_toplj;
  float dphi_x_toplj;
  float dr_xlep_tauh;
  float dphi_xlep_tauh;
  float dphi_met_toplj;
  float dphi_met_xlep;
  float dphi_met_tauh;
  // from both leg
  float dr_xlep_wj1;
  float dphi_xlep_wj1;
  float deta_xlep_wj1;
  float dr_xlep_wj2;
  float dphi_xlep_wj2;
  float deta_xlep_wj2;
  float dr_xlep_w;
  float dphi_xlep_w;
  float deta_xlep_w;
  float dr_xlep_leadbj;
  float dphi_xlep_leadbj; 
  float deta_xlep_leadbj; 
  float dr_tauh_wj1;
  float dphi_tauh_wj1;
  float deta_tauh_wj1;
  float dr_tauh_wj2;
  float dphi_tauh_wj2;
  float deta_tauh_wj2;
  float dr_tauh_leadbj;
  float dphi_tauh_leadbj;
  float deta_tauh_leadbj;
  float dr_leadbj_toplj;
  float dphi_leadbj_toplj;
  float deta_leadbj_toplj;
  float dr_w_toplj;
  float dphi_w_toplj;
  float deta_w_toplj;
  float dr_x_leadbj;
  float dphi_x_leadbj;
  float deta_x_leadbj;
  float dr_x_w;
  float dphi_x_w;
  float deta_x_w;
  float dr_t1_t2;
  float dphi_t1_t2;
  float deta_t1_t2;
  float total_mass;
  float total_vector_pt;
  float total_scalar_pt;
  float mt_xlep_met;
  // di-jets variables
  float dr_min_jets;
  float dr_max_jets;
  float ht_jets;
  float smin;
  float costhetaS_xlep_tau;
  float costhetaS_xlep_tauh;
  float minDr_XlepJets;
  float maxDr_XlepJets;
  
} TreeVariablesSL;

class MVASkim {
    
public:

  MVASkim(const std::string& filename, bool SL, bool DL);
  //virtual ~MVASkim() = default;
  virtual ~MVASkim() {};

  void fill(const TreeVariablesDL& varList);
  void fill(const TreeVariablesSL& varList);
  void close();
  void close(bool isDL, bool isSL);

  TFile* _mvaFile;
  TTree* _treeDL;
  TTree* _treeSL;

  TreeVariablesDL _varListDL;
  TreeVariablesSL _varListSL;
};
#endif
