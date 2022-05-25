#include <iostream>
#include <memory>
#include "TFile.h"
#include "TTree.h"
#include "MVASkim.h"

using std::string;
using std::cout;
using std::endl;

MVASkim::MVASkim(const string& filename, bool SL, bool DL) {
  _mvaFile = TFile::Open(filename.c_str(), "RECREATE", "Skimmed Tree");
  _mvaFile->cd();
  
  // ------------------------------------------------------------- //
  //                               DL Tree                         //
  // ------------------------------------------------------------- //
  if (DL) {
    _treeDL = new TTree("EventsDL", "Flat_Events_TTree_DL");
    
    _treeDL->Branch("event",          &_varListDL.event,         "event/I");
    _treeDL->Branch("event_wt",       &_varListDL.event_wt,      "event_wt/F");
    
    _treeDL->Branch("nleptons",       &_varListDL.nleptons,      "nleptons/I");
    _treeDL->Branch("njets",          &_varListDL.njets,         "njets/I");
    _treeDL->Branch("nbjets",         &_varListDL.nbjets,        "nbjets/I");
    _treeDL->Branch("nljets",         &_varListDL.nljets,        "nljets/I");
    _treeDL->Branch("ntauh",          &_varListDL.ntauh,         "ntauh/I");
    
    _treeDL->Branch("px_jet1",        &_varListDL.px_jet1,       "px_jet1/F");
    _treeDL->Branch("py_jet1",        &_varListDL.py_jet1,       "py_jet1/F");
    _treeDL->Branch("pz_jet1",        &_varListDL.pz_jet1,       "pz_jet1/F");
    _treeDL->Branch("energy_jet1",    &_varListDL.energy_jet1,   "energy_jet1/F");
    _treeDL->Branch("mass_jet1",      &_varListDL.mass_jet1,     "mass_jet1/F");
    _treeDL->Branch("pt_jet1",        &_varListDL.pt_jet1,       "pt_jet1/F");
    _treeDL->Branch("eta_jet1",       &_varListDL.eta_jet1,      "eta_jet1/F");

    _treeDL->Branch("px_jet2",        &_varListDL.px_jet2,       "px_jet2/F");
    _treeDL->Branch("py_jet2",        &_varListDL.py_jet2,       "py_jet2/F");
    _treeDL->Branch("pz_jet2",        &_varListDL.pz_jet2,       "pz_jet2/F");
    _treeDL->Branch("energy_jet2",    &_varListDL.energy_jet2,   "energy_jet2/F");
    _treeDL->Branch("mass_jet2",      &_varListDL.mass_jet2,     "mass_jet2/F");
    _treeDL->Branch("pt_jet2",        &_varListDL.pt_jet2,       "pt_jet2/F");
    _treeDL->Branch("eta_jet2",       &_varListDL.eta_jet2,      "eta_jet2/F");

    _treeDL->Branch("px_jet3",        &_varListDL.px_jet3,       "px_jet3/F");
    _treeDL->Branch("py_jet3",        &_varListDL.py_jet3,       "py_jet3/F");
    _treeDL->Branch("pz_jet3",        &_varListDL.pz_jet3,       "pz_jet3/F");
    _treeDL->Branch("energy_jet3",    &_varListDL.energy_jet3,   "energy_jet3/F");
    _treeDL->Branch("mass_jet3",      &_varListDL.mass_jet3,     "mass_jet3/F");
    _treeDL->Branch("pt_jet3",        &_varListDL.pt_jet3,       "pt_jet3/F");
    _treeDL->Branch("eta_jet3",       &_varListDL.eta_jet3,      "eta_jet3/F");

    _treeDL->Branch("px_jet4",        &_varListDL.px_jet4,       "px_jet4/F");
    _treeDL->Branch("py_jet4",        &_varListDL.py_jet4,       "py_jet4/F");
    _treeDL->Branch("pz_jet4",        &_varListDL.pz_jet4,       "pz_jet4/F");
    _treeDL->Branch("energy_jet4",    &_varListDL.energy_jet4,   "energy_jet4/F");
    _treeDL->Branch("mass_jet4",      &_varListDL.mass_jet4,     "mass_jet4/F");
    _treeDL->Branch("pt_jet4",        &_varListDL.pt_jet4,       "pt_jet4/F");
    _treeDL->Branch("eta_jet4",       &_varListDL.eta_jet4,      "eta_jet4/F");

    _treeDL->Branch("px_jet5",        &_varListDL.px_jet5,       "px_jet5/F");
    _treeDL->Branch("py_jet5",        &_varListDL.py_jet5,       "py_jet5/F");
    _treeDL->Branch("pz_jet5",        &_varListDL.pz_jet5,       "pz_jet5/F");
    _treeDL->Branch("energy_jet5",    &_varListDL.energy_jet5,   "energy_jet5/F");
    _treeDL->Branch("mass_jet5",      &_varListDL.mass_jet5,     "mass_jet5/F");
    _treeDL->Branch("pt_jet5",        &_varListDL.pt_jet5,       "pt_jet5/F");
    _treeDL->Branch("eta_jet5",       &_varListDL.eta_jet5,      "eta_jet5/F");

    _treeDL->Branch("px_jet6",        &_varListDL.px_jet6,       "px_jet6/F");
    _treeDL->Branch("py_jet6",        &_varListDL.py_jet6,       "py_jet6/F");
    _treeDL->Branch("pz_jet6",        &_varListDL.pz_jet6,       "pz_jet6/F");
    _treeDL->Branch("energy_jet6",    &_varListDL.energy_jet6,   "energy_jet6/F");
    _treeDL->Branch("mass_jet6",      &_varListDL.mass_jet6,     "mass_jet6/F");
    _treeDL->Branch("pt_jet6",        &_varListDL.pt_jet6,       "pt_jet6/F");
    _treeDL->Branch("eta_jet6",       &_varListDL.eta_jet6,      "eta_jet6/F");

    _treeDL->Branch("px_tauh",        &_varListDL.px_tauh,       "px_tauh/F");
    _treeDL->Branch("py_tauh",        &_varListDL.py_tauh,       "py_tauh/F");
    _treeDL->Branch("pz_tauh",        &_varListDL.pz_tauh,       "pz_tauh/F");
    _treeDL->Branch("energy_tauh",    &_varListDL.energy_tauh,   "energy_tauh/F");
    _treeDL->Branch("pt_tauh",        &_varListDL.pt_tauh,       "pt_tauh/F");
    _treeDL->Branch("eta_tauh",       &_varListDL.eta_tauh,      "eta_tauh/F");
    _treeDL->Branch("phi_tauh",       &_varListDL.phi_tauh,      "phi_tauh/F");

    _treeDL->Branch("px_met",         &_varListDL.px_met,        "px_met/F");
    _treeDL->Branch("py_met",         &_varListDL.py_met,        "py_met/F");
    _treeDL->Branch("pz_met",         &_varListDL.pz_met,        "pz_met/F");
    _treeDL->Branch("energy_met",     &_varListDL.energy_met,    "energy_met/F");
    _treeDL->Branch("met",            &_varListDL.met,           "met/F");
    _treeDL->Branch("phi_met",        &_varListDL.phi_met,       "phi_met/F");

    _treeDL->Branch("px_leadbj",      &_varListDL.px_leadbj,     "px_leadbj/F");
    _treeDL->Branch("py_leadbj",      &_varListDL.py_leadbj,     "py_leadbj/F");
    _treeDL->Branch("pz_leadbj",      &_varListDL.pz_leadbj,     "pz_leadbj/F");
    _treeDL->Branch("energy_leadbj",  &_varListDL.energy_leadbj, "energy_leadbj/F");
    _treeDL->Branch("pt_leadbj",      &_varListDL.pt_leadbj,     "pt_leadbj/F");
    _treeDL->Branch("eta_leadbj",     &_varListDL.eta_leadbj,    "eta_leadbj/F");
    _treeDL->Branch("phi_leadbj",     &_varListDL.phi_leadbj,    "phi_leadbj/F");

    _treeDL->Branch("px_leadlj",      &_varListDL.px_leadlj,     "px_leadlj/F");
    _treeDL->Branch("py_leadlj",      &_varListDL.py_leadlj,     "py_leadlj/F");
    _treeDL->Branch("pz_leadlj",      &_varListDL.pz_leadlj,     "pz_leadlj/F");
    _treeDL->Branch("energy_leadlj",  &_varListDL.energy_leadlj, "energy_leadlj/F");
    _treeDL->Branch("pt_leadlj",      &_varListDL.pt_leadlj,     "pt_leadlj/F");
    _treeDL->Branch("eta_leadlj",     &_varListDL.eta_leadlj,    "eta_leadlj/F");
    _treeDL->Branch("phi_leadlj",     &_varListDL.phi_leadlj,    "phi_leadlj/F");

    _treeDL->Branch("px_xlep",        &_varListDL.px_xlep,       "px_xlep/F");
    _treeDL->Branch("py_xlep",        &_varListDL.py_xlep,       "py_xlep/F");
    _treeDL->Branch("pz_xlep",        &_varListDL.pz_xlep,       "pz_xlep/F");
    _treeDL->Branch("energy_xlep",    &_varListDL.energy_xlep,   "energy_xlep/F");
    _treeDL->Branch("pt_xlep",        &_varListDL.pt_xlep,       "pt_xlep/F");
    _treeDL->Branch("eta_xlep",       &_varListDL.eta_xlep,      "eta_xlep/F");
    _treeDL->Branch("phi_xlep",       &_varListDL.phi_xlep,      "phi_xlep/F");

    _treeDL->Branch("px_wlep",        &_varListDL.px_wlep,       "px_wlep/F");
    _treeDL->Branch("py_wlep",        &_varListDL.py_wlep,       "py_wlep/F");
    _treeDL->Branch("pz_wlep",        &_varListDL.pz_wlep,       "pz_wlep/F");
    _treeDL->Branch("energy_wlep",    &_varListDL.energy_wlep,   "energy_wlep/F");
    _treeDL->Branch("pt_wlep",        &_varListDL.pt_wlep,       "pt_wlep/F");
    _treeDL->Branch("eta_wlep",       &_varListDL.eta_wlep,      "eta_wlep/F");
    _treeDL->Branch("phi_wlep",       &_varListDL.phi_wlep,      "phi_wlep/F");

    // --------- High level variables ---------- //                                                                                                                 
    // from one leg   
    _treeDL->Branch("met_along_tauh", &_varListDL.met_along_tauh,"met_along_tauh/F");
    _treeDL->Branch("frac_met_along_tauh", &_varListDL.frac_met_along_tauh, "frac_met_along_tauh/F");
    _treeDL->Branch("dr_xlep_tauh",    &_varListDL.dr_xlep_tauh,   "dr_xlep_tauh/F");
    _treeDL->Branch("dphi_xlep_tauh",  &_varListDL.dphi_xlep_tauh, "dphi_xlep_tauh/F");
    _treeDL->Branch("deta_xlep_tauh",  &_varListDL.deta_xlep_tauh, "deta_xlep_tauh/F");
  
    _treeDL->Branch("dr_xlep_leadlj",   &_varListDL.dr_xlep_leadlj,   "dr_xlep_leadlj/F");
    _treeDL->Branch("dphi_xlep_leadlj", &_varListDL.dphi_xlep_leadlj, "dphi_xlep_leadlj/F");
    _treeDL->Branch("deta_xlep_leadlj", &_varListDL.deta_xlep_leadlj, "deta_xlep_leadlj/F");

    _treeDL->Branch("dr_xlep_leadj",   &_varListDL.dr_xlep_leadj,   "dr_xlep_leadj/F");
    _treeDL->Branch("dphi_xlep_leadj", &_varListDL.dphi_xlep_leadj, "dphi_xlep_leadj/F");
    _treeDL->Branch("deta_xlep_leadj", &_varListDL.deta_xlep_leadj, "deta_xlep_leadj/F");

    _treeDL->Branch("dr_tauh_leadj",    &_varListDL.dr_tauh_leadj,    "dr_tauh_leadj/F");
    _treeDL->Branch("dphi_tauh_leadj",  &_varListDL.dphi_tauh_leadj,  "dphi_tauh_leadj/F");
    _treeDL->Branch("deta_tauh_leadj",  &_varListDL.deta_tauh_leadj,  "deta_tauh_leadj/F");

    _treeDL->Branch("dr_tauh_leadlj",    &_varListDL.dr_tauh_leadlj,   "dr_tauh_leadlj/F");
    _treeDL->Branch("dphi_tauh_leadlj",  &_varListDL.dphi_tauh_leadlj, "dphi_tauh_leadlj/F");
    _treeDL->Branch("deta_tauh_leadlj",  &_varListDL.deta_tauh_leadlj, "deta_tauh_leadlj/F");

    _treeDL->Branch("dphi_met_tauh",  &_varListDL.dphi_met_tauh, "dphi_met_tauh/F");
    _treeDL->Branch("dphi_met_leadlj",    &_varListDL.dphi_met_leadlj,   "dphi_met_leadlj/F");
    _treeDL->Branch("dphi_met_xlep",  &_varListDL.dphi_met_xlep, "dphi_met_xlep/F");

    _treeDL->Branch("dr_min_xlep_jets",   &_varListDL.dr_min_xlep_jets,   "dr_min_xlep_jets/F");
    _treeDL->Branch("dphi_min_xlep_jets", &_varListDL.dphi_min_xlep_jets, "dphi_min_xlep_jets/F");
    _treeDL->Branch("dr_max_xlep_jets",   &_varListDL.dr_max_xlep_jets,   "dr_max_xlep_jets/F");
    _treeDL->Branch("dphi_max_xlep_jets", &_varListDL.dphi_max_xlep_jets, "dphi_max_xlep_jets/F");

    _treeDL->Branch("mt_xlep_met",    &_varListDL.mt_xlep_met,     "mt_xlep_met/F");
    _treeDL->Branch("m_coll_xlep",    &_varListDL.m_coll_xlep,     "m_coll_xlep/F");
    _treeDL->Branch("effectiveMass_t1", &_varListDL.effectiveMass_t1, "effectiveMass_t1/F");
    _treeDL->Branch("invMass_x_test", &_varListDL.invMass_x_test,  "invMass_x_test/F");
    _treeDL->Branch("mt_x",           &_varListDL.mt_x,            "mt_x/F");

    // from other leg
    // dphi - met -others

    _treeDL->Branch("dphi_met_leadbj",  &_varListDL.dphi_met_leadbj, "dphi_met_leadbj/F");
    _treeDL->Branch("dphi_met_wlep",    &_varListDL.dphi_met_wlep,   "dphi_met_wlep/F");
  
    _treeDL->Branch("dr_wlep_leadbj",    &_varListDL.dr_wlep_leadbj,   "dr_wlep_leadbj/F");
    _treeDL->Branch("dphi_wlep_leadbj",  &_varListDL.dphi_wlep_leadbj, "dphi_wlep_leadbj/F");
    _treeDL->Branch("deta_wlep_leadbj", &_varListDL.deta_wlep_leadbj,"deta_wlep_leadbj/F");

    _treeDL->Branch("mt_wlep_met",    &_varListDL.mt_wlep_met,   "mt_wlep_met/F");
    _treeDL->Branch("m_coll_wlep",    &_varListDL.m_coll_wlep,   "m_coll_wlep/F");
    _treeDL->Branch("effectiveMass_t2", &_varListDL.effectiveMass_t2, "effectiveMass_t2/F");
  
    _treeDL->Branch("dr_tauh_leadbj",    &_varListDL.dr_tauh_leadbj,   "dr_tauh_leadbj/F");
    _treeDL->Branch("dphi_tauh_leadbj",  &_varListDL.dphi_tauh_leadbj, "dphi_tauh_leadbj/F");
    _treeDL->Branch("deta_tauh_leadbj",  &_varListDL.deta_tauh_leadbj, "deta_tauh_leadbj/F");

    _treeDL->Branch("dr_xlep_leadbj",    &_varListDL.dr_xlep_leadbj,   "dr_xlep_leadbj/F");
    _treeDL->Branch("dphi_xlep_leadbj",  &_varListDL.dphi_xlep_leadbj, "dphi_xlep_leadbj/F");
    _treeDL->Branch("deta_xlep_leadbj",  &_varListDL.deta_xlep_leadbj, "deta_xlep_leadbj/F");

    _treeDL->Branch("dr_xlep_wlep",    &_varListDL.dr_xlep_wlep,   "dr_xlep_wlep/F");
    _treeDL->Branch("dphi_xlep_wlep",  &_varListDL.dphi_xlep_wlep, "dphi_xlep_wlep/F");
    _treeDL->Branch("deta_xlep_Wlep",  &_varListDL.deta_xlep_wlep, "deta_xlep_wlep/F");
    _treeDL->Branch("scalarsumpt_xlep_wlep",&_varListDL.scalarsumpt_xlep_wlep,"scalarsumpt_xlep_wlep/F");
    _treeDL->Branch("vectorsumpt_xlep_wlep",&_varListDL.vectorsumpt_xlep_wlep,"vectorsumpt_xlep_wlep/F");

    _treeDL->Branch("dr_wlep_leadlj",    &_varListDL.dr_wlep_leadlj,   "dr_wlep_leadlj/F");
    _treeDL->Branch("dphi_wlep_leadlj",  &_varListDL.dphi_wlep_leadlj, "dphi_wlep_leadlj/F");
    _treeDL->Branch("deta_wlep_leadlj",  &_varListDL.deta_wlep_leadlj, "deta_wlep_leadlj/F");

    _treeDL->Branch("dr_wlep_tauh",      &_varListDL.dr_wlep_tauh,   "dr_wlep_tauh/F");
    _treeDL->Branch("dphi_wlep_tauh",    &_varListDL.dphi_wlep_tauh, "dphi_wlep_tauh/F");
    _treeDL->Branch("deta_wlep_tauh",    &_varListDL.deta_wlep_tauh, "deta_wlep_tauh/F");

    _treeDL->Branch("dr_min_wlep_jets",   &_varListDL.dr_min_wlep_jets,   "dr_min_wlep_jets/F");
    _treeDL->Branch("dphi_min_wlep_jets", &_varListDL.dphi_min_wlep_jets, "dphi_min_wlep_jets/F");
    _treeDL->Branch("dr_max_wlep_jets",   &_varListDL.dr_max_wlep_jets,   "dr_max_wlep_jets/F");
    _treeDL->Branch("dphi_max_wlep_jets", &_varListDL.dphi_max_wlep_jets, "dphi_max_wlep_jets/F");

    _treeDL->Branch("dr_min_jets",         &_varListDL.dr_min_jets,        "dr_min_jets/F");
    _treeDL->Branch("dr_max_jets",         &_varListDL.dr_max_jets,        "dr_max_jets/F");
    _treeDL->Branch("dr_leadbj_leadlj",    &_varListDL.dr_leadbj_leadlj,   "dr_leadbj_leadlj/F");
    _treeDL->Branch("dphi_leadbj_leadlj",  &_varListDL.dphi_leadbj_leadlj, "dphi_leadbj_leadlj/F");
    _treeDL->Branch("deta_leadbj_leadlj",  &_varListDL.deta_leadbj_leadlj, "deta_leadbj_leadlj/F");

    _treeDL->Branch("scalarsumpt_vis",  &_varListDL.scalarsumpt_vis,   "scalarsumpt_vis/F");
    _treeDL->Branch("scalarsumpt",      &_varListDL.scalarsumpt,       "scalarsumpt/F");
    _treeDL->Branch("effectivemass_vis",&_varListDL.effectivemass_vis, "effectivemass_vis/F");
    _treeDL->Branch("effectivemass",    &_varListDL.effectivemass,     "effectivemass/F");
    _treeDL->Branch("ht_jets",          &_varListDL.ht_jets,           "ht_jets/F");
    _treeDL->Branch("smin",             &_varListDL.smin,             "smin/F");
    _treeDL->Branch("costhetaS_xlep_tau",&_varListDL.costhetaS_xlep_tau,"costhetaS_xlep_tau/F");
    _treeDL->Branch("costhetaS_xlep_tauh",&_varListDL.costhetaS_xlep_tauh,"costhetaS_xlep_tauh/F");

  }  
  // ------------------------------------------------------------- //
  //                               SL Tree                         //
  // ------------------------------------------------------------- //
  else if (SL) {
    _treeSL = new TTree("EventsSL", "Flat_Events_TTree_SL");
    
    _treeSL->Branch("event",          &_varListSL.event,         "event/I");
    
    _treeSL->Branch("nleptons",       &_varListSL.nleptons,      "nleptons/I");
    _treeSL->Branch("njets",          &_varListSL.njets,         "njets/I");
    _treeSL->Branch("nbjets",         &_varListSL.nbjets,        "nbjets/I");
    _treeSL->Branch("nljets",         &_varListSL.nljets,        "nljets/I");
    _treeSL->Branch("ntauh",          &_varListSL.ntauh,         "ntauh/I");

    _treeSL->Branch("px_jet1",        &_varListSL.px_jet1,       "px_jet1/F");
    _treeSL->Branch("py_jet1",        &_varListSL.py_jet1,       "py_jet1/F");
    _treeSL->Branch("pz_jet1",        &_varListSL.pz_jet1,       "pz_jet1/F");
    _treeSL->Branch("energy_jet1",    &_varListSL.energy_jet1,   "energy_jet1/F");
    _treeSL->Branch("mass_jet1",      &_varListSL.mass_jet1,     "mass_jet1/F");
    _treeSL->Branch("pt_jet1",        &_varListSL.pt_jet1,       "pt_jet1/F");
    _treeSL->Branch("eta_jet1",       &_varListSL.eta_jet1,      "eta_jet1/F");

    _treeSL->Branch("px_jet2",        &_varListSL.px_jet2,       "px_jet2/F");
    _treeSL->Branch("py_jet2",        &_varListSL.py_jet2,       "py_jet2/F");
    _treeSL->Branch("pz_jet2",        &_varListSL.pz_jet2,       "pz_jet2/F");
    _treeSL->Branch("energy_jet2",    &_varListSL.energy_jet2,   "energy_jet2/F");
    _treeSL->Branch("mass_jet2",      &_varListSL.mass_jet2,     "mass_jet2/F");
    _treeSL->Branch("pt_jet2",        &_varListSL.pt_jet2,       "pt_jet2/F");
    _treeSL->Branch("eta_jet2",       &_varListSL.eta_jet2,      "eta_jet2/F");

    _treeSL->Branch("px_jet3",        &_varListSL.px_jet3,       "px_jet3/F");
    _treeSL->Branch("py_jet3",        &_varListSL.py_jet3,       "py_jet3/F");
    _treeSL->Branch("pz_jet3",        &_varListSL.pz_jet3,       "pz_jet3/F");
    _treeSL->Branch("energy_jet3",    &_varListSL.energy_jet3,   "energy_jet3/F");
    _treeSL->Branch("mass_jet3",      &_varListSL.mass_jet3,     "mass_jet3/F");
    _treeSL->Branch("pt_jet3",        &_varListSL.pt_jet3,       "pt_jet3/F");
    _treeSL->Branch("eta_jet3",       &_varListSL.eta_jet3,      "eta_jet3/F");

    _treeSL->Branch("px_jet4",        &_varListSL.px_jet4,       "px_jet4/F");
    _treeSL->Branch("py_jet4",        &_varListSL.py_jet4,       "py_jet4/F");
    _treeSL->Branch("pz_jet4",        &_varListSL.pz_jet4,       "pz_jet4/F");
    _treeSL->Branch("energy_jet4",    &_varListSL.energy_jet4,   "energy_jet4/F");
    _treeSL->Branch("mass_jet4",      &_varListSL.mass_jet4,     "mass_jet4/F");
    _treeSL->Branch("pt_jet4",        &_varListSL.pt_jet4,       "pt_jet4/F");
    _treeSL->Branch("eta_jet4",       &_varListSL.eta_jet4,      "eta_jet4/F");

    _treeSL->Branch("px_jet5",        &_varListSL.px_jet5,       "px_jet5/F");
    _treeSL->Branch("py_jet5",        &_varListSL.py_jet5,       "py_jet5/F");
    _treeSL->Branch("pz_jet5",        &_varListSL.pz_jet5,       "pz_jet5/F");
    _treeSL->Branch("energy_jet5",    &_varListSL.energy_jet5,   "energy_jet5/F");
    _treeSL->Branch("mass_jet5",      &_varListSL.mass_jet5,     "mass_jet5/F");
    _treeSL->Branch("pt_jet5",        &_varListSL.pt_jet5,       "pt_jet5/F");
    _treeSL->Branch("eta_jet5",       &_varListSL.eta_jet5,      "eta_jet5/F");

    _treeSL->Branch("px_jet6",        &_varListSL.px_jet6,       "px_jet6/F");
    _treeSL->Branch("py_jet6",        &_varListSL.py_jet6,       "py_jet6/F");
    _treeSL->Branch("pz_jet6",        &_varListSL.pz_jet6,       "pz_jet6/F");
    _treeSL->Branch("energy_jet6",    &_varListSL.energy_jet6,   "energy_jet6/F");
    _treeSL->Branch("mass_jet6",      &_varListSL.mass_jet6,     "mass_jet6/F");
    _treeSL->Branch("pt_jet6",        &_varListSL.pt_jet6,       "pt_jet6/F");
    _treeSL->Branch("eta_jet6",       &_varListSL.eta_jet6,      "eta_jet6/F");

    _treeSL->Branch("px_tauh",        &_varListSL.px_tauh,       "px_tauh/F");
    _treeSL->Branch("py_tauh",        &_varListSL.py_tauh,       "py_tauh/F");
    _treeSL->Branch("pz_tauh",        &_varListSL.pz_tauh,       "pz_tauh/F");
    _treeSL->Branch("energy_tauh",    &_varListSL.energy_tauh,   "energy_tauh/F");
    _treeSL->Branch("pt_tauh",        &_varListSL.pt_tauh,       "pt_tauh/F");
    _treeSL->Branch("eta_tauh",       &_varListSL.eta_tauh,      "eta_tauh/F");
    _treeSL->Branch("phi_tauh",       &_varListSL.phi_tauh,      "phi_tauh/F");

    _treeSL->Branch("px_met",        &_varListSL.px_met,       "px_met/F");
    _treeSL->Branch("py_met",        &_varListSL.py_met,       "py_met/F");
    _treeSL->Branch("pz_met",        &_varListSL.pz_met,       "pz_met/F");
    _treeSL->Branch("energy_met",    &_varListSL.energy_met,   "energy_met/F");
    _treeSL->Branch("pt_met",        &_varListSL.pt_met,       "pt_met/F");
    _treeSL->Branch("phi_met",       &_varListSL.phi_met,      "phi_met/F");

    _treeSL->Branch("px_leadbj",        &_varListSL.px_leadbj,       "px_leadbj/F");
    _treeSL->Branch("py_leadbj",        &_varListSL.py_leadbj,       "py_leadbj/F");
    _treeSL->Branch("pz_leadbj",        &_varListSL.pz_leadbj,       "pz_leadbj/F");
    _treeSL->Branch("energy_leadbj",    &_varListSL.energy_leadbj,   "energy_leadbj/F");
    _treeSL->Branch("pt_leadbj",        &_varListSL.pt_leadbj,       "pt_leadbj/F");
    _treeSL->Branch("eta_leadbj",       &_varListSL.eta_leadbj,      "eta_leadbj/F");
    _treeSL->Branch("phi_leadbj",       &_varListSL.phi_leadbj,      "phi_leadbj/F");

    _treeSL->Branch("px_toplj",        &_varListSL.px_toplj,       "px_toplj/F");
    _treeSL->Branch("py_toplj",        &_varListSL.py_toplj,       "py_toplj/F");
    _treeSL->Branch("pz_toplj",        &_varListSL.pz_toplj,       "pz_toplj/F");
    _treeSL->Branch("energy_toplj",    &_varListSL.energy_toplj,   "energy_toplj/F");
    _treeSL->Branch("pt_toplj",        &_varListSL.pt_toplj,       "pt_toplj/F");
    _treeSL->Branch("eta_toplj",       &_varListSL.eta_toplj,      "eta_toplj/F");
    _treeSL->Branch("phi_toplj",       &_varListSL.phi_toplj,      "phi_toplj/F");

    _treeSL->Branch("px_wjet1",        &_varListSL.px_wjet1,       "px_wjet1/F");
    _treeSL->Branch("py_wjet1",        &_varListSL.py_wjet1,       "py_wjet1/F");
    _treeSL->Branch("pz_wjet1",        &_varListSL.pz_wjet1,       "pz_wjet1/F");
    _treeSL->Branch("energy_wjet1",    &_varListSL.energy_wjet1,   "energy_wjet1/F");
    _treeSL->Branch("pt_wjet1",        &_varListSL.pt_wjet1,       "pt_wjet1/F");
    _treeSL->Branch("eta_wjet1",       &_varListSL.eta_wjet1,      "eta_wjet1/F");
    _treeSL->Branch("phi_wjet1",       &_varListSL.phi_wjet1,      "phi_wjet1/F");

    _treeSL->Branch("px_wjet2",        &_varListSL.px_wjet2,       "px_wjet2/F");
    _treeSL->Branch("py_wjet2",        &_varListSL.py_wjet2,       "py_wjet2/F");
    _treeSL->Branch("pz_wjet2",        &_varListSL.pz_wjet2,       "pz_wjet2/F");
    _treeSL->Branch("energy_wjet2",    &_varListSL.energy_wjet2,   "energy_wjet2/F");
    _treeSL->Branch("pt_wjet2",        &_varListSL.pt_wjet2,       "pt_wjet2/F");
    _treeSL->Branch("eta_wjet2",       &_varListSL.eta_wjet2,      "eta_wjet2/F");
    _treeSL->Branch("phi_wjet2",       &_varListSL.phi_wjet2,      "phi_wjet2/F");

    _treeSL->Branch("px_xlep",        &_varListSL.px_xlep,       "px_xlep/F");
    _treeSL->Branch("py_xlep",        &_varListSL.py_xlep,       "py_xlep/F");
    _treeSL->Branch("pz_xlep",        &_varListSL.pz_xlep,       "pz_xlep/F");
    _treeSL->Branch("energy_xlep",    &_varListSL.energy_xlep,   "energy_xlep/F");
    _treeSL->Branch("pt_xlep",        &_varListSL.pt_xlep,       "pt_xlep/F");
    _treeSL->Branch("eta_xlep",       &_varListSL.eta_xlep,      "eta_xlep/F");
    _treeSL->Branch("phi_xlep",       &_varListSL.phi_xlep,      "phi_xlep/F");

    _treeSL->Branch("px_tau",         &_varListSL.px_tau,        "px_tau/F");
    _treeSL->Branch("py_tau",         &_varListSL.py_tau,        "py_tau/F");
    _treeSL->Branch("pz_tau",         &_varListSL.pz_tau,        "pz_tau/F");
    _treeSL->Branch("energy_tau",     &_varListSL.energy_tau,    "energy_tau/F");
    _treeSL->Branch("pt_tau",         &_varListSL.pt_tau,        "pt_tau/F");
    _treeSL->Branch("eta_tau",        &_varListSL.eta_tau,       "eta_tau/F");
    _treeSL->Branch("phi_tau",        &_varListSL.phi_tau,       "phi_tau/F");

    _treeSL->Branch("px_nu",         &_varListSL.px_nu,        "px_nu/F");
    _treeSL->Branch("py_nu",         &_varListSL.py_nu,        "py_nu/F");
    _treeSL->Branch("pz_nu",         &_varListSL.pz_nu,        "pz_nu/F");
    _treeSL->Branch("energy_nu",     &_varListSL.energy_nu,    "energy_nu/F");
    _treeSL->Branch("pt_nu",         &_varListSL.pt_nu,        "pt_nu/F");
    _treeSL->Branch("eta_nu",        &_varListSL.eta_nu,       "eta_nu/F");
    _treeSL->Branch("phi_nu",        &_varListSL.phi_nu,       "phi_nu/F");

    _treeSL->Branch("px_w",         &_varListSL.px_w,        "px_w/F");
    _treeSL->Branch("py_w",         &_varListSL.py_w,        "py_w/F");
    _treeSL->Branch("pz_w",         &_varListSL.pz_w,        "pz_w/F");
    _treeSL->Branch("energy_w",     &_varListSL.energy_w,    "energy_w/F");
    _treeSL->Branch("pt_w",         &_varListSL.pt_w,        "pt_w/F");
    _treeSL->Branch("eta_w",        &_varListSL.eta_w,       "eta_w/F");
    _treeSL->Branch("phi_w",        &_varListSL.phi_w,       "phi_w/F");

    _treeSL->Branch("px_chi",         &_varListSL.px_chi,        "px_chi/F");
    _treeSL->Branch("py_chi",         &_varListSL.py_chi,        "py_chi/F");
    _treeSL->Branch("pz_chi",         &_varListSL.pz_chi,        "pz_chi/F");
    _treeSL->Branch("energy_chi",     &_varListSL.energy_chi,    "energy_chi/F");
    _treeSL->Branch("pt_chi",         &_varListSL.pt_chi,        "pt_chi/F");
    _treeSL->Branch("eta_chi",        &_varListSL.eta_chi,       "eta_chi/F");
    _treeSL->Branch("phi_chi",        &_varListSL.phi_chi,       "phi_chi/F");

    _treeSL->Branch("px_tSM",         &_varListSL.px_tSM,        "px_tSM/F");
    _treeSL->Branch("py_tSM",         &_varListSL.py_tSM,        "py_tSM/F");
    _treeSL->Branch("pz_tSM",         &_varListSL.pz_tSM,        "pz_tSM/F");
    _treeSL->Branch("energy_tSM",     &_varListSL.energy_tSM,    "energy_tSM/F");
    _treeSL->Branch("pt_tSM",         &_varListSL.pt_tSM,        "pt_tSM/F");
    _treeSL->Branch("eta_tSM",        &_varListSL.eta_tSM,       "eta_tSM/F");
    _treeSL->Branch("phi_tSM",        &_varListSL.phi_tSM,       "phi_tSM/F");

    _treeSL->Branch("px_tBSM",         &_varListSL.px_tBSM,        "px_tBSM/F");
    _treeSL->Branch("py_tBSM",         &_varListSL.py_tBSM,        "py_tBSM/F");
    _treeSL->Branch("pz_tBSM",         &_varListSL.pz_tBSM,        "pz_tBSM/F");
    _treeSL->Branch("energy_tBSM",     &_varListSL.energy_tBSM,    "energy_tBSM/F");
    _treeSL->Branch("pt_tBSM",         &_varListSL.pt_tBSM,        "pt_tBSM/F");
    _treeSL->Branch("eta_tBSM",        &_varListSL.eta_tBSM,       "eta_tBSM/F");
    _treeSL->Branch("phi_tBSM",        &_varListSL.phi_tBSM,       "phi_tBSM/F");

    _treeSL->Branch("invm_wj1_wj2",       &_varListSL.invm_wj1_wj2,      "invm_wj1_wj2/F");
    _treeSL->Branch("dr_wj1_wj2",         &_varListSL.dr_wj1_wj2,        "dr_wj1_wj2/F");
    _treeSL->Branch("dphi_wj1_wj2",       &_varListSL.dphi_wj1_wj2,      "dphi_wj1_wj2/F");
    _treeSL->Branch("deta_wj1_wj2",       &_varListSL.deta_wj1_wj2,      "deta_wj1_wj2/F");
    _treeSL->Branch("dr_w_leadbj",        &_varListSL.dr_w_leadbj,       "dr_w_leadbj/F");
    _treeSL->Branch("dphi_w_leadbj",      &_varListSL.dphi_w_leadbj,     "dphi_w_leadbj/F");
    _treeSL->Branch("dr_wj1_leadbj",      &_varListSL.dr_wj1_leadbj,     "dr_wj1_leadbj/F");
    _treeSL->Branch("dphi_wj1_leadbj",    &_varListSL.dphi_wj1_leadbj,   "dphi_wj1_leadbj/F");
    _treeSL->Branch("dr_wj2_leadbj",      &_varListSL.dr_wj2_leadbj,     "dr_wj2_leadbj/F");
    _treeSL->Branch("dphi_wj2_leadbj",    &_varListSL.dphi_wj2_leadbj,   "dphi_wj2_leadbj/F");
    _treeSL->Branch("invm_w_leadbj",      &_varListSL.invm_w_leadbj,     "invm_w_leadbj/F");
    _treeSL->Branch("dphi_met_leadbj",    &_varListSL.dphi_met_leadbj,   "dphi_met_leadbj/F");
    _treeSL->Branch("dphi_met_wj1",       &_varListSL.dphi_met_wj1,      "dphi_met_wj1/F");
    _treeSL->Branch("dphi_met_wj2",       &_varListSL.dphi_met_wj2,      "dphi_met_wj2/F");
    _treeSL->Branch("dphi_met_w",         &_varListSL.dphi_met_w,        "dphi_met_w/F");

    _treeSL->Branch("m_coll_x",           &_varListSL.m_coll_x,          "m_coll_x/F");
    _treeSL->Branch("invm_x_toplj",       &_varListSL.invm_x_toplj,      "invm_x_toplj/F");
    _treeSL->Branch("dr_xlep_toplj",      &_varListSL.dr_xlep_toplj,     "dr_xlep_toplj/F");
    _treeSL->Branch("dphi_xlep_toplj",    &_varListSL.dphi_xlep_toplj,   "dphi_xlep_toplj/F");
    _treeSL->Branch("dr_tauh_toplj",      &_varListSL.dr_tauh_toplj,     "dr_tauh_toplj/F");
    _treeSL->Branch("dphi_tauh_toplj",    &_varListSL.dphi_tauh_toplj,   "dphi_tauh_toplj/F");
    _treeSL->Branch("dr_x_toplj",         &_varListSL.dr_x_toplj,        "dr_x_toplj/F");
    _treeSL->Branch("dphi_x_toplj",       &_varListSL.dphi_x_toplj,      "dphi_x_toplj/F");
    _treeSL->Branch("dr_xlep_tauh",       &_varListSL.dr_xlep_tauh,      "dr_xlep_tauh/F");
    _treeSL->Branch("dphi_xlep_tauh",     &_varListSL.dphi_xlep_tauh,    "dphi_xlep_tauh/F");
    _treeSL->Branch("dphi_met_toplj",     &_varListSL.dphi_met_toplj,    "dphi_met_toplj/F");
    _treeSL->Branch("dphi_met_xlep",      &_varListSL.dphi_met_xlep ,    "dphi_met_xlep /F");
    _treeSL->Branch("dphi_met_tauh",      &_varListSL.dphi_met_tauh,     "dphi_met_tauh/F");

    _treeSL->Branch("dr_xlep_wj1",       &_varListSL.dr_xlep_wj1,      "dr_xlep_wj1/F");
    _treeSL->Branch("dphi_xlep_wj1",     &_varListSL.dphi_xlep_wj1,    "dphi_xlep_wj1/F");
    _treeSL->Branch("deta_xlep_wj1",     &_varListSL.deta_xlep_wj1,    "deta_xlep_wj1/F");
    _treeSL->Branch("dr_xlep_wj2",       &_varListSL.dr_xlep_wj2,      "dr_xlep_wj2/F");
    _treeSL->Branch("dphi_xlep_wj2",     &_varListSL.dphi_xlep_wj2,    "dphi_xlep_wj2/F");
    _treeSL->Branch("deta_xlep_wj2",     &_varListSL.deta_xlep_wj2,    "deta_xlep_wj2/F");
    _treeSL->Branch("dr_xlep_w",         &_varListSL.dr_xlep_w,        "dr_xlep_w/F");
    _treeSL->Branch("dphi_xlep_w",       &_varListSL.dphi_xlep_w,      "dphi_xlep_w/F");
    _treeSL->Branch("deta_xlep_w",       &_varListSL.deta_xlep_w,      "deta_xlep_w/F");
    _treeSL->Branch("dr_xlep_leadbj",    &_varListSL.dr_xlep_leadbj,   "dr_xlep_leadbj/F");
    _treeSL->Branch("dphi_xlep_leadbj",  &_varListSL.dphi_xlep_leadbj, "dphi_xlep_leadbj/F");
    _treeSL->Branch("deta_xlep_leadbj",  &_varListSL.deta_xlep_leadbj, "deta_xlep_leadbj/F");
    _treeSL->Branch("dr_tauh_wj1",       &_varListSL.dr_tauh_wj1,      "dr_tauh_wj1/F");
    _treeSL->Branch("dphi_tauh_wj1",     &_varListSL.dphi_tauh_wj1,    "dphi_tauh_wj1/F");
    _treeSL->Branch("deta_tauh_wj1",     &_varListSL.deta_tauh_wj1,    "deta_tauh_wj1/F");
    _treeSL->Branch("dr_tauh_wj2",       &_varListSL.dr_tauh_wj2,      "dr_tauh_wj2/F");
    _treeSL->Branch("dphi_tauh_wj2",     &_varListSL.dphi_tauh_wj2,    "dphi_tauh_wj2/F");
    _treeSL->Branch("deta_tauh_wj2",     &_varListSL.deta_tauh_wj2,    "deta_tauh_wj2/F");
    _treeSL->Branch("dr_tauh_leadbj",    &_varListSL.dr_tauh_leadbj,   "dr_tauh_leadbj/F");
    _treeSL->Branch("dphi_tauh_leadbj",  &_varListSL.dphi_tauh_leadbj, "dphi_tauh_leadbj/F");
    _treeSL->Branch("deta_tauh_leadbj",  &_varListSL.deta_tauh_leadbj, "deta_tauh_leadbj/F");
    _treeSL->Branch("dr_leadbj_toplj",  &_varListSL.dr_leadbj_toplj,  "dr_leadbj_toplj/F");
    _treeSL->Branch("dphi_leadbj_toplj",&_varListSL.dphi_leadbj_toplj,"dphi_leadbj_toplj/F");
    _treeSL->Branch("deta_leadbj_toplj",&_varListSL.deta_leadbj_toplj,"deta_leadbj_toplj/F");
    _treeSL->Branch("dr_w_toplj",       &_varListSL.dr_w_toplj,       "dr_w_toplj/F");
    _treeSL->Branch("dphi_w_toplj",     &_varListSL.dphi_w_toplj,    "dphi_w_toplj/F");
    _treeSL->Branch("deta_w_toplj",     &_varListSL.deta_w_toplj,    "deta_w_toplj/F");
    _treeSL->Branch("dr_x_leadbj",       &_varListSL.dr_x_leadbj,      "dr_x_leadbj/F");
    _treeSL->Branch("dphi_x_leadbj",     &_varListSL.dphi_x_leadbj,    "dphi_x_leadbj/F");
    _treeSL->Branch("deta_x_leadbj",     &_varListSL.deta_x_leadbj,    "deta_x_leadbj/F");
    _treeSL->Branch("dr_x_w",            &_varListSL.dr_x_w,           "dr_x_w/F");
    _treeSL->Branch("dphi_x_w",          &_varListSL.dphi_x_w,         "dphi_x_w/F");
    _treeSL->Branch("deta_x_w",          &_varListSL.deta_x_w,         "deta_x_w/F");
    _treeSL->Branch("dr_t1_t2",          &_varListSL.dr_t1_t2,         "dr_t1_t2/F");
    _treeSL->Branch("dphi_t1_t2",        &_varListSL.dphi_t1_t2,       "dphi_t1_t2/F");
    _treeSL->Branch("deta_t1_t2",        &_varListSL.deta_t1_t2,       "deta_t1_t2/F");
    _treeSL->Branch("total_mass",        &_varListSL.total_mass,       "total_mass/F");
    _treeSL->Branch("total_vector_pt",   &_varListSL.total_vector_pt,  "total_vector_pt/F");
    _treeSL->Branch("total_scalar_pt",   &_varListSL.total_scalar_pt,  "total_scalar_pt/F");
    _treeSL->Branch("mt_xlep_met",       &_varListSL.mt_xlep_met,      "mt_xlep_met/F");
    _treeSL->Branch("dr_min_jets",       &_varListSL.dr_min_jets,      "dr_min_jets/F");
    _treeSL->Branch("dr_max_jets",       &_varListSL.dr_max_jets,      "dr_max_jets/F");
    _treeSL->Branch("ht_jets",           &_varListSL.ht_jets,          "ht_jets/F");
    _treeSL->Branch("smin",              &_varListSL.smin,             "smin/F");
    _treeSL->Branch("costhetaS_xlep_tau",&_varListSL.costhetaS_xlep_tau,"costhetaS_xlep_tau/F");
    _treeSL->Branch("costhetaS_xlep_tauh",&_varListSL.costhetaS_xlep_tauh,"costhetaS_xlep_tauh/F");

    _treeSL->Branch("minDr_XlepJets",    &_varListSL.minDr_XlepJets,   "minDr_XlepJets/F");
    _treeSL->Branch("maxDr_XlepJets",    &_varListSL.maxDr_XlepJets,   "maxDr_XlepJets/F");
  }
  _mvaFile->ls();
}

//MVASkim::~MVASkim() {
//  if (_treeDL) delete _treeDL;  
//  if (_treeSL) delete _treeSL;  
//  if (_mvaFile) delete _mvaFile;
//}

void MVASkim::fill(const TreeVariablesDL& varList) {
  memcpy(&_varListDL, &varList, sizeof(varList));
  _mvaFile->cd();
  _treeDL->Fill();
}

void MVASkim::fill(const TreeVariablesSL& varList) {
  memcpy(&_varListSL, &varList, sizeof(varList));
  _mvaFile->cd();
  _treeSL->Fill();
}

void MVASkim::close() {
  _mvaFile->cd();
  _treeDL->Write();
  _treeSL->Write();
  _mvaFile->Save();
}

void MVASkim::close(bool isDL, bool isSL) {
  _mvaFile->cd();
  if (isDL) _treeDL->Write();
  else if (isSL) _treeSL->Write();
  _mvaFile->Save();
}
