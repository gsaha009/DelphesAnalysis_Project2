#include <iostream>
#include <memory>
#include "TFile.h"
#include "TTree.h"
#include "MVASkim.h"

using std::string;
using std::cout;
using std::endl;

MVASkim::MVASkim(const string& filename) {
  _mvaFile = TFile::Open(filename.c_str(), "RECREATE", "Skimmed Tree");
  _mvaFile->cd();
  // ------------------------------------------------------------- //
  //                               DL Tree                         //
  // ------------------------------------------------------------- //
  _treeDL = new TTree("EventsDL", "Flat_Events_TTree_DL");

  _treeDL->Branch("event",          &_varListDL.event,         "event/I");
  _treeDL->Branch("event_wt",       &_varListDL.event_wt,      "event_wt/F");

  _treeDL->Branch("nleptons",       &_varListDL.nleptons,      "nleptons/I");
  _treeDL->Branch("njets",          &_varListDL.njets,         "njets/I");
  _treeDL->Branch("nbjets",         &_varListDL.nbjets,        "nbjets/I");
  _treeDL->Branch("nljets",         &_varListDL.nljets,        "nljets/I");
  _treeDL->Branch("ntauh",          &_varListDL.ntauh,         "ntauh/I");

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

  //_treeDL->Branch("dr_met_tauh",    &_varListDL.dr_met_tauh,   "dr_met_tauh/F");
  _treeDL->Branch("dphi_met_tauh",  &_varListDL.dphi_met_tauh, "dphi_met_tauh/F");
  //_treeDL->Branch("deta_met_tauh",  &_varListDL.deta_met_tauh, "deta_met_tauh/F");

  _treeDL->Branch("dr_min_xlep_jets",   &_varListDL.dr_min_xlep_jets,   "dr_min_xlep_jets/F");
  _treeDL->Branch("dphi_min_xlep_jets", &_varListDL.dphi_min_xlep_jets, "dphi_min_xlep_jets/F");
  _treeDL->Branch("dr_max_xlep_jets",   &_varListDL.dr_max_xlep_jets,   "dr_max_xlep_jets/F");
  _treeDL->Branch("dphi_max_xlep_jets", &_varListDL.dphi_max_xlep_jets, "dphi_max_xlep_jets/F");

  _treeDL->Branch("mt_xlep_met",    &_varListDL.mt_xlep_met,     "mt_xlep_met/F");
  _treeDL->Branch("m_coll_xlep",    &_varListDL.m_coll_xlep,     "m_coll_xlep/F");
  _treeDL->Branch("effectiveMass_t1", &_varListDL.effectiveMass_t1, "effectiveMass_t1/F");
  _treeDL->Branch("invMass_x_test", &_varListDL.invMass_x_test,  "invMass_x_test/F");
  _treeDL->Branch("mt_x",           &_varListDL.mt_x,            "mt_x/F");

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
  
  // ------------------------------------------------------------- //
  //                               SL Tree                         //
  // ------------------------------------------------------------- //
  _treeSL = new TTree("EventsSL", "Flat_Events_TTree_SL");

  _treeSL->Branch("event",          &_varListSL.event,         "event/I");

  _treeSL->Branch("nleptons",       &_varListSL.nleptons,      "nleptons/I");
  _treeSL->Branch("njets",          &_varListSL.njets,         "njets/I");
  _treeSL->Branch("nbjets",         &_varListSL.nbjets,        "nbjets/I");
  _treeSL->Branch("nljets",         &_varListSL.nljets,        "nljets/I");
  _treeSL->Branch("ntauh",          &_varListSL.ntauh,         "ntauh/I");

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

  _treeSL->Branch("px_leadlj",        &_varListSL.px_leadlj,       "px_leadlj/F");
  _treeSL->Branch("py_leadlj",        &_varListSL.py_leadlj,       "py_leadlj/F");
  _treeSL->Branch("pz_leadlj",        &_varListSL.pz_leadlj,       "pz_leadlj/F");
  _treeSL->Branch("energy_leadlj",    &_varListSL.energy_leadlj,   "energy_leadlj/F");
  _treeSL->Branch("pt_leadlj",        &_varListSL.pt_leadlj,       "pt_leadlj/F");
  _treeSL->Branch("eta_leadlj",       &_varListSL.eta_leadlj,      "eta_leadlj/F");
  _treeSL->Branch("phi_leadlj",       &_varListSL.phi_leadlj,      "phi_leadlj/F");

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

  _treeSL->Branch("m_coll_x",           &_varListSL.m_coll_x,          "m_coll_x/F");
  _treeSL->Branch("m_coll_xtest",       &_varListSL.m_coll_xtest,      "m_coll_xtest/F");
  _treeSL->Branch("invm_x_leadlj",      &_varListSL.invm_x_leadlj,     "invm_x_leadlj/F");
  _treeSL->Branch("dr_xlep_leadlj",     &_varListSL.dr_xlep_leadlj,    "dr_xlep_leadlj/F");
  _treeSL->Branch("dphi_xlep_leadlj",   &_varListSL.dphi_xlep_leadlj,  "dphi_xlep_leadlj/F");
  _treeSL->Branch("dr_tauh_leadlj",     &_varListSL.dr_tauh_leadlj,    "dr_tauh_leadlj/F");
  _treeSL->Branch("dphi_tauh_leadlj",   &_varListSL.dphi_tauh_leadlj,  "dphi_tauh_leadlj/F");
  _treeSL->Branch("dr_x_leadlj",        &_varListSL.dr_x_leadlj,       "dr_x_leadlj/F");
  _treeSL->Branch("dphi_x_leadlj",      &_varListSL.dphi_x_leadlj,     "dphi_x_leadlj/F");
  _treeSL->Branch("dr_xlep_tauh",       &_varListSL.dr_xlep_tauh,      "dr_xlep_tauh/F");
  _treeSL->Branch("dphi_xlep_tauh",     &_varListSL.dphi_xlep_tauh,    "dphi_xlep_tauh/F");

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
  _treeSL->Branch("dr_leadbj_leadlj",  &_varListSL.dr_leadbj_leadlj, "dr_leadbj_leadlj/F");
  _treeSL->Branch("dphi_leadbj_leadlj",&_varListSL.dphi_leadbj_leadlj,"dphi_leadbj_leadlj/F");
  _treeSL->Branch("deta_leadbj_leadlj",&_varListSL.deta_leadbj_leadlj,"deta_leadbj_leadlj/F");
  _treeSL->Branch("dr_w_leadlj",       &_varListSL.dr_w_leadlj,      "dr_w_leadlj/F");
  _treeSL->Branch("dphi_w_leadlj",     &_varListSL.dphi_w_leadlj,    "dphi_w_leadlj/F");
  _treeSL->Branch("deta_w_leadlj",     &_varListSL.deta_w_leadlj,    "deta_w_leadlj/F");
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
  _treeSL->Branch("dphi_xlep_met",     &_varListSL.dphi_xlep_met,    "dphi_xlep_met/F");
  _treeSL->Branch("mt_xlep_met",       &_varListSL.mt_xlep_met,      "mt_xlep_met/F");
  _treeSL->Branch("dr_min_jets",       &_varListSL.dr_min_jets,      "dr_min_jets/F");
  _treeSL->Branch("dr_max_jets",       &_varListSL.dr_max_jets,      "dr_max_jets/F");
  _treeSL->Branch("ht_jets",           &_varListSL.ht_jets,          "ht_jets/F");

  _mvaFile->ls();
}

MVASkim::~MVASkim() {
  if (_treeDL) delete _treeDL;  
  if (_treeSL) delete _treeSL;  
  if (_mvaFile) delete _mvaFile;
}

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
