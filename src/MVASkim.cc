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
  _tree = new TTree("Events", "Flat_Events_TTree");

  _tree->Branch("event",          &_varList.event,         "event/I");

  _tree->Branch("nLeptons",       &_varList.nLeptons,      "nLeptons/I");
  _tree->Branch("nJets",          &_varList.nJets,         "nJets/I");
  _tree->Branch("nbJets",         &_varList.nbJets,        "nbJets/I");
  _tree->Branch("nlJets",         &_varList.nlJets,        "nlJets/I");
  _tree->Branch("nTauh",          &_varList.nTauh,         "nTauh/I");

  _tree->Branch("px_tauh1",       &_varList.px_tauh1,      "px_tauh1/F");
  _tree->Branch("py_tauh1",       &_varList.py_tauh1,      "py_tauh1/F");
  _tree->Branch("pz_tauh1",       &_varList.pz_tauh1,      "pz_tauh1/F");
  _tree->Branch("energy_tauh1",   &_varList.energy_tauh1,  "energy_tauh1/F");
  _tree->Branch("pt_tauh1",       &_varList.pt_tauh1,      "pt_tauh1/F");
  _tree->Branch("eta_tauh1",      &_varList.eta_tauh1,     "eta_tauh1/F");
  _tree->Branch("phi_tauh1",      &_varList.phi_tauh1,     "phi_tauh1/F");

  _tree->Branch("px_met",         &_varList.px_met,        "px_met/F");
  _tree->Branch("py_met",         &_varList.py_met,        "py_met/F");
  _tree->Branch("pz_met",         &_varList.pz_met,        "pz_met/F");
  _tree->Branch("energy_met",     &_varList.energy_met,    "energy_met/F");
  _tree->Branch("met",            &_varList.met,           "met/F");
  _tree->Branch("phi_met",        &_varList.phi_met,       "phi_met/F");

  _tree->Branch("px_bjet1",       &_varList.px_bjet1,      "px_bjet1/F");
  _tree->Branch("py_bjet1",       &_varList.py_bjet1,      "py_bjet1/F");
  _tree->Branch("pz_bjet1",       &_varList.pz_bjet1,      "pz_bjet1/F");
  _tree->Branch("energy_bjet1",   &_varList.energy_bjet1,  "energy_bjet1/F");
  _tree->Branch("pt_bjet1",       &_varList.pt_bjet1,      "pt_bjet1/F");
  _tree->Branch("eta_bjet1",      &_varList.eta_bjet1,     "eta_bjet1/F");
  _tree->Branch("phi_bjet1",      &_varList.phi_bjet1,     "phi_bjet1/F");

  _tree->Branch("px_ljet1",       &_varList.px_ljet1,      "px_ljet1/F");
  _tree->Branch("py_ljet1",       &_varList.py_ljet1,      "py_ljet1/F");
  _tree->Branch("pz_ljet1",       &_varList.pz_ljet1,      "pz_ljet1/F");
  _tree->Branch("energy_ljet1",   &_varList.energy_ljet1,  "energy_ljet1/F");
  _tree->Branch("pt_ljet1",       &_varList.pt_ljet1,      "pt_ljet1/F");
  _tree->Branch("eta_ljet1",      &_varList.eta_ljet1,     "eta_ljet1/F");
  _tree->Branch("phi_ljet1",      &_varList.phi_ljet1,     "phi_ljet1/F");

  _tree->Branch("px_Xlep",        &_varList.px_Xlep,       "px_Xlep/F");
  _tree->Branch("py_Xlep",        &_varList.py_Xlep,       "py_Xlep/F");
  _tree->Branch("pz_Xlep",        &_varList.pz_Xlep,       "pz_Xlep/F");
  _tree->Branch("energy_Xlep",    &_varList.energy_Xlep,   "energy_Xlep/F");
  _tree->Branch("pt_Xlep",        &_varList.pt_Xlep,       "pt_Xlep/F");
  _tree->Branch("eta_Xlep",       &_varList.eta_Xlep,      "eta_Xlep/F");
  _tree->Branch("phi_Xlep",       &_varList.phi_Xlep,      "phi_Xlep/F");

  _tree->Branch("px_Wlep",        &_varList.px_Wlep,       "px_Wlep/F");
  _tree->Branch("py_Wlep",        &_varList.py_Wlep,       "py_Wlep/F");
  _tree->Branch("pz_Wlep",        &_varList.pz_Wlep,       "pz_Wlep/F");
  _tree->Branch("energy_Wlep",    &_varList.energy_Wlep,   "energy_Wlep/F");
  _tree->Branch("pt_Wlep",        &_varList.pt_Wlep,       "pt_Wlep/F");
  _tree->Branch("eta_Wlep",       &_varList.eta_Wlep,      "eta_Wlep/F");
  _tree->Branch("phi_Wlep",       &_varList.phi_Wlep,      "phi_Wlep/F");

  _tree->Branch("scalarSumpt_XlwpWlep",    &_varList.scalarSumpt_XlwpWlep,   "scalarSumpt_XlwpWlep/F");
  _tree->Branch("vectorSumpt_XlwpWlep",    &_varList.vectorSumpt_XlwpWlep,   "vectorSumpt_XlwpWlep/F");
  _tree->Branch("dr_XlepWlep",    &_varList.dr_XlepWlep,   "dr_XlepWlep/F");
  _tree->Branch("dphi_XlepWlep",  &_varList.dphi_XlepWlep, "dphi_XlepWlep/F");
  _tree->Branch("deta_XlepWlep",  &_varList.deta_XlepWlep, "deta_XlepWlep/F");

  _tree->Branch("dr_Xleptauh",    &_varList.dr_Xleptauh,   "dr_Xleptauh/F");
  _tree->Branch("dphi_Xleptauh",  &_varList.dphi_Xleptauh, "dphi_Xleptauh/F");
  _tree->Branch("deta_Xleptauh",  &_varList.deta_Xleptauh, "deta_Xleptauh/F");
  _tree->Branch("dr_Wleptauh",    &_varList.dr_Wleptauh,   "dr_Wleptauh/F");
  _tree->Branch("dphi_Wleptauh",  &_varList.dphi_Wleptauh, "dphi_Wleptauh/F");
  _tree->Branch("deta_Wleptauh",  &_varList.deta_Wleptauh, "deta_Wleptauh/F");

  _tree->Branch("dr_tauhjet",     &_varList.dr_tauhjet,    "dr_tauhjet/F");
  _tree->Branch("dphi_tauhjet",   &_varList.dphi_tauhjet,  "dphi_tauhjet/F");
  _tree->Branch("deta_tauhjet",   &_varList.deta_tauhjet,  "deta_tauhjet/F");
  _tree->Branch("dr_tauhljet",    &_varList.dr_tauhljet,   "dr_tauhljet/F");
  _tree->Branch("dphi_tauhljet",  &_varList.dphi_tauhljet, "dphi_tauhljet/F");
  _tree->Branch("deta_tauhljet",  &_varList.deta_tauhljet, "deta_tauhljet/F");
  _tree->Branch("dr_tauhbjet",    &_varList.dr_tauhbjet,   "dr_tauhbjet/F");
  _tree->Branch("dphi_tauhbjet",  &_varList.dphi_tauhbjet, "dphi_tauhbjet/F");
  _tree->Branch("deta_tauhbjet",  &_varList.deta_tauhbjet, "deta_tauhbjet/F");

  _tree->Branch("dr_min_Xlepjets",   &_varList.dr_min_Xlepjets,   "dr_min_Xlepjets/F");
  _tree->Branch("dphi_min_Xlepjets", &_varList.dphi_min_Xlepjets, "dphi_min_Xlepjets/F");
  _tree->Branch("dr_max_Xlepjets",   &_varList.dr_max_Xlepjets,   "dr_max_Xlepjets/F");
  _tree->Branch("dphi_max_Xlepjets", &_varList.dphi_max_Xlepjets, "dphi_max_Xlepjets/F");
  _tree->Branch("dr_min_Wlepjets",   &_varList.dr_min_Wlepjets,   "dr_min_Wlepjets/F");
  _tree->Branch("dphi_min_Wlepjets", &_varList.dphi_min_Wlepjets, "dphi_min_Wlepjets/F");
  _tree->Branch("dr_max_Wlepjets",   &_varList.dr_max_Wlepjets,   "dr_max_Wlepjets/F");
  _tree->Branch("dphi_max_Wlepjets", &_varList.dphi_max_Wlepjets, "dphi_max_Wlepjets/F");

  _tree->Branch("dr_Xlepbjet",    &_varList.dr_Xlepbjet,   "dr_Xlepbjet/F");
  _tree->Branch("dphi_Xlepbjet",  &_varList.dphi_Xlepbjet, "dphi_Xlepbjet/F");
  _tree->Branch("deta_Xlepbjet",  &_varList.deta_Xlepbjet, "deta_Xlepbjet/F");
  _tree->Branch("dr_Wlepbjet",    &_varList.dr_Wlepbjet,   "dr_Wlepbjet/F");
  _tree->Branch("dphi_Wlepbjet",  &_varList.dphi_Wlepbjet, "dphi_Wlepbjet/F");
  _tree->Branch("deta_Wlepbjet",  &_varList.deta_Wlepbjet, "deta_Wlepbjet/F");

  _tree->Branch("dr_Xlepljet",    &_varList.dr_Xlepljet,   "dr_Xlepljet/F");
  _tree->Branch("dphi_Xlepljet",  &_varList.dphi_Xlepljet, "dphi_Xlepljet/F");
  _tree->Branch("deta_Xlepljet",  &_varList.deta_Xlepljet, "deta_Xlepljet/F");
  _tree->Branch("dr_Wlepljet",    &_varList.dr_Wlepljet,   "dr_Wlepljet/F");
  _tree->Branch("dphi_Wlepljet",  &_varList.dphi_Wlepljet, "dphi_Wlepljet/F");
  _tree->Branch("deta_Wlepljet",  &_varList.deta_Wlepljet, "deta_Wlepljet/F");

  _tree->Branch("dphi_Xlepmet",   &_varList.dphi_Xlepmet,  "dphi_Xlepmet/F");
  _tree->Branch("dphi_Wlepmet",   &_varList.dphi_Wlepmet,  "dphi_Wlepmet/F");
  _tree->Branch("mt_Xlepmet",     &_varList.mt_Xlepmet,    "mt_Xlepmet/F");
  _tree->Branch("mt_Wlepmet",     &_varList.mt_Wlepmet,    "mt_Wlepmet/F");

  _tree->Branch("dr_min_jets",    &_varList.dr_min_jets,   "dr_min_jets/F");
  _tree->Branch("dr_max_jets",    &_varList.dr_max_jets,   "dr_max_jets/F");
  _tree->Branch("dr_bjetljet",    &_varList.dr_bjetljet,   "dr_bjetljet/F");
  _tree->Branch("dphi_bjetljet",  &_varList.dphi_bjetljet, "dphi_bjetljet/F");
  _tree->Branch("deta_bjetljet",  &_varList.deta_bjetljet, "deta_bjetljet/F");

  _tree->Branch("scalarSumPtVis",  &_varList.scalarSumPtVis,   "scalarSumPtVis/F");
  _tree->Branch("scalarSumPt",     &_varList.scalarSumPt,      "scalarSumPt/F");
  _tree->Branch("effectiveMassVis",&_varList.effectiveMassVis, "effectiveMassVis/F");
  _tree->Branch("effectiveMass",   &_varList.effectiveMass,    "effectiveMassVis/F");
  _tree->Branch("M_coll_IC_Xlep",  &_varList.M_coll_IC_Xlep,   "M_coll_IC_Xlep/F");
  _tree->Branch("M_coll_IC_Wlep",  &_varList.M_coll_IC_Wlep,   "M_coll_IC_Wlep/F");
  _tree->Branch("HT_Jets",         &_varList.HT_Jets,          "HT_Jets/F");
  
  _mvaFile->ls();
}
MVASkim::~MVASkim() {
  if (_tree) delete _tree;  
  if (_mvaFile) delete _mvaFile;
}
void MVASkim::fill(const TreeVariables& varList) {
  memcpy(&_varList, &varList, sizeof(varList));
  _mvaFile->cd();
  _tree->Fill();
}
void MVASkim::close() {
  //_mvaFile = TFile::Open(filename.c_str(), "RECREATE", "Skimmed Tree");
  _mvaFile->cd();
  //_tree->Print();
  _tree->Write();
  _mvaFile->Save();
  //_mvaFile->Write();
  //_mvaFile->Close();
}
