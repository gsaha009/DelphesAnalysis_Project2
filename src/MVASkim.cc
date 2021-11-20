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
  _tree->Branch("nTauh",          &_varList.nTauh,         "nTauh/I");
  _tree->Branch("pt_tauh1",       &_varList.pt_tauh1,      "pt_tauh1/F");
  _tree->Branch("eta_tauh1",      &_varList.eta_tauh1,     "eta_tauh1/F");
  _tree->Branch("phi_tauh1",      &_varList.phi_tauh1,     "phi_tauh1/F");
  _tree->Branch("met",            &_varList.met,           "met/F");
  _tree->Branch("pt_bjet1",       &_varList.pt_bjet1,      "pt_bjet1/F");
  _tree->Branch("eta_bjet1",      &_varList.eta_bjet1,     "eta_bjet1/F");
  _tree->Branch("phi_bjet1",      &_varList.phi_bjet1,     "phi_bjet1/F");
  _tree->Branch("pt_lep1",        &_varList.pt_lep1,       "pt_lep1/F");
  _tree->Branch("pt_lep2",        &_varList.pt_lep2,       "pt_lep2/F");
  _tree->Branch("eta_lep1",       &_varList.eta_lep1,      "eta_lep1/F");
  _tree->Branch("eta_lep2",       &_varList.eta_lep2,      "eta_lep2/F");
  _tree->Branch("phi_lep1",       &_varList.phi_lep1,      "phi_lep1/F");
  _tree->Branch("phi_lep2",       &_varList.phi_lep2,      "phi_lep2/F");
  _tree->Branch("dr_lep1lep2",    &_varList.dr_lep1lep2,   "dr_lep1lep2/F");
  _tree->Branch("dphi_lep1lep2",  &_varList.dphi_lep1lep2, "dphi_lep1lep2/F");
  _tree->Branch("deta_lep1lep2",  &_varList.deta_lep1lep2, "deta_lep1lep2/F");
  _tree->Branch("dr_lep1tauh",    &_varList.dr_lep1tauh,   "dr_lep1tauh/F");
  _tree->Branch("dphi_lep1tauh",  &_varList.dphi_lep1tauh, "dphi_lep1tauh/F");
  _tree->Branch("deta_lep1tauh",  &_varList.deta_lep1tauh, "deta_lep1tauh/F");
  _tree->Branch("dr_lep2tauh",    &_varList.dr_lep2tauh,   "dr_lep2tauh/F");
  _tree->Branch("dphi_lep2tauh",  &_varList.dphi_lep2tauh, "dphi_lep2tauh/F");
  _tree->Branch("deta_lep2tauh",  &_varList.deta_lep2tauh, "deta_lep2tauh/F");
  _tree->Branch("dr_tauhjet",     &_varList.dr_tauhjet,    "dr_tauhjet/F");
  _tree->Branch("dphi_tauhjet",   &_varList.dphi_tauhjet,  "dphi_tauhjet/F");
  _tree->Branch("deta_tauhjet",   &_varList.deta_tauhjet,  "deta_tauhjet/F");
  _tree->Branch("dr_tauhljet",    &_varList.dr_tauhljet,   "dr_tauhljet/F");
  _tree->Branch("dphi_tauhljet",  &_varList.dphi_tauhljet, "dphi_tauhljet/F");
  _tree->Branch("deta_tauhljet",  &_varList.deta_tauhljet, "deta_tauhljet/F");
  _tree->Branch("dr_tauhbjet",    &_varList.dr_tauhbjet,   "dr_tauhbjet/F");
  _tree->Branch("dphi_tauhbjet",  &_varList.dphi_tauhbjet, "dphi_tauhbjet/F");
  _tree->Branch("deta_tauhbjet",  &_varList.deta_tauhbjet, "deta_tauhbjet/F");
  _tree->Branch("dr_lep1bjet",    &_varList.dr_lep1bjet,   "dr_lep1bjet/F");
  _tree->Branch("dphi_lep1bjet",  &_varList.dphi_lep1bjet, "dphi_lep1bjet/F");
  _tree->Branch("deta_lep1bjet",  &_varList.deta_lep1bjet, "deta_lep1bjet/F");
  _tree->Branch("dr_lep2bjet",    &_varList.dr_lep2bjet,   "dr_lep2bjet/F");
  _tree->Branch("dphi_lep2bjet",  &_varList.dphi_lep2bjet, "dphi_lep2bjet/F");
  _tree->Branch("deta_lep2bjet",  &_varList.deta_lep2bjet, "deta_lep2bjet/F");
  _tree->Branch("dr_lep1ljet",    &_varList.dr_lep1ljet,   "dr_lep1ljet/F");
  _tree->Branch("dphi_lep1ljet",  &_varList.dphi_lep1ljet, "dphi_lep1ljet/F");
  _tree->Branch("deta_lep1ljet",  &_varList.deta_lep1ljet, "deta_lep1ljet/F");
  _tree->Branch("dr_lep2ljet",    &_varList.dr_lep2ljet,   "dr_lep2ljet/F");
  _tree->Branch("dphi_lep2ljet",  &_varList.dphi_lep2ljet, "dphi_lep2ljet/F");
  _tree->Branch("deta_lep2ljet",  &_varList.deta_lep2ljet, "deta_lep2ljet/F");
  _tree->Branch("dphi_lep1met",   &_varList.dphi_lep1met,  "dphi_lep1met/F");
  _tree->Branch("dphi_lep2met",   &_varList.dphi_lep2met,  "dphi_lep2met/F");
  _tree->Branch("mt_lep1met",     &_varList.mt_lep1met,    "mt_lep1met/F");
  _tree->Branch("mt_lep2met",     &_varList.mt_lep2met,    "mt_lep2met/F");
  _tree->Branch("mt_muon1met",    &_varList.mt_muon1met,   "mt_muon1met/F");
  _tree->Branch("dr_min_jets",    &_varList.dr_min_jets,   "dr_min_jets/F");
  _tree->Branch("dr_max_jets",    &_varList.dr_max_jets,   "dr_max_jets/F");
  _tree->Branch("dr_bjetljet",    &_varList.dr_bjetljet,   "dr_bjetljet/F");
  _tree->Branch("dphi_bjetljet",  &_varList.dphi_bjetljet, "dphi_bjetljet/F");
  _tree->Branch("deta_bjetljet",  &_varList.deta_bjetljet, "deta_bjetljet/F");
  
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
