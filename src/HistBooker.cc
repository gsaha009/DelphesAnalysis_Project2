#include <memory>
#include <iostream>

#include "HistBooker.h"

void HistBooker::bookHist1D(const char* hname, const char* htitle, 
			    int nbins, float xlow, float xhigh) {
  new TH1D(hname, htitle, nbins, xlow, xhigh);
}

void HistBooker::bookHist2D(const char* hname, const char* htitle, 
			    int nbinsx, float xlow, float xhigh,
			    int nbinsy, float ylow, float yhigh) {
  new TH2D(hname, htitle, nbinsx, xlow, xhigh, nbinsy, ylow, yhigh);
}

// ------------------------------------------------------------------------------ //
//                                Common histograms                               //
// ------------------------------------------------------------------------------ //
void HistBooker::bookBasicHistograms() {
  // ------------------------------------------------------------------------------ //
  //                                Common histograms                               //
  // ------------------------------------------------------------------------------ //

  // Raw Muon 
  bookHist1D("Raw_mu_pt",                   "Raw #mu p_{T} (GeV)",      100, 0, 300);
  bookHist1D("Raw_mu_eta",                  "Raw #mu #eta",               50, -5, 5);
  bookHist1D("Raw_mu_phi",                  "Raw #mu #phi",           32, -3.2, 3.2);
  bookHist1D("Raw_mu_charge",               "Raw #mu charge",          3, -1.5, 1.5);
  bookHist1D("Raw_mu_isolationVar",         "Raw #mu Iso",                100, 0, 5);
  bookHist1D("Raw_mu_isolationVarRhoCorr",  "Raw #mu Iso #rho Corr",      100, 0, 5);
  bookHist1D("Raw_mu_D0",                   "Raw #mu D0",                 40, -2, 2);
  bookHist1D("Raw_mu_DZ",                   "Raw #mu DZ",                 40, -2, 2);
  // Raw Electron
  bookHist1D("Raw_ele_pt",                  "Raw #ele p_{T} (GeV)",     100, 0, 300);
  bookHist1D("Raw_ele_eta",                 "Raw #ele #eta",              50, -5, 5);
  bookHist1D("Raw_ele_phi",                 "Raw #ele #phi",          32, -3.2, 3.2);
  bookHist1D("Raw_ele_charge",              "Raw #ele charge",         3, -1.5, 1.5);
  bookHist1D("Raw_ele_isolationVar",        "Raw #ele Iso",               100, 0, 5);
  bookHist1D("Raw_ele_isolationVarRhoCorr", "Raw #ele Iso #rho Corr",     100, 0, 5);
  bookHist1D("Raw_ele_D0",                  "Raw #ele D0",                40, -2, 2);
  bookHist1D("Raw_ele_DZ",                  "Raw #ele DZ",                40, -2, 2);
  // Raw Jet
  bookHist1D("Raw_jet_pt",                  "Raw #jet p_{T} (GeV)",     100, 0, 300);
  bookHist1D("Raw_jet_eta",                 "Raw #jet #eta",              50, -5, 5);
  bookHist1D("Raw_jet_phi",                 "Raw #jet #phi",          32, -3.2, 3.2);
  bookHist1D("Raw_jet_mass",                "Raw #jet mass",             100, 0, 50);
  bookHist1D("Raw_jet_DeltaEta",            "Raw #jet radius in pseudorapidity",                     50, -5, 5);
  bookHist1D("Raw_jet_DeltaPhi",            "Raw #jet radius in azimuthal angle",                32, -3.2, 3.2);
  bookHist1D("Raw_jet_flavor",              "Raw #jet flavor",        10, -0.5, 9.5);
  bookHist1D("Raw_jet_tautag",              "Raw #jet tautag",         3, -1.5, 1.5);
  bookHist1D("Raw_jet_tauwt",               "Raw #jet probability for jet to be identified as #tau",  50, 0, 1);
  bookHist1D("Raw_jet_EhadOverEem",         "Raw #jet ratio of the hadronic versus electromagnetic energy deposited in the calorimeter", 100, 0, 5);
  bookHist1D("Raw_jet_BTag",                "Raw #jet btag",           2, -0.5, 1.5);
  bookHist1D("Raw_jet_BetaStar",            "Raw #jet BetaStar (sum pt of charged constituents coming from hard interaction)/(sum pt of charged constituents)", 100, 0, 5);
  // object and event cut flow
  bookHist1D("MuonCutFlow", "Muon CutFlow", 3, -0.5, 2.5);
  bookHist1D("ElectronCutFlow", "Electron CutFlow", 3, -0.5, 2.5);
  bookHist1D("TauhCutFlow", "Tauh CutFlow", 4, -0.5, 3.5);
  bookHist1D("JetCutFlow", "Jet CutFlow", 5, -0.5, 4.5);
  bookHist1D("PhotonCutFlow", "Photon CutFlow", 3, -0.5, 2.5);
  // basic level histogramming
  // n-objects
  bookHist1D("met", "Missing Tranverse Energy", 40, 0, 200);
  bookHist1D("nGoodMuon", "Number of Good Muons", 10, -0.5, 9.5);
  bookHist1D("nGoodElectron", "Number of Good Electrons", 10, -0.5, 9.5);
  bookHist1D("nGoodPhoton", "Number of Good Photon", 10, -0.5, 9.5);
  bookHist1D("nGoodLepton", "Number of Good Lepton", 10, -0.5, 9.5);
  bookHist1D("nGoodTauh", "Number of Good Tauhs", 10, -0.5, 9.5);
  bookHist1D("nGoodJet", "Number of Good Jets", 10, -0.5, 9.5);
  bookHist1D("nGoodlJet", "Number of Good Light Jets", 10, -0.5, 9.5);
  bookHist1D("nGoodBJet", "Number of Good BJets", 10, -0.5, 9.5);
  // muon
  bookHist1D("Muon1Pt", "pT of 1st Muon", 40, 0.0, 500.0);
  bookHist1D("Muon2Pt", "pT of 2nd Muon", 40, 0.0, 500.0);
  bookHist1D("Muon1Eta","Eta of 1st Muon", 40, 0.0, 5.0);
  bookHist1D("Muon2Eta","Eta of 2nd Muon", 40, 0.0, 5.0);
  bookHist1D("Muon1Phi","Phi of 1st Muon", 40, 0.0, 3.3);
  bookHist1D("Muon2Phi","Phi of 2nd Muon", 40, 0.0, 3.3);
  // electron
  bookHist1D("Electron1Pt", "pT of 1st Electron", 40, 0.0, 500.0);
  bookHist1D("Electron2Pt", "pT of 2nd Electron", 40, 0.0, 500.0);
  bookHist1D("Electron1Eta","Eta of 1st Electron", 40, 0.0, 5.0);
  bookHist1D("Electron2Eta","Eta of 2nd Electron", 40, 0.0, 5.0);
  bookHist1D("Electron1Phi","Phi of 1st Electron", 40, 0.0, 3.3);
  bookHist1D("Electron2Phi","Phi of 2nd Electron", 40, 0.0, 3.3);
  // tauh
  bookHist1D("Tauh_Pt",  "pT of 1st Tauh",  40, 0.0, 500.0);
  bookHist1D("Tauh_Eta", "Eta of 1st Tauh", 40, 0.0, 5.0);
  bookHist1D("Tauh_Phi", "Phi of 1st Tauh", 40, 0.0, 3.3);
  // bjet
  bookHist1D("BJet1Pt",  "pT of 1st bJet",  40, 0.0, 500.0);
  bookHist1D("BJet1Eta", "Eta of 1st bJet", 40, 0.0, 5.0);
  bookHist1D("BJet1Phi", "Phi of 1st bJet", 40, 0.0, 3.3);
  bookHist1D("BJet2Pt",  "pT of 2nd bJet",  40, 0.0, 500.0);
  bookHist1D("BJet2Eta", "Eta of 2nd bJet", 40, 0.0, 5.0);
  bookHist1D("BJet2Phi", "Phi of 2nd bJet", 40, 0.0, 3.3);
  // jet
  bookHist1D("Jet1Pt", "pT of 1st Jet", 40, 0.0, 500.0);
  bookHist1D("Jet2Pt", "pT of 2nd Jet", 40, 0.0, 500.0);
  bookHist1D("Jet3Pt", "pT of 3rd Jet", 40, 0.0, 500.0);
  bookHist1D("Jet4Pt", "pT of 4th Jet", 40, 0.0, 500.0);
  bookHist1D("Jet5Pt", "pT of 5th Jet", 40, 0.0, 500.0); 
  bookHist1D("Jet6Pt", "pT of 6th Jet", 40, 0.0, 500.0); 

  bookHist1D("lightJet1Pt", "pT of 1st lJet", 40, 0.0, 500.0);
  bookHist1D("lightJet2Pt", "pT of 2nd lJet", 40, 0.0, 500.0);
  bookHist1D("lightJet3Pt", "pT of 3rd lJet", 40, 0.0, 500.0);
  bookHist1D("lightJet4Pt", "pT of 4th lJet", 40, 0.0, 500.0);

  bookHist1D("lightJet1Eta", "#eta (jet_{1})", 40, 0.0, 5.0);
  bookHist1D("lightJet2Eta", "#eta (jet_{2})", 40, 0.0, 5.0);
  bookHist1D("lightJet3Eta", "#eta (jet_{3})", 40, 0.0, 5.0);
  bookHist1D("lightJet4Eta", "#eta (jet_{4})", 40, 0.0, 5.0);

  //DR and DPhi among jets and Muon
  bookHist1D("DR_Jet1Muon1", "#DeltaR(Jet_{1}, Muon_{1})", 40, 0.0, 5.0);
  bookHist1D("DR_Jet1Muon2", "#DeltaR(Jet_{1}, Muon_{2})", 40, 0.0, 5.0);
  bookHist1D("DR_Jet2Muon1", "#DeltaR(Jet_{2}, Muon_{1})", 40, 0.0, 5.0);
  bookHist1D("DR_Jet2Muon2", "#DeltaR(Jet_{2}, Muon_{2})", 40, 0.0, 5.0);
  bookHist1D("DR_Jet3Muon1", "#DeltaR(Jet_{3}, Muon_{1})", 40, 0.0, 5.0);
  bookHist1D("DR_Jet3Muon2", "#DeltaR(Jet_{3}, Muon_{2})", 40, 0.0, 5.0);
  bookHist1D("DR_Jet4Muon1", "#DeltaR(Jet_{4}, Muon_{1})", 40, 0.0, 5.0);
  bookHist1D("DR_Jet4Muon2", "#DeltaR(Jet_{4}, Muon_{2})", 40, 0.0, 5.0);

  bookHist1D("DR_bJet1Muon1", "#DeltaR(Jet_{1}, Muon_{1})", 40, 0.0, 5.0);
  bookHist1D("DR_bJet1Muon2", "#DeltaR(Jet_{1}, Muon_{2})", 40, 0.0, 5.0);
  bookHist1D("DR_bJet2Muon1", "#DeltaR(Jet_{2}, Muon_{1})", 40, 0.0, 5.0);
  bookHist1D("DR_bJet2Muon2", "#DeltaR(Jet_{2}, Muon_{2})", 40, 0.0, 5.0);

  bookHist1D("DPhi_Jet1Muon1", "#Delta#Phi(Jet_{1}, Muon_{1})", 40, 0.0, 3.3);
  bookHist1D("DPhi_Jet1Muon2", "#Delta#Phi(Jet_{1}, Muon_{2})", 40, 0.0, 3.3);
  bookHist1D("DPhi_Jet2Muon1", "#Delta#Phi(Jet_{2}, Muon_{1})", 40, 0.0, 3.3);
  bookHist1D("DPhi_Jet2Muon2", "#Delta#Phi(Jet_{2}, Muon_{2})", 40, 0.0, 3.3);
  bookHist1D("DPhi_Jet3Muon1", "#Delta#Phi(Jet_{3}, Muon_{1})", 40, 0.0, 3.3);
  bookHist1D("DPhi_Jet3Muon2", "#Delta#Phi(Jet_{3}, Muon_{2})", 40, 0.0, 3.3);
  bookHist1D("DPhi_Jet4Muon1", "#Delta#Phi(Jet_{4}, Muon_{1})", 40, 0.0, 3.3);
  bookHist1D("DPhi_Jet4Muon2", "#Delta#Phi(Jet_{4}, Muon_{2})", 40, 0.0, 3.3);

  bookHist1D("DPhi_bJet1Muon1", "#Delta#Phi(Jet_{1}, Muon_{1})", 40, 0.0, 3.3);
  bookHist1D("DPhi_bJet1Muon2", "#Delta#Phi(Jet_{1}, Muon_{2})", 40, 0.0, 3.3);
  bookHist1D("DPhi_bJet2Muon1", "#Delta#Phi(Jet_{2}, Muon_{1})", 40, 0.0, 3.3);
  bookHist1D("DPhi_bJet2Muon2", "#Delta#Phi(Jet_{2}, Muon_{2})", 40, 0.0, 3.3);

  //DR and DPhi among jets and Electron
  bookHist1D("DR_Jet1Electron1", "#DeltaR(Jet_{1}, Electron_{1})", 40, 0.0, 5.0);
  bookHist1D("DR_Jet1Electron2", "#DeltaR(Jet_{1}, Electron_{2})", 40, 0.0, 5.0);
  bookHist1D("DR_Jet2Electron1", "#DeltaR(Jet_{2}, Electron_{1})", 40, 0.0, 5.0);
  bookHist1D("DR_Jet2Electron2", "#DeltaR(Jet_{2}, Electron_{2})", 40, 0.0, 5.0);
  bookHist1D("DR_Jet3Electron1", "#DeltaR(Jet_{3}, Electron_{1})", 40, 0.0, 5.0);
  bookHist1D("DR_Jet3Electron2", "#DeltaR(Jet_{3}, Electron_{2})", 40, 0.0, 5.0);
  bookHist1D("DR_Jet4Electron1", "#DeltaR(Jet_{4}, Electron_{1})", 40, 0.0, 5.0);
  bookHist1D("DR_Jet4Electron2", "#DeltaR(Jet_{4}, Electron_{2})", 40, 0.0, 5.0);

  bookHist1D("DR_bJet1Electron1", "#DeltaR(Jet_{1}, Electron_{1})", 40, 0.0, 5.0);
  bookHist1D("DR_bJet1Electron2", "#DeltaR(Jet_{1}, Electron_{2})", 40, 0.0, 5.0);
  bookHist1D("DR_bJet2Electron1", "#DeltaR(Jet_{2}, Electron_{1})", 40, 0.0, 5.0);
  bookHist1D("DR_bJet2Electron2", "#DeltaR(Jet_{2}, Electron_{2})", 40, 0.0, 5.0);

  bookHist1D("DPhi_Jet1Electron1", "#Delta#Phi(Jet_{1}, Electron_{1})", 40, 0.0, 3.3);
  bookHist1D("DPhi_Jet1Electron2", "#Delta#Phi(Jet_{1}, Electron_{2})", 40, 0.0, 3.3);
  bookHist1D("DPhi_Jet2Electron1", "#Delta#Phi(Jet_{2}, Electron_{1})", 40, 0.0, 3.3);
  bookHist1D("DPhi_Jet2Electron2", "#Delta#Phi(Jet_{2}, Electron_{2})", 40, 0.0, 3.3);
  bookHist1D("DPhi_Jet3Electron1", "#Delta#Phi(Jet_{3}, Electron_{1})", 40, 0.0, 3.3);
  bookHist1D("DPhi_Jet3Electron2", "#Delta#Phi(Jet_{3}, Electron_{2})", 40, 0.0, 3.3);
  bookHist1D("DPhi_Jet4Electron1", "#Delta#Phi(Jet_{4}, Electron_{1})", 40, 0.0, 3.3);
  bookHist1D("DPhi_Jet4Electron2", "#Delta#Phi(Jet_{4}, Electron_{2})", 40, 0.0, 3.3);

  bookHist1D("DPhi_bJet1Electron1", "#Delta#Phi(Jet_{1}, Electron_{1})", 40, 0.0, 3.3);
  bookHist1D("DPhi_bJet1Electron2", "#Delta#Phi(Jet_{1}, Electron_{2})", 40, 0.0, 3.3);
  bookHist1D("DPhi_bJet2Electron1", "#Delta#Phi(Jet_{2}, Electron_{1})", 40, 0.0, 3.3);
  bookHist1D("DPhi_bJet2Electron2", "#Delta#Phi(Jet_{2}, Electron_{2})", 40, 0.0, 3.3);

  //PT among light Jets
  bookHist1D("PT_lightJet1Jet2", "PT (jet_{1}, jet_{2})", 40, 0, 500);
  bookHist1D("PT_lightJet1Jet3", "PT (jet_{1}, jet_{3})", 40, 0, 500);
  bookHist1D("PT_lightJet1Jet4", "PT (jet_{1}, jet_{4})", 40, 0, 500);
  bookHist1D("PT_lightJet2Jet3", "PT (jet_{2}, jet_{3})", 40, 0, 500);
  bookHist1D("PT_lightJet2Jet4", "PT (jet_{2}, jet_{4})", 40, 0, 500);
  bookHist1D("PT_lightJet3Jet4", "PT (jet_{3}, jet_{4})", 40, 0, 500);

  //IM among light Jets
  bookHist1D("IM_lightJet1Jet2", "M_{inv} (jet_{1}, jet_{2})", 40, 0, 500);
  bookHist1D("IM_lightJet1Jet3", "M_{inv} (jet_{1}, jet_{3})", 40, 0, 500);
  bookHist1D("IM_lightJet1Jet4", "M_{inv} (jet_{1}, jet_{4})", 40, 0, 500);
  bookHist1D("IM_lightJet2Jet3", "M_{inv} (jet_{2}, jet_{3})", 40, 0, 500);
  bookHist1D("IM_lightJet2Jet4", "M_{inv} (jet_{2}, jet_{4})", 40, 0, 500);
  bookHist1D("IM_lightJet3Jet4", "M_{inv} (jet_{3}, jet_{4})", 40, 0, 500);

  //DR among light Jets
  bookHist1D("DR_lightJet1Jet2", "#DeltaR(jet_{1}, jet_{2})", 40, 0.0, 5.0);
  bookHist1D("DR_lightJet1Jet3", "#DeltaR(jet_{1}, jet_{3})", 40, 0.0, 5.0);
  bookHist1D("DR_lightJet1Jet4", "#DeltaR(jet_{1}, jet_{4})", 40, 0.0, 5.0);
  bookHist1D("DR_lightJet2Jet3", "#DeltaR(jet_{2}, jet_{3})", 40, 0.0, 5.0);
  bookHist1D("DR_lightJet2Jet4", "#DeltaR(jet_{2}, jet_{4})", 40, 0.0, 5.0);
  bookHist1D("DR_lightJet3Jet4", "#DeltaR(jet_{3}, jet_{4})", 40, 0.0, 5.0);

  //DeltaPhi among light jets
  bookHist1D("DPhi_lightJet1Jet2", "#Delta #phi(jet_{1}, jet_{2})", 40, 0.0, 3.3);
  bookHist1D("DPhi_lightJet1Jet3", "#Delta #phi(jet_{1}, jet_{3})", 40, 0.0, 3.3);
  bookHist1D("DPhi_lightJet1Jet4", "#Delta #phi(jet_{1}, jet_{4})", 40, 0.0, 3.3);
  bookHist1D("DPhi_lightJet2Jet3", "#Delta #phi(jet_{2}, jet_{3})", 40, 0.0, 3.3);
  bookHist1D("DPhi_lightJet2Jet4", "#Delta #phi(jet_{2}, jet_{4})", 40, 0.0, 3.3);
  bookHist1D("DPhi_lightJet3Jet4", "#Delta #phi(jet_{3}, jet_{4})", 40, 0.0, 3.3);

  //DeltaEta among light jets
  bookHist1D("DEta_lightJet1Jet2", "#Delta #eta(jet_{1}, jet_{2})", 40, 0.0, 5.0);
  bookHist1D("DEta_lightJet1Jet3", "#Delta #eta(jet_{1}, jet_{3})", 40, 0.0, 5.0);
  bookHist1D("DEta_lightJet1Jet4", "#Delta #eta(jet_{1}, jet_{4})", 40, 0.0, 5.0);
  bookHist1D("DEta_lightJet2Jet3", "#Delta #eta(jet_{2}, jet_{3})", 40, 0.0, 5.0);
  bookHist1D("DEta_lightJet2Jet4", "#Delta #eta(jet_{2}, jet_{4})", 40, 0.0, 5.0);
  bookHist1D("DEta_lightJet3Jet4", "#Delta #eta(jet_{3}, jet_{4})", 40, 0.0, 5.0);

  bookHist1D("DR_minimum_allJets",    "minimum #DeltaR_{all jets}",   40, 0, 5.0);
  bookHist1D("DR_maximum_allJets",    "maximum #DeltaR_{all jets}",   40, 0, 5.0);
  bookHist1D("IM_minimum_lightJets",  "minimum InvMass_{light jets}", 40, 0, 500.0);
  // basic histogrammig ends here
  // tauh-leadbj 
  bookHist1D("DR_tauh_leadbj",   "#DeltaR (#tau_{h}, lead b-Jet)",    25, 0, 5);
  bookHist1D("DPhi_tauh_leadbj", "#Delta#phi (#tau_{h}, lead b-Jet)", 32, -3.2, 3.2);
  bookHist1D("DEta_tauh_leadbj", "#Delta#eta (#tau_{h}, lead b-Jet)", 50, -5, 5);
  // tauh-leadlj
  bookHist1D("DR_tauh_leadlj",   "#DeltaR (#tau_{h}, lead light-Jet)",    25, 0, 5);
  bookHist1D("DPhi_tauh_leadlj", "#Delta#phi (#tau_{h}, lead light-Jet)", 32, -3.2, 3.2);
  bookHist1D("DEta_tauh_leadlj", "#Delta#eta (#tau_{h}, lead light-Jet)", 50, -5, 5);
  // tauh-leadj
  bookHist1D("DR_tauh_leadj",    "#DeltaR (#tau_{h}, lead Jet)",    25, 0, 5);
  bookHist1D("DPhi_tauh_leadj",  "#Delta#phi (#tau_{h}, lead Jet)", 32, -3.2, 3.2);
  bookHist1D("DEta_tauh_leadj",  "#Delta#eta (#tau_{h}, lead Jet)", 50, -5, 5);
  // met-tauh,leadbj,leadlj 
  bookHist1D("DPhi_met_tauh",    "#Delta#phi (#slash{E}_{T}, #tau_{h})", 32, -3.2, 3.2);
  bookHist1D("DPhi_met_leadbj",  "#Delta#phi (#slash{E}_{T}, lead b-Jet)", 32, -3.2, 3.2);
  bookHist1D("DPhi_met_leadlj",  "#Delta#phi (#slash{E}_{T}, lead light-Jet)", 32, -3.2, 3.2);
  // leadlj-leadbj 
  bookHist1D("DR_leadbj_leadlj",   "#DeltaR (lead b-Jet, lead light-Jet)", 25, 0, 5);
  bookHist1D("DPhi_leadbj_leadlj", "#Delta#phi (lead b-Jet, lead light-Jet)", 32, -3.2, 3.2);
  bookHist1D("DEta_leadbj_leadlj", "#Delta#eta (lead b-Jet, lead light-Jet)", 50, -5, 5);

  bookHist1D("PT_RecoTau",  "Reconstructed #tau p_{T} (GeV) : Collinear approx", 50, 0, 500);
  bookHist1D("Eta_RecoTau", "Reconstructed #tau #eta : Collinear approx",        25, -2.5, 2.5);
  bookHist1D("PT_RecoNu",   "Reconstructed #nu p_{T} (GeV) : Collinear approx",  50, 0, 500);
  bookHist1D("Eta_RecoNu",  "Reconstructed #nu #eta : Collinear approx",         25, -2.5, 2.5);

  bookHist1D("BDT_Score",         "BDT Response", 40, -0.5, 0.5);
  bookHist1D("BDT_Score_fineBin", "BDT Response", 2000, -1.0, 1.0);  
}

void HistBooker::bookDLHistograms () {
  bookHist1D("evtCutFlow",        "Event Weight Sum",    10, -0.5, 9.5);
  bookHist1D("evtCutFlowWt",      "HepMCEvent Weight",   10, -0.5, 9.5);

  bookHist1D("XlepPt",  "#chi lepton p_{T} (GeV)", 50, 0., 200.);
  bookHist1D("XlepEta", "#chi lepton #eta",        25, -2.5, 2.5);
  bookHist1D("MET_along_tauh", "#slash{E}_{T} along #tau_{h}",                 50, 0., 300.);
  bookHist1D("FracMetAlongTauh", "#slash{E}_{T} along #tau_{h}/#slash{E}_{T}", 50, 0, 1);
  bookHist1D("DR_xlep_tauh",     "#DeltaR (#mu from #chi, #tau_{h})",          25, 0, 5);
  bookHist1D("DPhi_xlep_tauh",   "#Delta#phi (#mu from #chi, #tau_{h})",       32, -3.2, 3.2);
  bookHist1D("DEta_xlep_tauh",   "#Delta#eta (#mu from #chi, #tau_{h})",       50, -5, 5);
  bookHist1D("DR_xlep_leadlj",   "#DeltaR (#mu from #chi, lead light-Jet)",         25, 0, 5);
  bookHist1D("DPhi_xlep_leadlj", "#Delta#phi (#mu from #chi, lead light-Jet)",      32, -3.2, 3.2);
  bookHist1D("DEta_xlep_leadlj", "#Delta#eta (#mu from #chi, lead light-Jet)",      50, -5, 5);
  bookHist1D("DR_xlep_leadj",    "#DeltaR (#mu from #chi, lead Jet)",         25, 0, 5);
  bookHist1D("DPhi_xlep_leadj",  "#Delta#phi (#mu from #chi, lead Jet)",      32, -3.2, 3.2);
  bookHist1D("DEta_xlep_leadj",  "#Delta#eta (#mu from #chi, lead Jet)",      50, -5, 5);
  bookHist1D("DPhi_xlep_met",    "#Delta#phi (#mu from #chi, #slash{E}_{T})", 32, -3.2, 3.2);
  bookHist1D("EffectiveMass_t1", "Effective mass of top (GeV)",               50, 0, 500);
  bookHist1D("InvMass_X_test",   "Invariant mass of #chi (GeV) : improper",   50, 0, 200);
  bookHist1D("MT_X",             "Transverse mass of #chi (GeV)",             50, 0, 500);
  bookHist1D("InvM_coln_XlepTauh", "Collinear mass (GeV)",                    50, 0, 200);
  bookHist1D("InvM_coln_XlepTauh_MuMu", "Collinear mass (GeV) : #mu#mu",      50, 0, 200);
  bookHist1D("InvM_coln_XlepTauh_EleMu", "Collinear mass (GeV) : ele#mu",     50, 0, 200);
  bookHist1D("DR_XlepTauh_MuMu",      "DR_XlepTauh_MuMu",           25, 0, 5);
  bookHist1D("DR_XlepTauh_EleMu",      "DR_XlepTauh_EleMu",           25, 0, 5);
  bookHist1D("DR_WlepTauh_MuMu",      "DR_WlepTauh_MuMu",           25, 0, 5);
  bookHist1D("DR_WlepTauh_EleMu",      "DR_WlepTauh_EleMu",           25, 0, 5);

  bookHist1D("DPhi_XlepTauh_MuMu",      "DPhi_XlepTauh_MuMu",           16, 0, 3.2);
  bookHist1D("DPhi_XlepTauh_EleMu",     "DPhi_XlepTauh_EleMu",          16, 0, 3.2);
  bookHist1D("DPhi_WlepTauh_MuMu",      "DPhi_WlepTauh_MuMu",           16, 0, 3.2);
  bookHist1D("DPhi_WlepTauh_EleMu",     "DPhi_WlepTauh_EleMu",          16, 0, 3.2);

  bookHist2D("DR_XlepTauh_Vs_WlepTauh_MuMu", "DR_XlepTauh_Vs_WlepTauh_MuMu", 25, 0, 5, 25, 0, 5);
  bookHist2D("DR_XlepTauh_Vs_WlepTauh_EleMu", "DR_XlepTauh_Vs_WlepTauh_EleMu", 25, 0, 5, 25, 0, 5);
  bookHist1D("IsOppCharge_MuMu", "", 2, -0.5, 1.5);
  bookHist1D("IsOppCharge_EleMu", "", 2, -0.5, 1.5);
  bookHist2D("IsOppCharge_DR_XlepTauh_MuMu", "", 2, -0.5, 1.5, 25, 0, 5);
  bookHist2D("IsOppCharge_DR_WlepTauh_MuMu", "", 2, -0.5, 1.5, 25, 0, 5);
  bookHist2D("IsOppCharge_DR_XlepTauh_EleMu", "", 2, -0.5, 1.5, 25, 0, 5);
  bookHist2D("IsOppCharge_DR_WlepTauh_EleMu", "", 2, -0.5, 1.5, 25, 0, 5);

  bookHist1D("MT_xlep_met", "Transverse mass (#mu from #chi, #slash{E}_{T}) (GeV)", 50, 0, 500);

  bookHist1D("WlepPt",  "W#pm lepton p_{T} (GeV)", 50, 0., 200.);
  bookHist1D("WlepEta", "W#pm lepton #eta",        25, -2.5, 2.5);
  bookHist1D("DR_wlep_leadbj",    "#DeltaR (lepton from W, lead b-Jet)",    25, 0, 5);
  bookHist1D("DPhi_wlep_leadbj",  "#Delta#phi (lepton from W, lead b-Jet)", 32, -3.2, 3.2);
  bookHist1D("DEta_wlep_leadbj",  "#Delta#eta (lepton from W, lead b-Jet)", 50, -5, 5);
  bookHist1D("MT_wlep_met",       "M_{T} (lepton from W, #slash{E}_{T}) (GeV)",50, 0, 200);
  bookHist1D("InvM_coln_WlepTauh","Collinear mass of lepton from W and #tau (GeV) : Test", 100, 0, 1000);      

  bookHist1D("Scalarsumpt_xlep_wlep", "Scalar sum p_{T} (GeV)", 50, 0, 400);
  bookHist1D("Vectorsumpt_xlep_wlep", "Vector sum p_{T} (GeV)", 50, 0, 400);
  bookHist1D("DR_xlep_wlep",          "#DeltaR", 25, 0, 5);
  bookHist1D("DPhi_xlep_wlep",        "#Delta#phi", 32, -3.2, 3.2);
  bookHist1D("DEta_xlep_wlep",        "#Delta#eta", 50, -5, 5);
  bookHist1D("DR_xlep_leadbj",        "#DeltaR", 25, 0, 5);
  bookHist1D("DPhi_xlep_leadbj",      "#Delta#phi", 32, -3.2, 3.2);
  bookHist1D("DEta_xlep_leadbj",      "#Delta#eta", 50, -5, 5);
  bookHist1D("DR_wlep_leadlj",        "#DeltaR", 25, 0, 5);
  bookHist1D("DPhi_wlep_leadlj",      "#Delta#phi", 32, -3.2, 3.2);
  bookHist1D("DEta_wlep_leadlj",      "#Delta#eta", 50, -5, 5);
  bookHist1D("DR_wlep_tauh",          "#DeltaR", 25, 0, 5);
  bookHist1D("DPhi_wlep_tauh",        "#Delta#phi", 32, -3.2, 3.2);
  bookHist1D("DEta_wlep_tauh",        "#Delta#eta", 50, -5, 5);
  bookHist1D("Scalarsumpt_vis",       "Scalar sum p_{T} (GeV) visible", 100, 0, 1000);
  bookHist1D("Scalarsumpt",           "Scalar sum p_{T} (GeV)", 100, 0, 1000);
  bookHist1D("Effectivemass_vis",     "Effectivemass_vis (GeV)", 100, 0, 1000);
  bookHist1D("Effectivemass",         "Effectivemass (GeV)", 100, 0, 1000);

  bookHist1D("HT_Jets",           "H_{T} (GeV)",   60, 0, 600);
  bookHist1D("minDr_XlepJets",    "minimum #Delta R [#mu from #chi, all jets]",    25, 0, 5.0);
  bookHist1D("maxDr_XlepJets",    "maximum #Delta R [#mu from #chi, all jets]",    25, 0, 5.0);
  bookHist1D("minDphi_XlepJets",  "minimum #Delta #Phi [#mu from #chi, all jets]", 32, 0, 3.2);
  bookHist1D("maxDphi_XlepJets",  "maximum #Delta #Phi [#mu from #chi, all jets]", 32, 0, 3.2);
  bookHist1D("minDr_WlepJets",    "minimum #Delta R [lepton from W, all jets]",    25, 0, 5.0);
  bookHist1D("maxDr_WlepJets",    "maximum #Delta R [lepton from W, all jets]",    25, 0, 5.0);
  bookHist1D("minDphi_WlepJets",  "minimum #Delta #Phi [lepton from W, all jets]", 32, 0, 3.2);
  bookHist1D("maxDphi_WlepJets",  "maximum #Delta #Phi [lepton from W, all jets]", 32, 0, 3.2);
  bookHist1D("smin",              "#sqrt(#sedge(s)_{min})",                      100, 0, 1000);
}

void HistBooker::bookSLHistograms () {
  bookHist1D("evtCutFlow", "Event CutFlow", 9, -0.5, 8.5);
  bookHist1D("evtCutFlowWt", "HepMCEvent Weight", 9, -0.5, 8.5);

  bookHist1D("w_massdiff",    "mass diff to check w jets",             25, -100, 100);
  // book histograms for one leg
  bookHist1D("Wj1Pt",         "p_{T} (GeV) of leading-jet from W",     20, 0, 300);
  bookHist1D("Wj2Pt",         "p_{T} (GeV) of subleading-jet from W",  20, 0, 300);
  bookHist1D("WPt",           "p_{T} (GeV) of W",                      20, 0, 300);
  bookHist1D("InvM_Wj1Wj2",   "Invariant mass [W-jet1, W-jet2] (GeV)", 50, 0, 200);
  bookHist1D("DR_Wj1Wj2",     "#Delta R [W_{j1}, W_{j2}]",             25, 0, 5);
  bookHist1D("DPhi_Wj1Wj2",   "#Delta #Phi [W_{j1}, W_{j2}]",          32, -3.2, 3.2);
  bookHist1D("DEta_Wj1Wj2",   "#Delta #eta [W_{j1}, W_{j2}]",          25, -5, 5);
  bookHist1D("DR_w_leadb",    "#Delta R [W, leading b-jet]",           25, 0, 5);
  bookHist1D("DPhi_w_leadb",  "#Delta #Phi [W, leading b-jet]",        32, -3.2, 3.2);
  bookHist1D("DR_wj1_leadb",  "#Delta R [W_{j1}, leading b-jet]",     25, 0, 5);
  bookHist1D("DPhi_wj1_leadb","#Delta #Phi [W_{j1}, leading b-jet]",  32, -3.2, 3.2);
  bookHist1D("DR_wj2_leadb",  "#Delta R [W_{j2}, leading b-jet]",     25, 0, 5);
  bookHist1D("DPhi_wj2_leadb","#Delta #Phi [W_{j2}, leading b-jet]",  32, -3.2, 3.2);
  bookHist1D("InvM_w_leadbj", "Invariant Mass [W, leading b-jet] (GeV)", 50, 0, 500);
  // book histograms for another leg
  bookHist1D("XlepPt",          "p_{T} of #mu from #chi",         50, 0, 200);
  bookHist1D("XlepEta",         "#eta of #mu from #chi",          32, -3.2, 3.2);
  bookHist1D("M_coll_Xlep",     "M_{coll} (GeV)",                 200, 0, 1000);
  bookHist1D("M_coll_test",     "M_{coll} (GeV)",                 200, 0, 1000);
  bookHist1D("InvM_Xleadlj",    "Invariant mass [#chi, leading light-jet] (GeV)",  100, 0, 1000);
  bookHist1D("DR_Xlep_leadlj",  "#Delta R [#mu from #chi, leading light-jet]",     25, 0, 5);
  bookHist1D("DPhi_Xlep_leadlj","#Delta #Phi [#mu from #chi, leading light-jet]",  32, -3.2, 3.2);
  bookHist1D("DR_Xtauh_leadlj", "#Delta R [#tau_{h} from #chi, leading light-jet]", 25, 0, 5);
  bookHist1D("DPhi_Xtauh_leadlj","#Delta #Phi [#tau_{h} from #chi, leading light-jet]", 32, -3.2, 3.2);
  bookHist1D("DR_X_leadlj",     "#Delta R [#chi, leading light-jet]",        25, 0, 5);
  bookHist1D("DPhi_X_leadlj",   "#Delta #Phi [#chi, leading light-jet]",     32, -3.2, 3.2);
  bookHist1D("DR_Xlep_Xtau",    "#Delta R [(#mu, #tau) form #chi]",          25, 0, 5);
  bookHist1D("DPhi_Xlep_Xtau",  "#Delta #Phi [(#mu, #tau) form #chi]",       32, -3.2, 3.2);
  // communication between 2 legs
  bookHist1D("DR_Xlep_wj1",     "#Delta R",             25, 0, 5);     
  bookHist1D("DPhi_Xlep_wj1",   "#Delta #Phi",          32, -3.2, 3.2);
  bookHist1D("DEta_Xlep_wj1",   "#Delta #eta",          32, -3.2, 3.2);    
  bookHist1D("DR_Xlep_wj2",     "#Delta R",             25, 0, 5);
  bookHist1D("DPhi_Xlep_wj2",   "#Delta #Phi",          32, -3.2, 3.2);
  bookHist1D("DEta_Xlep_wj2",   "#Delta #eta",          32, -3.2, 3.2);    
  bookHist1D("DR_Xlep_w",       "#Delta R",             25, 0, 5); 
  bookHist1D("DPhi_Xlep_w",     "#Delta #Phi",          32, -3.2, 3.2);    
  bookHist1D("DEta_Xlep_w",     "#Delta #eta",          32, -3.2, 3.2);    
  bookHist1D("DR_Xlep_leadbj",       "#Delta R",             25, 0, 5);
  bookHist1D("DPhi_Xlep_leadbj",     "#Delta #Phi",          32, -3.2, 3.2);    
  bookHist1D("DEta_Xlep_leadbj",     "#Delta #eta",          32, -3.2, 3.2); 
  bookHist1D("DR_Xtau_wj1",     "#Delta R",             25, 0, 5);
  bookHist1D("DPhi_Xtau_wj1",   "#Delta #Phi",          32, -3.2, 3.2); 
  bookHist1D("DEta_Xtau_wj1",   "#Delta #eta",          32, -3.2, 3.2); 
  bookHist1D("DR_Xtau_wj2",     "#Delta R",             25, 0, 5);
  bookHist1D("DPhi_Xtau_wj2",   "#Delta #Phi",          32, -3.2, 3.2); 
  bookHist1D("DEta_Xtau_wj2",   "#Delta #eta",          32, -3.2, 3.2); 
  bookHist1D("DR_Xtau_w",       "#Delta R",             25, 0, 5);
  bookHist1D("DPhi_Xtau_w",     "#Delta #Phi",          32, -3.2, 3.2); 
  bookHist1D("DEta_Xtau_w",     "#Delta #eta",          32, -3.2, 3.2); 
  bookHist1D("DR_Xtau_leadbj",       "#Delta R",             25, 0, 5);
  bookHist1D("DPhi_Xtau_leadbj",     "#Delta #Phi",          32, -3.2, 3.2); 
  bookHist1D("DEta_Xtau_leadbj",     "#Delta #eta",          32, -3.2, 3.2); 
  //bookHist1D("DR_leadbj_leadlj",       "#Delta R",             25, 0, 5);
  //bookHist1D("DPhi_leadbj_leadlj",     "#Delta #Phi",          32, -3.2, 3.2); 
  //bookHist1D("DEta_leadbj_leadlj",     "#Delta #eta",          32, -3.2, 3.2); 
  bookHist1D("DR_w_leadlj",       "#Delta R",             25, 0, 5);
  bookHist1D("DPhi_w_leadlj",     "#Delta #Phi",          32, -3.2, 3.2); 
  bookHist1D("DEta_w_leadlj",     "#Delta #eta",          32, -3.2, 3.2); 
  bookHist1D("DR_X_leadbj",          "#Delta R",             25, 0, 5);
  bookHist1D("DPhi_X_leadbj",        "#Delta #Phi",          32, -3.2, 3.2); 
  bookHist1D("DEta_X_leadbj",        "#Delta #eta",          32, -3.2, 3.2); 
  bookHist1D("DR_X_w",          "#Delta R",             25, 0, 5);
  bookHist1D("DPhi_X_w",        "#Delta #Phi",          32, -3.2, 3.2); 
  bookHist1D("DEta_X_w",        "#Delta #eta",          32, -3.2, 3.2); 
  bookHist1D("DR_t1_t2",        "#Delta R",             25, 0, 5);
  bookHist1D("DPhi_t1_t2",      "#Delta #Phi",          32, -3.2, 3.2); 
  bookHist1D("DEta_t1_t2",      "#Delta #eta",          32, -3.2, 3.2); 
  bookHist1D("Total_mass",      "mass (GeV)",           100, 0, 2000); 
  bookHist1D("Total_vector_pt", "vec sum p_{T} (GeV)",  40, 0, 500); 
  bookHist1D("Total_scalar_pt", "scalar sum p_{T} (GeV)", 100, 0, 2000); 
  // some more high level variables
  bookHist1D("HT_Jets",         "H_{T} (GeV)",                100, 0, 2000);
  bookHist1D("minDr_XlepJets",  "min #Delta R [Xlep, jets]",  25, 0, 5);
  bookHist1D("maxDr_XlepJets",  "max #Delta R [Xlep, jets]",  25, 0, 5);
  bookHist1D("minDphi_XlepJets","minDphi_XlepJets",           32, 0, 3.2);
  bookHist1D("maxDphi_XlepJets", "maxDphi_XlepJets",          32, 0, 3.2);

  bookHist1D("smin",            "#sqrt(#sedge(s)_{min})",   100, 0, 1000);
}
