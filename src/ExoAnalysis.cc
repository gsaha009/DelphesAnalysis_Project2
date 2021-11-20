#include <iostream>
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <functional>
#include <numeric>
#include <string>
#include <climits>
#include <cassert>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <utility>
#include <typeinfo>
#include <memory>

#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TRandom.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TH1K.h"

#include "AnaUtil.h"
#include "ZCandidate.h"
#include "ExoAnalysis.h"


using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::ostringstream;
using std::vector;
using std::map;
using std::pair;
using std::setprecision;
using std::setw;
using std::ios;
using std::setiosflags;
using std::resetiosflags;

using namespace ZSpace;
//--------------------
//Constructor
//---------------------
ExoAnalysis::ExoAnalysis()
{
  fileList_.clear();
  hmap_.clear();
}

// ------------------
// Destructor
// -----------------
ExoAnalysis::~ExoAnalysis()
{}

// ------------------------------------------------------------------
// Prepare for the run, do necessary initialization etc.
// -----------------------------------------------------------------
bool ExoAnalysis::beginJob()
{
  histf = new TFile (histFile_.c_str(), "RECREATE");
  if (!histf) return false;
  histf->cd();
  //histf->mkdir("ObjectSelection");
  histf->mkdir("LepPairSelection");
  //histf->mkdir("Analysis");

  cout << "********** Booking Histograms **************" << endl;
  bookHistograms();

  if(_createMVATree) {
    #ifdef TRUE_CPP14
    skimObj_ = std::make_unique <MVASkim> (_mvaInputFile);
    #else
    skimObj_ = std::unique_ptr <MVASkim>(new MVASkim (_mvaInputFile));
    #endif
    if (!skimObj_) return false;
  }

   else if (_readMVA) {
     #ifdef TRUE_CPP14
     _mvaObj = std::make_unique<MVAnalysis>(_MVAnetwork, _MVAxmlFile);
     #else
     _mvaObj = std::unique_ptr <MVAnalysis>(new MVAnalysis (_MVAnetwork, _MVAxmlFile));
     #endif
     if (!_mvaObj) return false;
   }

   eventLoop(treeReader);

   return true;
}


// -----------------
// Book histograms
// -----------------
void ExoAnalysis::bookHistograms()
{
  histf->cd();
  //histf->cd("ObjectSelection");
  new TH1D("MuonCutFlow", "Muon CutFlow", 3, -0.5, 2.5);
  new TH1D("ElectronCutFlow", "Electron CutFlow", 3, -0.5, 2.5);
  new TH1D("TauhCutFlow", "Tauh CutFlow", 3, -0.5, 2.5);
  new TH1D("JetCutFlow", "Jet CutFlow", 5, -0.5, 4.5);
  new TH1D("PhotonCutFlow", "Photon CutFlow", 3, -0.5, 2.5);

  //histf->cd();
  //histf->cd("Analysis");
  new TH1D("evtCutFlow", "Event CutFlow", 8, -0.5, 7.5);
  new TH1D("evtCutFlowWt", "HepMCEvent Weight", 8, -0.5, 7.5);
  //new TH1D("nvtx", "Number of Primary vertices", 60, -0.5, 59.5);

  new TH1D("met", "Missing Tranverse Energy", 40, 0, 200);
  new TH1D("nGoodMuon", "Number of Good Muons", 10, -0.5, 9.5);
  new TH1D("nGoodElectron", "Number of Good Electrons", 10, -0.5, 9.5);
  new TH1D("nGoodPhoton", "Number of Good Photon", 10, -0.5, 9.5);
  new TH1D("nGoodLepton", "Number of Good Lepton", 10, -0.5, 9.5);
  new TH1D("nGoodTauh", "Number of Good Tauhs", 10, -0.5, 9.5);
  new TH1D("nGoodJet", "Number of Good Jets", 10, -0.5, 9.5);
  new TH1D("nGoodlJet", "Number of Good Light Jets", 10, -0.5, 9.5);
  new TH1D("nGoodBJet", "Number of Good BJets", 10, -0.5, 9.5);
  
  new TH1D("Muon1Pt", "pT of 1st Muon", 40, 0.0, 500.0);
  new TH1D("Muon2Pt", "pT of 2nd Muon", 40, 0.0, 500.0);
  new TH1D("Muon1Eta","Eta of 1st Muon", 40, 0.0, 5.0);
  new TH1D("Muon2Eta","Eta of 2nd Muon", 40, 0.0, 5.0);
  new TH1D("Muon1Phi","Phi of 1st Muon", 40, 0.0, 3.3);
  new TH1D("Muon2Phi","Phi of 2nd Muon", 40, 0.0, 3.3);

  new TH1D("Electron1Pt", "pT of 1st Electron", 40, 0.0, 500.0);
  new TH1D("Electron2Pt", "pT of 2nd Electron", 40, 0.0, 500.0);
  new TH1D("Electron1Eta","Eta of 1st Electron", 40, 0.0, 5.0);
  new TH1D("Electron2Eta","Eta of 2nd Electron", 40, 0.0, 5.0);
  new TH1D("Electron1Phi","Phi of 1st Electron", 40, 0.0, 3.3);
  new TH1D("Electron2Phi","Phi of 2nd Electron", 40, 0.0, 3.3);

  new TH1D("lep1Pt", "Leading lepton p_{T} (GeV)", 50, 0., 200.);
  new TH1D("lep2Pt", "Subleading lepton p_{T} (GeV)", 50, 0., 200.);
  new TH1D("lep1Eta", "Leading lepton #eta", 25, -2.5, 2.5);
  new TH1D("lep2Eta", "Subleading lepton #eta", 25, -2.5, 2.5);
  
  new TH1D("Tauh1Pt",  "pT of 1st Tauh",  40, 0.0, 500.0);
  new TH1D("Tauh1Eta", "Eta of 1st Tauh", 40, 0.0, 5.0);
  new TH1D("Tauh1Phi", "Phi of 1st Tauh", 40, 0.0, 3.3);

  new TH1D("BJet1Pt",  "pT of 1st bJet",  40, 0.0, 500.0);
  new TH1D("BJet1Eta", "Eta of 1st bJet", 40, 0.0, 5.0);
  new TH1D("BJet1Phi", "Phi of 1st bJet", 40, 0.0, 3.3);
  new TH1D("BJet2Pt",  "pT of 2nd bJet",  40, 0.0, 500.0);
  new TH1D("BJet2Eta", "Eta of 2nd bJet", 40, 0.0, 5.0);
  new TH1D("BJet2Phi", "Phi of 2nd bJet", 40, 0.0, 3.3);

  new TH1D("Jet1Pt", "pT of 1st Jet", 40, 0.0, 500.0);
  new TH1D("Jet2Pt", "pT of 2nd Jet", 40, 0.0, 500.0);
  new TH1D("Jet3Pt", "pT of 3rd Jet", 40, 0.0, 500.0);
  new TH1D("Jet4Pt", "pT of 4th Jet", 40, 0.0, 500.0);
  new TH1D("Jet5Pt", "pT of 5th Jet", 40, 0.0, 500.0); 
  new TH1D("Jet6Pt", "pT of 6th Jet", 40, 0.0, 500.0); 

  new TH1D("lightJet1Pt", "pT of 1st lJet", 40, 0.0, 500.0);
  new TH1D("lightJet2Pt", "pT of 2nd lJet", 40, 0.0, 500.0);
  new TH1D("lightJet3Pt", "pT of 3rd lJet", 40, 0.0, 500.0);
  new TH1D("lightJet4Pt", "pT of 4th lJet", 40, 0.0, 500.0);

  new TH1D("lightJet1Eta", "#eta (jet_{1})", 40, 0.0, 5.0);
  new TH1D("lightJet2Eta", "#eta (jet_{2})", 40, 0.0, 5.0);
  new TH1D("lightJet3Eta", "#eta (jet_{3})", 40, 0.0, 5.0);
  new TH1D("lightJet4Eta", "#eta (jet_{4})", 40, 0.0, 5.0);

  //DR and DPhi among light jets and Muon
  new TH1D("DR_Jet1Muon1", "#DeltaR(Jet_{1}, Muon_{1})", 40, 0.0, 5.0);
  new TH1D("DR_Jet1Muon2", "#DeltaR(Jet_{1}, Muon_{2})", 40, 0.0, 5.0);
  new TH1D("DR_Jet2Muon1", "#DeltaR(Jet_{2}, Muon_{1})", 40, 0.0, 5.0);
  new TH1D("DR_Jet2Muon2", "#DeltaR(Jet_{2}, Muon_{2})", 40, 0.0, 5.0);
  new TH1D("DR_Jet3Muon1", "#DeltaR(Jet_{3}, Muon_{1})", 40, 0.0, 5.0);
  new TH1D("DR_Jet3Muon2", "#DeltaR(Jet_{3}, Muon_{2})", 40, 0.0, 5.0);
  new TH1D("DR_Jet4Muon1", "#DeltaR(Jet_{4}, Muon_{1})", 40, 0.0, 5.0);
  new TH1D("DR_Jet4Muon2", "#DeltaR(Jet_{4}, Muon_{2})", 40, 0.0, 5.0);

  new TH1D("DR_bJet1Muon1", "#DeltaR(Jet_{1}, Muon_{1})", 40, 0.0, 5.0);
  new TH1D("DR_bJet1Muon2", "#DeltaR(Jet_{1}, Muon_{2})", 40, 0.0, 5.0);
  new TH1D("DR_bJet2Muon1", "#DeltaR(Jet_{2}, Muon_{1})", 40, 0.0, 5.0);
  new TH1D("DR_bJet2Muon2", "#DeltaR(Jet_{2}, Muon_{2})", 40, 0.0, 5.0);

  new TH1D("DPhi_Jet1Muon1", "#Delta#Phi(Jet_{1}, Muon_{1})", 40, 0.0, 3.3);
  new TH1D("DPhi_Jet1Muon2", "#Delta#Phi(Jet_{1}, Muon_{2})", 40, 0.0, 3.3);
  new TH1D("DPhi_Jet2Muon1", "#Delta#Phi(Jet_{2}, Muon_{1})", 40, 0.0, 3.3);
  new TH1D("DPhi_Jet2Muon2", "#Delta#Phi(Jet_{2}, Muon_{2})", 40, 0.0, 3.3);
  new TH1D("DPhi_Jet3Muon1", "#Delta#Phi(Jet_{3}, Muon_{1})", 40, 0.0, 3.3);
  new TH1D("DPhi_Jet3Muon2", "#Delta#Phi(Jet_{3}, Muon_{2})", 40, 0.0, 3.3);
  new TH1D("DPhi_Jet4Muon1", "#Delta#Phi(Jet_{4}, Muon_{1})", 40, 0.0, 3.3);
  new TH1D("DPhi_Jet4Muon2", "#Delta#Phi(Jet_{4}, Muon_{2})", 40, 0.0, 3.3);

  new TH1D("DPhi_bJet1Muon1", "#Delta#Phi(Jet_{1}, Muon_{1})", 40, 0.0, 3.3);
  new TH1D("DPhi_bJet1Muon2", "#Delta#Phi(Jet_{1}, Muon_{2})", 40, 0.0, 3.3);
  new TH1D("DPhi_bJet2Muon1", "#Delta#Phi(Jet_{2}, Muon_{1})", 40, 0.0, 3.3);
  new TH1D("DPhi_bJet2Muon2", "#Delta#Phi(Jet_{2}, Muon_{2})", 40, 0.0, 3.3);

  //DR and DPhi among light jets and Electron
  new TH1D("DR_Jet1Electron1", "#DeltaR(Jet_{1}, Electron_{1})", 40, 0.0, 5.0);
  new TH1D("DR_Jet1Electron2", "#DeltaR(Jet_{1}, Electron_{2})", 40, 0.0, 5.0);
  new TH1D("DR_Jet2Electron1", "#DeltaR(Jet_{2}, Electron_{1})", 40, 0.0, 5.0);
  new TH1D("DR_Jet2Electron2", "#DeltaR(Jet_{2}, Electron_{2})", 40, 0.0, 5.0);
  new TH1D("DR_Jet3Electron1", "#DeltaR(Jet_{3}, Electron_{1})", 40, 0.0, 5.0);
  new TH1D("DR_Jet3Electron2", "#DeltaR(Jet_{3}, Electron_{2})", 40, 0.0, 5.0);
  new TH1D("DR_Jet4Electron1", "#DeltaR(Jet_{4}, Electron_{1})", 40, 0.0, 5.0);
  new TH1D("DR_Jet4Electron2", "#DeltaR(Jet_{4}, Electron_{2})", 40, 0.0, 5.0);

  new TH1D("DR_bJet1Electron1", "#DeltaR(Jet_{1}, Electron_{1})", 40, 0.0, 5.0);
  new TH1D("DR_bJet1Electron2", "#DeltaR(Jet_{1}, Electron_{2})", 40, 0.0, 5.0);
  new TH1D("DR_bJet2Electron1", "#DeltaR(Jet_{2}, Electron_{1})", 40, 0.0, 5.0);
  new TH1D("DR_bJet2Electron2", "#DeltaR(Jet_{2}, Electron_{2})", 40, 0.0, 5.0);

  new TH1D("DPhi_Jet1Electron1", "#Delta#Phi(Jet_{1}, Electron_{1})", 40, 0.0, 3.3);
  new TH1D("DPhi_Jet1Electron2", "#Delta#Phi(Jet_{1}, Electron_{2})", 40, 0.0, 3.3);
  new TH1D("DPhi_Jet2Electron1", "#Delta#Phi(Jet_{2}, Electron_{1})", 40, 0.0, 3.3);
  new TH1D("DPhi_Jet2Electron2", "#Delta#Phi(Jet_{2}, Electron_{2})", 40, 0.0, 3.3);
  new TH1D("DPhi_Jet3Electron1", "#Delta#Phi(Jet_{3}, Electron_{1})", 40, 0.0, 3.3);
  new TH1D("DPhi_Jet3Electron2", "#Delta#Phi(Jet_{3}, Electron_{2})", 40, 0.0, 3.3);
  new TH1D("DPhi_Jet4Electron1", "#Delta#Phi(Jet_{4}, Electron_{1})", 40, 0.0, 3.3);
  new TH1D("DPhi_Jet4Electron2", "#Delta#Phi(Jet_{4}, Electron_{2})", 40, 0.0, 3.3);

  new TH1D("DPhi_bJet1Electron1", "#Delta#Phi(Jet_{1}, Electron_{1})", 40, 0.0, 3.3);
  new TH1D("DPhi_bJet1Electron2", "#Delta#Phi(Jet_{1}, Electron_{2})", 40, 0.0, 3.3);
  new TH1D("DPhi_bJet2Electron1", "#Delta#Phi(Jet_{2}, Electron_{1})", 40, 0.0, 3.3);
  new TH1D("DPhi_bJet2Electron2", "#Delta#Phi(Jet_{2}, Electron_{2})", 40, 0.0, 3.3);

  //PT among light Jets
  new TH1D("PT_lightJet1Jet2", "PT (jet_{1}, jet_{2})", 40, 0, 500);
  new TH1D("PT_lightJet1Jet3", "PT (jet_{1}, jet_{3})", 40, 0, 500);
  new TH1D("PT_lightJet1Jet4", "PT (jet_{1}, jet_{4})", 40, 0, 500);
  new TH1D("PT_lightJet2Jet3", "PT (jet_{2}, jet_{3})", 40, 0, 500);
  new TH1D("PT_lightJet2Jet4", "PT (jet_{2}, jet_{4})", 40, 0, 500);
  new TH1D("PT_lightJet3Jet4", "PT (jet_{3}, jet_{4})", 40, 0, 500);


  //IM among light Jets
  new TH1D("IM_lightJet1Jet2", "M_{inv} (jet_{1}, jet_{2})", 40, 0, 500);
  new TH1D("IM_lightJet1Jet3", "M_{inv} (jet_{1}, jet_{3})", 40, 0, 500);
  new TH1D("IM_lightJet1Jet4", "M_{inv} (jet_{1}, jet_{4})", 40, 0, 500);
  new TH1D("IM_lightJet2Jet3", "M_{inv} (jet_{2}, jet_{3})", 40, 0, 500);
  new TH1D("IM_lightJet2Jet4", "M_{inv} (jet_{2}, jet_{4})", 40, 0, 500);
  new TH1D("IM_lightJet3Jet4", "M_{inv} (jet_{3}, jet_{4})", 40, 0, 500);


  //DR among light Jets
  new TH1D("DR_lightJet1Jet2", "#DeltaR(jet_{1}, jet_{2})", 40, 0.0, 5.0);
  new TH1D("DR_lightJet1Jet3", "#DeltaR(jet_{1}, jet_{3})", 40, 0.0, 5.0);
  new TH1D("DR_lightJet1Jet4", "#DeltaR(jet_{1}, jet_{4})", 40, 0.0, 5.0);
  new TH1D("DR_lightJet2Jet3", "#DeltaR(jet_{2}, jet_{3})", 40, 0.0, 5.0);
  new TH1D("DR_lightJet2Jet4", "#DeltaR(jet_{2}, jet_{4})", 40, 0.0, 5.0);
  new TH1D("DR_lightJet3Jet4", "#DeltaR(jet_{3}, jet_{4})", 40, 0.0, 5.0);

  //DeltaPhi among light jets
  new TH1D("DPhi_lightJet1Jet2", "#Delta #phi(jet_{1}, jet_{2})", 40, 0.0, 3.3);
  new TH1D("DPhi_lightJet1Jet3", "#Delta #phi(jet_{1}, jet_{3})", 40, 0.0, 3.3);
  new TH1D("DPhi_lightJet1Jet4", "#Delta #phi(jet_{1}, jet_{4})", 40, 0.0, 3.3);
  new TH1D("DPhi_lightJet2Jet3", "#Delta #phi(jet_{2}, jet_{3})", 40, 0.0, 3.3);
  new TH1D("DPhi_lightJet2Jet4", "#Delta #phi(jet_{2}, jet_{4})", 40, 0.0, 3.3);
  new TH1D("DPhi_lightJet3Jet4", "#Delta #phi(jet_{3}, jet_{4})", 40, 0.0, 3.3);

  //DeltaEta among light jets
  new TH1D("DEta_lightJet1Jet2", "#Delta #eta(jet_{1}, jet_{2})", 40, 0.0, 5.0);
  new TH1D("DEta_lightJet1Jet3", "#Delta #eta(jet_{1}, jet_{3})", 40, 0.0, 5.0);
  new TH1D("DEta_lightJet1Jet4", "#Delta #eta(jet_{1}, jet_{4})", 40, 0.0, 5.0);
  new TH1D("DEta_lightJet2Jet3", "#Delta #eta(jet_{2}, jet_{3})", 40, 0.0, 5.0);
  new TH1D("DEta_lightJet2Jet4", "#Delta #eta(jet_{2}, jet_{4})", 40, 0.0, 5.0);
  new TH1D("DEta_lightJet3Jet4", "#Delta #eta(jet_{3}, jet_{4})", 40, 0.0, 5.0);

  new TH1D("MT", "Transverse mass (GeV)", 40, 0.0, 200.0);  
  new TH1D("DR_leadLightJet_Chi", "DR leadLightJet_Chi", 40, 0., 5.0);
  new TH1D("DPhi_leadLightJet_Chi", "DPhi leadLightJet_Chi", 40, 0.0, 3.3);
  new TH1D("DR_minimum_allJets", "minimum #DeltaR_{all jets}", 40, 0, 5.0);
  new TH1D("DR_maximum_allJets", "maximum #DeltaR_{all jets}", 40, 0, 5.0);
  new TH1D("IM_minimum_lightJets", "minimum InvMass_{light jets}", 40, 0, 500.0);

  
  new TH1D("DR_mu1tauh1", "#DeltaR(#mu_{1}, #tauh_{1})", 40, 0.0, 5.0);
  new TH1D("Deltaphi_mu1tauh1", "#Delta #phi(#mu_{1}, #tauh_{1})", 40, 0.0, 3.3);
  new TH1D("DR_mu2tauh1", "#DeltaR(#mu_{2}, #tauh_{1})", 40, 0.0, 5.0);
  new TH1D("Deltaphi_mu2tauh1", "#Delta #phi(#mu_{2}, #tauh_{1})", 40, 0.0, 3.3);

  new TH1D ("MT_X", "m_{T} (GeV)", 100, 0, 200);
  new TH1D ("InvM_coln_muTauh_GS", "Collinear mass of #mu and #tau_{h} (GeV)", 100, 0, 200);
  new TH1D ("InvM_coln_lep1Tauh_GS", "Collinear mass of leading-lepton and #tau_{h} (GeV)", 100, 0, 200);
  new TH1D ("InvM_coln_lep2Tauh_GS", "Collinear mass of subleading-lepton and #tau_{h} (GeV)", 100, 0, 200);
  new TH1D ("InvM_coln_muTauh_IC", "Collinear mass of #mu and #tau_{h} (GeV)", 100, 0, 200);
  new TH1D ("InvM_coln_lep1Tauh_IC", "Collinear mass of leading-lepton and #tau_{h} (GeV)", 100, 0, 200);
  new TH1D ("InvM_coln_lep2Tauh_IC", "Collinear mass of subleading-lepton and #tau_{h} (GeV)", 100, 0, 200);
  
  histf->cd();
  histf->ls();
}


// ---------------------------------------
// Clear vectors before event loop
// ---------------------------------------


void ExoAnalysis::clearLists() {
}


// -----------------
// The main event loop
// -----------------

void ExoAnalysis::eventLoop(ExRootTreeReader *treeReader)
{
  /* ************************* Objects for analysis **************** */
  //TClonesArray *BrGen      = treeReader->UseBranch("Particle");
  TClonesArray *BrElec     = treeReader->UseBranch("Electron");
  TClonesArray *BrMuon     = treeReader->UseBranch("Muon");
  TClonesArray *BrJet      = treeReader->UseBranch("Jet");
  TClonesArray *BrPhoton   = treeReader->UseBranch("Photon");
  TClonesArray *BrMet      = treeReader->UseBranch("MissingET");
  //TClonesArray *BrEvent    = treeReader->UseBranch("Event");
  //TClonesArray *BrWeight   = treeReader->UseBranch("Weight");
  
  
  size_t nEntries = treeReader->GetEntries();
  cout << "** Chain contains " << nEntries << " events" << endl;
  
  
  size_t nEvents = (maxEvt_ < 0) ? nEntries : maxEvt_;
  cout << "** Staring Analysis with " << nEvents << " events" << endl;
  
  
  double lumiFac = lumiWt(nEvents, maxEvt_);
  cout << endl
       << "evtWeightSum: " << setw(10) << setprecision(0) << nEvents << endl
       << "      lumiWt: " << setw(10) << setprecision(5) << lumiFac
       << endl;
  
  
  /* ************************ Event Loop **************************************** */
  for (size_t iEntry = 0; iEntry < nEvents; ++iEntry) {
    bool verbose {false};
    if ((iEntry/100000) > 0 && iEntry%100000 == 0) verbose = true;
    if (verbose) cout << " Events processed: " << "\t" << iEntry << endl;
    
    
    clearLists(); // reset analysis related lists for each event
    treeReader->ReadEntry(iEntry); // Load selected branches with data from spacified event
    
    
    /* ************************ EventWeight and EventWeightSum *********** */
    double evWt = 1.0;
    double lumiWt = evWt*lumiFac;  // Lumi Scaling
    
    
    /* ****************************************************************** */     
    histf->cd();
    //histf->cd("Analysis");
    
    
    AnaUtil::fillHist1D("evtCutFlow", 0, 1.0);   // events processed
    AnaUtil::fillHist1D("evtCutFlowWt", 0, lumiWt);  // events processed
    
    
    int nelec   = BrElec->GetEntriesFast();
    AnaUtil::fillHist1D("nRawEle", nelec, evWt);
    int nmuon   = BrMuon->GetEntriesFast();
    AnaUtil::fillHist1D("nRawMuon", nmuon, evWt);
    int nphoton = BrPhoton->GetEntriesFast();
    AnaUtil::fillHist1D("nRawPhoton", nphoton, evWt);
    int njet    = BrJet->GetEntriesFast();
    AnaUtil::fillHist1D("nRawJet", njet, evWt);
    MissingET* met_ = (MissingET*) BrMet->At(0);
    AnaUtil::fillHist1D("met", (met_ -> MET), evWt);
    
    
    //Object containers
    std::vector <Muon>      MuonColl;
    std::vector <Electron>  ElectronColl;
    std::vector <Photon>    PhotonColl;
    std::vector <Jet>       TauhColl;
    std::vector <Jet>       JetColl;
    std::vector <Jet>       BJetColl;
    std::vector <Jet>       lightJetColl;
    
    
    histf->cd();
    //histf->cd("ObjectSelection");
    
    
    /* **************************** Muon Selection **************************** */     
    for (TObject *_mu: *BrMuon){
       const Muon& mu = dynamic_cast<const Muon&> (*_mu);
       AnaUtil::fillHist1D("MuonCutFlow", 0, evWt);       
       if (mu.PT <= 15) continue;
       AnaUtil::fillHist1D("MuonCutFlow", 1, evWt);
       if (std::abs(mu.Eta) >= 2.4) continue;
       AnaUtil::fillHist1D("MuonCutFlow", 2, evWt);
       TLorentzVector mup4 = mu.P4();
       MuonColl.push_back(mu);
     }
    std::sort(std::begin(MuonColl), std::end(MuonColl), PtComparator<Muon>()); 
     
         
     /* ********************************** Electron Selection **************************** */     
     for (TObject *_el: *BrElec){
       const Electron& ele = dynamic_cast<const Electron&> (*_el);
       AnaUtil::fillHist1D("ElectronCutFlow", 0, evWt);
       
       if (ele.PT <= 20) continue;
       AnaUtil::fillHist1D("ElectronCutFlow", 1, evWt);
       if (std::abs(ele.Eta) >= 2.5) continue;
       AnaUtil::fillHist1D("ElectronCutFlow", 2, evWt);
       ElectronColl.push_back(ele);
     }
     std::sort(std::begin(ElectronColl), std::end(ElectronColl), PtComparator<Electron>());
     
     /* ********************************** Photon Selection *********************************** */
     for (TObject *_ph: *BrPhoton){
       const Photon& ph = dynamic_cast<const Photon&> (*_ph);
       AnaUtil::fillHist1D("PhotonCutFlow", 0, evWt);
       
       if (ph.PT <= 20) continue;
       AnaUtil::fillHist1D("PhotonCutFlow", 1, evWt);
       if (std::abs(ph.Eta) >= 2.4) continue;
       AnaUtil::fillHist1D("PhotonCutFlow", 2, evWt);
     }
     std::sort(std::begin(PhotonColl), std::end(PhotonColl), PtComparator<Photon>());
          
     /* ********************************** Jet Selection ************************************ */
     for (TObject *_j: *BrJet){
       const Jet& jet = dynamic_cast<const Jet&> (*_j);
       AnaUtil::fillHist1D("JetCutFlow", 0, evWt);

       if (jet.PT <= 30) continue;
       AnaUtil::fillHist1D("JetCutFlow", 1, evWt);
       if (std::abs(jet.Eta) >= 2.5) continue;
       AnaUtil::fillHist1D("JetCutFlow", 2, evWt);
       if (!(jetLeptonCleaning(jet, MuonColl, ElectronColl, 0.4))) continue;
       AnaUtil::fillHist1D("JetCutFlow", 3, evWt);
       if (!(jetTauhCleaning(jet, TauhColl, 0.4))) continue;
       AnaUtil::fillHist1D("JetCutFlow", 4, evWt);
       JetColl.push_back(jet);
     }
     std::sort(std::begin(JetColl), std::end(JetColl), PtComparator<Jet>());
     
     /* ********************************** Tau Selection ************************************ */
     //int iTauh = 0;
     for (TObject *_tau: *BrJet){
       const Jet& tauj = dynamic_cast<const Jet&> (*_tau);
       if (tauj.TauTag == 0) continue;
       AnaUtil::fillHist1D("TauhCutFlow", 0, evWt);
       if (tauj.PT <= 20) continue;
       AnaUtil::fillHist1D("TauhCutFlow", 1, evWt);
       if (std::abs(tauj.Eta) >= 2.4) continue;
       AnaUtil::fillHist1D("TauhCutFlow", 2, evWt);
       TauhColl.push_back(tauj);
     }
     std::sort(std::begin(TauhColl), std::end(TauhColl), PtComparator<Jet>());
     
     /* ************************************************************************************ */
     for (auto &jet : JetColl)
       (jet.BTag == 1 && std::fabs(jet.Eta) < 2.5) ? BJetColl.push_back(jet) : lightJetColl.push_back(jet);

     int nGoodMuon     = MuonColl.size();
     int nGoodEle      = ElectronColl.size();
     int nGoodPhoton   = PhotonColl.size();
     int nGoodTauh     = TauhColl.size();
     int nGoodJet      = JetColl.size();
     int nGoodlJet     = lightJetColl.size();
     int nGoodBJet     = BJetColl.size();

     // -------------- Initial Cuts ---------------- //
     if (nGoodMuon + nGoodEle + nGoodTauh + nGoodJet == 0) continue;
     AnaUtil::fillHist1D("evtCutFlow", 1, evWt);
     AnaUtil::fillHist1D("evtCutFlowWt", 1, lumiWt);
     
     if (nGoodMuon + nGoodEle != 2) continue;
     AnaUtil::fillHist1D("evtCutFlow", 2, evWt);
     AnaUtil::fillHist1D("evtCutFlowWt", 2, lumiWt);

     if (nGoodMuon < 1) continue;
     AnaUtil::fillHist1D("evtCutFlow", 3, evWt);
     AnaUtil::fillHist1D("evtCutFlowWt", 3, lumiWt);
          
     if (nGoodTauh != 1) continue;
     AnaUtil::fillHist1D("evtCutFlow", 4, evWt);
     AnaUtil::fillHist1D("evtCutFlowWt", 4, lumiWt);     

     bool has2ndBjet = (nGoodBJet > 1 && BJetColl[1].PT > 30) ? true : false;
     
     if (nGoodlJet < 1) continue;
     AnaUtil::fillHist1D("evtCutFlow", 5, evWt);
     AnaUtil::fillHist1D("evtCutFlowWt", 5, lumiWt);
     
     if (nGoodBJet < 1 || has2ndBjet) continue;
     AnaUtil::fillHist1D("evtCutFlow", 6, evWt);
     AnaUtil::fillHist1D("evtCutFlowWt", 6, lumiWt);

     std::vector<ZCandidate> ZCandList;
     if (nGoodMuon >= 2) ZSelector(MuonColl, ZCandList);
     if (nGoodEle  >= 2) ZSelector(ElectronColl, ZCandList);
     if (ZCandList.size() > 0 && ZCandList[0].massDiff < 15) continue;
     AnaUtil::fillHist1D("evtCutFlow", 7, evWt);
     AnaUtil::fillHist1D("evtCutFlowWt", 7, lumiWt);
     // ---------------------------------------------- //
     
     AnaUtil::fillHist1D("nGoodMuon", nGoodMuon);
     AnaUtil::fillHist1D("nGoodElectron", nGoodEle);
     AnaUtil::fillHist1D("nGoodLepton", (nGoodMuon + nGoodEle)); 
     AnaUtil::fillHist1D("nGoodPhoton", nGoodPhoton);
     AnaUtil::fillHist1D("nGoodTauh", nGoodTauh);
     AnaUtil::fillHist1D("nGoodJet", nGoodJet);
     AnaUtil::fillHist1D("nGoodBJet", nGoodBJet); 
     AnaUtil::fillHist1D("nGoodlJet", nGoodlJet);
     
     AnaUtil::fillHist1D("Tauh1Pt", TauhColl[0].PT, 1);
     AnaUtil::fillHist1D("Tauh1Eta", TauhColl[0].Eta, 1);
     AnaUtil::fillHist1D("Tauh1Phi", TauhColl[0].Phi, 1);

     // ---- Leading and sub-leading lepton pT ------ //
     TLorentzVector hLepP4 = MuonColl[0].P4();
     bool ismu = true;
     bool isele = false;
     if (nGoodEle > 0 && ElectronColl[0].PT > hLepP4.Pt()){
       isele   = true;
       ismu    = false;
       hLepP4  = ElectronColl[0].P4();
     }

     TLorentzVector h2LepP4;
     h2LepP4.SetPtEtaPhiM(0.,0.,0.,0.);
     if (isele) {
       if (nGoodEle > 1) h2LepP4 = (ElectronColl[1].PT > MuonColl[0].PT) ? ElectronColl[1].P4() : MuonColl[0].P4();
       else h2LepP4 = MuonColl[0].P4();
     }
     else if (ismu) {
       if (nGoodEle == 0) h2LepP4 = MuonColl[1].P4();
       else if (nGoodMuon == 1) h2LepP4 = ElectronColl[0].P4();
       else if (nGoodMuon > 1 && nGoodEle > 0) h2LepP4 = (MuonColl[1].PT > ElectronColl[0].PT) ? MuonColl[1].P4() : ElectronColl[0].P4();
     }
     // ----------------------------------------------- //

     AnaUtil::fillHist1D("lep1Pt", hLepP4.Pt(), 1);
     AnaUtil::fillHist1D("lep1Eta", hLepP4.Eta(), 1);
     AnaUtil::fillHist1D("lep2Pt", h2LepP4.Pt(), 1);
     AnaUtil::fillHist1D("lep2Eta", h2LepP4.Eta(), 1);
          
     size_t nLoop2Ele  = (nGoodEle >= 2)  ? 2 : nGoodEle;
     size_t nLoop2Mu   = (nGoodMuon >= 2) ? 2 : nGoodMuon;
     size_t nLoop4Jet  = (nGoodJet >= 4)  ? 4 : nGoodJet;
     size_t nLoop4lJet = (nGoodlJet >= 4) ? 4 : nGoodlJet;
     size_t nLoop6Jet  = (nGoodJet >= 6)  ? 6 : nGoodJet;
     size_t nLoopBJet  = (nGoodBJet >= 2) ? 2 : nGoodBJet;
     size_t nLoopTauh  = (nGoodTauh > 1)  ? 1 : nGoodTauh;
     
     // b-jet plots (for max 2 b-jets)
     for (size_t i = 0; i < nLoopBJet; ++i) {
       TLorentzVector bjP4 = BJetColl[i].P4();
       std::string pt_name  = "BJet"+std::to_string(i+1)+"Pt";
       std::string phi_name = "BJet"+std::to_string(i+1)+"Phi";
       std::string eta_name = "BJet"+std::to_string(i+1)+"Eta";
       AnaUtil::fillHist1D (pt_name.c_str(), bjP4.Pt(), 1);
       AnaUtil::fillHist1D (phi_name.c_str(), bjP4.Phi(), 1);
       AnaUtil::fillHist1D (eta_name.c_str(), bjP4.Eta(), 1);
     }

     // Loop over max 2 electrons
     for (size_t iele = 0; iele < nLoop2Ele; ++iele) {
       TLorentzVector eleP4 = ElectronColl[iele].P4();
       std::string elePt_name  = "Electron"+std::to_string(iele+1)+"Pt";
       std::string eleEta_name = "Electron"+std::to_string(iele+1)+"Eta";
       std::string elePhi_name = "Electron"+std::to_string(iele+1)+"Phi";
       AnaUtil::fillHist1D (elePt_name.c_str(), eleP4.Pt(), 1);
       AnaUtil::fillHist1D (eleEta_name.c_str(), eleP4.Eta(), 1);
       AnaUtil::fillHist1D (elePhi_name.c_str(), eleP4.Phi(), 1);
       for (size_t ijet = 0; ijet < nLoop4Jet; ++ijet) {
	 TLorentzVector jP4 = JetColl[ijet].P4();
	 std::string DR_jetEle_name   = "DR_Jet"+std::to_string(ijet+1)+"Electron"+std::to_string(iele+1);
	 std::string DPhi_jetEle_name = "DPhi_Jet"+std::to_string(ijet+1)+"Electron"+std::to_string(iele+1);
	 AnaUtil::fillHist1D (DR_jetEle_name.c_str(), jP4.DeltaR(eleP4), 1);
	 AnaUtil::fillHist1D (DPhi_jetEle_name.c_str(), std::fabs(TVector2::Phi_mpi_pi(jP4.Phi()-eleP4.Phi())), 1);
       }
       for (size_t ib = 0; ib < nLoopBJet; ++ib) {
	 TLorentzVector bjP4 = BJetColl[ib].P4();
	 std::string DR_bjetEle_name   = "DR_bJet"+std::to_string(ib+1)+"Electron"+std::to_string(iele+1);
	 std::string DPhi_bjetEle_name = "DPhi_bJet"+std::to_string(ib+1)+"Electron"+std::to_string(iele+1);
	 AnaUtil::fillHist1D (DR_bjetEle_name.c_str(), bjP4.DeltaR(eleP4), 1);
	 AnaUtil::fillHist1D (DPhi_bjetEle_name.c_str(), std::fabs(TVector2::Phi_mpi_pi(bjP4.Phi() - eleP4.Phi())), 1);	 
       }
     }

     // Loop over max 2 muons
     for (size_t imu = 0; imu < nLoop2Mu; ++imu) {
       TLorentzVector muP4 = MuonColl[imu].P4();
       std::string muPt_name  = "Muon"+std::to_string(imu+1)+"Pt";
       std::string muEta_name = "Muon"+std::to_string(imu+1)+"Eta";
       std::string muPhi_name = "Muon"+std::to_string(imu+1)+"Phi";
       AnaUtil::fillHist1D (muPt_name.c_str(), muP4.Pt(), 1);
       AnaUtil::fillHist1D (muEta_name.c_str(), muP4.Eta(), 1);
       AnaUtil::fillHist1D (muPhi_name.c_str(), muP4.Phi(), 1);
       for (size_t iTauh = 0; iTauh < nLoopTauh; ++iTauh) {
         TLorentzVector tauhP4 = TauhColl[iTauh].P4();
         std::string DR_tauhMu_name = "DR_mu"+std::to_string(imu+1)+"tauh"+std::to_string(iTauh+1);
	 std::string DPhi_tauhMu_name = "Deltaphi_mu"+std::to_string(imu+1)+"tauh"+std::to_string(iTauh+1);
	 AnaUtil::fillHist1D (DR_tauhMu_name.c_str(), tauhP4.DeltaR(muP4), 1);
	 AnaUtil::fillHist1D (DPhi_tauhMu_name.c_str(), std::fabs(TVector2::Phi_mpi_pi(tauhP4.Phi() - muP4.Phi())), 1);
       }
       for (size_t ijet = 0; ijet < nLoop4lJet; ++ijet) {
	 TLorentzVector jP4 = lightJetColl[ijet].P4();
	 std::string DR_jetMu_name = "DR_Jet"+std::to_string(ijet+1)+"Muon"+std::to_string(imu+1);
	 std::string DPhi_jetMu_name = "DPhi_Jet"+std::to_string(ijet+1)+"Muon"+std::to_string(imu+1);
	 AnaUtil::fillHist1D (DR_jetMu_name.c_str(), jP4.DeltaR(muP4), 1);
	 AnaUtil::fillHist1D (DPhi_jetMu_name.c_str(), std::fabs(TVector2::Phi_mpi_pi(jP4.Phi() - muP4.Phi())), 1);
       }
       for (size_t ib = 0; ib < nLoopBJet; ++ib) {
	 TLorentzVector bjP4 = BJetColl[ib].P4();
	 std::string DR_bjetMu_name   = "DR_bJet"+std::to_string(ib+1)+"Muon"+std::to_string(imu+1);
	 std::string DPhi_bjetMu_name = "DPhi_bJet"+std::to_string(ib+1)+"Muon"+std::to_string(imu+1);
	 AnaUtil::fillHist1D (DR_bjetMu_name.c_str(), bjP4.DeltaR(muP4), 1);
	 AnaUtil::fillHist1D (DPhi_bjetMu_name.c_str(), std::fabs(TVector2::Phi_mpi_pi(bjP4.Phi() - muP4.Phi())), 1);	 
       }
     }
     
     // Loop over max 6 good jets
     float minDR = 999.9;
     float maxDR = -999.9;
     for (size_t i = 0; i < nLoop6Jet; ++i) {
       std::string pt_name  = "Jet"+std::to_string(i+1)+"Pt";
       TLorentzVector jP4   = JetColl[i].P4();
       AnaUtil::fillHist1D (pt_name.c_str(), JetColl[i].PT, 1);
       for (size_t j = i+1; j < nLoop6Jet; ++j) {
	 TLorentzVector j2P4   = JetColl[j].P4();
	 float dR = jP4.DeltaR(j2P4);
	 if (dR < minDR) minDR = dR;
	 if (dR > maxDR) maxDR = dR;
       }
     }
     AnaUtil::fillHist1D("DR_minimum_allJets",minDR,1.0);
     AnaUtil::fillHist1D("DR_maximum_allJets",maxDR,1.0);
     // Loop over max 4 light jets
     float min_invM = 9999.9;
     for (size_t i = 0; i < nLoop4lJet; ++i) {
       std::string jetPtName  = "lightJet"+std::to_string(i+1)+"Pt";
       std::string jetEtaName = "lightJet"+std::to_string(i+1)+"Eta";
       TLorentzVector ljetiP4 = lightJetColl[i].P4();
       AnaUtil::fillHist1D (jetPtName.c_str(), ljetiP4.Pt(), 1);
       AnaUtil::fillHist1D (jetEtaName.c_str(), ljetiP4.Eta(), 1);
       for (size_t j = i+1; j < nLoop4lJet; ++j) {
	 TLorentzVector ljetjP4 =lightJetColl[j].P4();
	 std::string PT_histName   = "PT_lightJet"+std::to_string(i+1)+"Jet"+std::to_string(j+1);
	 std::string IM_histName   = "IM_lightJet"+std::to_string(i+1)+"Jet"+std::to_string(j+1);
	 std::string DR_histName   = "DR_lightJet"+std::to_string(i+1)+"Jet"+std::to_string(j+1);
	 std::string DEta_histName = "DEta_lightJet"+std::to_string(i+1)+"Jet"+std::to_string(j+1);
	 std::string DPhi_histName = "DPhi_lightJet"+std::to_string(i+1)+"Jet"+std::to_string(j+1);
	 float invM = (ljetiP4+ljetjP4).M();
	 AnaUtil::fillHist1D (PT_histName.c_str(), (ljetiP4+ljetjP4).Pt(), 1);
	 AnaUtil::fillHist1D (IM_histName.c_str(), invM, 1);
	 AnaUtil::fillHist1D (DR_histName.c_str(), ljetiP4.DeltaR(ljetjP4), 1);
	 AnaUtil::fillHist1D (DEta_histName.c_str(), std::fabs(ljetiP4.Eta() - ljetjP4.Eta()), 1);
	 AnaUtil::fillHist1D (DPhi_histName.c_str(), std::fabs(TVector2::Phi_mpi_pi(ljetiP4.Phi() - ljetjP4.Phi())), 1);
	 if (invM < min_invM) min_invM = invM;
       }
     }

     AnaUtil::fillHist1D("IM_minimum_lightJets",min_invM, 1.0);
     // -------------------- NEW --------------------- //
     //mT of leading lepton and met
     TLorentzVector metp4;
     metp4.SetPtEtaPhiE(met_->MET, 0.0, met_->Phi, met_->MET);

     float mT = Calculate_MT(hLepP4, metp4);
     AnaUtil::fillHist1D("met", metp4.Pt(), 1);
     AnaUtil::fillHist1D("MT", mT, 1);
     // ---------------------------------------------- //
     double XMT = Calculate_TotalMT(MuonColl[0].P4(), TauhColl[0].P4(), metp4);
     AnaUtil::fillHist1D ("MT_X", XMT);

     // Collinear mass
     // GS
     double corr = 1.0/(1 + metp4.Px()/TauhColl[0].P4().Px());
     TLorentzVector Ntap4;
     Ntap4.SetPtEtaPhiE((1.0/corr)*TauhColl[0].P4().Pt(), TauhColl[0].P4().Eta(), TauhColl[0].P4().Phi(), TauhColl[0].P4().E());
     double colnMass = (Ntap4 + MuonColl[0].P4()).M();
     double colnMass_lep1 = (Ntap4 + hLepP4).M();
     double colnMass_lep2 = (Ntap4 + h2LepP4).M();

     AnaUtil::fillHist1D("InvM_coln_muTauh_GS", colnMass);
     AnaUtil::fillHist1D("InvM_coln_lep1Tauh_GS", colnMass_lep1);
     AnaUtil::fillHist1D("InvM_coln_lep2Tauh_GS", colnMass_lep2);

     // IC
     double IM_mutau   = (MuonColl[0].P4()+TauhColl[0].P4()).M();
     double IM_lep1tau = (hLepP4 + TauhColl[0].P4()).M();
     double IM_lep2tau = (h2LepP4 + TauhColl[0].P4()).M();

     double pTvis_tau =  TauhColl[0].PT;
     double pTnu = (metp4.Px()*TauhColl[0].P4().Px() + metp4.Py()*TauhColl[0].P4().Py())/std::fabs(pTvis_tau);                         
     double x_vis = std::fabs(pTvis_tau)/(std::fabs(pTvis_tau) + std::fabs(pTnu));                  
     double M_coll_IC = IM_mutau/TMath::Sqrt(x_vis);  
     double M_coll_IC_lep1 = IM_lep1tau/TMath::Sqrt(x_vis);  
     double M_coll_IC_lep2 = IM_lep2tau/TMath::Sqrt(x_vis);  

     AnaUtil::fillHist1D("InvM_coln_muTauh_IC", M_coll_IC);
     AnaUtil::fillHist1D("InvM_coln_lep1Tauh_IC", M_coll_IC_lep1);
     AnaUtil::fillHist1D("InvM_coln_lep2Tauh_IC", M_coll_IC_lep2);


     // Filling the branches of a flat ntuple for skimming
     if (skimObj_) {
       TreeVariables varList;

       varList.event         = iEntry;
       // numbers of objects in final state
       varList.nLeptons      = nGoodMuon+nGoodEle;
       varList.nJets         = nGoodJet;
       varList.nbJets        = nGoodBJet;
       varList.nlJets        = nGoodlJet;
       varList.nTauh         = nGoodTauh;
       // properties of individual object
       varList.pt_tauh1      = TauhColl[0].PT;
       varList.eta_tauh1     = std::fabs(TauhColl[0].Eta);
       varList.phi_tauh1     = std::fabs(TauhColl[0].Phi);
       varList.met           = met_->MET;
       varList.pt_bjet1      = BJetColl[0].PT;
       varList.eta_bjet1     = std::fabs(BJetColl[0].Eta);
       varList.phi_bjet1     = std::fabs(BJetColl[0].Phi);
       varList.pt_ljet1      = lightJetColl[0].PT;
       varList.eta_ljet1     = std::fabs(lightJetColl[0].Eta);
       varList.phi_ljet1     = std::fabs(lightJetColl[0].Phi);
       varList.pt_lep1       = hLepP4.Pt();
       varList.pt_lep2       = h2LepP4.Pt();
       varList.eta_lep1      = std::fabs(hLepP4.Eta());
       varList.eta_lep2      = std::fabs(h2LepP4.Eta());
       varList.phi_lep1      = std::fabs(hLepP4.Phi());
       varList.phi_lep2      = std::fabs(h2LepP4.Phi());
       // di-lepton variables
       varList.dr_lep1lep2   = hLepP4.DeltaR(h2LepP4);
       varList.dphi_lep1lep2 = std::fabs(TVector2::Phi_mpi_pi(hLepP4.Phi() - h2LepP4.Phi()));
       varList.deta_lep1lep2 = std::fabs(hLepP4.Eta() - h2LepP4.Eta());
       // lepton-tauh variables
       varList.dr_lep1tauh   = hLepP4.DeltaR(TauhColl[0].P4());
       varList.dphi_lep1tauh = std::fabs(TVector2::Phi_mpi_pi(hLepP4.Phi() - TauhColl[0].Phi));
       varList.deta_lep1tauh = std::fabs(hLepP4.Eta() - TauhColl[0].Eta);
       varList.dr_lep2tauh   = h2LepP4.DeltaR(TauhColl[0].P4());
       varList.dphi_lep2tauh = std::fabs(TVector2::Phi_mpi_pi(h2LepP4.Phi() - TauhColl[0].Phi));
       varList.deta_lep2tauh = std::fabs(h2LepP4.Eta() - TauhColl[0].Eta);
       // tauh-jet variables
       varList.dr_tauhjet   = (TauhColl[0].P4()).DeltaR(JetColl[0].P4());
       varList.dphi_tauhjet = std::fabs(TVector2::Phi_mpi_pi(TauhColl[0].Phi - JetColl[0].Phi));
       varList.deta_tauhjet = std::fabs(JetColl[0].Eta - TauhColl[0].Eta);
       varList.dr_tauhljet   = (TauhColl[0].P4()).DeltaR(lightJetColl[0].P4());
       varList.dphi_tauhljet = std::fabs(TVector2::Phi_mpi_pi(TauhColl[0].Phi - lightJetColl[0].Phi));
       varList.deta_tauhljet = std::fabs(lightJetColl[0].Eta - TauhColl[0].Eta);
       varList.dr_tauhbjet   = (TauhColl[0].P4()).DeltaR(BJetColl[0].P4());
       varList.dphi_tauhbjet = std::fabs(TVector2::Phi_mpi_pi(TauhColl[0].Phi - BJetColl[0].Phi));
       varList.deta_tauhbjet = std::fabs(BJetColl[0].Eta - TauhColl[0].Eta);
       // lepton-jet variables
       // -- lep - bJets
       varList.dr_lep1bjet   = hLepP4.DeltaR(BJetColl[0].P4());
       varList.dr_lep2bjet   = h2LepP4.DeltaR(BJetColl[0].P4());
       varList.dphi_lep1bjet = std::fabs(TVector2::Phi_mpi_pi(hLepP4.Phi() - BJetColl[0].Phi));
       varList.dphi_lep2bjet = std::fabs(TVector2::Phi_mpi_pi(h2LepP4.Phi() - BJetColl[0].Phi));
       varList.deta_lep1bjet = std::fabs(BJetColl[0].Eta - hLepP4.Eta());
       varList.deta_lep2bjet = std::fabs(BJetColl[0].Eta - h2LepP4.Eta());
       // -- lep - lightJets
       varList.dr_lep1ljet   = hLepP4.DeltaR(lightJetColl[0].P4());
       varList.dr_lep2ljet   = h2LepP4.DeltaR(lightJetColl[0].P4());
       varList.dphi_lep1ljet = std::fabs(TVector2::Phi_mpi_pi(hLepP4.Phi() - lightJetColl[0].Phi));
       varList.dphi_lep2ljet = std::fabs(TVector2::Phi_mpi_pi(h2LepP4.Phi() - lightJetColl[0].Phi));
       varList.deta_lep1ljet = std::fabs(lightJetColl[0].Eta - hLepP4.Eta());
       varList.deta_lep2ljet = std::fabs(lightJetColl[0].Eta - h2LepP4.Eta());
       // lepton - met variables
       varList.dphi_lep1met  = std::fabs(TVector2::Phi_mpi_pi(hLepP4.Phi() - metp4.Phi()));
       varList.dphi_lep2met  = std::fabs(TVector2::Phi_mpi_pi(h2LepP4.Phi() - metp4.Phi()));
       varList.mt_lep1met    = Calculate_MT(hLepP4, metp4);
       varList.mt_lep2met    = Calculate_MT(h2LepP4, metp4);
       varList.mt_muon1met   = Calculate_MT(MuonColl[0].P4(), metp4);
       // di-jets variables
       varList.dr_min_jets   = minDR;
       varList.dr_max_jets   = maxDR;
       varList.dr_bjetljet   = (lightJetColl[0].P4()).DeltaR(BJetColl[0].P4());
       varList.dphi_bjetljet = std::fabs(TVector2::Phi_mpi_pi(BJetColl[0].Phi - lightJetColl[0].Phi));
       varList.deta_bjetljet = std::fabs(BJetColl[0].Eta - lightJetColl[0].Eta);
       
       skimObj_->fill(varList);
     }       
     histf->cd();
  } 
}


void ExoAnalysis::endJob() {
  
  histf->cd();
  //histf->cd("Analysis");
  vector<string> evLabels {
    "Events processed",
      "pass obj selection",
      "nMu + nEle == 2",
      "nMu >= 1",
      "nTau >= 1",
      "NLightJet >= 1",
      "nBJets == 1",
      "no Z"
      };
  AnaUtil::SetEvtCutFlowBinLabels("evtCutFlow",evLabels);
  AnaUtil::showEfficiency("evtCutFlow", evLabels, "Event Selection");
  AnaUtil::showEfficiency("evtCutFlowWt", evLabels, "Event Selection : Lumi Scaled");  
}

void ExoAnalysis::closeHistFiles(){
  histf->cd();
  histf->Write();
  histf->Close();
}
void ExoAnalysis::closeFiles(){
  closeHistFiles();
  if (skimObj_ != nullptr) skimObj_->close();
  delete treeReader;
  delete chain;
}
bool ExoAnalysis::readJob(const string& jobFile, int& nFiles)
{
  // Open the file containing the datacards
  ifstream fin(jobFile.c_str(), ios::in);    
  if (!fin) {
    cerr << "==> Input File: <<" << jobFile << ">> could not be opened!" << endl;
    return false;
  }

  static constexpr int BUF_SIZE = 256;
  char buf[BUF_SIZE];
  while (fin.getline(buf, BUF_SIZE, '\n')) {  // Pops off the newline character
    string line(buf);
    if (line.empty() || line == "START") continue;   
    
    // enable '#' and '//' style comments
    if (line.substr(0,1) == "#" || line.substr(0,2) == "//") continue;
    if (line == "END") break;
    
    // Split the line into words
    vector<string> tokens;
    AnaUtil::tokenize(line, tokens);
    int vsize = tokens.size();
    assert(vsize > 1);
    
    const string& key   = tokens.at(0);
    const string& value = tokens.at(1);
    if (key == "dataType") {
      string vtmp(value);
      std::transform(vtmp.begin(), vtmp.end(), vtmp.begin(), ::toupper);
      vector<string> dt;
      AnaUtil::tokenize(vtmp, dt, "#");
      if (dt.size()) {
        isMC_ = (dt.at(0) == "MC") ? true : false;
        if (isMC_ && dt.size() > 1) {
          isSignal_ = (dt.at(1) == "SIGNAL") ? true : false;
        }
      }
    }
    else if (key == "readGenInfo") 
      readGenInfo_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "inputFile")
      AnaUtil::buildList(tokens, fileList_);
    else if (key == "maxEvent") 
      maxEvt_ = std::stoi(value.c_str());
    else if (key == "histFile")
      histFile_ = value;
    else if (key == "readMVA")
      _readMVA = (atoi(value.c_str()) > 0) ? true : false;
    else if (key == "MVAnetwork")
      _MVAnetwork = value;
    else if (key == "MVAxmlFile")
      _MVAxmlFile = value;
   // else if (key == "MVAxmlFile2")
     // _MVAxmlFile2 = value;
    else if (key == "createMVATree")
      _createMVATree = (atoi(value.c_str()) > 0) ? true : false;
    else if (key == "mvaInputFile")
      _mvaInputFile = value;
    //else if (key == "mvaOutputFile")
    //  _mvaOutputFile = value;
   // else if (key == "isXZ")
   //   isXZ_ = (atoi(value.c_str()) > 0) ? true : false;
   // else if (key == "isH2Z")
     // isH2Z_ = (atoi(value.c_str()) > 0) ? true : false;
  //  else if (key == "CutOnMVAScorePlane")
    //  _CutOnMVAScorePlane = (atoi(value.c_str()) > 0) ? true : false;    
    else if (key == "eventId" && tokens.size() == 4)
      AnaUtil::buildMap(tokens, eventIdMap_);
    else {
      if (0) cout << "==> " << line << endl;
    AnaUtil::storeCuts(tokens, hmap_);
    }						
  }
  // Close the file
  fin.close();
  //  if (!isSignal_) readGenInfo_ = false;
  chain = new TChain("Delphes");
  // Build the chain of root files
  for (const auto& fname: fileList_) {
    cout << ">>> INFO. Adding input file " << fname << " to TChain " << endl;
    ++nFiles;
    chain->Add(fname.c_str());
    int nevt = static_cast<int>(chain->GetEntries());
    if (maxEvt_ > 0 && nevt >= maxEvt_) break;
  }
  if (!nFiles) {
    cerr << ">>> WARN. Input Root file list is empty! exiting ..." << endl;
    return false;
  } 
  //  cout<<nFiles<<endl;
  treeReader = new ExRootTreeReader(chain);
  if(!treeReader) return false;

  printJob();  
  return true;
}

void ExoAnalysis::printJob(ostream& os) const{
  os << endl;
  AnaUtil::showCuts(hmap_, os);
}

//******************************************Some More Functions (Make AnaBase and move them)***********************************************//

double ExoAnalysis::lumiWt(double evtWeightSum, int maxEvents) const{
  //  double nevt = (evtWeightSum > -1) ? evtWeightSum : AnaUtil::cutValue(lumiWtMap(), "nevents");
  double nevt = (evtWeightSum > -1) ? evtWeightSum : maxEvents;
  std::cout << "-- intLumi: " << AnaUtil::cutValue(lumiWtMap(), "intLumi")
	    <<setprecision(5)
	    << " xsec: " << AnaUtil::cutValue(lumiWtMap(), "xsec")
	    <<setprecision(5)
	    << " nevt: " << nevt << std::endl;
  return (AnaUtil::cutValue(lumiWtMap(), "intLumi") * AnaUtil::cutValue(lumiWtMap(), "xsec") / nevt);
}

int ExoAnalysis::getMotherId(const GenParticle& gp, TClonesArray *gen, int& mmid) const {
  int pdgid = gp.PID;
  int indx = gp.M1;
  if (indx < 0) return -1;
  GenParticle& mgp = dynamic_cast<GenParticle&> (*(gen->At(indx)));
  mmid = mgp.PID;
  if (std::abs(mmid) == 2212) return -1;
  while (mmid == pdgid) {
    indx = mgp.M1;
    if (indx < 0) continue;
    mgp  = dynamic_cast<GenParticle&>(*(gen->At(indx)));
    mmid = mgp.PID;
  }
  return indx;
}

bool ExoAnalysis::jetLeptonCleaning(const Jet& jet, std::vector <Muon> muList_, std::vector <Electron> eleList_, double dR) {
  TLorentzVector jetP4 = jet.P4();
  for (auto& m: muList_){
    AnaUtil::fillHist1D ("jetMuDR", jetP4.DeltaR(m.P4()), 1.0);
    if (jetP4.DeltaR(m.P4()) <= dR) return false;
  }
  for (auto& e: eleList_){
    AnaUtil::fillHist1D ("jetEleDR", jetP4.DeltaR(e.P4()), 1.0);
    if (jetP4.DeltaR(e.P4()) <= dR) return false;
  }
  return true;
}
      

bool ExoAnalysis::jetTauhCleaning(const Jet& jet, std::vector <Jet> TauhList_, double dR) {
   TLorentzVector jetP4 = jet.P4();
   for (auto& Th: TauhList_){
     if (Th.TauTag == 0) continue;
     AnaUtil::fillHist1D ("jetThDR", jetP4.DeltaR(Th.P4()), 1.0);
     if (jetP4.DeltaR(Th.P4()) <= dR) return false;
  }
  return true;
}

void ExoAnalysis::dumpGenInfo(TClonesArray *gen, ostream& os) {
  int ngen = gen->GetEntries();
  if (ngen == 0) return;

  os << setprecision(2);
  os << " -- # GenParticle: " << ngen << endl;
  os << "indx    status    pdgId     eta      phi      pt     energy  moIndx"
     << "      moID                   daughterID"
     << endl;
  int indx = 0;
  for (TObject* obj: *gen) {
    //const GenParticle& gp = *(dynamic_cast<GenParticle*> (obj));                                                                                                      
    const GenParticle& gp = dynamic_cast<const GenParticle&> (*obj);
    //int abs_pid = std::abs(gp.PID);                                                                                                                                   
    std::ostringstream mID;
    vector<int> m;
    m.push_back(gp.M1);
    m.push_back(gp.M2);
    for (int mi: m) {
      if (mi < 0 || mi >= ngen) continue;
      const GenParticle& mgp = dynamic_cast<const GenParticle&>(*(gen->At(mi)));
      mID << " " << mgp.PID;
    }
    string ms = mID.str();
    if (!ms.length()) ms = " -";

    std::ostringstream dID;
    vector<int> d;
    d.push_back(gp.D1);
    d.push_back(gp.D2);
    for (int di: d) {
      if (di < 0 || di >= ngen) continue;
      const GenParticle& dgp = dynamic_cast<const GenParticle&>(*(gen->At(di)));
      //  const GenParticle& dgp = (GenParticle*)gen->At(di);                                                                                                           
      double energy = dgp.E;
      int pdgid = dgp.PID;
      if (std::abs(pdgid) == 21 && energy <= 10) continue;
      dID << " " << dgp.PID;
    }
    string ds = dID.str();
    if (!ds.length()) ds = " -";
    os << setw(4)  << indx++
       << setw(8)  << gp.Status
       << setw(10) << gp.PID
       << setw(10) << gp.Eta
       << setw(9)  << gp.Phi
       << setw(9)  << gp.PT
       << setw(9)  << gp.E
       << setw(8)  << gp.M1
       << setw(10) << ms
       << setw(28) << ds
       << endl;
  }
}
//******************************************Some More Functions (Make AnaBase and move these)***********************************************//


































