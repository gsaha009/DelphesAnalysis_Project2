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
  cout << "********** Booking Histograms **************" << endl;
  bookHistograms();
  
  if(_createMVATree) {
#ifdef TRUE_CPP14
    skimObj_ = std::make_unique <MVASkim> (_mvaInputFile, isSL(), isDL());
#else
    skimObj_ = std::unique_ptr <MVASkim>(new MVASkim (_mvaInputFile, isSL(), isDL()));
#endif
    if (!skimObj_) return false;
  }
  
  else if (_readMVA) {
#ifdef TRUE_CPP14
    _mvaObj = std::make_unique<MVAnalysis>(_MVAnetwork, _MVAxmlFile, isSL(), isDL());
#else
    _mvaObj = std::unique_ptr <MVAnalysis>(new MVAnalysis (_MVAnetwork, _MVAxmlFile, isSL(), isDL()));
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
  HistBooker hbooker;
  histf->cd();
  hbooker.bookBasicHistograms();
  if (isDL())       hbooker.bookDLHistograms();
  else if (isSL())  hbooker.bookSLHistograms();
  else std::cerr << "No channel is mentioned " 
		 << std::endl;
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
  //TClonesArray *BrWeight   = treeReader->UseBranch("Weight");
  TClonesArray *BrEvent    = treeReader->UseBranch("Event");
  TClonesArray *BrElec     = treeReader->UseBranch("Electron");
  TClonesArray *BrMuon     = treeReader->UseBranch("Muon");
  TClonesArray *BrJet      = treeReader->UseBranch("Jet");
  TClonesArray *BrPhoton   = treeReader->UseBranch("Photon");
  TClonesArray *BrMet      = treeReader->UseBranch("MissingET");

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

    HepMCEvent *event   = (HepMCEvent*) BrEvent -> At(0); 
    //int evt_no          = event -> Number;
    float event_weight  = event -> Weight; // original event weight 
    //float cross_section = event -> CrossSection; // in pb
    //float cross_section_err = event -> CrossSectionError; // in pb

    /* ************************ EventWeight = 1 and EventWeightSum = 1 x nEvents ********************* */
    double evWt = 1.0;
    double lumiWt = evWt*lumiFac;  // Lumi Scaling
    
    /* ****************************************************************** */     
    histf->cd();
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
    
    // -------------------------------------------------------------------------------------------------------------------- //
    //                                                 Object Selection [start]                                             //
    // -------------------------------------------------------------------------------------------------------------------- //    
    std::vector <Muon>      MuonColl;
    std::vector <Electron>  ElectronColl;
    std::vector <Photon>    PhotonColl;
    std::vector <Jet>       TauhColl;
    std::vector <Jet>       JetColl;
    std::vector <Jet>       BJetColl;
    std::vector <Jet>       lightJetColl;
    
    histf->cd();
    // **** Muon Selection     
    for (TObject *_mu: *BrMuon){
      const Muon& mu = dynamic_cast<const Muon&> (*_mu);
      AnaUtil::fillHist1D("Raw_mu_pt",       mu.PT);
      AnaUtil::fillHist1D("Raw_mu_eta",      mu.Eta);
      AnaUtil::fillHist1D("Raw_mu_phi",      mu.Phi);
      AnaUtil::fillHist1D("Raw_mu_charge",   mu.Charge);
      AnaUtil::fillHist1D("Raw_mu_isolationVar",    mu.IsolationVar);
      AnaUtil::fillHist1D("Raw_mu_isolationVarRhoCorr",    mu.IsolationVarRhoCorr);
      AnaUtil::fillHist1D("Raw_mu_D0",       mu.D0);
      AnaUtil::fillHist1D("Raw_mu_DZ",       mu.DZ);

      AnaUtil::fillHist1D("MuonCutFlow", 0, evWt);       
      if (mu.PT <= 10) continue;
      AnaUtil::fillHist1D("MuonCutFlow", 1, evWt);
      if (std::abs(mu.Eta) >= 2.4) continue;
      AnaUtil::fillHist1D("MuonCutFlow", 2, evWt);
      TLorentzVector mup4 = mu.P4();
      MuonColl.push_back(mu);
    }
    std::sort(std::begin(MuonColl), std::end(MuonColl), PtComparator<Muon>()); 
    
    // **** Electron Selection
    for (TObject *_el: *BrElec){
      const Electron& ele = dynamic_cast<const Electron&> (*_el);
      AnaUtil::fillHist1D("Raw_ele_pt",       ele.PT);
      AnaUtil::fillHist1D("Raw_ele_eta",      ele.Eta);
      AnaUtil::fillHist1D("Raw_ele_phi",      ele.Phi);
      AnaUtil::fillHist1D("Raw_ele_charge",   ele.Charge);
      AnaUtil::fillHist1D("Raw_ele_isolationVar",    ele.IsolationVar);
      AnaUtil::fillHist1D("Raw_ele_isolationVarRhoCorr",    ele.IsolationVarRhoCorr);
      AnaUtil::fillHist1D("Raw_ele_D0",       ele.D0);
      AnaUtil::fillHist1D("Raw_ele_DZ",       ele.DZ);

      AnaUtil::fillHist1D("ElectronCutFlow", 0, evWt);      
      if (ele.PT <= 10) continue;
      AnaUtil::fillHist1D("ElectronCutFlow", 1, evWt);
      if (std::abs(ele.Eta) >= 2.5) continue;
      AnaUtil::fillHist1D("ElectronCutFlow", 2, evWt);
      ElectronColl.push_back(ele);
    }
    std::sort(std::begin(ElectronColl), std::end(ElectronColl), PtComparator<Electron>());
    
    // **** Photon Selection
    for (TObject *_ph: *BrPhoton){
      const Photon& ph = dynamic_cast<const Photon&> (*_ph);
      AnaUtil::fillHist1D("PhotonCutFlow", 0, evWt);
      if (ph.PT <= 20) continue;
      AnaUtil::fillHist1D("PhotonCutFlow", 1, evWt);
      if (std::abs(ph.Eta) >= 2.4) continue;
      AnaUtil::fillHist1D("PhotonCutFlow", 2, evWt);
    }
    std::sort(std::begin(PhotonColl), std::end(PhotonColl), PtComparator<Photon>());
    
    // **** Tau Selection
    for (TObject *_tau: *BrJet){
      const Jet& tauj = dynamic_cast<const Jet&> (*_tau);
      if (tauj.TauTag == 0) continue;
      AnaUtil::fillHist1D("TauhCutFlow", 0, evWt);
      if (tauj.PT <= 20) continue;
      AnaUtil::fillHist1D("TauhCutFlow", 1, evWt);
      if (std::abs(tauj.Eta) >= 2.4) continue;
      AnaUtil::fillHist1D("TauhCutFlow", 2, evWt);
      if (!(jetLeptonCleaning(tauj, MuonColl, ElectronColl, 0.4))) continue;
      AnaUtil::fillHist1D("TauhCutFlow", 3, evWt);
      TauhColl.push_back(tauj);
    }
    std::sort(std::begin(TauhColl), std::end(TauhColl), PtComparator<Jet>());

    // **** Jet Selection 
    for (TObject *_j: *BrJet){
      const Jet& jet = dynamic_cast<const Jet&> (*_j);
      AnaUtil::fillHist1D("Raw_jet_pt",       jet.PT);
      AnaUtil::fillHist1D("Raw_jet_eta",      jet.Eta);
      AnaUtil::fillHist1D("Raw_jet_phi",      jet.Phi);
      AnaUtil::fillHist1D("Raw_jet_mass",     jet.Mass);
      AnaUtil::fillHist1D("Raw_jet_DeltaEta", jet.DeltaEta);  // jet radius in pseudorapidity
      AnaUtil::fillHist1D("Raw_jet_DeltaPhi", jet.DeltaPhi);  // jet radius in azimuthal angle 
      AnaUtil::fillHist1D("Raw_jet_flavor",   jet.Flavor);
      AnaUtil::fillHist1D("Raw_jet_tautag",   jet.TauTag);    // 0 or 1 for a jet that has been tagged as a tau   
      AnaUtil::fillHist1D("Raw_jet_tauwt",    jet.TauWeight); // probability for jet to be identified as tau 
      AnaUtil::fillHist1D("Raw_jet_charge",   jet.Charge);    // only for tau
      AnaUtil::fillHist1D("Raw_jet_EhadOverEem", jet.EhadOverEem); // ratio of the hadronic versus electromagnetic energy deposited in the calorimeter  
      AnaUtil::fillHist1D("Raw_jet_BTag",     jet.BTag);
      AnaUtil::fillHist1D("Raw_jet_BetaStar", jet.BetaStar);  // (sum pt of charged constituents coming from hard interaction)/(sum pt of charged constituents)

      AnaUtil::fillHist1D("JetCutFlow", 0, evWt);
      if (jet.PT <= 20) continue;
      AnaUtil::fillHist1D("JetCutFlow", 1, evWt);
      if (std::abs(jet.Eta) >= 4.7) continue;
      AnaUtil::fillHist1D("JetCutFlow", 2, evWt);
      if (!(jetLeptonCleaning(jet, MuonColl, ElectronColl, 0.4))) continue;
      AnaUtil::fillHist1D("JetCutFlow", 3, evWt);
      if (!(jetTauhCleaning(jet, TauhColl, 0.4))) continue;
      AnaUtil::fillHist1D("JetCutFlow", 4, evWt);
      JetColl.push_back(jet);
    }
    std::sort(std::begin(JetColl), std::end(JetColl), PtComparator<Jet>());
        
    // **** b-jet selection
    for (auto &jet : JetColl)
      (jet.BTag == 1 && std::fabs(jet.Eta) < 2.5) ? BJetColl.push_back(jet) : lightJetColl.push_back(jet);
    // -------------------------------------------------------------------------------------------------------------------- //
    //                                                 Object Selection [end]                                               //
    // -------------------------------------------------------------------------------------------------------------------- //    

    int nGoodMuon     = MuonColl.size();
    int nGoodEle      = ElectronColl.size();
    int nGoodPhoton   = PhotonColl.size();
    int nGoodTauh     = TauhColl.size();
    int nGoodJet      = JetColl.size();
    int nGoodlJet     = lightJetColl.size();
    int nGoodBJet     = BJetColl.size();

    // -------------------------------------------------------------------------------------------------------------------- //
    //                                                   Basic histogramming [start]                                        //
    // -------------------------------------------------------------------------------------------------------------------- //
    AnaUtil::fillHist1D("nGoodMuon",     nGoodMuon);
    AnaUtil::fillHist1D("nGoodElectron", nGoodEle);
    AnaUtil::fillHist1D("nGoodLepton",   nGoodMuon + nGoodEle); 
    AnaUtil::fillHist1D("nGoodPhoton",   nGoodPhoton);
    AnaUtil::fillHist1D("nGoodTauh",     nGoodTauh);
    AnaUtil::fillHist1D("nGoodJet",      nGoodJet);
    AnaUtil::fillHist1D("nGoodBJet",     nGoodBJet); 
    AnaUtil::fillHist1D("nGoodlJet",     nGoodlJet);
     
    // Choose range of n-objects for plotting several variables
    size_t nLoop2Ele  = (nGoodEle >= 2)  ? 2 : nGoodEle;
    size_t nLoop2Mu   = (nGoodMuon >= 2) ? 2 : nGoodMuon;
    size_t nLoop4Jet  = (nGoodJet >= 4)  ? 4 : nGoodJet;
    size_t nLoop4lJet = (nGoodlJet >= 4) ? 4 : nGoodlJet;
    //size_t nLoop6Jet  = (nGoodJet >= 6)  ? 6 : nGoodJet;
    size_t nLoopBJet  = (nGoodBJet >= 2) ? 2 : nGoodBJet;
    size_t nLoopTauh  = (nGoodTauh > 1)  ? 1 : nGoodTauh;
     
    // b-jet plots (for max 2 b-jets)
    for (size_t i = 0; i < nLoopBJet; ++i) {
      TLorentzVector bjP4 = BJetColl[i].P4();
      std::string pt_name  = "BJet"+std::to_string(i+1)+"Pt";
      std::string phi_name = "BJet"+std::to_string(i+1)+"Phi";
      std::string eta_name = "BJet"+std::to_string(i+1)+"Eta";
      AnaUtil::fillHist1D (pt_name.c_str(),  bjP4.Pt());
      AnaUtil::fillHist1D (phi_name.c_str(), bjP4.Phi());
      AnaUtil::fillHist1D (eta_name.c_str(), bjP4.Eta());
    }
    // Loop over max 2 electrons
    for (size_t iele = 0; iele < nLoop2Ele; ++iele) {
      TLorentzVector eleP4 = ElectronColl[iele].P4();
      std::string elePt_name  = "Electron"+std::to_string(iele+1)+"Pt";
      std::string eleEta_name = "Electron"+std::to_string(iele+1)+"Eta";
      std::string elePhi_name = "Electron"+std::to_string(iele+1)+"Phi";
      AnaUtil::fillHist1D (elePt_name.c_str(),  eleP4.Pt());
      AnaUtil::fillHist1D (eleEta_name.c_str(), eleP4.Eta());
      AnaUtil::fillHist1D (elePhi_name.c_str(), eleP4.Phi());
      // nested loop over max 4 jets
      for (size_t ijet = 0; ijet < nLoop4Jet; ++ijet) {
	TLorentzVector jP4 = JetColl[ijet].P4();
	std::string DR_jetEle_name   = "DR_Jet"+std::to_string(ijet+1)+"Electron"+std::to_string(iele+1);
	std::string DPhi_jetEle_name = "DPhi_Jet"+std::to_string(ijet+1)+"Electron"+std::to_string(iele+1);
	AnaUtil::fillHist1D (DR_jetEle_name.c_str(),   jP4.DeltaR(eleP4));
	AnaUtil::fillHist1D (DPhi_jetEle_name.c_str(), std::fabs(TVector2::Phi_mpi_pi(jP4.Phi()-eleP4.Phi())));
      }
      // nested loop over max 2 bjets
      for (size_t ib = 0; ib < nLoopBJet; ++ib) {
	TLorentzVector bjP4 = BJetColl[ib].P4();
	std::string DR_bjetEle_name   = "DR_bJet"+std::to_string(ib+1)+"Electron"+std::to_string(iele+1);
	std::string DPhi_bjetEle_name = "DPhi_bJet"+std::to_string(ib+1)+"Electron"+std::to_string(iele+1);
	AnaUtil::fillHist1D (DR_bjetEle_name.c_str(),   bjP4.DeltaR(eleP4));
	AnaUtil::fillHist1D (DPhi_bjetEle_name.c_str(), std::fabs(TVector2::Phi_mpi_pi(bjP4.Phi() - eleP4.Phi())));	 
      }
    }
    // Loop over max 2 muons
    for (size_t imu = 0; imu < nLoop2Mu; ++imu) {
      TLorentzVector muP4 = MuonColl[imu].P4();
      std::string muPt_name  = "Muon"+std::to_string(imu+1)+"Pt";
      std::string muEta_name = "Muon"+std::to_string(imu+1)+"Eta";
      std::string muPhi_name = "Muon"+std::to_string(imu+1)+"Phi";
      AnaUtil::fillHist1D (muPt_name.c_str(),  muP4.Pt());
      AnaUtil::fillHist1D (muEta_name.c_str(), muP4.Eta());
      AnaUtil::fillHist1D (muPhi_name.c_str(), muP4.Phi());
      // nested loop over max 1 tau
      for (size_t iTauh = 0; iTauh < nLoopTauh; ++iTauh) {
	TLorentzVector tauhP4 = TauhColl[iTauh].P4();
	std::string DR_tauhMu_name = "DR_mu"+std::to_string(imu+1)+"tauh"+std::to_string(iTauh+1);
	std::string DPhi_tauhMu_name = "Deltaphi_mu"+std::to_string(imu+1)+"tauh"+std::to_string(iTauh+1);
	AnaUtil::fillHist1D (DR_tauhMu_name.c_str(),   tauhP4.DeltaR(muP4));
	AnaUtil::fillHist1D (DPhi_tauhMu_name.c_str(), std::fabs(TVector2::Phi_mpi_pi(tauhP4.Phi() - muP4.Phi())));
      }
      // nested loop over max 4 jets
      for (size_t ijet = 0; ijet < nLoop4Jet; ++ijet) {
	TLorentzVector jP4 = JetColl[ijet].P4();
	std::string DR_jetMu_name = "DR_Jet"+std::to_string(ijet+1)+"Muon"+std::to_string(imu+1);
	std::string DPhi_jetMu_name = "DPhi_Jet"+std::to_string(ijet+1)+"Muon"+std::to_string(imu+1);
	AnaUtil::fillHist1D (DR_jetMu_name.c_str(),   jP4.DeltaR(muP4));
	AnaUtil::fillHist1D (DPhi_jetMu_name.c_str(), std::fabs(TVector2::Phi_mpi_pi(jP4.Phi() - muP4.Phi())));
      }
      // nested loop over max 2 b-jets
      for (size_t ib = 0; ib < nLoopBJet; ++ib) {
	TLorentzVector bjP4 = BJetColl[ib].P4();
	std::string DR_bjetMu_name   = "DR_bJet"+std::to_string(ib+1)+"Muon"+std::to_string(imu+1);
	std::string DPhi_bjetMu_name = "DPhi_bJet"+std::to_string(ib+1)+"Muon"+std::to_string(imu+1);
	AnaUtil::fillHist1D (DR_bjetMu_name.c_str(),   bjP4.DeltaR(muP4));
	AnaUtil::fillHist1D (DPhi_bjetMu_name.c_str(), std::fabs(TVector2::Phi_mpi_pi(bjP4.Phi() - muP4.Phi())));	 
      }
    }
    // Loop over all good jets
    float minDR = 999.9;
    float maxDR = -999.9;
    for (size_t i = 0; i < static_cast<size_t>(nGoodJet); ++i) {
      std::string pt_name  = "Jet"+std::to_string(i+1)+"Pt";
      TLorentzVector jP4   = JetColl[i].P4();
      AnaUtil::fillHist1D (pt_name.c_str(), JetColl[i].PT, 1);
      for (size_t j = i+1; j < static_cast<size_t>(nGoodJet); ++j) {
	TLorentzVector j2P4   = JetColl[j].P4();
	float dR = jP4.DeltaR(j2P4);
	if (dR < minDR) minDR = dR;
	if (dR > maxDR) maxDR = dR;
      }
    }
    AnaUtil::fillHist1D("DR_minimum_allJets", minDR);
    AnaUtil::fillHist1D("DR_maximum_allJets", maxDR);
    // Loop over max 4 light jets
    float min_invM = 9999.9;
    for (size_t i = 0; i < nLoop4lJet; ++i) {
      std::string jetPtName  = "lightJet"+std::to_string(i+1)+"Pt";
      std::string jetEtaName = "lightJet"+std::to_string(i+1)+"Eta";
      TLorentzVector ljetiP4 = lightJetColl[i].P4();
      AnaUtil::fillHist1D (jetPtName.c_str(),  ljetiP4.Pt());
      AnaUtil::fillHist1D (jetEtaName.c_str(), ljetiP4.Eta());
      // nested loop over max 4 light jets
      for (size_t j = i+1; j < nLoop4lJet; ++j) {
	TLorentzVector ljetjP4 =lightJetColl[j].P4();
	std::string PT_histName   = "PT_lightJet"+std::to_string(i+1)+"Jet"+std::to_string(j+1);
	std::string IM_histName   = "IM_lightJet"+std::to_string(i+1)+"Jet"+std::to_string(j+1);
	std::string DR_histName   = "DR_lightJet"+std::to_string(i+1)+"Jet"+std::to_string(j+1);
	std::string DEta_histName = "DEta_lightJet"+std::to_string(i+1)+"Jet"+std::to_string(j+1);
	std::string DPhi_histName = "DPhi_lightJet"+std::to_string(i+1)+"Jet"+std::to_string(j+1);
	float invM = (ljetiP4+ljetjP4).M();
	AnaUtil::fillHist1D (PT_histName.c_str(), (ljetiP4+ljetjP4).Pt());
	AnaUtil::fillHist1D (IM_histName.c_str(), invM);
	AnaUtil::fillHist1D (DR_histName.c_str(), ljetiP4.DeltaR(ljetjP4));
	AnaUtil::fillHist1D (DEta_histName.c_str(), std::fabs(ljetiP4.Eta() - ljetjP4.Eta()));
	AnaUtil::fillHist1D (DPhi_histName.c_str(), std::fabs(TVector2::Phi_mpi_pi(ljetiP4.Phi() - ljetjP4.Phi())));
	if (invM < min_invM) min_invM = invM;
      }
    }
    AnaUtil::fillHist1D("IM_minimum_lightJets",min_invM);
    // -------------------------------------------------------------------------------------------------------------------- //
    //                                                  Basic histogramming [end]                                           //
    // -------------------------------------------------------------------------------------------------------------------- //

    
    // -------------------------------------------------------------------------------------------------------------------- //
    //                                                Pre-selection cuts [start]                                            //
    // -------------------------------------------------------------------------------------------------------------------- //
    if (nGoodMuon + nGoodEle + nGoodTauh + nGoodJet == 0) continue;
    AnaUtil::fillHist1D("evtCutFlow", 1, evWt);
    AnaUtil::fillHist1D("evtCutFlowWt", 1, lumiWt);
    
    if (isDL()) {
      if (nGoodMuon + nGoodEle != 2) continue;
      AnaUtil::fillHist1D("evtCutFlow", 2, evWt);
      AnaUtil::fillHist1D("evtCutFlowWt", 2, lumiWt);
      
      if (nGoodMuon < 1) continue;
      AnaUtil::fillHist1D("evtCutFlow", 3, evWt);
      AnaUtil::fillHist1D("evtCutFlowWt", 3, lumiWt);
      
      if (nGoodTauh != 1) continue;
      AnaUtil::fillHist1D("evtCutFlow", 4, evWt);
      AnaUtil::fillHist1D("evtCutFlowWt", 4, lumiWt);
      
      if (nGoodlJet < 1) continue;
      AnaUtil::fillHist1D("evtCutFlow", 5, evWt);
      AnaUtil::fillHist1D("evtCutFlowWt", 5, lumiWt);
      
      bool has2ndBjet = (nGoodBJet > 1 && BJetColl[1].PT > 30) ? true : false;
      if (nGoodBJet < 1 || has2ndBjet) continue;
      //if (nGoodBJet != 1) continue;
      AnaUtil::fillHist1D("evtCutFlow", 6, evWt);
      AnaUtil::fillHist1D("evtCutFlowWt", 6, lumiWt);
      
      std::vector<ZCandidate> ZCandList;
      if (nGoodMuon >= 2) ZSelector(MuonColl, ZCandList);
      if (nGoodEle  >= 2) ZSelector(ElectronColl, ZCandList);
      if (ZCandList.size() > 0 && ZCandList[0].massDiff < 10) continue;
      AnaUtil::fillHist1D("evtCutFlow", 7, evWt);
      AnaUtil::fillHist1D("evtCutFlowWt", 7, lumiWt);
    }

    else if (isSL()) {
      if (nGoodMuon != 1) continue;
      AnaUtil::fillHist1D("evtCutFlow", 2, evWt);
      AnaUtil::fillHist1D("evtCutFlowWt", 2, lumiWt);
      
      if (nGoodEle != 0) continue;
      AnaUtil::fillHist1D("evtCutFlow", 3, evWt);
      AnaUtil::fillHist1D("evtCutFlowWt", 3, lumiWt);

      if (nGoodTauh != 1) continue;
      AnaUtil::fillHist1D("evtCutFlow", 4, evWt);
      AnaUtil::fillHist1D("evtCutFlowWt", 4, lumiWt);     

      if (nGoodlJet < 3) continue;
      AnaUtil::fillHist1D("evtCutFlow", 5, evWt);
      AnaUtil::fillHist1D("evtCutFlowWt", 5, lumiWt);

      bool has2ndBjet = (nGoodBJet > 1 && BJetColl[1].PT > 30) ? true : false;
      if (nGoodBJet < 1 || has2ndBjet) continue;
      AnaUtil::fillHist1D("evtCutFlow", 6, evWt);
      AnaUtil::fillHist1D("evtCutFlowWt", 6, lumiWt);
    }

    else 
      std::cerr 
	<< "Channel isnt mentioned properly in the jobcard! No preselection cuts are applied ..."
	<< std::endl;
 
    // -------------------------------------------------------------------------------------------------------------------- //
    //                                                 Pre-selection cuts [end]                                             //
    // -------------------------------------------------------------------------------------------------------------------- //
   

    // Channel specific

    // met, tauh, tau (reconstructed from collinear mass correction), lead-b jet, lead-light jet are fixed for both SL and DL channels
    // these will be used throughout the whole analysis
    TLorentzVector metp4;
    metp4.SetPtEtaPhiE(met_->MET, 0.0, met_->Phi, met_->MET);
    AnaUtil::fillHist1D("met", metp4.Pt());
    
    TLorentzVector tauhp4   = TauhColl[0].P4();
    TLorentzVector leadbjp4 = BJetColl[0].P4();
    TLorentzVector leadljp4 = lightJetColl[0].P4();

    AnaUtil::fillHist1D("Tauh_Pt",       tauhp4.Pt());
    AnaUtil::fillHist1D("Tauh_Eta",      tauhp4.Eta());
    AnaUtil::fillHist1D("Tauh_Phi",      tauhp4.Phi());
    
    // tauh-leadbj
    float DR_tauh_leadbj   = tauhp4.DeltaR(leadbjp4);
    float DPhi_tauh_leadbj = TVector2::Phi_mpi_pi(tauhp4.Phi() - leadbjp4.Phi());
    float DEta_tauh_leadbj = tauhp4.Eta() - leadbjp4.Eta();
    AnaUtil::fillHist1D("DR_tauh_leadbj",   DR_tauh_leadbj);
    AnaUtil::fillHist1D("DPhi_tauh_leadbj", DPhi_tauh_leadbj);
    AnaUtil::fillHist1D("DEta_tauh_leadbj", DEta_tauh_leadbj);
    // tauh-leadlj
    float DR_tauh_leadlj   = tauhp4.DeltaR(leadljp4);
    float DPhi_tauh_leadlj = TVector2::Phi_mpi_pi(tauhp4.Phi() - leadljp4.Phi());
    float DEta_tauh_leadlj = tauhp4.Eta() - leadljp4.Eta();
    AnaUtil::fillHist1D("DR_tauh_leadlj",   DR_tauh_leadlj);
    AnaUtil::fillHist1D("DPhi_tauh_leadlj", DPhi_tauh_leadlj);
    AnaUtil::fillHist1D("DEta_tauh_leadlj", DEta_tauh_leadlj);
    // leadlj-leadbj
    float DR_leadbj_leadlj   = leadbjp4.DeltaR(leadljp4);
    float DPhi_leadbj_leadlj = TVector2::Phi_mpi_pi(leadbjp4.Phi() - leadljp4.Phi());
    float DEta_leadbj_leadlj = leadbjp4.Eta() - leadljp4.Eta();
    AnaUtil::fillHist1D("DR_leadbj_leadlj",   DR_leadbj_leadlj);
    AnaUtil::fillHist1D("DPhi_leadbj_leadlj", DPhi_leadbj_leadlj);
    AnaUtil::fillHist1D("DEta_leadbj_leadlj", DEta_leadbj_leadlj);
    // met-tauh,leadbj,leadlj
    float DPhi_met_tauh   = TVector2::Phi_mpi_pi(metp4.Phi() - tauhp4.Phi());
    float DPhi_met_leadbj = TVector2::Phi_mpi_pi(metp4.Phi() - leadbjp4.Phi());
    float DPhi_met_leadlj = TVector2::Phi_mpi_pi(metp4.Phi() - leadljp4.Phi());
    AnaUtil::fillHist1D("DPhi_met_tauh",    DPhi_met_tauh);
    AnaUtil::fillHist1D("DPhi_met_leadbj",  DPhi_met_leadbj);
    AnaUtil::fillHist1D("DPhi_met_leadlj",  DPhi_met_leadlj);
    //tauh-leadj
    float DR_tauh_leadj   = (tauhp4).DeltaR(JetColl[0].P4());
    float DPhi_tauh_leadj = TVector2::Phi_mpi_pi(tauhp4.Phi() - JetColl[0].Phi);
    float DEta_tauh_leadj = JetColl[0].Eta - tauhp4.Eta();
    AnaUtil::fillHist1D("DR_tauh_leadj",    DR_tauh_leadj);
    AnaUtil::fillHist1D("DPhi_tauh_leadj",  DPhi_tauh_leadj);
    AnaUtil::fillHist1D("DEta_tauh_leadj",  DEta_tauh_leadj);

    float pTvis_tau      = tauhp4.Pt();
    float pTnu           = (metp4.Px()*tauhp4.Px() + metp4.Py()*tauhp4.Py())/std::fabs(pTvis_tau);
    float x_vis          = std::fabs(pTvis_tau)/(std::fabs(pTvis_tau) + std::fabs(pTnu));
    TLorentzVector taup4; 
    // [tauh_visible + tau_neutrino] P4
    // Not suitable for DL becuase W has a real neutrino as its decay product. So total met 
    // doesnt represent the contribution of tau neutrino (decay product of Chi) only.
    // It should be appropriate for SL.
    taup4.SetPtEtaPhiE((1.0/x_vis)*tauhp4.Pt(), tauhp4.Eta(), tauhp4.Phi(), (1.0/x_vis)*tauhp4.E());
    // 2202.04336
    TLorentzVector nup4;
    nup4.SetPtEtaPhiE((1-x_vis)*taup4.Pt(), taup4.Eta(), taup4.Phi(), (1-x_vis)*taup4.E());

    float tau_pt  = taup4.Pt();
    float tau_eta = taup4.Eta();
    AnaUtil::fillHist1D("PT_RecoTau",  tau_pt);
    AnaUtil::fillHist1D("Eta_RecoTau", tau_eta);
    AnaUtil::fillHist1D("PT_RecoNu",   nup4.Pt());
    AnaUtil::fillHist1D("Eta_RecoNu",  nup4.Eta());
    // ------------------------------------------------------------------------------ //
    //                                        DL                                      //
    //            p p > t t~, (t > b w, w > l nu ), (t~ > X j, X > mu ta)             //
    // ------------------------------------------------------------------------------ //
    if (isDL()) {
      TLorentzVector XLepP4, WLepP4;
      bool MuMu  = (MuonColl.size() == 2) ? true : false;
      bool EleMu = (MuMu) ? false : true;
      bool oppCharge {false};
      
      // ------------ MuMu channel ------------- //
      if (MuMu) { 
	bool leadMuCloseToTauh = (tauhp4.DeltaR(MuonColl[0].P4()) <= tauhp4.DeltaR(MuonColl[1].P4())) ? true : false;
	if (leadMuCloseToTauh) {	                            // if leading muon is the closest muon to tauh
	  if (MuonColl[0].Charge * TauhColl[0].Charge < 0.0) {      // lead-mu and tauh : opposite charge
	    oppCharge = true;                                       // So, oppCharge = true
	    XLepP4 = MuonColl[0].P4();                              // XLep = leading muon 
	    WLepP4 = MuonColl[1].P4();                              // Wlep = subleading muon 
	  }
	  else if (MuonColl[1].Charge * TauhColl[0].Charge < 0.0) { // if above condition is failed,
                                                                    // i.e. lead-mu and tauh : same charge
                                                                    // sublead-mu and tauh : opposite charge
	    oppCharge = true;                                       // So, oppCharge = true 
	    XLepP4 = MuonColl[1].P4();                              // XLep = subleading muon
	    WLepP4 = MuonColl[0].P4();                              // Wlep = leading muon
	  }
	  else oppCharge = false;                                   // both muons and tauh : same charge
	}

	else {                                                      // if sublead muon is the closest one to tauh
	  if (MuonColl[1].Charge * TauhColl[0].Charge < 0.0) {      // sublead-mu and tauh : opposite charge
	    oppCharge =true;                                        // So, oppCharge = true 
	    XLepP4 = MuonColl[1].P4();                              // Xlep = subleading muon
	    WLepP4 = MuonColl[0].P4();                              // Wlep = leading muon 
	  }
	  else if (MuonColl[0].Charge * TauhColl[0].Charge < 0.0) { // if above condition is failed, 
	                                                            // i.e. sublead-mu and tauh : same charge
                                                                    // lead-mu and tauh : opposite charge  
	    oppCharge = true;                                       // So, oppCharge = true   
	    XLepP4 = MuonColl[0].P4();                              // XLep = leading muon  
	    WLepP4 = MuonColl[1].P4();	                            // Wlep = subleading muon 
	  }
	  else oppCharge = false;                                   // both muons and tauh : same charge 
	}
      }
  
      /*	if (isLowMassX()) {
	  XLepP4 = MuonColl[1].P4();
	  WLepP4 = MuonColl[0].P4();
	  oppCharge = (MuonColl[1].Charge * TauhColl[0].Charge < 0.0) ? true : false; 
	}
	else {
	  XLepP4 = MuonColl[0].P4();
	  WLepP4 = MuonColl[1].P4();
	  oppCharge = (MuonColl[0].Charge * TauhColl[0].Charge < 0.0) ? true : false;
	}
	}
      */

      // ------------ EleMu channel ------------- //
      else if (EleMu) {
	XLepP4 = MuonColl[0].P4();
	WLepP4 = ElectronColl[0].P4();
	oppCharge = (MuonColl[0].Charge * TauhColl[0].Charge < 0.0) ? true : false;
      }

      else std::cerr << "Wrong lepton selection !!! :( \n";

      //if (XLepP4.DeltaR(WLepP4) < 0.4) continue;
      AnaUtil::fillHist1D("evtCutFlow", 8, evWt);
      AnaUtil::fillHist1D("evtCutFlowWt", 8, lumiWt);
      
      if (!oppCharge) continue;
      AnaUtil::fillHist1D("evtCutFlow", 9, evWt);
      AnaUtil::fillHist1D("evtCutFlowWt", 9, lumiWt);

      //////////////////////////////////////////////////////
      ////////// variables from one leg (t/t~) /////////////
      //////////////////////////////////////////////////////

      float XlepPt           = XLepP4.Pt();
      float XlepEta          = XLepP4.Eta();
      float MET_along_tauh   = taup4.Pt();
      float FracMetAlongTauh = MET_along_tauh/metp4.Pt();
      // xlep - tauh
      float DR_xlep_tauh     = XLepP4.DeltaR(tauhp4);
      float DPhi_xlep_tauh   = TVector2::Phi_mpi_pi(XLepP4.Phi() - tauhp4.Phi());
      float DEta_xlep_tauh   = XLepP4.Eta() - tauhp4.Eta();
      // xlep - lead lj
      float DR_xlep_leadlj   = XLepP4.DeltaR(leadljp4);
      float DPhi_xlep_leadlj = TVector2::Phi_mpi_pi(XLepP4.Phi() - leadljp4.Phi());
      float DEta_xlep_leadlj = XLepP4.Eta() - leadljp4.Eta();
      // xlep - lead j
      float DR_xlep_leadj    = XLepP4.DeltaR(JetColl[0].P4());
      float DPhi_xlep_leadj  = TVector2::Phi_mpi_pi(XLepP4.Phi() - JetColl[0].Phi);
      float DEta_xlep_leadj  = XLepP4.Eta() - JetColl[0].Eta;
      // dphi - met - other
      float DPhi_xlep_met    = TVector2::Phi_mpi_pi(XLepP4.Phi() - metp4.Phi());  
      // masses
      float MT_xlep_met      = Calculate_MT(XLepP4, metp4);
      float M_coll_Xlep      = (XLepP4 + tauhp4).M()/TMath::Sqrt(x_vis);  
      float EffectiveMass_t1 = (taup4 + XLepP4 + leadljp4).M();
      float InvMass_X_test   = (taup4 + tauhp4).M();      
      float XMT              = Calculate_TotalMT(XLepP4, tauhp4, metp4);

      AnaUtil::fillHist1D("XlepPt",  XLepP4.Pt());
      AnaUtil::fillHist1D("XlepEta", XLepP4.Eta());
      AnaUtil::fillHist1D("MET_along_tauh", MET_along_tauh);
      AnaUtil::fillHist1D("FracMetAlongTauh", FracMetAlongTauh);
      AnaUtil::fillHist1D("DR_xlep_tauh", DR_xlep_tauh);
      AnaUtil::fillHist1D("DPhi_xlep_tauh", DPhi_xlep_tauh);
      AnaUtil::fillHist1D("DEta_xlep_tauh", DEta_xlep_tauh);
      AnaUtil::fillHist1D("DR_xlep_leadlj", DR_xlep_leadlj);
      AnaUtil::fillHist1D("DPhi_xlep_leadlj", DPhi_xlep_leadlj);
      AnaUtil::fillHist1D("DEta_xlep_leadlj", DEta_xlep_leadlj);
      AnaUtil::fillHist1D("DR_xlep_leadj", DR_xlep_leadj);
      AnaUtil::fillHist1D("DPhi_xlep_leadj", DPhi_xlep_leadj);
      AnaUtil::fillHist1D("DEta_xlep_leadj", DEta_xlep_leadj);
      AnaUtil::fillHist1D("DPhi_xlep_met", DPhi_xlep_met);
      AnaUtil::fillHist1D("EffectiveMass_t1", EffectiveMass_t1);
      AnaUtil::fillHist1D("InvMass_X_test", InvMass_X_test);
      AnaUtil::fillHist1D("MT_X", XMT);
      AnaUtil::fillHist1D("InvM_coln_XlepTauh", M_coll_Xlep);
      AnaUtil::fillHist1D("MT_xlep_met", MT_xlep_met);

      //////////////////////////////////////////////////////
      ////////// variables from other leg (t~/t) ///////////
      //////////////////////////////////////////////////////
      float WlepPt  = WLepP4.Pt();
      float WlepEta = WLepP4.Eta();
      float DR_wlep_leadbj = WLepP4.DeltaR(leadbjp4);
      float DPhi_wlep_leadbj = TVector2::Phi_mpi_pi(WLepP4.Phi() - leadbjp4.Phi());
      float DEta_wlep_leadbj = WLepP4.Eta() - leadbjp4.Eta();
      float DPhi_wlep_met    = TVector2::Phi_mpi_pi(WLepP4.Phi() - metp4.Phi());
      float MT_wlep_met = Calculate_MT(WLepP4, metp4);
      float M_coll_Wlep = (WLepP4 + tauhp4).M()/TMath::Sqrt(x_vis);
      float EffectiveMass_t2 = (WLepP4 + metp4 + leadbjp4).M();      

      AnaUtil::fillHist1D("WlepPt",  WLepP4.Pt());
      AnaUtil::fillHist1D("WlepEta", WLepP4.Eta());      
      AnaUtil::fillHist1D("DR_wlep_leadbj", DR_wlep_leadbj);
      AnaUtil::fillHist1D("DPhi_wlep_leadbj", DPhi_wlep_leadbj);
      AnaUtil::fillHist1D("DEta_wlep_leadbj", DEta_wlep_leadbj);
      AnaUtil::fillHist1D("MT_wlep_met",      MT_wlep_met);
      AnaUtil::fillHist1D("InvM_coln_WlepTauh", M_coll_Wlep);      

      //////////////////////////////////////////////////////
      ///////////// communicate between 2 legs /////////////
      //////////////////////////////////////////////////////

      float Scalarsumpt_xlep_wlep = XLepP4.Pt()+WLepP4.Pt();
      float Vectorsumpt_xlep_wlep = (XLepP4+WLepP4).Pt();
      float DR_xlep_wlep          = XLepP4.DeltaR(WLepP4);
      float DPhi_xlep_wlep        = TVector2::Phi_mpi_pi(XLepP4.Phi() - WLepP4.Phi());
      float DEta_xlep_wlep        = XLepP4.Eta() - WLepP4.Eta(); 
      float DR_xlep_leadbj        = XLepP4.DeltaR(leadbjp4);
      float DPhi_xlep_leadbj      = TVector2::Phi_mpi_pi(XLepP4.Phi() - leadbjp4.Phi()); 
      float DEta_xlep_leadbj      = leadbjp4.Eta() - XLepP4.Eta();
      float DR_wlep_leadlj        = WLepP4.DeltaR(leadljp4);
      float DPhi_wlep_leadlj      = TVector2::Phi_mpi_pi(WLepP4.Phi() - leadljp4.Phi()); 
      float DEta_wlep_leadlj      = leadljp4.Eta() - WLepP4.Eta();
      float DR_wlep_tauh          = WLepP4.DeltaR(tauhp4);
      float DPhi_wlep_tauh        = TVector2::Phi_mpi_pi(WLepP4.Phi() - tauhp4.Phi()); 
      float DEta_wlep_tauh        = tauhp4.Eta() - WLepP4.Eta();

      float Scalarsumpt_vis       = XLepP4.Pt() + WLepP4.Pt() + leadbjp4.Pt() + leadljp4.Pt();
      float Scalarsumpt           = XLepP4.Pt() + WLepP4.Pt() + leadbjp4.Pt() + leadljp4.Pt() + metp4.Pt();
      float Effectivemass_vis     = (XLepP4 + WLepP4 + leadbjp4 + leadljp4).M();
      float Effectivemass         = (XLepP4 + WLepP4 + leadbjp4 + leadljp4 + metp4).M();

      if (MuMu) {
	AnaUtil::fillHist1D("InvM_coln_XlepTauh_MuMu", M_coll_Xlep);
	AnaUtil::fillHist1D("DR_XlepTauh_MuMu", DR_xlep_tauh);
	AnaUtil::fillHist1D("DR_WlepTauh_MuMu", DR_wlep_tauh);
	AnaUtil::fillHist1D("DPhi_XlepTauh_MuMu", DPhi_xlep_tauh);
	AnaUtil::fillHist1D("DPhi_WlepTauh_MuMu", DPhi_wlep_tauh);
	AnaUtil::fillHist2D("DR_XlepTauh_Vs_WlepTauh_MuMu", DR_xlep_tauh, DR_wlep_tauh);
	AnaUtil::fillHist1D("IsOppCharge_MuMu", oppCharge ? 1 : 0);
	AnaUtil::fillHist2D("IsOppCharge_DR_XlepTauh_MuMu", oppCharge ? 1 : 0, DR_xlep_tauh);
	AnaUtil::fillHist2D("IsOppCharge_DR_WlepTauh_MuMu", oppCharge ? 1 : 0, DR_wlep_tauh);
      }
      else if (EleMu) {
	AnaUtil::fillHist1D("InvM_coln_XlepTauh_EleMu", M_coll_Xlep);
	AnaUtil::fillHist1D("DR_XlepTauh_EleMu", DR_xlep_tauh);
	AnaUtil::fillHist1D("DR_WlepTauh_EleMu", DR_wlep_tauh);
	AnaUtil::fillHist1D("DPhi_XlepTauh_EleMu", DPhi_xlep_tauh);
	AnaUtil::fillHist1D("DPhi_WlepTauh_EleMu", DPhi_wlep_tauh);
	AnaUtil::fillHist2D("DR_XlepTauh_Vs_WlepTauh_EleMu", DR_xlep_tauh, DR_wlep_tauh);
	AnaUtil::fillHist1D("IsOppCharge_EleMu", oppCharge ? 1 : 0);
	AnaUtil::fillHist2D("IsOppCharge_DR_XlepTauh_EleMu", oppCharge ? 1 : 0, DR_xlep_tauh);
	AnaUtil::fillHist2D("IsOppCharge_DR_WlepTauh_EleMu", oppCharge ? 1 : 0, DR_wlep_tauh);
      }
      else     
	std::cerr << "Wrong category ... not MuMu or EleMu \n";


      AnaUtil::fillHist1D("Scalarsumpt_xlep_wlep", Scalarsumpt_xlep_wlep);
      AnaUtil::fillHist1D("Vectorsumpt_xlep_wlep", Vectorsumpt_xlep_wlep);
      AnaUtil::fillHist1D("DR_xlep_wlep",          DR_xlep_wlep);
      AnaUtil::fillHist1D("DPhi_xlep_wlep",        DPhi_xlep_wlep);
      AnaUtil::fillHist1D("DEta_xlep_wlep",        DEta_xlep_wlep);
      AnaUtil::fillHist1D("DR_xlep_leadbj",        DR_xlep_leadbj);
      AnaUtil::fillHist1D("DPhi_xlep_leadbj",      DPhi_xlep_leadbj);
      AnaUtil::fillHist1D("DEta_xlep_leadbj",      DEta_xlep_leadbj);
      AnaUtil::fillHist1D("DR_wlep_leadlj",        DR_wlep_leadlj);
      AnaUtil::fillHist1D("DPhi_wlep_leadlj",      DPhi_wlep_leadlj);
      AnaUtil::fillHist1D("DEta_wlep_leadlj",      DEta_wlep_leadlj);
      AnaUtil::fillHist1D("DR_wlep_tauh",          DR_wlep_tauh);
      AnaUtil::fillHist1D("DPhi_wlep_tauh",        DPhi_wlep_tauh);
      AnaUtil::fillHist1D("DEta_wlep_tauh",        DEta_wlep_tauh);
      AnaUtil::fillHist1D("Scalarsumpt_vis",       Scalarsumpt_vis);
      AnaUtil::fillHist1D("Scalarsumpt",           Scalarsumpt);
      AnaUtil::fillHist1D("Effectivemass_vis",     Effectivemass_vis);
      AnaUtil::fillHist1D("Effectivemass",         Effectivemass);


      float HT_jets = 0.0;
      double minDr_XlepJets = 999.9;
      double minDr_WlepJets = 999.9;
      double maxDr_XlepJets = -999.9;
      double maxDr_WlepJets = -999.9;
      double minDphi_XlepJets = 999.9;
      double minDphi_WlepJets = 999.9;
      double maxDphi_XlepJets = -999.9;
      double maxDphi_WlepJets = -999.9;
      for (auto &jet : JetColl) {
	HT_jets += jet.PT;
	double dr_jetXlep   = XLepP4.DeltaR(jet.P4());
	double dphi_jetXlep = std::fabs(TVector2::Phi_mpi_pi(jet.P4().Phi() - XLepP4.Phi())); 
	double dr_jetWlep   = WLepP4.DeltaR(jet.P4());
	double dphi_jetWlep = std::fabs(TVector2::Phi_mpi_pi(jet.P4().Phi() - WLepP4.Phi()));
	if (dr_jetXlep <= minDr_XlepJets) minDr_XlepJets = dr_jetXlep;
	if (dr_jetWlep <= minDr_WlepJets) minDr_WlepJets = dr_jetWlep;
	if (dr_jetXlep > maxDr_XlepJets)  maxDr_XlepJets = dr_jetXlep;
	if (dr_jetWlep > maxDr_WlepJets)  maxDr_WlepJets = dr_jetWlep;
	if (dphi_jetXlep <= minDphi_XlepJets) minDphi_XlepJets = dphi_jetXlep;
	if (dphi_jetWlep <= minDphi_WlepJets) minDphi_WlepJets = dphi_jetWlep;
	if (dphi_jetXlep > maxDphi_XlepJets)  maxDphi_XlepJets = dphi_jetXlep;
	if (dphi_jetWlep > maxDphi_WlepJets)  maxDphi_WlepJets = dphi_jetWlep;
      }
      AnaUtil::fillHist1D ("HT_Jets",          HT_jets);
      AnaUtil::fillHist1D ("minDr_XlepJets",   minDr_XlepJets);
      AnaUtil::fillHist1D ("maxDr_XlepJets",   maxDr_XlepJets);
      AnaUtil::fillHist1D ("minDphi_XlepJets", minDphi_XlepJets);
      AnaUtil::fillHist1D ("maxDphi_XlepJets", maxDphi_XlepJets);
      AnaUtil::fillHist1D ("minDr_WlepJets",   minDr_WlepJets);
      AnaUtil::fillHist1D ("maxDr_WlepJets",   maxDr_WlepJets);
      AnaUtil::fillHist1D ("minDphi_WlepJets", minDphi_WlepJets);
      AnaUtil::fillHist1D ("maxDphi_WlepJets", maxDphi_WlepJets);
      
      float smin = comp_smin(XLepP4, WLepP4, tauhp4, leadbjp4, leadljp4, metp4);
      AnaUtil::fillHist1D ("smin", smin);

      // Filling the branches of a flat ntuple for skimming
      if (skimObj_) {
	TreeVariablesDL varList;
	
	varList.event         = iEntry;
	varList.event_wt      = event_weight;
	// numbers of objects in final state
	varList.nleptons      = nGoodMuon + nGoodEle;
	varList.njets         = nGoodJet;
	varList.nbjets        = nGoodBJet;
	varList.nljets        = nGoodlJet;
	varList.ntauh         = nGoodTauh;
	// --------- Low level variables ---------- //
	// ==> Tauh
	varList.px_tauh       = tauhp4.Px();
	varList.py_tauh       = tauhp4.Py();
	varList.pz_tauh       = tauhp4.Pz();
	varList.energy_tauh   = tauhp4.E();
	varList.pt_tauh       = tauhp4.Pt();
	varList.eta_tauh      = std::fabs(tauhp4.Eta());
	varList.phi_tauh      = std::fabs(tauhp4.Phi());
	// ==> MET
	varList.px_met        = metp4.Px();
	varList.py_met        = metp4.Py();
	varList.pz_met        = metp4.Pz();
	varList.energy_met    = metp4.E();
	varList.met           = met_->MET;
	varList.phi_met       = metp4.Phi();
	// ==> bjet
	varList.px_leadbj     = leadbjp4.Px();
	varList.py_leadbj     = leadbjp4.Py();
	varList.pz_leadbj     = leadbjp4.Pz();
	varList.energy_leadbj = leadbjp4.E();
	varList.pt_leadbj     = leadbjp4.Pt();
	varList.eta_leadbj    = std::fabs(leadbjp4.Eta());
	varList.phi_leadbj    = std::fabs(leadbjp4.Phi());
	// ==> light jet
	varList.px_leadlj     = leadljp4.Px();
	varList.py_leadlj     = leadljp4.Py();
	varList.pz_leadlj     = leadljp4.Pz();
	varList.energy_leadlj = leadljp4.E();
	varList.pt_leadlj     = leadljp4.Pt();
	varList.eta_leadlj    = std::fabs(leadljp4.Eta());
	varList.phi_leadlj    = std::fabs(leadljp4.Phi());
	// ==> lepton from Chi
	varList.px_xlep       = XLepP4.Px();
	varList.py_xlep       = XLepP4.Py();
	varList.pz_xlep       = XLepP4.Pz();
	varList.energy_xlep   = XLepP4.E();
	varList.pt_xlep       = XlepPt;
	varList.eta_xlep      = std::fabs(XlepEta);
	varList.phi_xlep      = std::fabs(XLepP4.Phi());
	// ==> lepton from W
	varList.px_wlep       = WLepP4.Px();
	varList.py_wlep       = WLepP4.Py();
	varList.pz_wlep       = WLepP4.Pz();
	varList.energy_xlep   = XLepP4.E();
	varList.pt_wlep       = WlepPt;
	varList.eta_wlep      = std::fabs(WlepEta);
	varList.phi_wlep      = std::fabs(WLepP4.Phi());

	// --------- High level variables ---------- //
	// from one leg
	varList.met_along_tauh       = MET_along_tauh;
	varList.frac_met_along_tauh  = FracMetAlongTauh;
	// xlep - tauh
	varList.dr_xlep_tauh         = DR_xlep_tauh;
	varList.dphi_xlep_tauh       = std::fabs(DPhi_xlep_tauh);
	varList.deta_xlep_tauh       = std::fabs(DEta_xlep_tauh);
	// xlep - lead lj
	varList.dr_xlep_leadlj       = DR_xlep_leadlj;
	varList.dphi_xlep_leadlj     = std::fabs(DPhi_xlep_leadlj);
	varList.deta_xlep_leadlj     = std::fabs(DEta_xlep_leadlj);
	// xlep - lead j
	varList.dr_xlep_leadj        = DR_xlep_leadj;
	varList.dphi_xlep_leadj      = std::fabs(DPhi_xlep_leadj);
	varList.deta_xlep_leadj      = std::fabs(DEta_xlep_leadj);
	// tauh - lead j
	varList.dr_tauh_leadj        = DR_tauh_leadj;
	varList.dphi_tauh_leadj      = std::fabs(DPhi_tauh_leadj);
	varList.deta_tauh_leadj      = std::fabs(DEta_tauh_leadj);
	// tauh - lead lj
	varList.dr_tauh_leadlj       = DR_tauh_leadlj;
	varList.dphi_tauh_leadlj     = std::fabs(DPhi_tauh_leadlj);
	varList.deta_tauh_leadlj     = std::fabs(DEta_tauh_leadlj);
	// dphi - met - other
	varList.dphi_met_tauh        = std::fabs(DPhi_met_tauh);
	varList.dphi_met_leadlj      = std::fabs(DPhi_met_leadlj);
	varList.dphi_met_xlep        = std::fabs(DPhi_xlep_met);
	// dr dphi xlep-jets min/max
	varList.dr_min_xlep_jets     = minDr_XlepJets;
	varList.dphi_min_xlep_jets   = minDphi_XlepJets;
	varList.dr_max_xlep_jets     = maxDr_XlepJets;
	varList.dphi_max_xlep_jets   = maxDphi_XlepJets;
	// masses
	varList.mt_xlep_met          = MT_xlep_met;
	varList.m_coll_xlep          = M_coll_Xlep;
	varList.effectiveMass_t1     = EffectiveMass_t1;
	varList.invMass_x_test       = InvMass_X_test;
	varList.mt_x                 = XMT;

	// from other leg
	// dphi - met -others
	varList.dphi_met_leadbj       = std::fabs(DPhi_met_leadbj);
	varList.dphi_met_wlep         = std::fabs(DPhi_wlep_met);
	// wlep - lead bj
	varList.dr_wlep_leadbj        = DR_wlep_leadbj;
	varList.dphi_wlep_leadbj      = std::fabs(DPhi_wlep_leadbj);
	varList.deta_wlep_leadbj      = std::fabs(DEta_wlep_leadbj);
	// masses
	varList.mt_wlep_met           = MT_wlep_met;
	varList.m_coll_wlep           = M_coll_Wlep;
	varList.effectiveMass_t2      = EffectiveMass_t2; 

	// communicate between 2 legs
	// tauh-bjet variables
	varList.dr_tauh_leadbj        = DR_tauh_leadbj;
	varList.dphi_tauh_leadbj      = std::fabs(DPhi_tauh_leadbj);
	varList.deta_tauh_leadbj      = std::fabs(DEta_tauh_leadbj);
	// xlep-bjet
	varList.dr_xlep_leadbj        = DR_xlep_leadbj;
	varList.dphi_xlep_leadbj      = std::fabs(DPhi_xlep_leadbj);
	varList.deta_xlep_leadbj      = std::fabs(DEta_xlep_leadbj);
	// xlep - wlep
	varList.dr_xlep_wlep          = DR_xlep_wlep;
	varList.dphi_xlep_wlep        = std::fabs(DPhi_xlep_wlep);
	varList.deta_xlep_wlep        = std::fabs(DEta_xlep_wlep);
	varList.scalarsumpt_xlep_wlep = Scalarsumpt_xlep_wlep;
	varList.vectorsumpt_xlep_wlep = Vectorsumpt_xlep_wlep;
	// wlep - lead lj
	varList.dr_wlep_leadlj        = DR_wlep_leadlj;
	varList.dphi_wlep_leadlj      = std::fabs(DPhi_wlep_leadlj);
	varList.deta_wlep_leadlj      = std::fabs(DEta_wlep_leadlj);
	// wlep - tauh
	varList.dr_wlep_tauh          = DR_wlep_tauh;
	varList.dphi_wlep_tauh        = std::fabs(DPhi_wlep_tauh);
	varList.deta_wlep_tauh        = std::fabs(DEta_wlep_tauh);
	// lepton-jet variables
	varList.dr_min_wlep_jets      = minDr_WlepJets;
	varList.dphi_min_wlep_jets    = minDphi_WlepJets;
	varList.dr_max_wlep_jets      = maxDr_WlepJets;
	varList.dphi_max_wlep_jets    = maxDphi_WlepJets;
	// di-jets variables
	varList.dr_min_jets           = minDR;
	varList.dr_max_jets           = maxDR;
	varList.dr_leadbj_leadlj      = DR_leadbj_leadlj;
	varList.dphi_leadbj_leadlj    = std::fabs(DPhi_leadbj_leadlj);
	varList.deta_leadbj_leadlj    = std::fabs(DEta_leadbj_leadlj);
	// Other observables
	varList.scalarsumpt_vis       = Scalarsumpt_vis;
	varList.scalarsumpt           = Scalarsumpt;
	varList.effectivemass_vis     = Effectivemass_vis;
	varList.effectivemass         = Effectivemass;
	varList.ht_jets               = HT_jets;
	varList.smin                  = smin;

	skimObj_->fill(varList);
      }       
      histf->cd();

      double mvaOut = -999.9;
      
      if (_readMVA) {
	InputVariablesDL varList;
	/*
	varList.dphi_met_tauh        = std::fabs(DPhi_met_tauh);
	varList.dphi_wlep_leadbj     = std::fabs(DPhi_wlep_leadbj);
	varList.dphi_xlep_wlep       = std::fabs(DPhi_xlep_wlep);
	varList.dr_xlep_tauh         = DR_xlep_tauh;
	varList.effectivemass        = Effectivemass;
	varList.ht_jets              = HT_jets;
	varList.mt_wlep_met          = MT_wlep_met;
	varList.mt_x                 = XMT;
	varList.met                  = met_->MET;
	varList.scalarsumpt          = Scalarsumpt;
	varList.dphi_leadbj_leadlj   = std::fabs(DPhi_leadbj_leadlj);
	varList.m_coll_xlep          = M_coll_Xlep;
	varList.m_coll_wlep          = M_coll_Wlep;
	*/
	varList.pt_tauh=                 tauhp4.Pt();
	varList.met=                     met_->MET; 
	varList.pt_leadbj=               leadbjp4.Pt();
	varList.pt_xlep=                 XlepPt;
	varList.pt_wlep=                 WlepPt;
	varList.vectorsumpt_xlep_wlep=   Vectorsumpt_xlep_wlep;
	varList.dr_xlep_wlep=            DR_xlep_wlep;
	varList.dphi_wlep_tauh=          std::fabs(DPhi_wlep_tauh);
	varList.dphi_wlep_leadbj=        std::fabs(DPhi_wlep_leadbj);
	varList.ht_jets=                 HT_jets;
	varList.dr_min_xlep_jets=        minDr_XlepJets;
	varList.dr_min_wlep_jets=        minDr_WlepJets;
	varList.dphi_xlep_tauh=          std::fabs(DPhi_xlep_tauh);
	varList.effectivemass=           Effectivemass;
	varList.dr_min_jets=             minDR;
	varList.mt_wlep_met=             MT_wlep_met;
	varList.dphi_tauh_leadbj=        std::fabs(DPhi_tauh_leadbj);
	varList.dphi_met_wlep=           std::fabs(DPhi_wlep_met);
	varList.dphi_leadbj_leadlj=      std::fabs(DPhi_leadbj_leadlj);
	

	mvaOut = _mvaObj->evaluate(_MVAnetwork, varList);
      }
      histf->cd();
      
      AnaUtil::fillHist1D("BDT_Score", mvaOut);
      AnaUtil::fillHist1D("BDT_Score_fineBin", mvaOut);
    }

    // ------------------------------------------------------------------------------ //
    //                                        SL                                      //
    //            p p > t t~, (t > b w, w > j j ), (t~ > X j, X > mu ta)              //
    // ------------------------------------------------------------------------------ //
    else if (isSL()) {
      //////////////////////////////////////////////////////
      ////////// variables from one leg (t/t~) /////////////
      //////////////////////////////////////////////////////
      float W_mass_nom = 80.379;
      std::vector<Jet>WJetsPair;
      wjetsSelector(lightJetColl, WJetsPair);

      if (WJetsPair.size() == 0) continue;
      AnaUtil::fillHist1D("evtCutFlow", 7, evWt);
      AnaUtil::fillHist1D("evtCutFlowWt", 7, lumiWt);
      TLorentzVector wj1p4 = WJetsPair[0].P4();
      TLorentzVector wj2p4 = WJetsPair[1].P4();
      TLorentzVector wp4   = wj1p4+wj2p4;
      float massdiff = wp4.M() - W_mass_nom;
      AnaUtil::fillHist1D("w_massdiff", massdiff);

      if (std::fabs(massdiff) > 30.0) continue;
      AnaUtil::fillHist1D("evtCutFlow", 8, evWt);
      AnaUtil::fillHist1D("evtCutFlowWt", 8, lumiWt);

      TLorentzVector t1p4   = wp4 + leadbjp4;

      float InvM_Wj1Wj2     = wp4.M();
      float DR_Wj1Wj2       = wj1p4.DeltaR(wj2p4);
      float DPhi_Wj1Wj2     = TVector2::Phi_mpi_pi(wj1p4.Phi() - wj2p4.Phi());
      float DEta_Wj1Wj2     = wj1p4.Eta() - wj2p4.Eta();

      float DR_w_leadbj     = leadbjp4.DeltaR(wp4);
      float DPhi_w_leadbj   = TVector2::Phi_mpi_pi(wp4.Phi() - leadbjp4.Phi()); 
      float InvM_w_leadbj   = t1p4.M();

      float DR_wj1_leadbj   = leadbjp4.DeltaR(wj1p4);
      float DPhi_wj1_leadbj = TVector2::Phi_mpi_pi(wj1p4.Phi() - leadbjp4.Phi()); 

      float DR_wj2_leadbj   = leadbjp4.DeltaR(wj2p4);
      float DPhi_wj2_leadbj = TVector2::Phi_mpi_pi(wj2p4.Phi() - leadbjp4.Phi()); 

      AnaUtil::fillHist1D("Wj1Pt",          wj1p4.Pt());
      AnaUtil::fillHist1D("Wj2Pt",          wj2p4.Pt());
      AnaUtil::fillHist1D("WPt",            wp4.Pt());
      AnaUtil::fillHist1D("InvM_Wj1Wj2",    InvM_Wj1Wj2);
      AnaUtil::fillHist1D("DR_Wj1Wj2",      DR_Wj1Wj2);
      AnaUtil::fillHist1D("DPhi_Wj1Wj2",    DPhi_Wj1Wj2);
      AnaUtil::fillHist1D("DEta_Wj1Wj2",    DEta_Wj1Wj2);
      AnaUtil::fillHist1D("DR_w_leadb",     DR_w_leadbj);
      AnaUtil::fillHist1D("DPhi_w_leadb",   DPhi_w_leadbj);
      AnaUtil::fillHist1D("DR_wj1_leadb",   DR_wj1_leadbj);
      AnaUtil::fillHist1D("DPhi_wj1_leadb", DPhi_wj1_leadbj);
      AnaUtil::fillHist1D("DR_wj2_leadb",   DR_wj2_leadbj);
      AnaUtil::fillHist1D("DPhi_wj2_leadb", DPhi_wj2_leadbj);
      AnaUtil::fillHist1D("InvM_w_leadbj",  InvM_w_leadbj);
      
      //////////////////////////////////////////////////////
      ////////// variables from other leg (t~/t) ///////////
      //////////////////////////////////////////////////////
      TLorentzVector XLepP4  = MuonColl[0].P4();
      float XlepPt           = XLepP4.Pt();
      float XlepEta          = XLepP4.Eta();
      float M_coll_Xlep      = (XLepP4 + tauhp4).M()/TMath::Sqrt(x_vis); // Collinear mass formation
      TLorentzVector Xp4     = taup4 + XLepP4;
      TLorentzVector t2p4    = Xp4 + leadljp4;
      float M_coll_test      = Xp4.M();
      float InvM_Xleadlj     = t2p4.M(); 

      float DR_Xlep_leadlj   = XLepP4.DeltaR(leadljp4);
      float DPhi_Xlep_leadlj = TVector2::Phi_mpi_pi(XLepP4.Phi() - leadljp4.Phi()); 

      float DR_Xtauh_leadlj  = tauhp4.DeltaR(leadljp4);
      float DPhi_Xtauh_leadlj= TVector2::Phi_mpi_pi(tauhp4.Phi() - leadljp4.Phi()); 

      float DR_X_leadlj      = Xp4.DeltaR(leadljp4);
      float DPhi_X_leadlj    = TVector2::Phi_mpi_pi(Xp4.Phi() - leadljp4.Phi());

      float DR_Xlep_Xtau     = XLepP4.DeltaR(tauhp4);
      float DPhi_Xlep_Xtau   = TVector2::Phi_mpi_pi(XLepP4.Phi() - tauhp4.Phi());

      AnaUtil::fillHist1D("XlepPt",        XlepPt);
      AnaUtil::fillHist1D("XlepEta",       XlepEta);
      AnaUtil::fillHist1D("M_coll_Xlep",   M_coll_Xlep);
      AnaUtil::fillHist1D("M_coll_test",   M_coll_test);
      AnaUtil::fillHist1D("InvM_Xleadlj",  InvM_Xleadlj);
      AnaUtil::fillHist1D("DR_Xlep_leadlj",    DR_Xlep_leadlj);
      AnaUtil::fillHist1D("DPhi_Xlep_leadlj",  DPhi_Xlep_leadlj);
      AnaUtil::fillHist1D("DR_Xtauh_leadlj",   DR_Xtauh_leadlj);
      AnaUtil::fillHist1D("DPhi_Xtauh_leadlj", DPhi_Xtauh_leadlj);
      AnaUtil::fillHist1D("DR_X_leadlj",    DR_X_leadlj);
      AnaUtil::fillHist1D("DPhi_X_leadlj",  DPhi_X_leadlj);
      AnaUtil::fillHist1D("DR_Xlep_Xtau",   DR_Xlep_Xtau);
      AnaUtil::fillHist1D("DPhi_Xlep_Xtau", DPhi_Xlep_Xtau);


      //////////////////////////////////////////////////////
      /////////// communication between 2 legs /////////////
      //////////////////////////////////////////////////////

      // Xlep-w
      float DR_Xlep_wj1       = XLepP4.DeltaR(wj1p4);
      float DPhi_Xlep_wj1     = TVector2::Phi_mpi_pi(XLepP4.Phi() - wj1p4.Phi());
      float DEta_Xlep_wj1     = XLepP4.Eta() - wj1p4.Eta();
      float DR_Xlep_wj2       = XLepP4.DeltaR(wj2p4);
      float DPhi_Xlep_wj2     = TVector2::Phi_mpi_pi(XLepP4.Phi() - wj2p4.Phi());
      float DEta_Xlep_wj2     = XLepP4.Eta() - wj2p4.Eta();
      float DR_Xlep_w         = XLepP4.DeltaR(wp4);
      float DPhi_Xlep_w       = TVector2::Phi_mpi_pi(XLepP4.Phi() - wp4.Phi());
      float DEta_Xlep_w       = XLepP4.Eta() - wp4.Eta();
      // Xlep-leadb
      float DR_Xlep_leadbj    = XLepP4.DeltaR(leadbjp4);
      float DPhi_Xlep_leadbj  = TVector2::Phi_mpi_pi(XLepP4.Phi() - leadbjp4.Phi());
      float DEta_Xlep_leadbj  = XLepP4.Eta() - leadbjp4.Eta();
      // Xtau-w
      float DR_Xtau_wj1       = tauhp4.DeltaR(wj1p4);
      float DPhi_Xtau_wj1     = TVector2::Phi_mpi_pi(tauhp4.Phi() - wj1p4.Phi());
      float DEta_Xtau_wj1     = tauhp4.Eta() - wj1p4.Eta();
      float DR_Xtau_wj2       = tauhp4.DeltaR(wj2p4);
      float DPhi_Xtau_wj2     = TVector2::Phi_mpi_pi(tauhp4.Phi() - wj2p4.Phi());
      float DEta_Xtau_wj2     = tauhp4.Eta() - wj2p4.Eta();
      float DR_Xtau_w         = tauhp4.DeltaR(wp4);
      float DPhi_Xtau_w       = TVector2::Phi_mpi_pi(tauhp4.Phi() - wp4.Phi());
      float DEta_Xtau_w       = tauhp4.Eta() - wp4.Eta();
      // Xtau-leadb
      float DR_Xtau_leadbj         = tauhp4.DeltaR(leadbjp4);
      float DPhi_Xtau_leadbj       = TVector2::Phi_mpi_pi(tauhp4.Phi() - leadbjp4.Phi());
      float DEta_Xtau_leadbj       = tauhp4.Eta() - leadbjp4.Eta();
      // leadb-leadl
      //float DR_leadbj_leadlj       = leadbjp4.DeltaR(leadljp4);
      //float DPhi_leadbj_leadlj     = TVector2::Phi_mpi_pi(leadbjp4.Phi() - leadljp4.Phi());
      //float DEta_leadbj_leadlj     = leadbjp4.Eta() - leadljp4.Eta();
      // W-leadl
      float DR_w_leadlj         = wp4.DeltaR(leadljp4);
      float DPhi_w_leadlj       = TVector2::Phi_mpi_pi(wp4.Phi() - leadljp4.Phi());
      float DEta_w_leadlj       = wp4.Eta() - leadljp4.Eta();
      // Chi-leadb
      float DR_X_leadbj            = Xp4.DeltaR(leadbjp4);
      float DPhi_X_leadbj          = TVector2::Phi_mpi_pi(Xp4.Phi() - leadbjp4.Phi());
      float DEta_X_leadbj          = Xp4.Eta() - leadbjp4.Eta();
      // Chi-W
      float DR_X_w            = Xp4.DeltaR(wp4);
      float DPhi_X_w          = TVector2::Phi_mpi_pi(Xp4.Phi() - wp4.Phi());
      float DEta_X_w          = Xp4.Eta() - wp4.Eta();
      // t-t~
      float DR_t1_t2          = t1p4.DeltaR(t2p4);
      float DPhi_t1_t2        = TVector2::Phi_mpi_pi(t1p4.Phi() - t2p4.Phi());
      float DEta_t1_t2        = t1p4.Eta() - t2p4.Eta();
      // total
      float Total_mass        = (t1p4 + t2p4).M();
      float Total_vector_pt   = (t1p4 + t2p4).Pt();
      float Total_scalar_pt   = XlepPt + taup4.Pt() + leadljp4.Pt() + leadbjp4.Pt() + wj1p4.Pt() + wj2p4.Pt();

      AnaUtil::fillHist1D("DR_Xlep_wj1",             DR_Xlep_wj1);       
      AnaUtil::fillHist1D("DPhi_Xlep_wj1",           DPhi_Xlep_wj1);     
      AnaUtil::fillHist1D("DEta_Xlep_wj1",           DEta_Xlep_wj1);     
      AnaUtil::fillHist1D("DR_Xlep_wj2",             DR_Xlep_wj2);
      AnaUtil::fillHist1D("DPhi_Xlep_wj2",           DPhi_Xlep_wj2);     
      AnaUtil::fillHist1D("DEta_Xlep_wj2",           DEta_Xlep_wj2);     
      AnaUtil::fillHist1D("DR_Xlep_w",               DR_Xlep_w);
      AnaUtil::fillHist1D("DPhi_Xlep_w",             DPhi_Xlep_w);       
      AnaUtil::fillHist1D("DEta_Xlep_w",             DEta_Xlep_w);       
      AnaUtil::fillHist1D("DR_Xlep_leadbj",          DR_Xlep_leadbj);
      AnaUtil::fillHist1D("DPhi_Xlep_leadbj",        DPhi_Xlep_leadbj);       
      AnaUtil::fillHist1D("DEta_Xlep_leadbj",        DEta_Xlep_leadbj);
      AnaUtil::fillHist1D("DR_Xtau_wj1",             DR_Xtau_wj1);
      AnaUtil::fillHist1D("DPhi_Xtau_wj1",           DPhi_Xtau_wj1);
      AnaUtil::fillHist1D("DEta_Xtau_wj1",           DEta_Xtau_wj1);
      AnaUtil::fillHist1D("DR_Xtau_wj2",             DR_Xtau_wj2);
      AnaUtil::fillHist1D("DPhi_Xtau_wj2",           DPhi_Xtau_wj2);
      AnaUtil::fillHist1D("DEta_Xtau_wj2",           DEta_Xtau_wj2);
      AnaUtil::fillHist1D("DR_Xtau_w",               DR_Xtau_w);
      AnaUtil::fillHist1D("DPhi_Xtau_w",             DPhi_Xtau_w);
      AnaUtil::fillHist1D("DEta_Xtau_w",             DEta_Xtau_w);
      AnaUtil::fillHist1D("DR_Xtau_leadbj",          DR_Xtau_leadbj);
      AnaUtil::fillHist1D("DPhi_Xtau_leadbj",        DPhi_Xtau_leadbj);
      AnaUtil::fillHist1D("DEta_Xtau_leadbj",        DEta_Xtau_leadbj);
      AnaUtil::fillHist1D("DR_leadbj_leadlj",        DR_leadbj_leadlj);
      AnaUtil::fillHist1D("DPhi_leadbj_leadlj",      DPhi_leadbj_leadlj);
      AnaUtil::fillHist1D("DEta_leadbj_leadlj",      DEta_leadbj_leadlj);
      AnaUtil::fillHist1D("DR_w_leadlj",             DR_w_leadlj);
      AnaUtil::fillHist1D("DPhi_w_leadlj",           DPhi_w_leadlj);
      AnaUtil::fillHist1D("DEta_w_leadlj",           DEta_w_leadlj);
      AnaUtil::fillHist1D("DR_X_leadbj",             DR_X_leadbj);
      AnaUtil::fillHist1D("DPhi_X_leadbj",           DPhi_X_leadbj); 
      AnaUtil::fillHist1D("DEta_X_leadbj",           DEta_X_leadbj);
      AnaUtil::fillHist1D("DR_X_w",                  DR_X_w);
      AnaUtil::fillHist1D("DPhi_X_w",                DPhi_X_w); 
      AnaUtil::fillHist1D("DEta_X_w",                DEta_X_w);
      AnaUtil::fillHist1D("DR_t1_t2",                DR_t1_t2);
      AnaUtil::fillHist1D("DPhi_t1_t2",              DPhi_t1_t2);
      AnaUtil::fillHist1D("DEta_t1_t2",              DEta_t1_t2);
      AnaUtil::fillHist1D("Total_mass",              Total_mass);
      AnaUtil::fillHist1D("Total_vector_pt",         Total_vector_pt);
      AnaUtil::fillHist1D("Total_scalar_pt",         Total_scalar_pt);

      //////////////////////////////////////////////////////
      ////////// A few more high level variables ///////////
      //////////////////////////////////////////////////////
      float HT_jets = 0.0;
      double minDr_XlepJets   = 999.9;
      double maxDr_XlepJets   = -999.9;
      double minDphi_XlepJets = 999.9;
      double maxDphi_XlepJets = -999.9;
      for (auto &jet : JetColl) {
	HT_jets += jet.PT;
	double dr_jetXlep   = XLepP4.DeltaR(jet.P4());
	double dphi_jetXlep = std::fabs(TVector2::Phi_mpi_pi(jet.P4().Phi() - XLepP4.Phi())); 
	if (dr_jetXlep <= minDr_XlepJets) minDr_XlepJets = dr_jetXlep;
	if (dr_jetXlep > maxDr_XlepJets)  maxDr_XlepJets = dr_jetXlep;
	if (dphi_jetXlep <= minDphi_XlepJets) minDphi_XlepJets = dphi_jetXlep;
	if (dphi_jetXlep > maxDphi_XlepJets)  maxDphi_XlepJets = dphi_jetXlep;
      }
      AnaUtil::fillHist1D ("HT_Jets",          HT_jets);
      AnaUtil::fillHist1D ("minDr_XlepJets",   minDr_XlepJets);
      AnaUtil::fillHist1D ("maxDr_XlepJets",   maxDr_XlepJets);
      AnaUtil::fillHist1D ("minDphi_XlepJets", minDphi_XlepJets);
      AnaUtil::fillHist1D ("maxDphi_XlepJets", maxDphi_XlepJets);

      float smin = comp_smin(XLepP4, tauhp4, leadbjp4, leadljp4, wj1p4, wj2p4, metp4);
      AnaUtil::fillHist1D ("smin", smin);

      // needs to be done
      if (skimObj_) {
	TreeVariablesSL varList;
	
	varList.event         = iEntry;
	varList.event_wt      = event_weight;
	// numbers of objects in final state
	varList.nleptons      = nGoodMuon+nGoodEle;
	varList.njets         = nGoodJet;
	varList.nbjets        = nGoodBJet;
	varList.nljets        = nGoodlJet;
	varList.ntauh         = nGoodTauh;
	// --------- Low level variables ---------- //
	// ==> Tauh vars
	varList.px_tauh       = tauhp4.Px();
	varList.py_tauh       = tauhp4.Py();
	varList.pz_tauh       = tauhp4.Pz();
	varList.energy_tauh   = tauhp4.E();
	varList.pt_tauh       = tauhp4.Pt();
	varList.eta_tauh      = std::fabs(tauhp4.Eta());
	varList.phi_tauh      = std::fabs(tauhp4.Phi());
	// ==> MET vars
	varList.px_met        = metp4.Px();
	varList.py_met        = metp4.Py();
	varList.pz_met        = metp4.Pz();
	varList.energy_met    = metp4.E();
	varList.pt_met        = met_->MET;
	varList.phi_met       = metp4.Phi();
	// ==> bjet vars
	varList.px_leadbj     = leadbjp4.Px();
	varList.py_leadbj     = leadbjp4.Py();
	varList.pz_leadbj     = leadbjp4.Pz();
	varList.energy_leadbj = leadbjp4.E();
	varList.pt_leadbj     = leadbjp4.Pt();
	varList.eta_leadbj    = std::fabs(leadbjp4.Eta());
	varList.phi_leadbj    = std::fabs(leadbjp4.Phi());
	// ==> leading light jet vars
	varList.px_leadlj     = leadljp4.Px();
	varList.py_leadlj     = leadljp4.Py();
	varList.pz_leadlj     = leadljp4.Pz();
	varList.energy_leadlj = leadljp4.E();
	varList.pt_leadlj     = leadljp4.Pt();
	varList.eta_leadlj    = std::fabs(leadljp4.Eta());
	varList.phi_leadlj    = std::fabs(leadljp4.Phi());
	// ==> w jet1 vars
	varList.px_wjet1      = wj1p4.Px();
	varList.py_wjet1      = wj1p4.Py();
	varList.pz_wjet1      = wj1p4.Pz();
	varList.energy_wjet1  = wj1p4.E();
	varList.pt_wjet1      = wj1p4.Pt();
	varList.eta_wjet1     = std::fabs(wj1p4.Eta());
	varList.phi_wjet1     = std::fabs(wj1p4.Phi());
	// ==> w jet2 vars
	varList.px_wjet2      = wj2p4.Px();
	varList.py_wjet2      = wj2p4.Py();
	varList.pz_wjet2      = wj2p4.Pz();
	varList.energy_wjet2  = wj2p4.E();
	varList.pt_wjet2      = wj2p4.Pt();
	varList.eta_wjet2     = std::fabs(wj2p4.Eta());
	varList.phi_wjet2     = std::fabs(wj2p4.Phi());
	// ==> lepton from Chi
	varList.px_xlep       = XLepP4.Px();
	varList.py_xlep       = XLepP4.Py();
	varList.pz_xlep       = XLepP4.Pz();
	varList.energy_xlep   = XLepP4.E();
	varList.pt_xlep       = XlepPt;
	varList.eta_xlep      = std::fabs(XlepEta);
	varList.phi_xlep      = std::fabs(XLepP4.Phi());

	// --------- High level variables ---------- //
	// from one leg
	varList.invm_wj1_wj2    = InvM_Wj1Wj2;
	varList.dr_wj1_wj2      = DR_Wj1Wj2;
	varList.dphi_wj1_wj2    = std::fabs(DPhi_Wj1Wj2);
	varList.deta_wj1_wj2    = std::fabs(DEta_Wj1Wj2);
	varList.dr_w_leadbj     = DR_w_leadbj;
	varList.dphi_w_leadbj   = std::fabs(DPhi_w_leadbj);
	varList.dr_wj1_leadbj   = DR_wj1_leadbj;
	varList.dphi_wj1_leadbj = std::fabs(DPhi_wj1_leadbj);
	varList.dr_wj2_leadbj   = DR_wj2_leadbj;
	varList.dphi_wj2_leadbj = std::fabs(DPhi_wj2_leadbj);
	varList.invm_w_leadbj   = InvM_w_leadbj;
	// from other leg
	varList.m_coll_x          = M_coll_Xlep;
	varList.m_coll_xtest      = M_coll_test;
	varList.invm_x_leadlj     = InvM_Xleadlj;
	varList.dr_xlep_leadlj    = DR_Xlep_leadlj;
	varList.dphi_xlep_leadlj  = std::fabs(DPhi_Xlep_leadlj);
	varList.dr_tauh_leadlj    = DR_Xtauh_leadlj;
	varList.dphi_tauh_leadlj  = std::fabs(DPhi_Xtauh_leadlj);
	varList.dr_x_leadlj       = DR_X_leadlj;
	varList.dphi_x_leadlj     = std::fabs(DPhi_X_leadlj);
	varList.dr_xlep_tauh      = DR_Xlep_Xtau;
	varList.dphi_xlep_tauh    = std::fabs(DPhi_Xlep_Xtau);
	// from both leg
	varList.dr_xlep_wj1      = DR_Xlep_wj1;
	varList.dphi_xlep_wj1    = std::fabs(DPhi_Xlep_wj1);
	varList.deta_xlep_wj1    = std::fabs(DEta_Xlep_wj1);
	varList.dr_xlep_wj2      = DR_Xlep_wj2;
	varList.dphi_xlep_wj2    = std::fabs(DPhi_Xlep_wj2);
	varList.deta_xlep_wj2    = std::fabs(DEta_Xlep_wj2);
	varList.dr_xlep_w        = DR_Xlep_w;
	varList.dphi_xlep_w      = std::fabs(DPhi_Xlep_w);
	varList.deta_xlep_w      = std::fabs(DEta_Xlep_w);
	varList.dr_xlep_leadbj   = DR_Xlep_leadbj;
	varList.dphi_xlep_leadbj = std::fabs(DPhi_Xlep_leadbj);
	varList.deta_xlep_leadbj = std::fabs(DEta_Xlep_leadbj);
	varList.dr_tauh_wj1      = DR_Xtau_wj1;
	varList.dphi_tauh_wj1    = std::fabs(DPhi_Xtau_wj1);
	varList.deta_tauh_wj1    = std::fabs(DEta_Xtau_wj1);
	varList.dr_tauh_wj2      = DR_Xtau_wj2;
	varList.dphi_tauh_wj2    = std::fabs(DPhi_Xtau_wj2);
	varList.deta_tauh_wj2    = std::fabs(DEta_Xtau_wj2);
	varList.dr_tauh_leadbj    = DR_Xtau_leadbj;
	varList.dphi_tauh_leadbj  = std::fabs(DPhi_Xtau_leadbj);
	varList.deta_tauh_leadbj  = std::fabs(DEta_Xtau_leadbj);
	varList.dr_leadbj_leadlj  = DR_leadbj_leadlj;
	varList.dphi_leadbj_leadlj  = std::fabs(DPhi_leadbj_leadlj);
	varList.deta_leadbj_leadlj  = std::fabs(DEta_leadbj_leadlj);
	varList.dr_w_leadlj       = DR_w_leadlj;
	varList.dphi_w_leadlj     = std::fabs(DPhi_w_leadlj);
	varList.deta_w_leadlj     = std::fabs(DEta_w_leadlj);
	varList.dr_x_leadbj       = DR_X_leadbj;
	varList.dphi_x_leadbj     = std::fabs(DPhi_X_leadbj);
	varList.deta_x_leadbj     = std::fabs(DEta_X_leadbj);
	varList.dr_x_w           = DR_X_w; 
	varList.dphi_x_w         = std::fabs(DPhi_X_w);
	varList.deta_x_w         = std::fabs(DEta_X_w);
	varList.dr_t1_t2         = DR_t1_t2; 
	varList.dphi_t1_t2       = std::fabs(DPhi_t1_t2);
	varList.deta_t1_t2       = std::fabs(DEta_t1_t2);
	varList.total_mass       = Total_mass; 
	varList.total_vector_pt  = Total_vector_pt;
	varList.total_scalar_pt  = Total_scalar_pt;
	varList.dphi_xlep_met    = std::fabs(TVector2::Phi_mpi_pi(XLepP4.Phi() - metp4.Phi()));
	varList.mt_xlep_met      = Calculate_MT(XLepP4, metp4);
	// di-jets variables
	varList.dr_min_jets      = minDR;
	varList.dr_max_jets      = maxDR;
	varList.ht_jets          = HT_jets;
	varList.smin             = smin;
	
	skimObj_->fill(varList);
      }       
      histf->cd();

      double mvaOut = -999.9;
      
      if (_readMVA) {
	InputVariablesSL varList;
	
	varList.pt_tauh              = TauhColl[0].PT;
	varList.met                  = met_->MET;
	
	mvaOut = _mvaObj->evaluate(_MVAnetwork, varList);
      }
      histf->cd();
      
      AnaUtil::fillHist1D("BDT_Score", mvaOut);
      AnaUtil::fillHist1D("BDT_Score_fineBin", mvaOut);
    }

    else 
      std::cerr << "Wrong channel. Please mention either SL or DL through jobcard"
		<< std::endl;
  
  } // Event loop ends
} // function scope ends

void ExoAnalysis::endJob() {
  
  histf->cd();
  //histf->cd("Analysis");
  vector<string> evLabels;
  if (isDL()) 
    evLabels = 
      {
	"Events processed",
	"pass obj selection",
	"nMu + nEle == 2",
	"nMu >= 1",
	"nTau >= 1",
	"NLightJet >= 1",
	"nBJets == 1",
	"no Z",
	"Xlep-Wlep clean [NA]",
	"Xlep-tauh opp charge"
      };
  else if (isSL())
    evLabels =
      {
        "Events processed",
	"pass obj selection",
	"nMu = 1",
	"nEle = 0",
	"nTau = 1",
	"NLightJet >= 3",
	"nBJets == 1",
	"has probable 2 W jets",
	"|w reco mass - mW| < 40 GeV"
      };
  
  AnaUtil::SetEvtCutFlowBinLabels("evtCutFlow",evLabels);
  AnaUtil::showEfficiency("evtCutFlow", evLabels, "Event Selection");
  AnaUtil::showEfficiency("evtCutFlowWt", evLabels, "Event Selection : Lumi Scaled");  
}

void ExoAnalysis::wjetsSelector(const std::vector <Jet>& alljetList_, std::vector <Jet>& wjetspair) {
  float Wmass = 80.379;
  float minMassDiff = 9999.9;
  unsigned int j1_idx = 0;
  unsigned int j2_idx = 0;
  // 0th light jet comes from top
  for (unsigned int i = 1; i < alljetList_.size(); ++i) {
    const auto& ip = alljetList_[i];
    TLorentzVector j1p4 = ip.P4();
    for (unsigned int j = i+1; j < alljetList_.size(); ++j) {
      const auto& jp = alljetList_[j];
      TLorentzVector j2p4 = jp.P4();      
      float invM = (j1p4 + j2p4).M();
      float diff = std::fabs(Wmass - invM);
      if (diff < minMassDiff) {
	minMassDiff = diff;
	j1_idx = i;
	j2_idx = j;
      }
    }
  }
  wjetspair.push_back(alljetList_[j1_idx]);
  wjetspair.push_back(alljetList_[j2_idx]);
}

float ExoAnalysis::comp_smin(const TLorentzVector &lep1p4, 
			     const TLorentzVector &tauh1p4,
			     const TLorentzVector &bj1p4,
			     const TLorentzVector &lj1p4,
			     const TLorentzVector &wj1p4, 
			     const TLorentzVector &wj2p4,
			     const TLorentzVector &metp4) {
  TLorentzVector visp4 = lep1p4 + tauh1p4 + bj1p4 + wj1p4 + wj2p4 + lj1p4;
  float vispt = visp4.Pt();
  float vism  = visp4.M();
  float viset = TMath::Sqrt(std::pow(vism, 2) + std::pow(vispt, 2));
  float metet = metp4.E();
  return TMath::Sqrt(std::pow(vism, 2) + 2 * (viset * metet - (visp4.Px() * metp4.Px() + visp4.Py() * metp4.Py())));
}

float ExoAnalysis::comp_smin(const TLorentzVector &lep1p4,
			     const TLorentzVector &lep2p4,
			     const TLorentzVector &tauh1p4,
			     const TLorentzVector &bj1p4, 
			     const TLorentzVector &lj1p4,
			     const TLorentzVector &metp4) {
  TLorentzVector visp4 = lep1p4 + lep2p4 + tauh1p4 + bj1p4 + lj1p4;
  float vispt = visp4.Pt();
  float vism  = visp4.M();
  float viset = TMath::Sqrt(std::pow(vism, 2) + std::pow(vispt, 2));
  float metet = metp4.E();
  return TMath::Sqrt(std::pow(vism, 2) + 2 * (viset * metet - (visp4.Px() * metp4.Px() + visp4.Py() * metp4.Py())));
}

void ExoAnalysis::closeHistFiles(){
  histf->cd();
  histf->Write();
  histf->Close();
}
void ExoAnalysis::closeFiles(){
  closeHistFiles();
  //if (skimObj_ != nullptr) skimObj_->close();
  if (skimObj_ != nullptr) skimObj_->close(isDL(), isSL());
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
   else if (key == "isDL")
      _isDL = (atoi(value.c_str()) > 0) ? true : false;
   else if (key == "isSL")
      _isSL = (atoi(value.c_str()) > 0) ? true : false;
   else if (key == "highXmass")
      _highXmass = (atoi(value.c_str()) > 0) ? true : false;
   else if (key == "lowXmass")
      _lowXmass = (atoi(value.c_str()) > 0) ? true : false;
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
     //if (Th.TauTag == 0) continue;
     //AnaUtil::fillHist1D ("jetThDR", jetP4.DeltaR(Th.P4()), 1.0);
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
