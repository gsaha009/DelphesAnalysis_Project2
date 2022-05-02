#include <iostream>
#include <memory>
#include "MVAnalysis.h"

using std::string;
using std::cout;
using std::endl;

MVAnalysis::MVAnalysis(const string& mva_algo, const string& xmlfile, bool SL, bool DL) 
{
  //  reader_ = std::make_unique<TMVA::Reader>(new TMVA::Reader("!Silent"));
  if (DL) {
    readerDL_ = new TMVA::Reader("!Silent");
    /*
    readerDL_->AddVariable("dphi_met_tauh",       &varListDL_.dphi_met_tauh);
    readerDL_->AddVariable("dphi_wlep_leadbj",    &varListDL_.dphi_wlep_leadbj);
    readerDL_->AddVariable("dphi_xlep_wlep",      &varListDL_.dphi_xlep_wlep);
    readerDL_->AddVariable("dr_xlep_tauh",        &varListDL_.dr_xlep_tauh);
    readerDL_->AddVariable("effectivemass",       &varListDL_.effectivemass);
    readerDL_->AddVariable("ht_jets",             &varListDL_.ht_jets);
    readerDL_->AddVariable("mt_wlep_met",         &varListDL_.mt_wlep_met);
    readerDL_->AddVariable("mt_x",                &varListDL_.mt_x);
    readerDL_->AddVariable("met",                 &varListDL_.met);
    readerDL_->AddVariable("scalarsumpt",         &varListDL_.scalarsumpt);
    readerDL_->AddVariable("dphi_leadbj_leadlj",  &varListDL_.dphi_leadbj_leadlj);
    readerDL_->AddVariable("m_coll_xlep",         &varListDL_.m_coll_xlep);
    readerDL_->AddVariable("m_coll_wlep",         &varListDL_.m_coll_wlep);
    */

    readerDL_->AddVariable("pt_tauh",                 &varListDL_.pt_tauh);
    readerDL_->AddVariable("met",                     &varListDL_.met); 
    readerDL_->AddVariable("pt_leadbj",               &varListDL_.pt_leadbj);
    readerDL_->AddVariable("pt_xlep",                 &varListDL_.pt_xlep);
    readerDL_->AddVariable("pt_wlep",                 &varListDL_.pt_wlep);
    readerDL_->AddVariable("vectorsumpt_xlep_wlep",   &varListDL_.vectorsumpt_xlep_wlep);
    readerDL_->AddVariable("dr_xlep_wlep",            &varListDL_.dr_xlep_wlep);
    readerDL_->AddVariable("dphi_wlep_tauh",          &varListDL_.dphi_wlep_tauh);
    readerDL_->AddVariable("dphi_wlep_leadbj",        &varListDL_.dphi_wlep_leadbj);
    readerDL_->AddVariable("ht_jets",                 &varListDL_.ht_jets);
    readerDL_->AddVariable("dr_min_xlep_jets",        &varListDL_.dr_min_xlep_jets);
    readerDL_->AddVariable("dr_min_wlep_jets",        &varListDL_.dr_min_wlep_jets);
    readerDL_->AddVariable("dphi_xlep_tauh",          &varListDL_.dphi_xlep_tauh);
    readerDL_->AddVariable("effectivemass",           &varListDL_.effectivemass);
    readerDL_->AddVariable("dr_min_jets",             &varListDL_.dr_min_jets);
    readerDL_->AddVariable("mt_wlep_met",             &varListDL_.mt_wlep_met);
    readerDL_->AddVariable("dphi_tauh_leadbj",        &varListDL_.dphi_tauh_leadbj);
    readerDL_->AddVariable("dphi_met_wlep",           &varListDL_.dphi_met_wlep);
    readerDL_->AddVariable("dphi_leadbj_leadlj",      &varListDL_.dphi_leadbj_leadlj);
    
    
    std::cout<<"XML: "<<xmlfile<<"\n";
    std::cout<<"algo: "<<mva_algo<<"\n";
    readerDL_->BookMVA(mva_algo.c_str(), xmlfile);
  }

  else if (SL) {
    readerSL_ = new TMVA::Reader("!Silent");
    
    readerSL_->AddVariable("pt_tauh",             &varListSL_.pt_tauh);
    readerSL_->AddVariable("met",                 &varListSL_.met);
    
    readerSL_->BookMVA(mva_algo.c_str(), xmlfile);
  }
}

double MVAnalysis::evaluate(const string& mva_algo, const InputVariablesDL& varListDL) {
  memcpy(&varListDL_, &varListDL, sizeof(varListDL)); // use move syntax here
  return readerDL_->EvaluateMVA(mva_algo.c_str());
}
double MVAnalysis::evaluate(const string& mva_algo, const InputVariablesSL& varListSL) {
  memcpy(&varListSL_, &varListSL, sizeof(varListSL)); // use move syntax here
  return readerSL_->EvaluateMVA(mva_algo.c_str());
}
