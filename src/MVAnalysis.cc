#include <iostream>
#include <memory>
#include "MVAnalysis.h"

using std::string;
using std::cout;
using std::endl;

MVAnalysis::MVAnalysis(const string& mva_algo, const string& xmlfile) 
{
  //  reader_ = std::make_unique<TMVA::Reader>(new TMVA::Reader("!Silent"));
  reader_ = new TMVA::Reader("!Silent");

  reader_->AddVariable("pt_tauh1",             &varList_.pt_tauh1);
  reader_->AddVariable("met",                  &varList_.met);
  reader_->AddVariable("mt_Wlepmet",           &varList_.mt_Wlepmet);
  reader_->AddVariable("pt_bjet1",             &varList_.pt_bjet1);
  reader_->AddVariable("dphi_Wleptauh",        &varList_.dphi_Wleptauh);
  reader_->AddVariable("dphi_Xleptauh",        &varList_.dphi_Xleptauh);
  reader_->AddVariable("pt_Xlep",              &varList_.pt_Xlep);
  reader_->AddVariable("pt_Wlep",              &varList_.pt_Wlep);
  reader_->AddVariable("vectorSumpt_XlwpWlep", &varList_.vectorSumpt_XlwpWlep);
  reader_->AddVariable("dr_XlepWlep",          &varList_.dr_XlepWlep);
  reader_->AddVariable("dphi_tauhbjet",        &varList_.dphi_tauhbjet);
  reader_->AddVariable("dr_min_Xlepjets",      &varList_.dr_min_Xlepjets);
  reader_->AddVariable("dr_min_Wlepjets",      &varList_.dr_min_Wlepjets);
  reader_->AddVariable("deta_Wlepbjet",        &varList_.deta_Wlepbjet);
  reader_->AddVariable("dphi_Xlepmet",         &varList_.dphi_Xlepmet);
  reader_->AddVariable("dphi_Wlepmet",         &varList_.dphi_Wlepmet);
  reader_->AddVariable("dr_min_jets",          &varList_.dr_min_jets);
  reader_->AddVariable("dphi_bjetljet",        &varList_.dphi_bjetljet);
  reader_->AddVariable("effectiveMass",        &varList_.effectiveMass);
  reader_->AddVariable("HT_Jets",              &varList_.HT_Jets);
  
  
  std::cout<<"XML: "<<xmlfile<<"\n";
  std::cout<<"algo: "<<mva_algo<<"\n";
  reader_->BookMVA(mva_algo.c_str(), xmlfile);
}
double MVAnalysis::evaluate(const string& mva_algo, const InputVariables& varList) {
  memcpy(&varList_, &varList, sizeof(varList)); // use move syntax here
  return reader_->EvaluateMVA(mva_algo.c_str());
}
