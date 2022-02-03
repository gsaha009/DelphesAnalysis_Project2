#include <iostream>
#include <memory>
#include "MVAnalysis.h"

using std::string;
using std::cout;
using std::endl;

MVAnalysis::MVAnalysis(const string& mva_algo, const string& xmlfile) 
{
  //  reader_ = std::make_unique<TMVA::Reader>(new TMVA::Reader("!Silent"));
  readerDL_ = new TMVA::Reader("!Silent");

  readerDL_->AddVariable("pt_tauh",             &varListDL_.pt_tauh);
  readerDL_->AddVariable("met",                 &varListDL_.met);
  


  readerSL_ = new TMVA::Reader("!Silent");

  readerSL_->AddVariable("pt_tauh",             &varListSL_.pt_tauh);
  readerSL_->AddVariable("met",                 &varListSL_.met);
  
  std::cout<<"XML: "<<xmlfile<<"\n";
  std::cout<<"algo: "<<mva_algo<<"\n";
  readerDL_->BookMVA(mva_algo.c_str(), xmlfile);
  readerSL_->BookMVA(mva_algo.c_str(), xmlfile);
}

double MVAnalysis::evaluate(const string& mva_algo, const InputVariablesDL& varListDL) {
  memcpy(&varListDL_, &varListDL, sizeof(varListDL)); // use move syntax here
  return readerDL_->EvaluateMVA(mva_algo.c_str());
}
double MVAnalysis::evaluate(const string& mva_algo, const InputVariablesSL& varListSL) {
  memcpy(&varListSL_, &varListSL, sizeof(varListSL)); // use move syntax here
  return readerSL_->EvaluateMVA(mva_algo.c_str());
}
