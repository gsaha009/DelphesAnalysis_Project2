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

  reader_->AddVariable("met", &varList_.met);
  reader_->AddVariable("hLepPt", &varList_.hLepPt);
  reader_->AddVariable("X_tauhPt", &varList_.X_tauhPt);
  reader_->AddVariable("DR_Z_l1l2", &varList_.DR_Z_l1l2);
  reader_->AddVariable("DR_Xmu_Xtau", &varList_.DR_Xmu_Xtau);
  reader_->AddVariable("X_muMetMT", &varList_.X_muMetMT);
  reader_->AddVariable("DPhi_Z_Xmu", &varList_.DPhi_Z_Xmu);
  reader_->AddVariable("XHelicity_Angle", &varList_.XHelicity_Angle);
  reader_->AddVariable("XZPlane_Angle", &varList_.XZPlane_Angle);

  std::cout<<"XML: "<<xmlfile<<"\n";
  std::cout<<"algo: "<<mva_algo<<"\n";
  reader_->BookMVA(mva_algo.c_str(), xmlfile);
}
double MVAnalysis::evaluate(const string& mva_algo, const InputVariables& varList) {
  memcpy(&varList_, &varList, sizeof(varList)); // use move syntax here
  return reader_->EvaluateMVA(mva_algo.c_str());
}