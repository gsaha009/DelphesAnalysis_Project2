#ifndef __MVAnalysis__h
#define __MVAnalysis__h

#include <fstream>
#include <string>

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

typedef struct  
{
  float pt_tauh;
  float met;
  
} InputVariablesDL;

typedef struct  
{
  float pt_tauh;
  float met;
  
} InputVariablesSL;

class MVAnalysis {
    
public:

  MVAnalysis(const std::string& mva_algo, const std::string& xmlfile);
  virtual ~MVAnalysis() {}

  double evaluate(const std::string& tag, const InputVariablesDL& varList);
  double evaluate(const std::string& tag, const InputVariablesSL& varList);

  InputVariablesDL varListDL_;
  InputVariablesSL varListSL_;
  //  std::unique_ptr<TMVA::Reader> reader_;
  TMVA::Reader* readerDL_;
  TMVA::Reader* readerSL_;
};
#endif
