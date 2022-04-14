#ifndef __HISTBOOKER__H
#define __HISTBOOKER__H

#include <string>
#include "AnaUtil.h"

class HistBooker {
 public: 
  HistBooker() {}
  virtual ~HistBooker() {}
  void bookHist1D(const char* hname, const char* htitle, 
		  int nbins, float xlow, float xhigh);
  void bookHist2D(const char* hname, const char* htitle, 
		  int nbinsx, float xlow, float xhigh, 
		  int nbinsy, float ylow, float yhigh);
  void bookBasicHistograms();
  void bookDLHistograms();
  void bookSLHistograms();
};
#endif
