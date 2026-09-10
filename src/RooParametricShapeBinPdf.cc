//---------------------------------------------------------------------------
#include "RooFit.h"

#include "Riostream.h"
#include <TMath.h>
#include <cassert>
#include <cmath>
#include <iostream>
#include <math.h>

#include "../interface/RooParametricShapeBinPdf.h"
#include "RooRealVar.h"
#include "RooRealVarSharedProperties.h"
#include "RooArgList.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"

using namespace std;
using namespace RooFit;

#if ROOT_VERSION_CODE < ROOT_VERSION(6,24,0)
// This (and its later uses) is a performance patch for older root versions
namespace {
  class RooRealVarSharedPropertiesSmart : public RooRealVarSharedProperties {
    public:
      void setHashTableSize(int size) { _altBinning.setHashTableSize(size); }
      int getHashTableSize() const { return _altBinning.getHashTableSize(); }
  };
  class RooRealVarSmart : public RooRealVar {
    public:
      void setHashTableSize(int size) { ((RooRealVarSharedPropertiesSmart*)sharedProp().get())->setHashTableSize(size); }
      int getHashTableSize() const { return ((RooRealVarSharedPropertiesSmart*)sharedProp().get())->getHashTableSize(); }
  };
}
#endif

ClassImp(RooParametricShapeBinPdf)
//---------------------------------------------------------------------------

RooParametricShapeBinPdf::RooParametricShapeBinPdf(const char *name, const char *title, RooAbsReal& _pdf, 
			       RooAbsReal& _x, const RooArgList& _pars, const TH1 &_shape ) : RooAbsPdf(name, title), 
  x("x", "x Observable", this, _x),
  pars("pars","pars",this),
  mypdf("mypdf","mypdf", this, _pdf),
  // underflowIntegral(nullptr),
  // overflowIntegral(nullptr),
  xBins(0),
  xMax(0),
  xMin(0)
  // xMinReal(x.min()),
  // xMaxReal(x.max())
{
  memset(&xArray, 0, sizeof(xArray));
  pars.add(_pars);
  setTH1Binning(_shape, /*resetRangeName=*/"");

}
//---------------------------------------------------------------------------
RooParametricShapeBinPdf::RooParametricShapeBinPdf(const RooParametricShapeBinPdf& other, const char* name) : RooAbsPdf(other, name), 
   x("x", this, other.x),
   pars("pars",this,RooListProxy()),
   mypdf("mypdf",this,other.mypdf),
   myintegrals(),
   // underflowIntegral(nullptr),
   // overflowIntegral(nullptr),
   xBins(other.xBins),
   xMax(other.xMax),
   xMin(other.xMin),
   // xMinReal(other.xMinReal),
   // xMaxReal(other.xMaxReal),
   extraRangeName(other.extraRangeName)
{
  //memset(&xArray, 0, sizeof(xArray));
  for (Int_t i=0; i<xBins+1; i++){
    xArray[i] = other.xArray[i];
  }
  
  pars.add(other.pars);
  // myintegrals = other.myintegrals;
}
//---------------------------------------------------------------------------
void RooParametricShapeBinPdf::setTH1Binning(const TH1 &_Hnominal, const char* resetRangeName) {
  myintegrals.clear();
  // underflowIntegral = nullptr;
  // overflowIntegral = nullptr;
  extraRangeName = std::string(resetRangeName);
  xBins = _Hnominal.GetXaxis()->GetNbins();
  xMin = _Hnominal.GetXaxis()->GetBinLowEdge(1);
  xMax = _Hnominal.GetXaxis()->GetBinUpEdge(xBins);
  memset(&xArray, 0, sizeof(xArray));
  for (Int_t i=0; i<xBins+1; i++){
    xArray[i] =  _Hnominal.GetXaxis()->GetBinLowEdge(i+1);
  }
  // const TAxis* axis = _Hnominal.GetXaxis();
  // bool foundMin = false;
  // int binNum = 0;
  // for (int i = 1; i <= axis->GetNbins(); ++i) {
  //   double low = axis->GetBinLowEdge(i);
  //   double up = axis->GetBinUpEdge(i);
  //   if (!foundMin) {
  //     if (low >= x.min()) {
  //       xMin = low;
  //       xArray[binNum++] = low;
  //       foundMin = true;
  //     } else {
  //       continue;
  //     }
  //   }
  //   if (up <= x.max()) {
  //     xArray[binNum++] = up;
  //     xMax = up;
  //   } else {
  //     break;
  //   }
  // }
  // xBins = binNum - 1;

  RooRealVar x_rrv = dynamic_cast<const RooRealVar &>(x.arg());
#if ROOT_VERSION_CODE < ROOT_VERSION(6, 24, 0)
    //modify x to use hash table for range lookup
    RooRealVarSmart* x_smart(static_cast<RooRealVarSmart*>(&x_rrv));
    x_smart->setHashTableSize(1);
#endif
  for (Int_t iBin = 0; iBin < xBins; iBin++) {
    std::string rangeName = Form("%s_%s_range_bin%d", GetName(), x.GetName(), iBin) + (extraRangeName.empty() ? "" : ("_" + extraRangeName));
    
    Double_t xLow = xArray[iBin];
    Double_t xHigh = xArray[iBin + 1];
    // never check here, because we want to reset the range every time we set the binning
    if (x.arg().hasRange(rangeName.c_str())) {
      double oldMin = x_rrv.getMin(rangeName.c_str());
      double oldMax = x_rrv.getMax(rangeName.c_str());
      std::cout << ">>> Range " << rangeName << " already exists, resetting it from [" << oldMin << ", " << oldMax << "] to [" << xLow << ", " << xHigh << "]" << std::endl;
    }
    x_rrv.setRange(rangeName.c_str(), xLow, xHigh);
  }
  std::string rangeName = Form("%s_%s_range_FULL", GetName(), x.GetName()) + (extraRangeName.empty() ? "" : ("_" + extraRangeName));
  if (x.arg().hasRange(rangeName.c_str())) {
    double oldMin = x_rrv.getMin(rangeName.c_str());
    double oldMax = x_rrv.getMax(rangeName.c_str());
    std::cout << ">>> Range " << rangeName << " already exists, resetting it from [" << oldMin << ", " << oldMax << "] to [" << xMin << ", " << xMax << "]" << std::endl;
  }
  x_rrv.setRange(rangeName.c_str(), xMin, xMax);
  // std::string overflowRangeName = Form("%s_%s_range_overflow", GetName(), x.GetName()) + (extraRangeName.empty() ? "" : ("_" + extraRangeName));
  // std::string underflowRangeName = Form("%s_%s_range_underflow", GetName(), x.GetName()) + (extraRangeName.empty() ? "" : ("_" + extraRangeName));
  // if (x.min() < xMin) {
  //   if (x.arg().hasRange(underflowRangeName.c_str())) {
  //     double oldMin = x_rrv.getMin(underflowRangeName.c_str());
  //     double oldMax = x_rrv.getMax(underflowRangeName.c_str());
  //     std::cout << ">>> Range " << underflowRangeName << " already exists, resetting it from [" << oldMin << ", " << oldMax << "] to [" << x.min() << ", " << xMin << "]" << std::endl;
  //   }
  //   x_rrv.setRange(underflowRangeName.c_str(), x.min(), xMin);
  // }
  // if (x.max() > xMax) {
  //   if (x.arg().hasRange(overflowRangeName.c_str())) {
  //     double oldMin = x_rrv.getMin(overflowRangeName.c_str());
  //     double oldMax = x_rrv.getMax(overflowRangeName.c_str());
  //     std::cout << ">>> Range " << overflowRangeName << " already exists, resetting it from [" << oldMin << ", " << oldMax << "] to [" << xMax << ", " << x.max() << "]" << std::endl;
  //   }
  //   x_rrv.setRange(overflowRangeName.c_str(), xMax, x.max());
  // }
}
//---------------------------------------------------------------------------
/// Return the parameteric p.d.f
RooAbsPdf* RooParametricShapeBinPdf::getPdf() const {
  auto* arg = mypdf.absArg();
  return arg ? static_cast<RooAbsPdf*>(arg) : nullptr;
}
//---------------------------------------------------------------------------
/// Return the bin-by-bin integrals
RooAbsReal* RooParametricShapeBinPdf::getIntegral(int index) const {
  if (myintegrals.empty()) {
    RooAbsReal* myintegral;
    RooRealVar x_rrv = dynamic_cast<const RooRealVar&>(x.arg());
#if ROOT_VERSION_CODE < ROOT_VERSION(6, 24, 0)
    //modify x to use hash table for range lookup
    RooRealVarSmart* x_smart(static_cast<RooRealVarSmart*>(&x_rrv));
    x_smart->setHashTableSize(1);
#endif

    RooArgSet obs;
    obs.add(x.arg());
    for (Int_t iBin = 0; iBin < xBins; iBin++) {
      std::string rangeName = Form("%s_%s_range_bin%d", GetName(), x.GetName(), iBin) + (extraRangeName.empty() ? "" : ("_" + extraRangeName));
      if (!x.arg().hasRange(rangeName.c_str())) {
        Double_t xLow = xArray[iBin];
        Double_t xHigh = xArray[iBin + 1];
        x_rrv.setRange(rangeName.c_str(), xLow, xHigh);
      }
      myintegral = getPdf()->createIntegral(obs, RooFit::NormSet(obs), Range(rangeName.c_str()));
      myintegrals.push_back(std::unique_ptr<RooAbsReal>(myintegral));
    }
    // cout << getPdf()->createIntegral(obs)->getVal() << std::endl;
    //last one is full integral (full range)
    std::string rangeName = Form("%s_%s_range_FULL", GetName(), x.GetName()) + (extraRangeName.empty() ? "" : ("_" + extraRangeName));
    if (!x.arg().hasRange(rangeName.c_str()) || x_rrv.getMin(rangeName.c_str()) != xMin || x_rrv.getMax(rangeName.c_str()) != xMax) {
      if (x.arg().hasRange(rangeName.c_str())) {
        std::cout << ">>> Range " << rangeName << " already exists, BUT MISMATCHED, resetting it from [" << x_rrv.getMin(rangeName.c_str()) << ", " << x_rrv.getMax(rangeName.c_str()) << "] to [" << xMin << ", " << xMax << "]" << std::endl;
      }
      x_rrv.setRange(rangeName.c_str(), xMin, xMax);
    }
    myintegral = getPdf()->createIntegral(obs, RooFit::NormSet(obs), Range(rangeName.c_str()));
    myintegrals.push_back(std::unique_ptr<RooAbsReal>(myintegral));
  }

  return myintegrals.at(index) ? myintegrals.at(index).get() : 0;
}

// RooAbsReal* RooParametricShapeBinPdf::getUnderflowIntegral() const {
//   double xmin_real = x.min();
//   if (xmin_real < xMin) {
//     RooRealVar x_rrv = dynamic_cast<const RooRealVar&>(x.arg());
//     std::string underflowRangeName = Form("%s_%s_range_underflow", GetName(), x.GetName()) + (extraRangeName.empty() ? "" : ("_" + extraRangeName));
//     if (xmin_real != xMinReal || !x.arg().hasRange(underflowRangeName.c_str())) {
//       if (xmin_real != xMinReal) {
//         std::cout << "WARNING: RooParametricShapeBinPdf::getUnderflowIntegral(): x minimum has been modified" << std::endl;
//         xMinReal = xmin_real;
//       }
//       x_rrv.setRange(underflowRangeName.c_str(), xMinReal, xMin);
//     }
//     if (!underflowIntegral) {
//       RooArgSet obs;
//       obs.add(x.arg());
//       underflowIntegral = std::unique_ptr<RooAbsReal>(getPdf()->createIntegral(obs, RooFit::NormSet(obs), Range(underflowRangeName.c_str())));
//     }
//     return underflowIntegral.get();
//   }
//   return 0;
// }

// RooAbsReal* RooParametricShapeBinPdf::getOverflowIntegral() const {
//   double xmax_real = x.max();
//   if (xmax_real > xMax) {
//     std::string overflowRangeName = Form("%s_%s_range_overflow", GetName(), x.GetName()) + (extraRangeName.empty() ? "" : ("_" + extraRangeName));
//     RooRealVar x_rrv = dynamic_cast<const RooRealVar&>(x.arg());
//     if (xmax_real != xMaxReal || !x.arg().hasRange(overflowRangeName.c_str())) {
//       if (xmax_real != xMaxReal) {
//         std::cout << "WARNING: RooParametricShapeBinPdf::getOverflowIntegral(): x maximum has been modified" << std::endl;
//         xMaxReal = xmax_real;
//       }
//       x_rrv.setRange(overflowRangeName.c_str(), xMax, xMaxReal);
//     }
//     if (!overflowIntegral) {
//       RooArgSet obs;
//       obs.add(x.arg());
//       overflowIntegral = std::unique_ptr<RooAbsReal>(getPdf()->createIntegral(obs, RooFit::NormSet(obs), Range(overflowRangeName.c_str())));
//     }
//     return overflowIntegral.get();
//   }
//   return 0;
// }
//---------------------------------------------------------------------------
Double_t RooParametricShapeBinPdf::evaluate() const
{
  Double_t integral = 0.0;
  Int_t iBin;

  // if (x < xMin && x.min() < xMin) {
  //   return getUnderflowIntegral()->getVal() / (xMin - x.min());
  // }
  // if (x >= xMax && x.max() > xMax) {
  //   return getOverflowIntegral()->getVal() / (x.max() - xMax);
  // }

  for(iBin=0; iBin<xBins; iBin++) {  
    if (x>=xArray[iBin] && x < xArray[iBin+1] ) break;
  }
  
  if(iBin < 0 || iBin >= xBins) {
    //cout << "in bin " << iBin << " which is outside of range" << endl;
    return 0.0;
  }

  Double_t xLow = xArray[iBin];
  Double_t xHigh = xArray[iBin+1];

  // check again if x variable has the right range already defined 
  // needed when combining multiple workspaces, and taking variable x from only one of them!
  std::string rangeName  = Form("%s_%s_range_bin%d", GetName(), x.GetName(), iBin) + (extraRangeName.empty() ? "" : ("_" + extraRangeName));
  if (!x.arg().hasRange(rangeName.c_str())) {
    RooRealVar x_rrv = dynamic_cast<const RooRealVar &>(x.arg());

#if ROOT_VERSION_CODE < ROOT_VERSION(6,24,0)
    //modify x to use hash table for range lookup
    RooRealVarSmart* x_smart(static_cast<RooRealVarSmart*>(&x_rrv));
    if(x_smart->getHashTableSize()==0) x_smart->setHashTableSize(1);
#endif

    Double_t xLow = xArray[iBin];
    Double_t xHigh = xArray[iBin+1];
    x_rrv.setRange(rangeName.c_str(),xLow,xHigh);
  }

  integral = getIntegral(iBin)->getVal() / (xHigh-xLow);
  
  if (integral>0.0) {
    return integral;
  } else return 0;

}

// //---------------------------------------------------------------------------
Int_t RooParametricShapeBinPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const{
  if (matchArgs(allVars, analVars, x)) return 1;
  return 0;
}

// //---------------------------------------------------------------------------
Double_t RooParametricShapeBinPdf::analyticalIntegral(Int_t code, const char* rangeName) const{

  Double_t xRangeMin = x.min(rangeName); Double_t xRangeMax = x.max(rangeName);
  
  Double_t integral = 0.0;
  
  RooArgSet obs;
  obs.add(x.arg());

  if (code==1 && xRangeMin == xMin && xRangeMax == xMax){
    integral = getIntegral(xBins)->getVal();
    // std::cout << ">>> RooParametricShapeBinPdf::analyticalIntegral(): using precomputed integral over full range and gives " << integral << std::endl;
    return integral;
  }
  else if(code==1) {
    RooAbsReal* myintegral = getPdf()->createIntegral(obs, RooFit::NormSet(obs), Range(rangeName));
    integral = myintegral->getVal();
    // std::cout << ">>> RooParametricShapeBinPdf::analyticalIntegral(): computing integral with analytical internal function for range " << rangeName << " and gives " << integral << std::endl;
    return integral;
   } else {
    cout << "WARNING IN RooParametricShapeBinPdf: integration code is not correct" << endl;
    cout << "                           what are you integrating on?" << endl;
    return 1.0;
  }
}

const std::vector<std::unique_ptr<RooAbsReal>>& RooParametricShapeBinPdf::checkIntegrals() const
{
  return myintegrals;
}
// //---------------------------------------------------------------------------

