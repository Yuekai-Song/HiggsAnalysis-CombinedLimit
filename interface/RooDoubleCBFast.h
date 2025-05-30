#ifndef ROODOUBLECB
#define ROODOUBLECB

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"

class RooDoubleCBFast : public RooAbsPdf {
public:
  RooDoubleCBFast();
  RooDoubleCBFast(const char *name, const char *title,
              RooAbsReal& _x,
              RooAbsReal& _mean,
              RooAbsReal& _width,
              RooAbsReal& _alpha1,
              RooAbsReal& _n1,
              RooAbsReal& _alpha2,
              RooAbsReal& _n2
           );
  RooDoubleCBFast(const RooDoubleCBFast& other, const char* name=0) ;
  TObject* clone(const char* newname) const override { return new RooDoubleCBFast(*this,newname); }
  inline ~RooDoubleCBFast() override { }
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const override ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const override ;

protected:

  RooRealProxy x ;
  RooRealProxy mean;
  RooRealProxy width;
  RooRealProxy alpha1;
  RooRealProxy n1;
  RooRealProxy alpha2;
  RooRealProxy n2;
  
  Double_t evaluate() const override ;

private:

  ClassDefOverride(RooDoubleCBFast,1)
};
#endif
