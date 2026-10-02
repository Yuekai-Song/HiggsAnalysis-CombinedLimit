/*****************************************************************************
 * RooBernsteinFastOld
 * Josh Bendavid (CERN)
 * Fast templated version of RooBernstein class using SMatrices
 * 
 * 
 *   
 *       
 *****************************************************************************/
#ifndef ROO_BERNSTEINFASTOLD
#define ROO_BERNSTEINFASTOLD

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooChangeTracker.h"
#include "TMath.h"
#include "Math/SMatrix.h"

class RooRealVar;
class RooArgList ;

// Use N parameters for degree-N polynomial, with the first coefficient c0=1 implicit. The coefficients in the RooArgList are c1..cN.
template<int N> class RooBernsteinFastOld : public RooAbsPdf {
public:
  
  RooBernsteinFastOld() {}
  
  RooBernsteinFastOld(const char *name, const char *title,
                RooAbsReal& x, const RooArgList& coefList) :
    RooAbsPdf(name, title),
    _x("x", "Dependent", this, x),
    _coefList("coefList","List of coefficients",this)
    {
      _xmin = _x.min();
      _xmax = _x.max();
      _coefList.add(coefList);
      
      
      //precompute coefficients for integral
      for (int ipow=0; ipow<=N; ++ipow) {
        _rvector(ipow) = 1.0/((double)ipow+1.0);
      }

      //precompute coefficients for conversion from bernstein basis to power basis
      for (int ibern=0; ibern<=N; ++ibern) {
        for (int ipow=0; ipow<ibern; ++ipow) {
          _cmatrix(ipow,ibern) = 0.;
        }
        for (int ipow=ibern; ipow<=N; ++ipow) {
          _cmatrix(ipow, ibern) = pow(-1.,ipow-ibern)*TMath::Binomial(N,ipow)*TMath::Binomial(ipow,ibern);
        }
      }
      
    }                
                

  RooBernsteinFastOld(const RooBernsteinFastOld& other, const char* name = 0) :
      RooAbsPdf(other, name), 
    _x("x", this, other._x), 
    _xmin(other._xmin),
    _xmax(other._xmax),
    _coefList("coefList",this,other._coefList),
    _cmatrix(other._cmatrix),
    _rvector(other._rvector),
    _bernvector(other._bernvector),
    _powvector(other._powvector),
    _xvector(other._xvector) {}
  
  
  TObject* clone(const char* newname) const override { return new RooBernsteinFastOld(*this, newname); }
  ~RooBernsteinFastOld() override { }
  
  
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const override
    {

      // No analytical calculation available (yet) of integrals over subranges (as for standard RooBernstein)
      // if (rangeName && strlen(rangeName)) {
      //   return 0 ;
      // }

      if (matchArgs(allVars, analVars, _x)) return 1;
      return 0;
    }
    
    
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const override
    {

      _bernvector[0] = 1.0;
      for (int ipow=1; ipow<=N; ++ipow) {
        _bernvector[ipow] = static_cast<RooAbsReal*>(_coefList.at(ipow-1))->getVal();
      }     
      
      _powvector = _cmatrix*_bernvector;
        
      // double xmin = _x.min();
      // double xmax = _x.max();    
      if (_xmin > _x.min() || _xmax < _x.max()) {
        std::cout << "WARNING: RooBernsteinFastOld analyticalIntegral(): x range is modified to be larger than the original range." << std::endl;
      }
      double xmin = _xmin > _x.min() ? _x.min() : _xmin;
      double xmax = _xmax < _x.max() ? _x.max() : _xmax;
      const double xlo = (_x.min(rangeName) - xmin) / (xmax - xmin);
      const double xhi = (_x.max(rangeName) - xmin) / (xmax - xmin);

      double norm = 0;
      for (int ipow=0; ipow<=N; ++ipow) {
        norm += (TMath::Power(xhi,ipow+1) - TMath::Power(xlo,ipow+1)) * _powvector[ipow] * _rvector[ipow];
      }
      return (xmax-xmin)*norm;

    }

protected:

  typedef ROOT::Math::SMatrix<double,N+1,N+1,ROOT::Math::MatRepStd<double,N+1,N+1> > MType;
  typedef ROOT::Math::SVector<double,N+1> VType;
  
  RooRealProxy _x;
  double _xmin, _xmax;
  RooListProxy _coefList ;
  MType _cmatrix;            //conversion matrix between bernstein and power bases
  VType _rvector;            //vector of integration coefficients
  mutable VType _bernvector; //coefficients in bernstein basis
  mutable VType _powvector;  //coefficients in power basis
  mutable VType _xvector;    //vector of powers of x variable
  
  Double_t evaluate() const override
    {

      _bernvector[0] = 1.0;
      bool changed = false;
      for (int ipow=1; ipow<=N; ++ipow) {
        double rval = static_cast<RooAbsReal*>(_coefList.at(ipow-1))->getVal();
        if (_bernvector[ipow] != rval) {
          _bernvector[ipow] = rval;
          changed = true;
        }
      }
      
      if (changed) {
        _powvector = _cmatrix*_bernvector;   
      }
      
      // double xmin = _x.min();
      // double xmax = _x.max();
      if (_xmin > _x.min() || _xmax < _x.max()) {
        std::cout << "WARNING: RooBernsteinFastOld evaluate(): x range is modified to be larger than the original range." << std::endl;
      }
      double xmin = _xmin > _x.min() ? _x.min() : _xmin;
      double xmax = _xmax < _x.max() ? _x.max() : _xmax;
      double x = (_x - xmin) / (xmax - xmin); // rescale to [0,1]
      _xvector[0] = 1.;
      for (int ipow=1; ipow<=N; ++ipow) {
        _xvector[ipow] = x*_xvector[ipow-1];
      }
      
      return ROOT::Math::Dot(_powvector,_xvector);
      
    }

  ClassDefOverride(RooBernsteinFastOld,1) // Polynomial PDF
};

#endif