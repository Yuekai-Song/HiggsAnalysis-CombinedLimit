#include "../interface/RooEFTScalingFunction.h"
#include <RooAbsArg.h>

ClassImp(RooEFTScalingFunction)

    RooEFTScalingFunction::RooEFTScalingFunction(const char *name,
                                                 const char *title,
                                                 const std::map<std::string, double> &coeffs,
                                                 const RooArgList &terms)
    : RooAbsReal(name, title),
      coeffs_(coeffs),
      terms_("!terms", "RooArgList of Wilson coefficients", this),
      offset_(1.0) {
  // Add Wilson coefficients to terms_ container
  for (RooAbsArg *a : terms) {
    RooAbsReal *rar = dynamic_cast<RooAbsReal *>(a);
    if (!rar) {
      throw std::invalid_argument(std::string("Term ") + a->GetName() + " of RooEFTScalingFunction is a " +
                                  a->ClassName());
    }
    terms_.add(*rar);
  }

  // Loop over elements in mapping: add components to vector depending on string
  for (auto const &x : coeffs) {
    TString term_name = x.first;
    double term_prefactor = x.second;

    if (term_name.Contains("_")) {
      TString first_term = dynamic_cast<TObjString *>(term_name.Tokenize("_")->At(0))->GetString();
      TString second_term = dynamic_cast<TObjString *>(term_name.Tokenize("_")->At(1))->GetString();

      // Squared-quadratic components
      int rar1 = terms_.index(first_term);
      int rar2 = rar1;
      if (second_term == "2") {
        if (rar1 != -1) {
          std::vector<int> vterms;
          // RooAbsReal *rar1 = dynamic_cast<RooAbsReal *>( terms.find(first_term) );
          // RooAbsReal *rar2 = dynamic_cast<RooAbsReal *>( terms.find(first_term) );
          vterms.push_back(rar1);
          vterms.push_back(rar2);
          vcomponents_.emplace(vterms, term_prefactor);
        }
      } else {
        rar2 = terms_.index(second_term);
        // Cross-quadratic components
        if (rar1 != -1 && rar2 != -1) {
          std::vector<int> vterms;
          // RooAbsReal *rar1 = dynamic_cast<RooAbsReal *>( terms.find(first_term) );
          // RooAbsReal *rar2 = dynamic_cast<RooAbsReal *>( terms.find(first_term) );
          vterms.push_back(rar1);
          vterms.push_back(rar2);
          vcomponents_.emplace(vterms, term_prefactor);
        }
      }
    } else {
      int rar = terms_.index(term_name);
      if (rar != -1) {
        std::vector<int> vterms;
        //   RooAbsReal *rar = dynamic_cast<RooAbsReal *>( terms.find(term_name) );
        vterms.push_back(rar);
        vcomponents_.emplace(vterms, term_prefactor);
      }
    }
  }
}

RooEFTScalingFunction::RooEFTScalingFunction(const RooEFTScalingFunction &other, const char *name)
    : RooAbsReal(other, name),
      coeffs_(other.coeffs_),
      terms_("!terms", this, other.terms_),
      vcomponents_(other.vcomponents_),
      offset_(other.offset_) {}

Double_t RooEFTScalingFunction::evaluate() const {
  if (vcomponents_.empty()) {
    return offset_;
  }
  double ret = offset_;
  for (auto const &x : vcomponents_) {
    double res = x.second;
    for (auto const &y : x.first) {
      res *= dynamic_cast<RooAbsReal *>(terms_.at(y))->getVal();
    }
    ret += res;
  }
  return ret;
}
