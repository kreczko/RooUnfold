//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldTUnfold.h 274 2011-02-04 22:07:32Z T.J.Adye $
//
// Description:
//      Unfolding class using TUnfold from ROOT to do the actual unfolding.
//
// Authors: Richard Claridge <richard.claridge@stfc.ac.uk> & Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

#ifndef ROOUNFOLDTUNFOLD_H_
#define ROOUNFOLDTUNFOLD_H_

#include "RooUnfold.h"
#include "TUnfold.h"

class RooUnfoldResponse;
class TH1;
class TH1D;
class TH2D;

class RooUnfoldTUnfold : public RooUnfold {

public:

  RooUnfoldTUnfold(); // default constructor
  RooUnfoldTUnfold (const char*    name, const char*    title); // named constructor
  RooUnfoldTUnfold (const TString& name, const TString& title); // named constructor
  RooUnfoldTUnfold (const RooUnfoldTUnfold& rhs); // copy constructor
  virtual ~RooUnfoldTUnfold(); // destructor
  RooUnfoldTUnfold& operator= (const RooUnfoldTUnfold& rhs); // assignment operator
  virtual RooUnfoldTUnfold* Clone (const char* newname= 0) const;
  RooUnfoldTUnfold (const RooUnfoldResponse* res, const TH1* meas,TUnfold::ERegMode reg=TUnfold::kRegModeDerivative,
  const char* name= 0, const char* title= 0);

  void Reset();
  TUnfold* Impl();
  void FixTau(double tau);
  void OptimiseTau();
  virtual void SetRegParm(double parm);
  double GetTau() const;
  virtual double GetRegParm() const;
  void SetRegMethod (TUnfold::ERegMode regmethod);
  TUnfold::ERegMode GetRegMethod() const;

protected:
  void Init();
  void Destroy();
  virtual void Unfold();
  virtual void GetCov();
  virtual void GetSettings();
  void Assign   (const RooUnfoldTUnfold& rhs); // implementation of assignment operator
  void CopyData (const RooUnfoldTUnfold& rhs);

private:
  TUnfold::ERegMode _reg_method; //Regularisation method
  TUnfold* _unf; //! Implementation in TUnfold object (no streamer)
  bool tau_set;
  double _tau;

public:

  ClassDef (RooUnfoldTUnfold, 1)   // Interface to TUnfold
};

// Inline method definitions

inline
RooUnfoldTUnfold::RooUnfoldTUnfold()
: RooUnfold()
{
  // Default constructor. Use Setup() to prepare for unfolding.
  Init();
}

inline
RooUnfoldTUnfold::RooUnfoldTUnfold (const char* name, const char* title)
: RooUnfold(name,title)
{
  // Basic named constructor. Use Setup() to prepare for unfolding.
  Init();
}

inline
RooUnfoldTUnfold::RooUnfoldTUnfold (const TString& name, const TString& title)
: RooUnfold(name,title)
{
  // Basic named constructor. Use Setup() to prepare for unfolding.
  Init();
}

inline
RooUnfoldTUnfold& RooUnfoldTUnfold::operator= (const RooUnfoldTUnfold& rhs)
{
  // Assignment operator for copying RooUnfoldTUnfold settings.
  Assign(rhs);
  return *this;
}

inline
RooUnfoldTUnfold::~RooUnfoldTUnfold()
{
  Destroy();
}


inline
void RooUnfoldTUnfold::SetRegParm(double parm)
{
  // Set regularisation parameter (tau)
  FixTau(parm);
}

inline
double RooUnfoldTUnfold::GetTau() const
{
  // Return regularisation parameter (tau)
  return _tau;
}

inline
double RooUnfoldTUnfold::GetRegParm() const
{
  // Return regularisation parameter (tau)
  return _tau;
}

inline
TUnfold::ERegMode RooUnfoldTUnfold::GetRegMethod() const
{
  return _reg_method;
}

#endif /*ROOUNFOLDTUNFOLD_H_*/
