#ifndef ROOUNFOLDPARMS_H_
#define ROOUNFOLDPARMS_H_

#include "TNamed.h"

class TH1;
class RooUnfold;
class TProfile;

class RooUnfoldParms : public TNamed {
	public:
	RooUnfoldParms(const RooUnfold* unfold_in=0,Int_t err=1,const TH1* truth=0);
	virtual ~RooUnfoldParms();
	TProfile* GetChi2();
	TProfile* GetRMSError();
	TProfile* GetMeanResiduals();
	TH1* GetRMSResiduals();
	const RooUnfold* unfold; // Input object from RooUnfold
	Int_t doerror; // Set error calculation method
	const TH1* hTrue; // Truth Distribution
	void SetMinParm(double min);
	void SetMaxParm(double max);
	void SetStepSizeParm(double size);
	
	private:
	bool _done_math;
	TH1* hrms; // Output plot
	TProfile* hch2; // Output plot
	TProfile* herr; // Output plot
	TProfile* hres; // Output plot
	void DoMath();
	void Init();
	Double_t _maxparm; //Maximum parameter
	Double_t _minparm; //Minimum parameter
	Double_t _stepsizeparm; //Step size
public:
	ClassDef (RooUnfoldParms, 0)
};
#endif /*ROOUNFOLDPARMS_H_*/
