/*
 * TauSVDUnfold.h
 *
 *  Created on: 31 Dec 2014
 *      Author: phxlk
 */

#ifndef TAUSVDUNFOLD_H_
#define TAUSVDUNFOLD_H_
#include "TSVDUnfold_local.h"

class TauSVDUnfold: public TSVDUnfold_local {
public:
	TauSVDUnfold(const TH1D* bdat, const TH1D* bini, const TH1D* xini, const TH2D* Adet);
	TauSVDUnfold(const TH1D* bdat, TH2D* Bcov, const TH1D* bini, const TH1D* xini, const TH2D* Adet);
	TauSVDUnfold(const TauSVDUnfold& other);
	virtual ~TauSVDUnfold();

	void SetTau(double tau);
	double GetTau() const;

	// Get Curvature of last unfolding run
	double GetCurv() const;

	// Get Histo of Weights of LAST Unfolding run
	TH1D* GetWeights() const;

	// Do the unfolding
	// "tau"   - number of singular values used (regularisation)
	TH1D* Unfold(double tau);

	double kToTau(int kreg) const;

protected:
	double fTau;
	double fCurv;
	TH1D* fWeights;
	TVectorD fASV;

	ClassDef( TauSVDUnfold, 0 );
};

#endif /* TAUSVDUNFOLD_H_ */
