/*
 * TauSVDUnfold.h
 *
 * original:
 * https://github.com/eschliec/TopAnalysis/tree/master/Configuration/analysis/unfolding
 *
 *  Modified on: 31 Dec 2014
 *      Author: phxlk
 */

#ifndef TAUSVDUNFOLD_H_
#define TAUSVDUNFOLD_H_
#include "TSVDUnfold_local.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"

class TauSVDUnfold: public TSVDUnfold_local {
public:
	/**
	 * For documentation on constructors see TSVDUnfold::TSVDUnfold
	 */
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

	TGraph* ScanLCurve(unsigned int n_point, double tau_min, double tau_max);
	double GetLcurveX() const;
	double GetLcurveY() const;
	virtual TH2D* GetUnfoldCovMatrix( const TH2D* cov, Int_t ntoys, Int_t seed );
	virtual TH2D* GetAdetCovMatrix( Int_t ntoys, Int_t seed );
	/**
	 * Added by Joern: adapted from TSVDUnfold::GetUnfoldCovMatrix to include
	 * the normalisation in the procedure
	 * Determine for given input error matrix covariance matrix of unfolded
	 * spectrum from toy simulation given the passed covariance matrix on measured spectrum
	 * "cov"    - covariance matrix on the measured spectrum, to be propagated
	 * "ntoys"  - number of pseudo experiments used for the propagation
	 * "seed"   - seed for pseudo experiments
	 * "normType" - set type of normalisation (1=extrinsic, 2=intrinsic)
	 * Note that this covariance matrix will contain effects of forced normalisation if spectrum is normalised to unit area.
	 */
//	TH2D* GetUnfoldCovMatrixNorm(const TH2D* cov, Int_t ntoys, Int_t seed = 1, Int_t normType = 2, Int_t verbose = 0);
//	static TH1D* IntNormalizeSVDDistribution(TH1D* inputHist);
	/**
	 * Calculates the covariance matrix (dispersion matrix)
	 * Takes a 1D histogram and returns a diagonal 2D histogram with entries
	 * cov(i,i) = error(i)^2
	 */
	static TH2D* get_data_covariance_matrix(const TH1D* data_histogram);

	static double get_global_correlation(const TH2D* stat_cov_hist, const TH1D* data_histogram);

	static TH1D* get_global_correlation_hist(const TH2D* covariance_hist, const TH1D* data_histogram);

	TVectorD getASV() const;
protected:
	double fTau;
	double fCurv;
	TH1D* fWeights;
	TVectorD fASV;

ClassDef( TauSVDUnfold, 0 )
	;
};

#endif /* TAUSVDUNFOLD_H_ */
