/*
 * TauSVDUnfold.cxx
 *
 *  Created on: 31 Dec 2014
 *      Author: kreczko
 */

#include "../include/TauSVDUnfold.h"
#include "TDecompSVD.h"
#include "TMath.h"
#include "TH1D.h"
#include "TRandom3.h"
#include <vector>
#include <iostream>

using namespace std;

TauSVDUnfold::TauSVDUnfold(const TH1D* bdat, const TH1D* bini, const TH1D* xini, const TH2D* Adet) :
				TSVDUnfold_local(bdat, bini, xini, Adet),
				fTau(0.),
				fCurv(-1.),
				fWeights(),
				fASV() {
}

TauSVDUnfold::TauSVDUnfold(const TH1D* bdat, TH2D* Bcov, const TH1D* bini, const TH1D* xini, const TH2D* Adet) :
				TSVDUnfold_local(bdat, Bcov, bini, xini, Adet),
				fTau(0.),
				fCurv(-1.),
				fWeights(),
				fASV() {

}
TauSVDUnfold::TauSVDUnfold(const TauSVDUnfold& other) :
				TSVDUnfold_local(other),
				fTau(other.GetTau()),
				fCurv(other.GetCurv()),
				fWeights(other.GetWeights()),
				fASV() {

}

TH1D* TauSVDUnfold::GetWeights() const {
	return fWeights;
}

double TauSVDUnfold::GetTau() const {
	return fTau;
}

void TauSVDUnfold::SetTau(double tau) {
	fTau = tau;
}

double TauSVDUnfold::GetCurv() const {
	return fCurv;
}

TauSVDUnfold::~TauSVDUnfold() {
}

TH1D* TauSVDUnfold::Unfold(double tau) {
	fTau = tau;
	// Make the histos
	if (!fToyMode && !fMatToyMode)
		InitHistos();

	// Create vectors and matrices
	TVectorD vb(fNdim), vbini(fNdim), vxini(fNdim), vberr(fNdim);
	TMatrixD mB(fNdim, fNdim), mA(fNdim, fNdim), mCurv(fNdim, fNdim), mC(fNdim, fNdim);

	Double_t eps = 1e-12;
	Double_t sreg;
	// Copy histogams entries into vector
	if (fToyMode) {
		H2V(fToyhisto, vb);
		H2Verr(fToyhisto, vberr);
	} else {
		H2V(fBdat, vb);
		H2Verr(fBdat, vberr);
	}
	H2M(fBcov, mB);
	H2V(fBini, vbini);
	H2V(fXini, vxini);
	if (fMatToyMode)
		H2M(fToymat, mA);
	else
		H2M(fAdet, mA);

	// Fill and invert the second derivative matrix
	FillCurvatureMatrix(mCurv, mC);
	// Inversion of mC by help of SVD
	TDecompSVD CSVD(mC);
	TMatrixD CUort = CSVD.GetU();
	TMatrixD CVort = CSVD.GetV();
	TVectorD CSV = CSVD.GetSig();

	TMatrixD CSVM(fNdim, fNdim);
	for (Int_t i = 0; i < fNdim; i++)
		CSVM(i, i) = 1 / CSV(i);

	CUort.Transpose(CUort);
	TMatrixD mCinv = (CVort * CSVM) * CUort;
	//    // Rescale matrix and vectors by error of data vector. Replaced by using full covmat now
	//    vbini = VecDiv   ( vbini, vberr );
	//    vb    = VecDiv   ( vb,    vberr, 1 );
	//    mA    = MatDivVec( mA,    vberr, 1 );
	//    vberr = VecDiv   ( vberr, vberr, 1 );

	//Rescale using the data covariance matrix
	TDecompSVD BSVD(mB);
	TMatrixD QT = BSVD.GetU();
	QT.Transpose(QT);
	TVectorD B2SV = BSVD.GetSig();
	TVectorD BSV(B2SV);
	for (int i = 0; i < fNdim; i++) {
		BSV(i) = TMath::Sqrt(B2SV(i));
	}
	TMatrixD mAtmp(fNdim, fNdim);
	TVectorD vbtmp(fNdim);
	mAtmp *= 0;
	vbtmp *= 0;
	for (int i = 0; i < fNdim; i++) {
		for (int j = 0; j < fNdim; j++) {
			if (BSV(i)) {
				vbtmp(i) += QT(i, j) * vb(j) / BSV(i);
			}
			for (int m = 0; m < fNdim; m++) {
				if (BSV(i)) {
					mAtmp(i, j) += QT(i, m) * mA(m, j) / BSV(i);
				}
			}
		}
	}
	mA = mAtmp;
	vb = vbtmp;
	// Singular value decomposition and matrix operations
	TDecompSVD ASVD(mA * mCinv);
	TMatrixD Uort = ASVD.GetU();
	TMatrixD Vort = ASVD.GetV();
	// TODO: This part is failing for tau unfolding
	TVectorD ASV = ASVD.GetSig();
//	fASV = TVectorD(ASV);
	if (!fToyMode && !fMatToyMode) {
		V2H(ASV, *fSVHist);
	}

	TMatrixD Vreg = mCinv * Vort;
	Uort.Transpose(Uort);
	TVectorD vd = Uort * vb;

	if (!fToyMode && !fMatToyMode) {
		V2H(vd, *fDHist);
	}

	Int_t k = GetKReg() - 1;
	TVectorD vx(fNdim); // Return variable

	// Damping factors
	TVectorD vdz(fNdim);
	TMatrixD Z(fNdim, fNdim);
	for (Int_t i = 0; i < fNdim; i++) {
		if (ASV(i) < ASV(0) * eps)
			sreg = ASV(0) * eps;
		else
			sreg = ASV(i);

		vdz(i) = sreg / (sreg * sreg + fTau * fTau);
		Z(i, i) = vdz(i) * vdz(i);
	}
	TVectorD vz = TSVDUnfold_local::CompProd(vd, vdz);

	TMatrixD VortT(Vort);
	VortT.Transpose(VortT);
	TMatrixD W = mCinv * Vort * Z * VortT * mCinv;

	TMatrixD Xtau(fNdim, fNdim);
	TMatrixD Xinv(fNdim, fNdim);
	Xtau *= 0;
	Xinv *= 0;
	for (Int_t i = 0; i < fNdim; i++) {
		for (Int_t j = 0; j < fNdim; j++) {
			Xtau(i, j) = vxini(i) * vxini(j) * W(i, j);

			double a = 0;
			for (Int_t m = 0; m < fNdim; m++) {
				a += mA(m, i) * mA(m, j);
			}
			if (vxini(i) * vxini(j))
				Xinv(i, j) = a / vxini(i) / vxini(j);
		}
	}
	// Compute the weights
	TVectorD vw = Vreg * vz;

	// DAVID
	vx = TSVDUnfold_local::CompProd(vw, vxini);

	if (fNormalize) { // Scale result to unit area
		Double_t scale = vx.Sum();
		if (scale > 0) {
			vx *= 1.0 / scale;
			Xtau *= 1. / scale / scale;
			Xinv *= scale * scale;
		}
	}

	if (!fToyMode && !fMatToyMode) {
		M2H(Xtau, *fXtau);
		M2H(Xinv, *fXinv);
	}

	// DAVID
	// Speichere die Kruemmung ab!
	fCurv = GetCurvature(vw, mCurv);
	// Get Curvature and also chi2 in case of MC unfolding
	if (!fToyMode && !fMatToyMode) {
		Info("Unfold", "Unfolding param: %i", k + 1);
		Info("Unfold", "Curvature of weight distribution: %f", fCurv);
	}

	TH1D* h = (TH1D*) fBdat->Clone("unfoldingresult");
	for (int i = 1; i <= fNdim; i++) {
		h->SetBinContent(i, 0.);
		h->SetBinError(i, 0.);
	}
	V2H(vx, *h);
	// DAVID
	// Save the weights
	// but only if this is not a "Toy"-Run
	if (!fToyMode && !fMatToyMode) {
		if (fWeights != NULL) {
			delete fWeights;
			fWeights = NULL;
		}
		fWeights = (TH1D*) fBdat->Clone("Weights");
		V2H(vw, *fWeights);
	}

	return h;
}

double TauSVDUnfold::kToTau(int kreg) const {
	double tau(0.);
	if (fASV.NonZeros() == 0 && kreg > 0)
		tau = fASV(kreg);
	return tau;
}

TH2D* TauSVDUnfold::get_data_covariance_matrix(const TH1D* data_histogram) {
	// get bins from data_histogram
	unsigned int n_bins = data_histogram->GetNbinsX();
	vector<double> bin_edges;
	for (unsigned int i = 1; i <= n_bins; ++i) {
		double low_edge = data_histogram->GetBinLowEdge(i);
		double bin_width = data_histogram->GetBinWidth(i);
		bin_edges.push_back(low_edge);
		bin_edges.push_back(low_edge + bin_width);
	}

	// create 2D histogram
	TH2D * data_covariance_matrix = new TH2D("data_covariance_matrix", "data_covariance_matrix", n_bins, &bin_edges[0],
			n_bins, &bin_edges[0]);

	// fill it
	for (unsigned int i = 1; i <= n_bins; ++i) {
		double variance = pow(data_histogram->GetBinError(i), 2.);
		data_covariance_matrix->SetBinContent(i, i, variance);
	}

	return data_covariance_matrix;
}

TH1D* TauSVDUnfold::get_global_correlation_hist(const TH1D* data_histogram) {
	TH2D* covariance_hist = (TH2D*) TauSVDUnfold::get_data_covariance_matrix(data_histogram);

	unsigned int n_bins = data_histogram->GetNbinsX();
	vector<int> bin_map;
	bin_map.resize(n_bins);

	unsigned int number_of_bins_to_skip(0);
	unsigned int bin_counter(0);

	// all bins except underflow (0) and overflow (nbins + 1)
	for (auto i = 1; i <= n_bins; ++i) {
		double data = data_histogram->GetBinContent(i);
		if (data <= 0.) {
			// Through out bin with no data
			bin_map[i - 1] = -1;
			++number_of_bins_to_skip;
		} else {
			// Search for  bins with empty rows/columns
			bool skip_bin = true;
			for (auto j = 1; j <= n_bins; ++j) {
				double value = covariance_hist->GetBinContent(i, j);
				if (value != 0.) {
					skip_bin = false;
				}
			}
			// Through out bins with empty rows/columns
			if (skip_bin == true) {
				bin_map[i - 1] = -1;
				++number_of_bins_to_skip;
			} else {
				bin_map[i - 1] = bin_counter;
				bin_counter++;
			}
		}
	}

	unsigned int matrix_dimension = n_bins - number_of_bins_to_skip;
	TMatrixDSym covariance_matrix(matrix_dimension);

	// New Matrix
	// Beware the side bins of the problem
	// AND the side bins of the TH2D object
	for (unsigned int i = 1; i <= n_bins; i++) {
		for (unsigned int j = 1; j <= n_bins; j++) {

			// Is this bin to be skipped?
			bool skip_bin = false;
			if (bin_map[i - 1] == -1)
				skip_bin = true;
			if (bin_map[j - 1] == -1)
				skip_bin = true;
			// Set Element
			if (skip_bin == false) {
				double value = covariance_hist->GetBinContent(i, j);
				int bin_nr_i = bin_map[i - 1];
				int bin_nr_j = bin_map[j - 1];
				covariance_matrix[bin_nr_i][bin_nr_j] = value;
			}
		}
	}

	// Determinant
	double *det_covariance_matrix = new double(0.);

	// Invert the whole thing
	TMatrixDSym covariance_matrix_invers = covariance_matrix;
	covariance_matrix_invers.Invert(det_covariance_matrix);

	// Check Invertibility
	bool is_invertible = *det_covariance_matrix != 0.;
	if (is_invertible == false) {
		cout << "Error in TauSVDUnfold::get_global_correlation_hist() " << endl;
		cout << "Covariance Matrix cannot be inverted." << endl;
		cout << "Check the reason for this now." << endl;
		exit(1);
	}

	// Create new Histo for global correlation
	TH1D* global_correlation_hist = new TH1D();
	global_correlation_hist->SetNameTitle("global_correlation_hist", "global_correlation_hist");
	if (covariance_hist->GetXaxis()->IsVariableBinSize() == true) {
		const TArrayD* xbins = covariance_hist->GetXaxis()->GetXbins();
		global_correlation_hist->SetBins(n_bins, xbins->GetArray());
	} else {
		double xmin = covariance_hist->GetXaxis()->GetXmin();
		double xmax = covariance_hist->GetXaxis()->GetXmax();
		global_correlation_hist->SetBins(n_bins, xmin, xmax);
	}

	// Fill the histo you just created
	for (unsigned int i = 1; i <= n_bins; i++) {
		double global_correlation = 0.;

		// Find out the "true" bin number
		int binnr = bin_map[i - 1];

		// Skip bad bins
		bool skipThis = false;
		if (binnr == -1)
			skipThis = true;

		// Run over good bins
		if (skipThis == false) {
			double cov = covariance_matrix[binnr][binnr];
			double covinv = covariance_matrix_invers[binnr][binnr];
			// The product cov*covinv should be greater than zero
			double cov_prod = cov * covinv;
			double var(0.);
			if (cov_prod > 0)
				var = 1. / cov_prod;
			double global_correlation_squared = 0.;
//			cout << "cov: " << cov << endl;
//			cout << "covinv: " << covinv << endl;
//			cout << "cov_prod: " << cov_prod << endl;
//			cout << "var: " << var << endl;
			if (var > 0.)
				global_correlation_squared = 1. - var;

			if (global_correlation_squared > 0) {
				global_correlation = 100. * sqrt(global_correlation_squared);
			} else {
				global_correlation = 0.;
			}
		} else {
			global_correlation = 0.;
		}

		// Set the value
		global_correlation_hist->SetBinContent(i, global_correlation);
		global_correlation_hist->SetBinError(i, 0.);
	}
	delete covariance_hist;

	return global_correlation_hist;
}

double TauSVDUnfold::get_global_correlation(const TH1D* data_histogram) {
	TH1D * global_correlation_hist = TauSVDUnfold::get_global_correlation_hist(data_histogram);

	double sum = 0.;
	double average = 0.;
	int bincounter = 0;
	unsigned int n_bins = global_correlation_hist->GetNbinsX();

	for (unsigned int i = 1; i <= n_bins; i++) {
		if (i == 1)
			continue;
		if (i == n_bins)
			continue;
		double globcorr = global_correlation_hist->GetBinContent(i);
		sum += globcorr;
		bincounter++;
	}

	// Averaging
	if (bincounter > 0) {
		average = sum / ((double) bincounter);
	}

	delete global_correlation_hist;

	return average;
}

TH2D* TauSVDUnfold::GetUnfoldCovMatrix(const TH2D* cov, Int_t ntoys, Int_t seed) {
	// Determine for given input error matrix covariance matrix of unfolded
	// spectrum from toy simulation given the passed covariance matrix on measured spectrum
	// "cov"    - covariance matrix on the measured spectrum, to be propagated
	// "ntoys"  - number of pseudo experiments used for the propagation
	// "seed"   - seed for pseudo experiments
	// Note that this covariance matrix will contain effects of forced normalisation if spectrum is normalised to unit area.

	fToyMode = true;
	TH1D* unfres = 0;
	TH2D* unfcov = (TH2D*) fAdet->Clone("unfcovmat");
	unfcov->SetTitle("Toy covariance matrix");
	for (int i = 1; i <= fNdim; i++)
		for (int j = 1; j <= fNdim; j++)
			unfcov->SetBinContent(i, j, 0.);

	// Code for generation of toys (taken from RooResult and modified)
	// Calculate the elements of the upper-triangular matrix L that
	// gives Lt*L = C, where Lt is the transpose of L (the "square-root method")
	TMatrixD L(fNdim, fNdim);
	L *= 0;

	for (Int_t iPar = 0; iPar < fNdim; iPar++) {

		// Calculate the diagonal term first
		L(iPar, iPar) = cov->GetBinContent(iPar + 1, iPar + 1);
		for (Int_t k = 0; k < iPar; k++)
			L(iPar, iPar) -= TMath::Power(L(k, iPar), 2);
		if (L(iPar, iPar) > 0.0)
			L(iPar, iPar) = TMath::Sqrt(L(iPar, iPar));
		else
			L(iPar, iPar) = 0.0;

		// ...then the off-diagonal terms
		for (Int_t jPar = iPar + 1; jPar < fNdim; jPar++) {
			L(iPar, jPar) = cov->GetBinContent(iPar + 1, jPar + 1);
			for (Int_t k = 0; k < iPar; k++)
				L(iPar, jPar) -= L(k, iPar) * L(k, jPar);
			if (L(iPar, iPar) != 0.)
				L(iPar, jPar) /= L(iPar, iPar);
			else
				L(iPar, jPar) = 0;
		}
	}

	// Remember it
	TMatrixD *Lt = new TMatrixD(TMatrixD::kTransposed, L);
	TRandom3 random(seed);

	fToyhisto = (TH1D*) fBdat->Clone("toyhisto");
	TH1D *toymean = (TH1D*) fBdat->Clone("toymean");
	for (Int_t j = 1; j <= fNdim; j++)
		toymean->SetBinContent(j, 0.);

	// Get the mean of the toys first
	for (int i = 1; i <= ntoys; i++) {

		// create a vector of unit Gaussian variables
		TVectorD g(fNdim);
		for (Int_t k = 0; k < fNdim; k++)
			g(k) = random.Gaus(0., 1.);

		// Multiply this vector by Lt to introduce the appropriate correlations
		g *= (*Lt);

		// Add the mean value offsets and store the results
		for (int j = 1; j <= fNdim; j++) {
			fToyhisto->SetBinContent(j, fBdat->GetBinContent(j) + g(j - 1));
			fToyhisto->SetBinError(j, fBdat->GetBinError(j));
		}
		unfres = Unfold(GetTau());

		for (Int_t j = 1; j <= fNdim; j++) {
			toymean->SetBinContent(j, toymean->GetBinContent(j) + unfres->GetBinContent(j) / ntoys);
		}
		delete unfres;
		unfres = 0;
	}

	// Reset the random seed
	random.SetSeed(seed);

	//Now the toys for the covariance matrix
	for (int i = 1; i <= ntoys; i++) {

		// Create a vector of unit Gaussian variables
		TVectorD g(fNdim);
		for (Int_t k = 0; k < fNdim; k++)
			g(k) = random.Gaus(0., 1.);

		// Multiply this vector by Lt to introduce the appropriate correlations
		g *= (*Lt);

		// Add the mean value offsets and store the results
		for (int j = 1; j <= fNdim; j++) {
			fToyhisto->SetBinContent(j, fBdat->GetBinContent(j) + g(j - 1));
			fToyhisto->SetBinError(j, fBdat->GetBinError(j));
		}
		unfres = Unfold(GetTau());

		for (Int_t j = 1; j <= fNdim; j++) {
			for (Int_t k = 1; k <= fNdim; k++) {
				unfcov->SetBinContent(j, k,
						unfcov->GetBinContent(j, k)
								+ ((unfres->GetBinContent(j) - toymean->GetBinContent(j))
										* (unfres->GetBinContent(k) - toymean->GetBinContent(k)) / (ntoys - 1)));
			}
		}
		delete unfres;
		unfres = 0;
	}
	delete Lt;
	delete toymean;
	fToyMode = kFALSE;

	return unfcov;
}

TH2D* TauSVDUnfold::GetAdetCovMatrix(Int_t ntoys, Int_t seed = 1) {
	// Determine covariance matrix of unfolded spectrum from finite statistics in
	// response matrix using pseudo experiments
	// "ntoys"  - number of pseudo experiments used for the propagation
	// "seed"   - seed for pseudo experiments

	fMatToyMode = true;
	TH1D* unfres = 0;
	TH2D* unfcov = (TH2D*) fAdet->Clone("unfcovmat");
	unfcov->SetTitle("Toy covariance matrix");
	for (int i = 1; i <= fNdim; i++)
		for (int j = 1; j <= fNdim; j++)
			unfcov->SetBinContent(i, j, 0.);

	//Now the toys for the detector response matrix
	TRandom3 random(seed);

	fToymat = (TH2D*) fAdet->Clone("toymat");
	TH1D *toymean = (TH1D*) fXini->Clone("toymean");
	for (Int_t j = 1; j <= fNdim; j++)
		toymean->SetBinContent(j, 0.);

	for (int i = 1; i <= ntoys; i++) {
		for (Int_t k = 1; k <= fNdim; k++) {
			for (Int_t m = 1; m <= fNdim; m++) {
				if (fAdet->GetBinContent(k, m)) {
					fToymat->SetBinContent(k, m, random.Poisson(fAdet->GetBinContent(k, m)));
				}
			}
		}

		unfres = Unfold(GetTau());

		for (Int_t j = 1; j <= fNdim; j++) {
			toymean->SetBinContent(j, toymean->GetBinContent(j) + unfres->GetBinContent(j) / ntoys);
		}
		delete unfres;
		unfres = 0;
	}

	// Reset the random seed
	random.SetSeed(seed);

	for (int i = 1; i <= ntoys; i++) {
		for (Int_t k = 1; k <= fNdim; k++) {
			for (Int_t m = 1; m <= fNdim; m++) {
				if (fAdet->GetBinContent(k, m))
					fToymat->SetBinContent(k, m, random.Poisson(fAdet->GetBinContent(k, m)));
			}
		}

		unfres = Unfold(GetTau());

		for (Int_t j = 1; j <= fNdim; j++) {
			for (Int_t k = 1; k <= fNdim; k++) {
				unfcov->SetBinContent(j, k,
						unfcov->GetBinContent(j, k)
								+ ((unfres->GetBinContent(j) - toymean->GetBinContent(j))
										* (unfres->GetBinContent(k) - toymean->GetBinContent(k)) / (ntoys - 1)));
			}
		}
		delete unfres;
		unfres = 0;
	}
	delete toymean;
	fMatToyMode = kFALSE;

	return unfcov;
}
