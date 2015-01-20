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
	TVectorD ASV = ASVD.GetSig();
	fASV = ASV;

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

	TVectorD vz = CompProd(vd, vdz);

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
	vx = CompProd(vw, vxini);

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
