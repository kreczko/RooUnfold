/*
 * File:   TauSVDUnfoldTests.cpp
 * Author: ale
 *
 * Created on August 13, 2014, 6:49 PM
 */
#include <boost/test/unit_test.hpp>
#include "../include/TauSVDUnfold.h"
#include "../include/RooUnfoldSvd.h"
#include "../include/RooUnfoldResponse.h"
#include <iostream>
#include "TROOT.h"

using namespace std;

struct TauSVDUnfoldSetup {
	TauSVDUnfoldSetup() :
					nbins(6),
					n_toy(100),
					kreg(3),
					taureg(1),
					data(new TH1D("data", "data", nbins, 0, nbins)),
					gen_var(new TH1D("gen_var", "gen_var", nbins, 0, nbins)),
					reco_var(new TH1D("reco_var", "reco_var", nbins, 0, nbins)),
					gen_vs_reco(new TH2D("gen_vs_reco", "gen_vs_reco", nbins, 0, nbins, nbins, 0, nbins)),
					roo_response(),
					roounfold_svd_k(),
					roounfold_svd_tau(){
		gROOT->SetBatch(1);
//		gROOT->ProcessLine("gErrorIgnoreLevel = 1001;");
		gROOT->ProcessLine("gErrorAbortLevel = 1001;");
		// from toy MC 1
		data->SetBinContent(1, 365);
		data->SetBinContent(2, 578);
		data->SetBinContent(3, 531);
		data->SetBinContent(4, 198);
		data->SetBinContent(5, 58);
		data->SetBinContent(6, 31);
		set_sqrt_N_error(data);

		// from toy MC 2
		reco_var->SetBinContent(1, 367);
		reco_var->SetBinContent(2, 586);
		reco_var->SetBinContent(3, 523);
		reco_var->SetBinContent(4, 195);
		reco_var->SetBinContent(5, 61);
		reco_var->SetBinContent(6, 29);
		set_sqrt_N_error (reco_var);

		gen_var->SetBinContent(1, 3441);
		gen_var->SetBinContent(2, 5181);
		gen_var->SetBinContent(3, 4050);
		gen_var->SetBinContent(4, 1546);
		gen_var->SetBinContent(5, 449);
		gen_var->SetBinContent(6, 149);
		set_sqrt_N_error(gen_var);

		gen_vs_reco->SetBinContent(1, 1, 236);
		gen_vs_reco->SetBinContent(1, 2, 92);
		gen_vs_reco->SetBinContent(1, 3, 4);

		gen_vs_reco->SetBinContent(2, 1, 163);
		gen_vs_reco->SetBinContent(2, 2, 297);
		gen_vs_reco->SetBinContent(2, 3, 61);
		gen_vs_reco->SetBinContent(2, 4, 1);

		gen_vs_reco->SetBinContent(3, 1, 24);
		gen_vs_reco->SetBinContent(3, 2, 167);
		gen_vs_reco->SetBinContent(3, 3, 216);
		gen_vs_reco->SetBinContent(3, 3, 23);

		gen_vs_reco->SetBinContent(4, 1, 2);
		gen_vs_reco->SetBinContent(4, 2, 4);
		gen_vs_reco->SetBinContent(4, 3, 61);
		gen_vs_reco->SetBinContent(4, 4, 74);
		gen_vs_reco->SetBinContent(4, 5, 5);

		gen_vs_reco->SetBinContent(5, 4, 12);
		gen_vs_reco->SetBinContent(5, 5, 26);
		gen_vs_reco->SetBinContent(5, 6, 3);

		gen_vs_reco->SetBinContent(6, 4, 1);
		gen_vs_reco->SetBinContent(6, 5, 5);
		gen_vs_reco->SetBinContent(6, 6, 11);
		set_sqrt_N_error(gen_vs_reco);

		roo_response = new RooUnfoldResponse(reco_var, gen_var, gen_vs_reco);
		roounfold_svd_k = new RooUnfoldSvd(roo_response, data, kreg, n_toy);
		roounfold_svd_tau = new RooUnfoldSvd(roo_response, data, taureg, n_toy);

	}
	~TauSVDUnfoldSetup() {
		delete data;
		delete gen_var;
		delete reco_var;
		delete gen_vs_reco;
		delete roo_response;
//		delete tau_svd_unfold;
		delete roounfold_svd_k;
		delete roounfold_svd_tau;
	}

	void set_sqrt_N_error(TH1D* hist) {
		for (auto i = 1; i <= hist->GetNbinsX(); ++i) {
			hist->SetBinError(i, sqrt(hist->GetBinContent(i)));
		}
	}

	void set_sqrt_N_error(TH2D* hist) {
		for (auto i = 1; i <= hist->GetNbinsX(); ++i) {
			for (auto j = 1; j <= hist->GetNbinsY(); ++j) {
				hist->SetBinError(i, j, sqrt(hist->GetBinContent(i, j)));
			}
		}
	}

	unsigned int nbins;
	unsigned int n_toy;
	int kreg;
	double taureg;
	TH1D* data, *gen_var, *reco_var;
	TH2D* gen_vs_reco;
	RooUnfoldResponse* roo_response;
//	TauSVDUnfold* tau_svd_unfold;
	RooUnfoldSvd* roounfold_svd_k;
	RooUnfoldSvd* roounfold_svd_tau;

};

BOOST_AUTO_TEST_SUITE (TauSVDUnfoldTestSuite)
BOOST_FIXTURE_TEST_CASE(test_get_data_covariance_matrix, TauSVDUnfoldSetup) {
	TH2D* cov = TauSVDUnfold::get_data_covariance_matrix(data);
	for (auto i = 1; i <= nbins; ++i) {
		double data_error = data->GetBinError(i);
		BOOST_CHECK_EQUAL(cov->GetBinContent(i, i), data_error * data_error);
	}
	delete cov;
}

BOOST_FIXTURE_TEST_CASE(test_get_global_correlation, TauSVDUnfoldSetup) {
	double corr = TauSVDUnfold::get_global_correlation(data);
	BOOST_CHECK_EQUAL(corr, -1.);
}

BOOST_FIXTURE_TEST_CASE(test_get_global_correlation_hist, TauSVDUnfoldSetup) {
	TH1D* corr_hist = TauSVDUnfold::get_global_correlation_hist(data);
	for (auto i = 1; i <= nbins; ++i) {
		BOOST_CHECK_EQUAL(corr_hist->GetBinContent(i), i * i);
	}
	delete corr_hist;
}

BOOST_AUTO_TEST_SUITE_END()
