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
					cov_hist(),
					roo_response(),
					roounfold_svd_k(),
					roounfold_svd_tau() {
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
		set_sqrt_N_error(reco_var);

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

		roounfold_svd_tau->Hreco(RooUnfold::kNoError);
		TH2D* cov = TauSVDUnfold::get_data_covariance_matrix(data);
		cov_hist = roounfold_svd_tau->Impl()->GetUnfoldCovMatrix(cov, n_toy);
		delete cov;
	}
	~TauSVDUnfoldSetup() {
		delete data;
		delete gen_var;
		delete reco_var;
		delete gen_vs_reco;
		delete cov_hist;
		delete roo_response;
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

	uint32_t nbins;
	unsigned int n_toy;
	int kreg;
	double taureg;
	TH1D* data, *gen_var, *reco_var;
	TH2D* gen_vs_reco;
	TH2D* cov_hist;
	RooUnfoldResponse* roo_response;
	RooUnfoldSvd* roounfold_svd_k;
	RooUnfoldSvd* roounfold_svd_tau;

};

BOOST_AUTO_TEST_SUITE (TauSVDUnfoldTestSuite)
BOOST_FIXTURE_TEST_CASE(test_get_data_covariance_matrix, TauSVDUnfoldSetup) {
	TH2D* cov = TauSVDUnfold::get_data_covariance_matrix(data);
	for (uint32_t i = 1; i <= nbins; ++i) {
		double data_error = data->GetBinError(i);
		BOOST_CHECK_EQUAL(cov->GetBinContent(i, i), data_error * data_error);
	}
	delete cov;
}

BOOST_FIXTURE_TEST_CASE(test_get_global_correlation, TauSVDUnfoldSetup) {
	double corr = TauSVDUnfold::get_global_correlation(cov_hist, data);
	BOOST_CHECK_CLOSE(corr, 99, 1);
}

BOOST_FIXTURE_TEST_CASE(test_get_global_correlation_hist, TauSVDUnfoldSetup) {
	TH1D* corr_hist = TauSVDUnfold::get_global_correlation_hist(cov_hist, data);
	// for the above example, the first bin should be around 95 %
	BOOST_CHECK_CLOSE(corr_hist->GetBinContent(1), 95, 1);
	// for all others it should be between 99 and 100 %
	for (uint32_t i = 2; i <= nbins; ++i) {
		BOOST_CHECK_CLOSE(corr_hist->GetBinContent(i), 99, 1);
	}
	delete corr_hist;
}

BOOST_FIXTURE_TEST_CASE(test_k_to_tau, TauSVDUnfoldSetup) {
	const TauSVDUnfold* tau_svd = (TauSVDUnfold*) roounfold_svd_tau->Impl();
	double tau = tau_svd->GetTau();
	BOOST_CHECK_EQUAL(tau, taureg);
	TVectorD ASV = tau_svd->getASV();
	uint32_t n_elements = ASV.GetNoElements();
	BOOST_CHECK_EQUAL(n_elements, nbins + 1);
	tau = tau_svd->kToTau(kreg);
	BOOST_CHECK_CLOSE(tau, 3.83, 0.1);
}

BOOST_FIXTURE_TEST_CASE(test_get_tau, TauSVDUnfoldSetup) {
	double tau = ((TauSVDUnfold*) roounfold_svd_tau->Impl())->GetTau();
	BOOST_CHECK_EQUAL(tau, taureg);
}

BOOST_FIXTURE_TEST_CASE(test_get_ASV, TauSVDUnfoldSetup) {
	const TauSVDUnfold* tau_svd = (TauSVDUnfold*) roounfold_svd_tau->Impl();
	TVectorD ASV = tau_svd->getASV();
	uint32_t n_elements = ASV.GetNoElements();
	BOOST_CHECK_EQUAL(n_elements, nbins + 1);
	BOOST_CHECK(ASV.NonZeros() > 0);
}
BOOST_AUTO_TEST_SUITE_END()
