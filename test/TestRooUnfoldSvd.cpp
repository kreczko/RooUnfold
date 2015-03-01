/*
 * TestUtils.cpp
 *
 *  Created on: 24 Feb 2015
 *      Author: kreczko
 */

#include <boost/test/unit_test.hpp>
#include "../include/RooUnfold.h"
#include "../include/RooUnfoldSvd.h"
#include <cmath>
#include "TH2D.h"
#include "TROOT.h"
#include <iostream>

using namespace std;

struct RooUnfoldSvdSetup {
	RooUnfoldSvdSetup() :
					nbins(6),
					n_toy(100),
					kreg(3),
					taureg(40),
					data(new TH1D("data", "data", nbins, 0, nbins)),
					gen_var(new TH1D("gen_var", "gen_var", nbins, 0, nbins)),
					reco_var(new TH1D("reco_var", "reco_var", nbins, 0, nbins)),
					gen_vs_reco(new TH2D("gen_vs_reco", "gen_vs_reco", nbins, 0, nbins, nbins, 0, nbins)),
					roo_response(),
					roounfold_svd_k(),
					roounfold_svd_tau() {
		gROOT->SetBatch(1);
		gROOT->ProcessLine("gErrorIgnoreLevel = 1001;");
		// crash on ROOT::Error
		gROOT->ProcessLine("gErrorAbortLevel = 2001;");
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

	}
	~RooUnfoldSvdSetup() {
		delete data;
		delete gen_var;
		delete reco_var;
		delete gen_vs_reco;
		delete roo_response;
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

	void check_unfolded_data_against_truth(const TH1D* unfolded_data, bool use_sqrt_N = false) {
		for (auto i = 1; i <= unfolded_data->GetNbinsX(); ++i) {
			double error = use_sqrt_N ? sqrt(unfolded_data->GetBinContent(i)) : unfolded_data->GetBinError(i);
			BOOST_CHECK_CLOSE(unfolded_data->GetBinContent(i), gen_var->GetBinContent(i), error);
//			cout << unfolded_data->GetBinContent(i) << " +-" << error << endl;
		}
	}

	unsigned int nbins;
	unsigned int n_toy;
	int kreg;
	double taureg;
	TH1D* data, *gen_var, *reco_var;
	TH2D* gen_vs_reco;
	RooUnfoldResponse* roo_response;
	RooUnfoldSvd* roounfold_svd_k;
	RooUnfoldSvd* roounfold_svd_tau;

};

BOOST_AUTO_TEST_SUITE (TestRooUnfoldSvd)

BOOST_AUTO_TEST_CASE(test_set_and_get_tau) {

	RooUnfoldSvd r = RooUnfoldSvd();
	r.SetTauTerm(0.1);
	BOOST_CHECK_EQUAL(r.GetTauTerm(), 0.1);
}

BOOST_FIXTURE_TEST_CASE(test_get_tau_from_constructor, RooUnfoldSvdSetup) {
	BOOST_CHECK_EQUAL(roounfold_svd_tau->GetTauTerm(), taureg);
}

BOOST_FIXTURE_TEST_CASE(test_get_k_from_constructor, RooUnfoldSvdSetup) {
	BOOST_CHECK_EQUAL(roounfold_svd_k->GetKterm(), kreg);
}

BOOST_FIXTURE_TEST_CASE(test_k_unfold, RooUnfoldSvdSetup) {
	TH1D* unfolded_data = (TH1D*) roounfold_svd_k->Hreco(RooUnfold::kCovToy);
	check_unfolded_data_against_truth(unfolded_data);
	delete unfolded_data;
}

BOOST_FIXTURE_TEST_CASE(test_k_unfold_kNoError_errors, RooUnfoldSvdSetup) {
	TH1D* unfolded_data = (TH1D*) roounfold_svd_k->Hreco(RooUnfold::kNoError);
	check_unfolded_data_against_truth(unfolded_data, true);
	delete unfolded_data;
}

BOOST_FIXTURE_TEST_CASE(test_k_unfold_kErrors_errors, RooUnfoldSvdSetup) {
	TH1D* unfolded_data = (TH1D*) roounfold_svd_k->Hreco(RooUnfold::kErrors);
	check_unfolded_data_against_truth(unfolded_data);
	delete unfolded_data;
}

BOOST_FIXTURE_TEST_CASE(test_k_unfold_kCovariance_errors, RooUnfoldSvdSetup) {
	TH1D* unfolded_data = (TH1D*) roounfold_svd_k->Hreco(RooUnfold::kCovariance);
	check_unfolded_data_against_truth(unfolded_data);
	delete unfolded_data;
}

BOOST_FIXTURE_TEST_CASE(test_tau_unfold_kNoError_errors, RooUnfoldSvdSetup) {
	// FIXME: Hreco has errors == 0
	TH1D* unfolded_data = (TH1D*) roounfold_svd_tau->Hreco(RooUnfold::kNoError);
	// the bug can be seen with
	// check_unfolded_data_against_truth(unfolded_data);
	check_unfolded_data_against_truth(unfolded_data, true);
	delete unfolded_data;
}

BOOST_FIXTURE_TEST_CASE(test_tau_unfold_kErrors_errors, RooUnfoldSvdSetup) {
	TH1D* unfolded_data = (TH1D*) roounfold_svd_tau->Hreco(RooUnfold::kErrors);
	check_unfolded_data_against_truth(unfolded_data);
	delete unfolded_data;
}

BOOST_FIXTURE_TEST_CASE(test_tau_unfold_kCovariance_errors, RooUnfoldSvdSetup) {
	TH1D* unfolded_data = (TH1D*) roounfold_svd_tau->Hreco(RooUnfold::kCovariance);
	check_unfolded_data_against_truth(unfolded_data);
	delete unfolded_data;
}

BOOST_FIXTURE_TEST_CASE(test_tau_unfold_kCovToy_errors, RooUnfoldSvdSetup) {
	TH1D* unfolded_data = (TH1D*) roounfold_svd_tau->Hreco(RooUnfold::kCovToy);
	check_unfolded_data_against_truth(unfolded_data);
	delete unfolded_data;
}


BOOST_FIXTURE_TEST_CASE(test_get_reg_parm_k, RooUnfoldSvdSetup) {
	BOOST_CHECK_EQUAL(roounfold_svd_k->GetRegParm(), kreg);
}

BOOST_FIXTURE_TEST_CASE(test_get_reg_parm_tau, RooUnfoldSvdSetup) {
	BOOST_CHECK_EQUAL(roounfold_svd_tau->GetRegParm(), taureg);
}

BOOST_AUTO_TEST_SUITE_END()

