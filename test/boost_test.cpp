#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE RooUnfoldTests
#include <boost/test/unit_test.hpp>
#include "TH1.h"
#include "TROOT.h"

struct GlobalSetup {
	GlobalSetup() {
		TH1::AddDirectory(false);
		gROOT->SetBatch(1);
		gROOT->ProcessLine("gErrorIgnoreLevel = 1001;");
		gROOT->ProcessLine("gErrorAbortLevel = 1001;");
	}
	~GlobalSetup() {
	}
};

BOOST_GLOBAL_FIXTURE(GlobalSetup);
