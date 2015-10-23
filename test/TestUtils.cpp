/*
 * TestUtils.cpp
 *
 *  Created on: 24 Feb 2015
 *      Author: kreczko
 */

#include <boost/test/unit_test.hpp>
#include "../include/RooUnfoldUtils.h"
#include <cmath>

BOOST_AUTO_TEST_SUITE (TestUtils)

BOOST_AUTO_TEST_CASE(test_calculate_tau_scan_points) {

	TVectorD* points = RooUnfoldUtils::calculate_tau_scan_points(0.1, 1000, 4);
	BOOST_CHECK_EQUAL(points->GetNrows(), 4); //right size
	double step_size =  (log10(1000) - log10(0.1)) / 3;

	for (auto i = 0; i < points->GetNrows(); ++i){
		BOOST_CHECK_EQUAL((*points)[i], pow(10, log10(0.1) + i* step_size));
	}
}

BOOST_AUTO_TEST_SUITE_END()


