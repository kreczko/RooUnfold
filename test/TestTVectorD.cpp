/*
 * File:   TauSVDUnfoldTests.cpp
 * Author: ale
 *
 * Created on August 13, 2014, 6:49 PM
 */
#include <boost/test/unit_test.hpp>
#include <iostream>
#include "TVectorD.h"
#include <exception>

using namespace std;

BOOST_AUTO_TEST_SUITE (TVectorDTestSuite)
BOOST_AUTO_TEST_CASE(test_assign) {
	uint32_t n(3);
	double elements[n] = { 1., 2., 3. };
	TVectorD vector = TVectorD(n, elements);
	BOOST_CHECK_EQUAL(vector.GetNoElements(), n);

	// initialise vector with default constructor
	TVectorD copy = TVectorD();
	copy.Clear();
	// does not work without resize
	copy.ResizeTo(n);
	// this returns "vectors not compatible" without copy.ResizeTo(n);
	copy = vector;
	BOOST_CHECK_EQUAL(copy.GetNoElements(), n);
	BOOST_CHECK_EQUAL(copy[0], 1.);
	BOOST_CHECK_EQUAL(copy[1], 2.);
	BOOST_CHECK_EQUAL(copy[2], 3.);
}

BOOST_AUTO_TEST_SUITE_END()
