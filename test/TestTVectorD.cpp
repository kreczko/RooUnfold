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
#define N 3
using namespace std;

BOOST_AUTO_TEST_SUITE (TVectorDTestSuite)
BOOST_AUTO_TEST_CASE(test_assign) {
	double elements[N] = { 1., 2., 3. };
	TVectorD vector = TVectorD(N, elements);
	BOOST_CHECK_EQUAL(vector.GetNoElements(), N);

	// initialise vector with default constructor
	TVectorD copy = TVectorD();
	copy.Clear();
	// does not work without resize
	copy.ResizeTo(N);
	// this returns "vectors not compatible" without copy.ResizeTo(n);
	copy = vector;
	BOOST_CHECK_EQUAL(copy.GetNoElements(), N);
	BOOST_CHECK_EQUAL(copy[0], 1.);
	BOOST_CHECK_EQUAL(copy[1], 2.);
	BOOST_CHECK_EQUAL(copy[2], 3.);
}

BOOST_AUTO_TEST_SUITE_END()
