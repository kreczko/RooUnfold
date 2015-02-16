#define BOOST_TEST_MODULE SimpleTest
#include <boost/test/included/unit_test.hpp>
int add(int i, int j) {
	return i + j;
}

BOOST_AUTO_TEST_SUITE(SimpleTest)
BOOST_AUTO_TEST_CASE( one_and_one ) {
	BOOST_REQUIRE(add(1, 1) == 2);
}
BOOST_AUTO_TEST_SUITE_END( )
