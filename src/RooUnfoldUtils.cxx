/*
 * RooUnfoldUtils.cxx
 *
 *  Created on: 24 Feb 2015
 *      Author: kreczko
 */
#include "../include/RooUnfoldUtils.h"
#include <cmath>
namespace RooUnfoldUtils {
TVectorD* calculate_tau_scan_points(double tau_min, double tau_max, unsigned int n_points) {
	// Use 3 scan points minimum
	if (n_points < 3)
		n_points = 3;

	// Setup Vector
	TVectorD* points = new TVectorD(n_points);

	// Find the scan points
	// Use equidistant steps on a logarithmic scale
	double steplog = (log10(tau_max) - log10(tau_min)) / ((double) (n_points - 1));
	for (auto i = 0; i < n_points; i++) {
		double posScanPoint = pow(10., (log10(tau_min) + ((double) i) * steplog));
		(*points)[i] = posScanPoint;
	}

	return points;
}
}
