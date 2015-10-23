/*
 * RooUnfoldUtils.h
 *
 *  Created on: 24 Feb 2015
 *      Author: kreczko
 */

#ifndef INCLUDE_ROOUNFOLDUTILS_H
#define INCLUDE_ROOUNFOLDUTILS_H

#include "TVectorD.h"
namespace RooUnfoldUtils {

/**
 * Calculates equidistant points on a logarithmic scale
 *
 * parameters:
 *  - tau_min: minimal tau value (> 0)
 *  - tau_max: maximal tau value
 *  - n_points: number of points (> 3)
 */
TVectorD* calculate_tau_scan_points(double tau_min, double tau_max, unsigned int n_points);
}
#endif /* INCLUDE_ROOUNFOLDUTILS_H */
