/*
 * Cerres (c) 2012-2014 EPFL, Lausanne
 */

package cerres.algorithms

import ceres.common.{QuadDouble => QD}


object FixedPointIteration extends ErrorEstimates {


  def iterate(phi: Double => Double, start: Double, tol: Double): Double = {
    var xn = start
    var xn_1 = start + 2 * tol
    var count = 0
    while (math.abs(xn - xn_1) > tol && count < 100) {
      xn_1 = xn
      xn = phi(xn)
      count += 1
    }
    if (count == 100) { println("Diverged!"); return Double.NaN }
    return xn
  }

  def iterateQD(phi: QD => QD, start: QD, tol: QD): QD = {
    var xn = start
    var xn_1 = start + 2 * tol
    var count = 0
    while (QD.abs(xn - xn_1) > tol && count < 100) {
      xn_1 = xn
      xn = phi(xn)
      count += 1
    }
    if (count == 100) { println("Diverged!"); return Double.NaN }
    return xn
  }


  def iterateWithEstimates(phi: Double => Double, start: Double, tol: Double, L: Double):
    (Double, Double, Double) = {
    var x0 = start
    var x1 = 0.0

    var xn = start
    var xn_1 = start + 2 * tol
    var count = 0
    while (math.abs(xn - xn_1) > tol && count < 100) {
      xn_1 = xn
      xn = phi(xn)
      count += 1
      if (count == 1) x1 = xn
    }

    val estimate1 = computeErrorEstimate1(L, xn, xn_1)
    val estimate2 = computeErrorEstimate2(L, x1, x0, count)

    if (count == 100) { println("Diverged!"); xn = Double.NaN }
    return (xn, estimate1, estimate2)
  }




}
