/*
 * Cerres (c) 2012-2014 EPFL, Lausanne
 */

package cerres.algorithms

import ceres.common.{QuadDouble => QD}
import QD._

object SecantMethod extends ErrorEstimates {

  def computeRoot(f: Double => Double, start1: Double, start2: Double,
    tol: Double): Double = {

    var err = 2 * tol
    var xn = start1
    var xn1 = start2
    var count = 0

    while(err > tol && count < 100) {
      val x_new = xn - f(xn) * (xn - xn1)/(f(xn) - f(xn1))
      xn1 = xn
      xn = x_new
      err = math.abs(f(xn))
      count += 1
    }
    if (count == 100) { println("Diverged!"); return Double.NaN }
    return xn
  }

  def computeRootQD(f: QD => QD, start1: QD, start2: QD,
    tol: QD): QD = {

    var err = 2 * tol
    var xn = start1
    var xn1 = start2
    var count = 0

    while(err > tol && count < 100) {
      val x_new = xn - f(xn) * (xn - xn1)/(f(xn) - f(xn1))
      xn1 = xn
      xn = x_new

      println("xn = " + xn + "   xn1 = " + xn1)
      err = abs(f(xn))
      count += 1
    }
    if (count == 100) { println("Diverged!"); return Double.NaN }
    return xn

  }
}
