/*
 * Cerres (c) 2012-2014 EPFL, Lausanne
 */
 
package cerres.algorithms


trait ErrorEstimates {

  def computeErrorEstimate1(L: Double, xn: Double, xn_1: Double): Double = {
    math.abs(xn - xn_1) * L/(1-L)
  }

  def computeErrorEstimate2(L: Double, x1: Double, x0: Double, k: Int): Double = {
    math.abs(x1 - x0) * math.pow(L, k)/(1-L)
  }

}
