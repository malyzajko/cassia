/*
 * Cerres (c) 2012-2014 EPFL, Lausanne
 */
 
package cerres.evaluation

import cerres.algorithms._
import cerres.macros.Macros

import ceres.affine._
import ceres.common._

import math._

trait ErrorComparison {
  // Comparison of our errors with the error estimates from Contraction Mapping Theorem
  import Macros._

  // Fix this, these are Fixedpoint computations...
  /*def verhulstModelEST(tol: Double): (Double, Double, Double) = {
    val r = 4.0; val K = 1.11; val x0 = 0.1
    val L = lipschitzConstant(), Interval(3.3, 3.4))
    val f = (x: Double) => (r*x) / (1 + (x/K))
    FixedPointIteration.iterateWithEstimates(f, x0, tol, L)
  }

  def predatorPreyModelEST(tol: Double): (Double, Double, Double) = {
    val r = 4.0; val K = 1.11; val x0 = 0.7
    val f = (x: Double) => (r*x*x) / (1 + pow(x/K, 2))
    val L = lipschitzConstant((x: AffineForm) => (r*x*x) / (1 + (x/K)*(x/K)), Interval(4.6, 4.7))
    FixedPointIteration.iterateWithEstimates(f, x0, tol, L)
  }*/

  def rodSystemEST(tol: Double): (Double, Double, Double) = {
    val a1 = 10.0; val a2 = 13.0; val a3 = 8.0; val a4 = 10.0
    val b = 0.4; val x0 = -0.1
    val L = lipschitzConstant((x: AffineForm) => (AffineForm(a1)/a4) * AffineForm.sin(x) - AffineForm.sin(b - x),
      Interval(-0.09, -0.08))
    val f = (x: Double) =>
      a1/a2 * cos(b) - a1/a4 * cos(x) - cos(b - x) + (a1*a1 + a2*a2 - a3*a3 + a4*a4)/(2*a2*a4)
    val fPrime = (x: Double) => (a1/a4) * sin(x) - sin(b - x)
    NewtonsMethod.computeRootWithEstimates(f, fPrime, x0, tol, L)
  }



  private def lipschitzConstant(der: AffineForm => AffineForm, range: Interval): Double = {
    import math._
    val derBound = der(AffineForm(range)).interval
    max(abs(derBound.xlo), abs(derBound.xhi))
  }
}
