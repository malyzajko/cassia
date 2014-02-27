/*
 * Cerres (c) 2012-2014 EPFL, Lausanne
 */

package cerres.algorithms

import ceres.common.{QuadDouble => QD}
import ceres.affine._
import ceres.common._
import AffineForm._

import cerres.macros.Macros._

object NewtonsMethod extends ErrorEstimates {

  /**
    Computes a root of the function f near the starting value using Newton's method.
    Stops when |f(xn)| becomes smaller than the given tolerance.
   */
  def computeRoot(f: Double => Double, fprime: Double => Double, start: Double,
    tol: Double): Double = {

    var err = 2 * tol
    var xn = start
    var count = 0
    while(err > tol && count < 100) {
      xn = xn - f(xn)/fprime(xn)
      err = math.abs(f(xn))
      count += 1
    }

    if (count == 100) { println("Diverged!"); return Double.NaN }
    return xn
  }

  /**
    Solves a system of n nonlinear equations with Newton's method.
    @param f functions to be solved
    @param J Jacobian
    @param x0 starting points

    @return solutions
  */
  def computeRoot(f: Array[(Double, Double) => Double], J: Array[Array[(Double, Double) => Double]],
    x0: Array[Double], tol: Double): Array[Double] = {
    assert(f.length == 2 && J.length == 2 && J(0).length == 2 && x0.length == 2, "Wrong dimension")
    assert(tol >= 0, "Tolerance must be given as a positive number")
    var xn = x0

    var error = 2 * tol
    var count = 0
    var f_xn = evaluateVector(f, xn)

    while (error > tol && count < 100) {
      val partialDer = evaluateMatrix(J, xn)
      val delta = solveSystem(partialDer, negateVector(f_xn))
      xn = addVectors(xn, delta)
      f_xn = evaluateVector(f, xn)
      error = norm(f_xn)
      count += 1
    }
    if (count == 100) { println("Diverged!"); return Array.fill(x0.length)(Double.NaN) }
    xn
  }

  /**
    Solves a system of n nonlinear equations with Newton's method.
    @param f functions to be solved
    @param J Jacobian
    @param x0 starting points

    @return solutions
  */
  def computeRoot3(f: Array[(Double, Double, Double) => Double], J: Array[Array[(Double, Double, Double) => Double]],
    x0: Array[Double], tol: Double): Array[Double] = {
    assert(f.length == 3 && J.length == 3 && J(0).length == 3 && x0.length == 3, "Wrong dimension")
    assert(tol >= 0, "Tolerance must be given as a positive number")

    var xn = new DoubleVector(3, x0)
    var error = 2 * tol
    var count = 0
    var f_xn = evalVector(f, xn.data)

    while (error > tol && count < 100) {
      val partialDer = evaluateMatrix(J, xn.data)
      val delta = new DoubleVector(3, solveSystem(partialDer.data, (-f_xn).data))
      xn = xn + delta
      f_xn = evaluateVector3(f, xn.data)
      error = f_xn.norm
      count += 1
    }
    if (count == 100) { println("Diverged!"); return Array.fill(x0.length)(Double.NaN) }
    xn.data
  }

 /**
    Solves a system of n nonlinear equations with Newton's method.
    @param f functions to be solved
    @param J Jacobian
    @param x0 starting points

    @return solutions
  */
  def computeRootQD(f: Array[(QD, QD) => QD], J: Array[Array[(QD, QD) => QD]],
    x0: Array[QD], tol: QD): Array[QD] = {

    var xn = x0

    var error = 2 * tol
    var count = 0
    var f_xn = evaluateVectorQD(f, xn)

    while (error > tol && count < 100) {
      val partialDer = evaluateMatrixQD(J, xn)
      val delta = solveSystemQD(partialDer, negateVectorQD(f_xn))
      xn = addVectorsQD(xn, delta)
      f_xn = evaluateVectorQD(f, xn)
      error = normQD(f_xn)
      count += 1
    }
    if (count == 100) { println("Diverged!"); return Array.fill(x0.length)(QD(Double.NaN)) }
    xn
  }

  def computeRoot3QD(f: Array[(QD, QD, QD) => QD], J: Array[Array[(QD, QD, QD) => QD]],
    x0: Array[QD], tol: QD): Array[QD] = {

    var xn = x0

    var error = 2 * tol
    var count = 0
    var f_xn = evaluateVectorQD(f, xn)

    while (error > tol && count < 100) {
      val partialDer = evaluateMatrixQD(J, xn)
      val delta = solveSystemQD(partialDer, negateVectorQD(f_xn))
      xn = addVectorsQD(xn, delta)
      f_xn = evaluateVectorQD(f, xn)
      error = normQD(f_xn)
      count += 1
    }
    if (count == 100) { println("Diverged!"); return Array.fill(x0.length)(QD(Double.NaN)) }
    xn
  }


  /**
    Computes the root of f near start, and uses the Lipschitz constant L to compute
    two error estimates given from the Contraction Mapping Theorem.
    @return (root, estimate based on last two values, estimate based on first two values)
  */
  def computeRootWithEstimates(f: Double => Double, fprime: Double => Double, start: Double,
    tol: Double, L: Double): (Double, Double, Double) = {

    var x0 = start
    var x1 = 0.0
    var xn_1 = 0.0

    var err = 2 * tol
    var xn = start
    var count = 0
    while(err > tol && count < 100) {
      xn_1 = xn
      xn = xn - f(xn)/fprime(xn)
      err = math.abs(f(xn))
      count += 1
      if (count == 1) x1 = xn
    }
    val estimate1 = computeErrorEstimate1(L, xn, xn_1)
    val estimate2 = computeErrorEstimate2(L, x1, x0, count)

    if (count == 100) { println("Diverged!"); xn = Double.NaN }
    return (xn, estimate1, estimate2)
  }

  def computeRootQD(f: QD => QD, fprime: QD => QD, start: QD,
    tol: QD): QD = {
    import QD._

    var err = 2 * tol
    var xn = start
    var count = 0
    while(err > tol && count < 100) {
      xn = xn - f(xn)/fprime(xn)
      err = abs(f(xn))
      count += 1
    }
    if (count == 100) { println("Diverged!"); return Double.NaN }
    return xn
  }

}
