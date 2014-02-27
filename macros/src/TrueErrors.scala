/*
 * Cerres (c) 2012-2014 EPFL, Lausanne
 */

package cerres.macros

import scala.Double.{PositiveInfinity => PlusInf}
import scala.Double.{NegativeInfinity => MinusInf}
import ceres.affine._
import ceres.common._

trait TrueErrors extends ASTTranslation {

  def computeTrueErrorAffine(f: AffineForm => AffineForm, fp: AffineForm => AffineForm,
    xn: Double, tol: Double): Interval = {
    import AffineForm._
    assert(tol > 0, "tolerance must be positive")
    if (xn != xn) return Double.NaN

    val xbar = Interval(xn - tol, xn + tol)
    val yDelta: AffineForm = f(new AForm(xn))
    val slope = fp(AffineForm(xbar))

    if (slope.interval.contains(0.0))
      println("Failed to compute bounded error! You may have a multiple root.")

    return (yDelta/slope).interval
  }

  def computeTrueErrorInterval(f: (Interval => Interval), fp: (Interval => Interval),
    xn: Double, tol: Double): Interval = {
    import Interval._
    assert(tol > 0, "tolerance must be positive")
    if (xn != xn) return Double.NaN

    val xbar = Interval(xn - tol, xn + tol)
    val yDelta: Interval = f(Interval(xn))
    val slope: Interval = fp(xbar)

    if (slope.contains(0.0))
      println("Failed to compute bounded error! You may have a multiple root.")

    return yDelta/slope
  }



  def front2DAffine(f: Array[(AffineForm, AffineForm) => AffineForm],
    J: Array[Array[(AffineForm, AffineForm) => AffineForm]], x: DoubleVector, tol: Double): Array[Interval] = {
    val dim = f.length
    val xbar: AffineVector = x +/- tol
    val A = new AffineMatrix(dim)
    for (i <- 0 until dim) {
      for (j <- 0 until dim) {
        A(i,j) = J(i)(j)(xbar(0), xbar(1))
      }
    }
    val b = new AffineVector(dim)
    for (i <- 0 until dim) {
      b(i) = - f(i)(new AForm(x(0)), new AForm(x(1)))
    }
    computeTrueErrorXAffine(b, A, x, tol)
  }

  def front2DInterval(f: Array[(Interval, Interval) => Interval],
    J: Array[Array[(Interval, Interval) => Interval]], x: DoubleVector, tol: Double): Array[Interval] = {
    val dim = f.length
    val xbar: IntervalVector = x plusMinus tol
    val A = new IntervalMatrix(dim)
    for (i <- 0 until dim) {
      for (j <- 0 until dim) {
        A(i,j) = J(i)(j)(xbar(0), xbar(1))
      }
    }
    val b = new IntervalVector(dim)
    for (i <- 0 until dim) {
      b(i) = - f(i)(Interval(x(0)), Interval(x(1)))
    }
    computeTrueErrorXInterval(b, A, x, tol)
  }

  def front3DAffine(f: Array[(AffineForm, AffineForm, AffineForm) => AffineForm],
    J: Array[Array[(AffineForm, AffineForm, AffineForm) => AffineForm]], x: DoubleVector,
    tol: Double): Array[Interval] = {
    val dim = f.length
    val xbar: AffineVector = x +/- tol
    val A = new AffineMatrix(dim)
    for (i <- 0 until dim) {
      for (j <- 0 until dim) {
        A(i,j) = J(i)(j)(xbar(0), xbar(1), xbar(2))
      }
    }
    val b = new AffineVector(dim)
    for (i <- 0 until dim) {
      b(i) = - f(i)(new AForm(x(0)), new AForm(x(1)), new AForm(x(2)))
    }
    computeTrueErrorXAffine(b, A, x, tol)
  }

  def front3DInterval(f: Array[(Interval, Interval, Interval) => Interval],
    J: Array[Array[(Interval, Interval, Interval) => Interval]], x: DoubleVector,
    tol: Double): Array[Interval] = {
    val dim = f.length
    val xbar: IntervalVector = x plusMinus tol
    val A = new IntervalMatrix(dim)
    for (i <- 0 until dim) {
      for (j <- 0 until dim) {
        A(i,j) = J(i)(j)(xbar(0), xbar(1), xbar(2))
      }
    }
    val b = new IntervalVector(dim)
    for (i <- 0 until dim) {
      b(i) = - f(i)(Interval(x(0)), Interval(x(1)), Interval(x(2)))
    }
    computeTrueErrorXInterval(b, A, x, tol)
  }

  def frontXDInterval(f: Array[Array[Interval] => Interval],
    J: Array[Array[Array[Interval] => Interval]], x: DoubleVector, tol: Double): Array[Interval] = {
    val dim = f.length
    val xbar: IntervalVector = x plusMinus tol
    val A = new IntervalMatrix(dim)
    for (i <- 0 until dim) {
      for (j <- 0 until dim) {
        A(i,j) = J(i)(j)(xbar.data)
      }
    }
    val b = new IntervalVector(dim)
    for (i <- 0 until dim) {
      b(i) = - f(i)(x.toIntervalArray)
    }
    computeTrueErrorXInterval(b, A, x, tol)
  }

  def frontXDAffine(f: Array[Array[AffineForm] => AffineForm],
    J: Array[Array[Array[AffineForm] => AffineForm]], x: DoubleVector, tol: Double): Array[Interval] = {
    val dim = f.length
    val xbar: AffineVector = x +/- tol
    val A = new AffineMatrix(dim)
    for (i <- 0 until dim) {
      for (j <- 0 until dim) {
        A(i,j) = J(i)(j)(xbar.data)
      }
    }
    val b = new AffineVector(dim)
    for (i <- 0 until dim) {
      b(i) = - f(i)(x.toAffineArray)
    }
    computeTrueErrorXAffine(b, A, x, tol)
  }
  def computeTrueErrorXAffine(b: AffineVector, A: AffineMatrix, x: DoubleVector, tol: Double): Array[Interval] = {
    val dim = x.data.length
    //val xbar: AffineVector = x +/- tol
    //val A = J.eval(xbar)
    //val b: AffineVector = - f.eval(x)
    val R: DoubleMatrix = A.mid.inverse
    val X = new DoubleVector(dim) +/- tol
    val I = identityMatrix(dim)
    val z = R*b + (I - R*A)*X

    for (aa <- z.data)
      if (aa == FullForm || aa == EmptyForm)
        println("Failed to compute bounded error! You may have a multiple root.")

    return z.toInterval.data
  }

  def computeTrueErrorXInterval(b: IntervalVector, A: IntervalMatrix, x: DoubleVector, tol: Double): Array[Interval] = {
    val dim = x.data.length
    //val xbar: AffineVector = x +/- tol
    //val A = J.eval(xbar)
    //val b: AffineVector = - f.eval(x)
    val R: DoubleMatrix = A.mid.inverse
    val X = new DoubleVector(dim) plusMinus tol
    val I = intervalIdentityMatrix(dim)
    val z = R*b + (I - R*A)*X

    for (aa <- z.data)
      if (aa == EmptyInterval || aa.xlo == MinusInf || aa.xhi == PlusInf)
        println("Failed to compute bounded error! You may have a multiple root.")

    return z.data
  }

}
