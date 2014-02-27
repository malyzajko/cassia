/*
 * Cerres (c) 2012-2014 EPFL, Lausanne
 */

package cerres.examples

import cerres.algorithms._
import math._
import cerres.macros._

import ceres.affine._
import ceres.common._
import AffineForm.{int2AffineForm, double2AffineForm}
import ceres.common.{QuadDouble => QD}

trait UnaryApplications {
  import Macros._

  def verhulstModel(tol: Double): (Double, Interval) = {
    val r = 4.0; val K = 1.11; val x0 = 0.1
    val x = FixedPointIteration.iterate((x: Double) => (r*x) / (1 + (x/K)), x0, tol)
    val err = errorBound((x: Double) => r / (1 + (x/K)) - 1, x, tol)
    //val err = Interval(0.0)
    return (x, err)
  }

  def predatorPreyModel(tol: Double): (Double, Interval) = {
    val r = 4.0; val K = 1.11; val x0 = 0.7
    val x = FixedPointIteration.iterate((x: Double) => (r*x*x) / (1 + pow(x/K, 2)), x0, tol)
    val error = errorBound((x: Double) => (r*x) / (1 + (x/K)*(x/K)) - 1, x, tol)
    //val error = Interval(0.0)
    return (x, error)
  }

  def rodSystem(tol: Double): (Double, Interval) = {
    val a1 = 10.0; val a2 = 13.0; val a3 = 8.0; val a4 = 10.0
    val b = 0.4; val x0 = -0.1
    val f = (x: Double) =>
      a1/a2 * cos(b) - a1/a4 * cos(x) - cos(b - x) + (a1*a1 + a2*a2 - a3*a3 + a4*a4)/(2*a2*a4)
    val fPrime = (x: Double) => (a1/a4) * sin(x) - sin(b - x)
    val alpha = NewtonsMethod.computeRoot(f, fPrime, x0, tol)
    val error = errorBound(f, alpha, tol)
    //val error = Interval(0.0)

    return (alpha, error)
  }

  def carbonGasState(tol: Double): (Double, Interval) = {
    val T = 300; val a = 0.401; val b = 42.7e-6; val N = 1000
    val p = 3.5e7; val k = 1.3806503e-23; val x0 = 0.1
    val f = (V: Double) => (p + a * (N / V) * (N / V)) * (V - N * b) - k * N * T
    val fPrime = (V: Double) => -2*a * (N*N)/(V*V*V) * (V - N*b) + p + a * (N/V) * (N/V)
    val V = NewtonsMethod.computeRoot(f, fPrime, x0, tol)
    val err = errorBound(f, V, tol)
    //val err = Interval(0.0)

    return (V, err)
  }

  def butlerVolmer(tol: Double): (Double, Interval) = {
    val alpha = 0.2; val beta = 2.0; val x0 = 2.0
    val f = (x: Double) => exp(alpha*x) - exp(-(1 - alpha)*x) - beta
    val fPrime = (x: Double) => alpha * exp(alpha*x) + (1 - alpha) * exp(-(1 - alpha)*x)
    val root = NewtonsMethod.computeRoot(f, fPrime, x0, tol)
    val error = errorBound(f, root, tol)
    //val error = Interval(0.0)
    return (root, error)
  }

  def sinXSquared(tol: Double): (Double, Interval) = {
    val x0 = 1.6
    val f = (x: Double) => (x/2)*(x/2) - sin(x)
    val fPrime = (x: Double) => x/2 - cos(x)
    val root = NewtonsMethod.computeRoot(f, fPrime, x0, tol)
    val error = errorBound(f, root, tol)
    //val error = Interval(0.0)
    return (root, error)
  }

  def exponents(tol: Double): (Double, Interval) = {
    val x0 = 1.7
    val f = (x: Double) => exp(x)*(x - 1) - exp(-x)*(x + 1)
    val fPrime = (x: Double) => exp(x)*x - exp(-x)*x
    val root = NewtonsMethod.computeRoot(f, fPrime, x0, tol)
    val error = errorBound(f, root, tol)
    //val error = Interval(0.0)
    return (root, error)
  }

  def polynomial1(tol: Double): (Double, Interval) = {
    val x0 = 1.2
    val f = (x: Double) => x*x*x/3.0 - 2 * x*x + 4.5
    val fPrime = (x: Double) => x*x - 4*x
    val root = NewtonsMethod.computeRoot(f, fPrime, x0, tol)
    val error = errorBound(f, root, tol)
    return (root, error)
  }

  def polynomial2(tol: Double): (Double, Interval) = {
    val x0 = 6.5
    val f = (x: Double) => pow(x,6) + 4.2*pow(x,5) -72.3*pow(x,4) -214.4*pow(x, 3) + 1127.1*pow(x, 2) + 1602.9*x - 5040.5
    val fPrime = (x: Double) => 6*pow(x,5) + 21.0*pow(x,4) -289.2*pow(x,3) -643.2*pow(x, 2) + 2254.2*x + 1602.9
    val root = NewtonsMethod.computeRoot(f, fPrime, x0, tol)
    val error = errorBound(f, root, tol)
    return (root, error)
  }

  /* ====================================
      VERSION WITH QuadDouble
   ====================================== */
  def verhulstModelQD(tol: QD): QD = {
    val r = QD(4.0); val K = QD(1.11); val x0 = QD(0.1)
    FixedPointIteration.iterateQD((x: QD) => (r*x) / (1 + (x/K)), x0, tol)
  }

  def predatorPreyModelQD(tol: QD): QD = {
    val r = QD(4.0); val K = QD(1.11); val x0 = QD(0.7)
    FixedPointIteration.iterateQD((x: QD) => (r*x*x) / (1 + QD.pow(x/K, 2)), x0, tol)
  }

  def rodSystemQD(tol: QD): QD = {
    import QD._
    val a1 = QD(10.0); val a2 = QD(13.0); val a3 = QD(8.0); val a4 = QD(10.0)
    val b = QD(0.4); val x0 = QD(-0.1)
    val f = (x: QD) =>
      a1/a2 * cos(b) - a1/a4 * cos(x) - cos(b - x) + (a1*a1 + a2*a2 - a3*a3 + a4*a4)/(2*a2*a4)
    val fPrime = (x: QD) => (a1/a4) * sin(x) - sin(b - x)
    NewtonsMethod.computeRootQD(f, fPrime, x0, tol)
  }

  def carbonGasStateQD(tol: QD): QD = {
    val T = QD(300); val a = QD(0.401); val b = QD(42.7e-6); val N = QD(1000)
    val p = QD(3.5e7); val k = QD(1.3806503e-23); val x0 = QD(0.1)
    val f = (V: QD) => (p + a * (N / V) * (N / V)) * (V - N * b) - k * N * T
    val fPrime = (V: QD) => -2*a * (N*N)/(V*V*V) * (V - N*b) + p + a * (N/V) * (N/V)
    NewtonsMethod.computeRootQD(f, fPrime, x0, tol)
  }

  def butlerVolmerQD(tol: QD): QD = {
    import QD._
    val alpha = QD(0.2); val beta = QD(2.0); val x0 = QD(2.0)
    val f = (x: QD) => exp(alpha*x) - exp(-(1 - alpha)*x) - beta
    val fPrime = (x: QD) => alpha * exp(alpha*x) + (1 - alpha) * exp(-(1 - alpha)*x)
    NewtonsMethod.computeRootQD(f, fPrime, x0, tol)
  }

  def sinXSquaredQD(tol: QD): QD = {
    val x0 = QD(1.6)
    val f = (x: QD) => (x/2)*(x/2) - QD.sin(x)
    val fPrime = (x: QD) => x/2 - QD.cos(x)
    NewtonsMethod.computeRootQD(f, fPrime, x0, tol)
  }

  def exponentsQD(tol: QD): QD = {
    import QD._
    val x0 = QD(1.7)
    val f = (x: QD) => exp(x)*(x - 1) - exp(-x)*(x + 1)
    val fPrime = (x: QD) => exp(x)*x - exp(-x)*x
    NewtonsMethod.computeRootQD(f, fPrime, x0, tol)
  }

  def polynomial1QD(tol: QD): QD = {
    import QD._
    val x0 = QD(1.2)
    val f = (x: QD) => x*x*x/3.0 - 2 * x*x + 4.5
    val fPrime = (x: QD) => x*x - 4*x
    NewtonsMethod.computeRootQD(f, fPrime, x0, tol)
  }

  def polynomial2QD(tol: QD): QD = {
    import QD._
    val x0 = QD(6.5)
    val f = (x: QD) => pow(x,6) + 4.2*pow(x,5) -72.3*pow(x,4) -214.4*pow(x, 3) + 1127.1*pow(x, 2) + 1602.9*x - 5040.5
    val fPrime = (x: QD) => 6*pow(x,5) + 21.0*pow(x,4) -289.2*pow(x,3) -643.2*pow(x, 2) + 2254.2*x + 1602.9
    NewtonsMethod.computeRootQD(f, fPrime, x0, tol)
  }

}
