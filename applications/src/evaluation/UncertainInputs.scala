/*
 * Cerres (c) 2012-2014 EPFL, Lausanne
 */

package cerres.evaluation

import cerres.algorithms.NewtonsMethod._
import ceres.common.{QuadDouble => QD}
import ceres.common._
import ceres.affine._
import cerres.macros._
import cerres.macros.Macros._
import ceres.smartfloat.SmartFloat

import scala.{Double => D}

object UncertainInputs extends App {
  import math._

  //example1
  /*val tol = 1e-9
  val T = 300; val a = 0.401; val b = 42.7e-6;
  val p = 3.5e7; val k = 1.3806503e-23; val x0 = 0.1
  var N = 1000
  val f = (V: Double) => (p + a * (N / V) * (N / V)) * (V - N * b) - k * N * T
  val fPrime = (V: Double) => -2*a * (N*N)/(V*V*V) * (V - N*b) + p + a * (N/V) * (N/V)
  var V = computeRoot(f, fPrime, x0, tol)
  println("N = " + N + ", V = " + V)
  V = computeRoot(f, derivative(f), x0, tol)
  println("N = " + N + ", V = " + V)


  N = 1005
  V = computeRoot(f, fPrime, x0, tol)
  println("\nN = " + N + ", V = " + V)
  V = computeRoot(f, derivative(f), x0, tol)
  println("N = " + N + ", V = " + V)


  N = 995
  V = computeRoot(f, fPrime, x0, tol)
  println("\nN = " + N + ", V = " + V)
  V = computeRoot(f, derivative(f), x0, tol)
  println("N = " + N + ", V = " + V)
*/

  Interval.doubleString = "%1.16e"

  //println()
  example2


  def example1 = {
    var tol = 1e-10
    val a1 = 10.0; val a2 = 13.0; val a3 = 8.0; val a4 = 10.0
    val b = 0.4 +/- 0.01;
    val x0 = -0.1
    val f = (x: Double) =>
      a1/a2 * cos(b.mid) - a1/a4 * cos(x) - cos(b.mid - x) + (a1*a1 + a2*a2 - a3*a3 + a4*a4)/(2*a2*a4)
    val fPrime = (x: Double) => (a1/a4) * sin(x) - sin(b.mid - x)
    val root = computeRoot(f, fPrime, x0, tol)
    println("root : " + root)

    val error = errorBound(f, root, 0.03)
    println(error)
  }

  def example2 = {
    val tol = 1e-9
    val T = 300; val a = 0.401; val b = 42.7e-6;
    val p = 3.5e7; val k = 1.3806503e-23; val x0 = 0.1
    val N = 1000 +/- 5
    val Nm = N.mid
    val f = (V: Double) => (p + a * (N.mid / V) * (N.mid / V)) * (V - N.mid * b) - k * N.mid * T
    val V = computeRoot(f, derivative(f), x0, tol)
    println("N = " + N + ", V = " + V)
    val Vcert = certify(V, assertBound(f, V, 0.0005))
    println("certified V: " + Vcert.interval)
  }

  def minPositionToWall: SmartFloat = new SmartFloat(3.4)
  val l = 0.3

  def pendulumExample = {
    import SmartFloat.{sin => _sin, certainly}
    val k = 0.3; val tol = 1e-13; val x0 = Array(0.18, 0.25)
    val f1 = (x1: D, x2: D) => tan(x1) - k * (2*sin(x1) + sin(x2))
    val f2 = (x1: D, x2: D) => tan(x2) - 2*k * (sin(x1) + sin(x2))
    val r = computeRoot(Array(f1, f2), jacobian(f1, f2), x0, tol)
    val roots = certify(r, errorBound(f1, f2, r(0), r(1), tol))

    val L = _sin(roots(0)) * l + _sin(roots(1)) * l

    if (certainly(L <= minPositionToWall))
      // continue with computation as it's safe
      println("We're good")
  }
}
