/*
 * Cerres (c) 2012-2014 EPFL, Lausanne
 */

import cerres.macros._
import ceres.common.{QuadDouble => QD}

import ceres.smartfloat.SmartFloat

import ceres.common._
import ceres.affine._

/**
 * Regression tests for macros.
 */
object MacroTest extends App {
  import Macros._
  import math._

  /*testDerivative1
  testDerivative2
  testDerivative3*/
  testDerivativeN

  /*
  testUnary
  testBinary
  testTernary
  testAnyN*/

  //testFunctionTransform

  def testFunctionTransform = {
    val f1 = (x: Double) => sin(x + 0.4) + log(3*x)/(x*x)
    val f2 = (x: Double, y: Double) => sin(x + y) / log(x*x) + 1

    val f1Integ = double2Interval(f1)
    println(f1Integ(Interval(0.3)).mid + "   " + f1(0.3))
    println(f1Integ(Interval(0.45)).mid + "   " + f1(0.45))
    val f2Integ = double2Interval(f2)
    println(f2Integ(Interval(0.3), Interval(0.75)).mid + "   " + f2(0.3, 0.75))
    println(f2Integ(Interval(0.45), Interval(0.98)).mid + "   " + f2(0.45, 0.98))

    val f1Aff = double2Affine(f1)
    println(f1Aff(AffineForm(0.3)).x0 + "   " + f1(0.3))
    println(f1Aff(AffineForm(0.45)).x0 + "   " + f1(0.45))
    val f2Aff = double2Affine(f2)
    println(f2Aff(AffineForm(0.3), AffineForm(0.75)).x0 + "   " + f2(0.3, 0.75))
    println(f2Aff(AffineForm(0.45), AffineForm(0.98)).x0 + "   " + f2(0.45, 0.98))
  }


  def testDerivative1 = {
    val a1 = 10.0; val a2 = 13.0; val a3 = 8.0; val a4 = 10.0; val b = 0.4
    val f = (x: Double) =>
      a1/a2 * cos(b) - a1/a4 * cos(x) - cos(b - x) +
      (a1*a1 + a2*a2 - a3*a3 + a4*a4)/(2*a2*a4)
    val fPrime = derivative(f)
    println("\ntestDerivative1")
    println(fPrime(0.356) == 0.30454203373225586)
  }

  def testDerivative2 = {
    val f1 = (x1: Double, x2: Double) => x1 * cos(x2) + x2 * cos(x1) - 0.9
    val f2 = (x1: Double, x2: Double) => x1 * sin(x2) + x2 * sin(x1) - 0.1

    val j = jacobian(f1, f2)
    println("\ntestDerivative2")
    println(j(0)(0)(0.33, 0.33) == 0.8391081441580804)
    println(j(0)(1)(0.25, 0.37) == 0.8785085637194042)
    println(j(1)(0)(0.33, 0.33) == 0.636237001759236)
    println(j(1)(1)(0.33, 0.77) == 0.5609535493664797)
  }

  def testDerivative3 = {
    import math._
    import scala.{Double => D}
    val f1 = (v: D, w: D, r: D) => 3 + 2/(r*r) - 0.125*(3-2*v)*(w*w*r*r)/(1-v) - 4.5
    val f2 = (v: D, w: D, r: D) => 6*v - 0.5 * v * (w*w*r*r) / (1-v) - 2.5
    val f3 = (v: D, w: D, r: D) => 3 - 2/(r*r) - 0.125 * (1+2*v) * (w*w*r*r) / (1-v) - 0.5
    val j = jacobian(f1, f2, f3)

    println("\ntestDerivative3")
    println(j(0)(0)(0.33, 0.33, 0.33)  == -0.00330229728224549)
    println(j(0)(1)(0.25, 0.37, 0.264) == -0.021489600000000005)
    println(j(0)(2)(0.25, 0.37, 0.001) == -4.000000000000114E9)

    println(j(1)(0)(0.33, 0.33, 0.33)  == 5.986790810871018)
    println(j(1)(1)(0.25, 0.37, 0.264)  == -0.00859584)
    println(j(1)(2)(0.25, 0.37, 0.001)  == -4.563333333333333E-5)

    println(j(2)(0)(0.33, 0.33, 0.33)  == -0.00990689184673647)
    println(j(2)(1)(0.25, 0.37, 0.264)  == -0.012893759999999999)
    println(j(2)(2)(0.25, 0.37, 0.001)  == 3.9999999999999313E9)
  }

  def testDerivativeN = {
    val f1 = (x1: Double, x2: Double, x3: Double, x4: Double) =>
      x1 * cos(x2) + x2 * cos(x1) - 0.9 + cos(2*x3) - cos(2*x4) - 0.4
    val f2 = (x1: Double, x2: Double, x3: Double, x4: Double) =>
      x1 * sin(x2) + x2 * sin(x1) - 0.1 + 2*(x4 - x3) + sin(2*x4) - sin(2*x3) - 1.2
    val f3 = (x1: Double, x2: Double, x3: Double, x4: Double) =>
      x3 * cos(x4) + x4 * cos(x3) - 0.9 + cos(2*x1) - cos(2*x2) - 0.4
    val f4 = (x1: Double, x2: Double, x3: Double, x4: Double) =>
      x3 * sin(x4) + x4 * sin(x3) - 0.1 + 2*(x2 - x1) + sin(2*x2) - sin(2*x1) - 1.2

    val j = jacobian(List(f1, f2, f3, f4))

    println("\ntestDerivative X")
    println(j(0)(0)(0.33, 0.33, 0.11, 0.33)  == 0.8391081441580804)
    println(j(0)(1)(0.25, 0.37, 0.11, 0.33)  == 0.8785085637194042)
    println(j(1)(0)(0.33, 0.11, 0.33, 0.33)  == 0.2138429586252974)
    println(j(1)(1)(0.33, 0.77, 0.25, 0.77)  == 0.5609535493664797)
  }



  def testUnary = {
    val r = 4.0; val K = 1.11; val x0 = 0.1; val tol = 1e-9
    val x = 3.329999999902176
    val err = errorBound((x: Double) => r / (1 + (x/K)) - 1, x, tol)
    val err2 = assertBound((x: Double) => r / (1 + (x/K)) - 1, x, tol)
    val err3 = certify((x: Double) => r / (1 + (x/K)) - 1, x, tol)


    println("\n Unary test")
    println(err.toString == "[-9.7827e-11,-9.7823e-11]")
    println(err2.toString == "[-9.7827e-11,-9.7823e-11]")
    println(err3.toString == "[3.329999999804349,3.330000000000003]  (1.110223024690388E-16)")
  }

  def testBinary = {
    import math._
    val tol = 1e-9
    val f1 = (x1: Double, x2: Double) => cos(2*x1) - cos(2*x2) - 0.4
    val f2 = (x1: Double, x2: Double) => 2*(x2 - x1) + sin(2*x2) - sin(2*x1) - 1.2
    //val x0 = Array(0.1, 0.5)
    val roots = Array(0.1565200696473004, 0.4933763741817787)
    val errors = errorBound(f1, f2, roots(0), roots(1), tol)
    val errors2 = assertBound(f1, f2, roots(0), roots(1), tol)
    //val errors3 = certify(f1, f2, roots(0), roots(1), tol)
    println("\n Binary test")
    println(errors.deep.mkString(", ") == "[3.5835e-11,3.5836e-11], [4.1466e-11,4.1467e-11]")
    println(errors2.deep.mkString(", ") == "[3.5835e-11,3.5836e-11], [4.1466e-11,4.1467e-11]")
    //println(errors3.deep.mkString(", "))
  }


  def testTernary = {
    import math._
    import scala.{Double => D}
    val tol = 1e-9
    val f1 = (v: D, w: D, r: D) =>
      3 + 2/(r*r) - 0.125*(3-2*v)*(w*w*r*r)/(1-v) - 4.5
    val f2 = (v: D, w: D, r: D) =>
      6*v - 0.5 * v * (w*w*r*r) / (1-v) - 2.5
    val f3 = (v: D, w: D, r: D) =>
      3 - 2/(r*r) - 0.125 * (1+2*v) * (w*w*r*r) / (1-v) - 0.5
    //val x0 = Array(0.75, 0.5, 0.5)
    val roots = Array(0.5, 1.0000000000018743, 0.9999999999970013)
    val errors = errorBound(f1, f2, f3, roots(0), roots(1), roots(2), tol)
    val errors2 = assertBound(f1, f2, f3, roots(0), roots(1), roots(2), tol)
    //val errors3 = certify(f1, f2, f3, roots(0), roots(1), roots(2), tol)
    println("\n Ternary test")
    println(errors.deep.mkString(", ") == "[-4.2403e-16,5.0730e-16], [-1.8757e-12,-1.8730e-12], [2.9983e-12,2.9991e-12]")
    println(errors2.deep.mkString(", ") == "[-4.2403e-16,5.0730e-16], [-1.8757e-12,-1.8730e-12], [2.9983e-12,2.9991e-12]")
  }

  def testAnyN = {
    import math._
    import scala.{Double => D}
    val tol = 1e-9
    val f1 = (v: D, w: D, r: D) =>
      3 + 2/(r*r) - 0.125*(3-2*v)*(w*w*r*r)/(1-v) - 4.5
    val f2 = (v: D, w: D, r: D) =>
      6*v - 0.5 * v * (w*w*r*r) / (1-v) - 2.5
    val f3 = (v: D, w: D, r: D) =>
      3 - 2/(r*r) - 0.125 * (1+2*v) * (w*w*r*r) / (1-v) - 0.5
    //val x0 = Array(0.75, 0.5, 0.5)
    val roots = Array(0.5, 1.0000000000018743, 0.9999999999970013)
    val errors = errorBoundX(List(f1, f2, f3), roots, tol)
    val errors2 = assertBoundX(List(f1, f2, f3), roots, tol)
    //val errors3 = certifyX(f1, f2, f3, roots(0), roots(1), roots(2), tol)
    println("\n Any N test")
    println(errors.deep.mkString(", ") == "[-4.2403e-16,5.0730e-16], [-1.8757e-12,-1.8730e-12], [2.9983e-12,2.9991e-12]")
    println(errors2.deep.mkString(", ") == "[-4.2403e-16,5.0730e-16], [-1.8757e-12,-1.8730e-12], [2.9983e-12,2.9991e-12]")
  }
}
