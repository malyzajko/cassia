/*
 * Cerres (c) 2012-2014 EPFL, Lausanne
 */

package cerres.examples

import math._
import cerres.algorithms.NewtonsMethod._
import cerres.macros._
import ceres.common._
import ceres.affine._
import ceres.common.{QuadDouble => QD}

trait MultivariateApplications {
  import Macros._

  def stressDistribution(tol: Double): (Array[Double], Array[Interval]) = {
    import math._
    val f1 = (x1: Double, x2: Double) => cos(2*x1) - cos(2*x2) - 0.4
    val f2 = (x1: Double, x2: Double) => 2*(x2 - x1) + sin(2*x2) - sin(2*x1) - 1.2
    val f = Array(f1, f2)
    val f11 = (x1: Double, x2: Double) => -2 * sin(2*x1)
    val f12 = (x1: Double, x2: Double) => 2 * sin(2*x2)
    val f21 = (x1: Double, x2: Double) => -2.0 - 2*cos(2*x1)
    val f22 = (x1: Double, x2: Double) => 2.0 + 2 * cos(2*x2)
    val J = Array(Array(f11, f12), Array(f21, f22))
    val x0 = Array(0.1, 0.5)
    val roots = computeRoot(f, J, x0, tol)
    val errors = assertBound(f1, f2, roots(0), roots(1), tol)
    (roots, errors)
  }

  def sinCosSystem(tol: Double): (Array[Double], Array[Interval]) = {
    import math._
    val f1 = (x1: Double, x2: Double) => x1 * cos(x2) + x2 * cos(x1) - 0.9
    val f2 = (x1: Double, x2: Double) => x1 * sin(x2) + x2 * sin(x1) - 0.1
    val f = Array(f1, f2)
    val f11 = (x1: Double, x2: Double) => cos(x2) - x2 * sin(x1)
    val f12 = (x1: Double, x2: Double) => - x1 * sin(x2) + cos(x1)
    val f21 = (x1: Double, x2: Double) => sin(x2) + x2 * cos(x1)
    val f22 = (x1: Double, x2: Double) => x1 * cos(x2) + sin(x1)
    val J = Array(Array(f11, f12), Array(f21, f22))
    val x0 = Array(0.0, 1.0)
    val roots = computeRoot(f, J, x0, tol)
    val errors = errorBound(f1, f2, roots(0), roots(1), tol)
    (roots, errors)
  }

  // Numerical Methods in Scientific Computing, Dahlquist, BJoerk
  def doublePendulum(tol: Double): (Array[Double], Array[Interval]) = {
    import math._
    val k = 0.3
    val f1 = (x1: Double, x2: Double) => tan(x1) - k * (2*sin(x1) + sin(x2))
    val f2 = (x1: Double, x2: Double) => tan(x2) - 2*k * (sin(x1) + sin(x2))
    val f = Array(f1, f2)
    /*val f11 = (x1: Double, x2: Double) => 1.0 + tan(x1)*tan(x1) -2*k*cos(x1)
    val f12 = (x1: Double, x2: Double) => -k*cos(x2)
    val f21 = (x1: Double, x2: Double) => -2*k*cos(x1)
    val f22 = (x1: Double, x2: Double) => 1.0 + tan(x2)*tan(x2) -2*k*cos(x2)
    val J = Array(Array(f11, f12), Array(f21, f22))*/
    /*val x1 = 0.179779407413519; val x2 = 0.250801371841681
    println("f1: " + f1(x1, x2))
    println("f2: " + f2(x1, x2))
    println("error: " + errorBound(f1, f2, x1, x2, 1e-13).toList)
    */
    val x0 = Array(0.18, 0.25)
    val roots = computeRoot(f, jacobian(f1, f2), x0, tol)
    val errors = errorBound(f1, f2, roots(0), roots(1), tol)
    (roots, errors)
    //(Array(x1, x2), errors)
  }


  def circleParabola(tol: Double): (Array[Double], Array[Interval]) = {
    import math._
    val f1 = (x: Double, y: Double) => x*x - 2*x - y + 1.0
    val f2 = (x: Double, y: Double) => x*x + y*y - 1.0
    val f = Array(f1, f2)
    val f11 = (x: Double, y: Double) => 2*x - 2.0
    val f12 = (x: Double, y: Double) => -1.0
    val f21 = (x: Double, y: Double) => 2*x
    val f22 = (x: Double, y: Double) => 2*y
    val J = Array(Array(f11, f12), Array(f21, f22))
    val x0 = Array(0.9, 0.2)
    val roots = computeRoot(f, J, x0, tol)
    val errors = assertBound(f1, f2, roots(0), roots(1), tol)
    (roots, errors)
  }

  def noNameQuadratic(tol: Double): (Array[Double], Array[Interval]) = {
    import math._
    val f1 = (x: Double, y: Double) => x*x + y*y - 4*x
    val f2 = (x: Double, y: Double) => y*y + 2*x - 2.0
    val f = Array(f1, f2)
    val f11 = (x: Double, y: Double) => 2*x - 4.0
    val f12 = (x: Double, y: Double) => 2*y
    val f21 = (x: Double, y: Double) => 2.0
    val f22 = (x: Double, y: Double) => 2*y
    val J = Array(Array(f11, f12), Array(f21, f22))
    val x0 = Array(0.35, 1.15)
    val roots = computeRoot(f, J, x0, tol)
    val errors = assertBound(f1, f2, roots(0), roots(1), tol)
    (roots, errors)
  }

  def turbine(tol: Double): (Array[Double], Array[Interval]) = {
    import math._
    import scala.{Double => D}
    val f1 = (v: D, w: D, r: D) =>
      3 + 2/(r*r) - 0.125*(3-2*v)*(w*w*r*r)/(1-v) - 4.5
    val f2 = (v: D, w: D, r: D) =>
      6*v - 0.5 * v * (w*w*r*r) / (1-v) - 2.5
    val f3 = (v: D, w: D, r: D) =>
      3 - 2/(r*r) - 0.125 * (1+2*v) * (w*w*r*r) / (1-v) - 0.5
    val f = Array(f1, f2, f3)
    /*
    val f11 = (v: Double, w: Double, r: Double) => - (w*w*r*r) / (8*pow(1-v, 2))
    val f12 = (v: Double, w: Double, r: Double) => - ((3-2*v)*w*r*r) / (4*(1-v))
    val f13 = (v: Double, w: Double, r: Double) => - 4/(r*r) - ((3-2*v)*w*w*r) / (4*(1-v))
    val f21 = (v: Double, w: Double, r: Double) => 6 + (w*w*r*r) / (2*pow(1-v, 2))
    val f22 = (v: Double, w: Double, r: Double) => - (v*w*r*r) / (1-v)
    val f23 = (v: Double, w: Double, r: Double) => -(v*w*w*r) / (1-v)
    val f31 = (v: Double, w: Double, r: Double) => - (3*w*w*r*r) / (8*pow(1-v, 2))
    val f32 = (v: Double, w: Double, r: Double) => - ((1+2*v)*w*r*r) / (4*(1-v))
    val f33 = (v: Double, w: Double, r: Double) => 4/(r*r) - ((1+2*v)*w*w*r) / (4*(1-v))
    val J = Array(Array(f11, f12, f13), Array(f21, f22, f23), Array(f31, f32, f33))
    */
    val x0 = Array(0.75, 0.5, 0.5)
    //val roots = computeRoot3(f, J, x0, tol)
    val roots = computeRoot3(f, jacobian(f1, f2, f3), x0, tol)
    val errors = assertBound(f1, f2, f3, roots(0), roots(1), roots(2), tol)
    (roots, errors)
  }

  def noNameQuadratic2(tol: Double): (Array[Double], Array[Interval]) = {
    import math._
    val f1 = (x1: Double, x2: Double, x3: Double) => 4*x1*x1 + 2*x2 - 4.0
    val f2 = (x1: Double, x2: Double, x3: Double) => x1 + 2*x2*x2 + 2*x3 - 4.0
    val f3 = (x1: Double, x2: Double, x3: Double) => x2 + 2*x3*x3 - 4.0
    val f = Array(f1, f2, f3)

    val f11 = (x1: Double, x2: Double, x3: Double) => 8*x1
    val f12 = (x1: Double, x2: Double, x3: Double) => 2.0
    val f13 = (x1: Double, x2: Double, x3: Double) => 0.0

    val f21 = (x1: Double, x2: Double, x3: Double) => 1.0
    val f22 = (x1: Double, x2: Double, x3: Double) => 4*x2
    val f23 = (x1: Double, x2: Double, x3: Double) => 2.0

    val f31 = (x1: Double, x2: Double, x3: Double) => 0.0
    val f32 = (x1: Double, x2: Double, x3: Double) => 1.0
    val f33 = (x1: Double, x2: Double, x3: Double) => 4*x3
    val J = Array(Array(f11, f12, f13), Array(f21, f22, f23), Array(f31, f32, f33))
    val x0 = Array(0.3, 0.5, 0.5)
    val roots = computeRoot3(f, J, x0, tol)
    val errors = errorBound(f1, f2, f3, roots(0), roots(1), roots(2), tol)
    (roots, errors)
  }


  /* ====================================
            Quad QD VERSIONs
   ====================================== */
  def stressDistributionQD(tol: QD): Array[QD] = {
    import QD._
    val f1 = (x1: QD, x2: QD) => cos(2*x1) - cos(2*x2) - 0.4
    val f2 = (x1: QD, x2: QD) => 2*(x2 - x1) + sin(2*x2) - sin(2*x1) - 1.2
    val f = Array(f1, f2)
    val f11 = (x1: QD, x2: QD) => -2 * sin(2*x1)
    val f12 = (x1: QD, x2: QD) => 2 * sin(2*x2)
    val f21 = (x1: QD, x2: QD) => -2.0 - 2*cos(2*x1)
    val f22 = (x1: QD, x2: QD) => 2.0 + 2 * cos(2*x2)
    val J = Array(Array(f11, f12), Array(f21, f22))
    val x0 = Array(QD(0.1), QD(0.5))
    computeRootQD(f, J, x0, tol)
  }

  // Numerical Methods with Worked Examples, C. Woodford, C. Phillips
  // @return (computed roots, certified errors)
  def sinCosSystemQD(tol: QD): Array[QD] = {
    import QD._
    val f1 = (x1: QD, x2: QD) => x1 * cos(x2) + x2 * cos(x1) - 0.9
    val f2 = (x1: QD, x2: QD) => x1 * sin(x2) + x2 * sin(x1) - 0.1
    val f = Array(f1, f2)
    val f11 = (x1: QD, x2: QD) => cos(x2) - x2 * sin(x1)
    val f12 = (x1: QD, x2: QD) => - x1 * sin(x2) + cos(x1)
    val f21 = (x1: QD, x2: QD) => sin(x2) + x2 * cos(x1)
    val f22 = (x1: QD, x2: QD) => x1 * cos(x2) + sin(x1)
    val J = Array(Array(f11, f12), Array(f21, f22))
    val x0 = Array(QD(0.0), QD(1.0))
    computeRootQD(f, J, x0, tol)
  }

  // Numerical Methods in Scientific Computing, Dahlquist, BJoerk
  def doublePendulumQD(tol: QD): Array[QD] = {
    import QD._
    val k = 0.3
    val f1 = (x1: QD, x2: QD) => tan(x1) - k * (2*sin(x1) + sin(x2))
    val f2 = (x1: QD, x2: QD) => tan(x2) - 2*k * (sin(x1) + sin(x2))
    val f = Array(f1, f2)
    val f11 = (x1: QD, x2: QD) => 1.0 + tan(x1)*tan(x1) -2*k*cos(x1)
    val f12 = (x1: QD, x2: QD) => -k*cos(x2)
    val f21 = (x1: QD, x2: QD) => -2*k*cos(x1)
    val f22 = (x1: QD, x2: QD) => 1.0 + tan(x2)*tan(x2) -2*k*cos(x2)
    val J = Array(Array(f11, f12), Array(f21, f22))
    val x0 = Array(QD(0.18), QD(0.25))
    computeRootQD(f, J, x0, tol)
  }


  def circleParabolaQD(tol: QD): Array[QD] = {
    import QD._
    val f1 = (x: QD, y: QD) => x*x - 2*x - y + 1.0
    val f2 = (x: QD, y: QD) => x*x + y*y - 1.0
    val f = Array(f1, f2)
    val f11 = (x: QD, y: QD) => 2*x - 2.0
    val f12 = (x: QD, y: QD) => QD(-1.0)
    val f21 = (x: QD, y: QD) => 2*x
    val f22 = (x: QD, y: QD) => 2*y
    val J = Array(Array(f11, f12), Array(f21, f22))
    val x0 = Array(QD(0.9), QD(0.2))
    computeRootQD(f, J, x0, tol)
  }

  def noNameQuadraticQD(tol: QD): Array[QD] = {
    import QD._
    val f1 = (x: QD, y: QD) => x*x + y*y - 4*x
    val f2 = (x: QD, y: QD) => y*y + 2*x - 2.0
    val f = Array(f1, f2)
    val f11 = (x: QD, y: QD) => 2*x - 4.0
    val f12 = (x: QD, y: QD) => 2*y
    val f21 = (x: QD, y: QD) => QD(2.0)
    val f22 = (x: QD, y: QD) => 2*y
    val J = Array(Array(f11, f12), Array(f21, f22))
    val x0 = Array(QD(0.35), QD(1.15))
    computeRootQD(f, J, x0, tol)
  }

  def turbineQD(tol: QD): Array[QD] = {
    import QD._
    val f1 = (v: QD, w: QD, r: QD) =>
      3 + 2/(r*r) - 0.125*(3-2*v)*(w*w*r*r)/(1-v) - 4.5
    val f2 = (v: QD, w: QD, r: QD) =>
      6*v - 0.5 * v * (w*w*r*r) / (1-v) - 2.5
    val f3 = (v: QD, w: QD, r: QD) =>
      3 - 2/(r*r) - 0.125 * (1+2*v) * (w*w*r*r) / (1-v) - 0.5
    val f = Array(f1, f2, f3)

    val f11 = (v: QD, w: QD, r: QD) => - (w*w*r*r) / (8*pow(1-v, 2))
    val f12 = (v: QD, w: QD, r: QD) => - ((3-2*v)*w*r*r) / (4*(1-v))
    val f13 = (v: QD, w: QD, r: QD) => - 4/(r*r) - ((3-2*v)*w*w*r) / (4*(1-v))

    val f21 = (v: QD, w: QD, r: QD) => 6 + (w*w*r*r) / (2*pow(1-v, 2))
    val f22 = (v: QD, w: QD, r: QD) => - (v*w*r*r) / (1-v)
    val f23 = (v: QD, w: QD, r: QD) => -(v*w*w*r) / (1-v)

    val f31 = (v: QD, w: QD, r: QD) => - (3*w*w*r*r) / (8*pow(1-v, 2))
    val f32 = (v: QD, w: QD, r: QD) => - ((1+2*v)*w*r*r) / (4*(1-v))
    val f33 = (v: QD, w: QD, r: QD) => 4/(r*r) - ((1+2*v)*w*w*r) / (4*(1-v))
    val J = Array(Array(f11, f12, f13), Array(f21, f22, f23), Array(f31, f32, f33))
    val x0 = Array(QD(0.75), QD(0.5), QD(0.5))
    computeRoot3QD(f, J, x0, tol)
  }

  def noNameQuadratic2QD(tol: QD): Array[QD] = {
    import QD._
    val f1 = (x1: QD, x2: QD, x3: QD) => 4*x1*x1 + 2*x2 - 4.0
    val f2 = (x1: QD, x2: QD, x3: QD) => x1 + 2*x2*x2 + 2*x3 - 4.0
    val f3 = (x1: QD, x2: QD, x3: QD) => x2 + 2*x3*x3 - 4.0
    val f = Array(f1, f2, f3)

    val f11 = (x1: QD, x2: QD, x3: QD) => 8*x1
    val f12 = (x1: QD, x2: QD, x3: QD) => QD(2.0)
    val f13 = (x1: QD, x2: QD, x3: QD) => QD(0.0)

    val f21 = (x1: QD, x2: QD, x3: QD) => QD(1.0)
    val f22 = (x1: QD, x2: QD, x3: QD) => 4*x2
    val f23 = (x1: QD, x2: QD, x3: QD) => QD(2.0)

    val f31 = (x1: QD, x2: QD, x3: QD) => QD(0.0)
    val f32 = (x1: QD, x2: QD, x3: QD) => QD(1.0)
    val f33 = (x1: QD, x2: QD, x3: QD) => 4*x3
    val J = Array(Array(f11, f12, f13), Array(f21, f22, f23), Array(f31, f32, f33))
    val x0 = Array(QD(0.3), QD(0.5), QD(0.5))
    computeRoot3QD(f, J, x0, tol)
  }



}
