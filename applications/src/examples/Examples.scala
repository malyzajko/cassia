/*
 * Cerres (c) 2012-2014 EPFL, Lausanne
 */

package cerres.examples

import ceres.common._
import ceres.affine._
import ceres.common.{QuadDouble => QD}

object Examples extends App with UnaryApplications with MultivariateApplications {
  Interval.doubleString = "%1.16e"
  //-----------------------------------------
  //    Unary examples
  //-----------------------------------------
  unary("rods", 1e-10, rodSystem)
  unary("verhulst", 1e-9, verhulstModel)
  unary("predator-prey", 1e-9, predatorPreyModel)
  unary("carbon gas state equation", 1e-12, carbonGasState)
  unary("Butler-Volmer equation", 1e-10, butlerVolmer)
  unary("x^2 - sin", 1e-10, sinXSquared)
  unary("exponents", 1e-8, exponents)
  unary("polynomial1", 1e-7, polynomial1)
  unary("polynomial2", 1e-5, polynomial2)

  def unary(name: String, tol: Double, func: Double => (Double, Interval)) = {
    print(name)
    val r = func(tol)
    print(", root: " + r._1)
    println(", error: " + r._2)
  }

  //-----------------------------------------
  //    Multivariate examples
  //-----------------------------------------
  multivariate("stress distribution", 1e-10, stressDistribution)
  multivariate("sin-cos system", 1e-7, sinCosSystem) // choosing 1e-8 or smaller results in one more iteration
  multivariate("double pendulum", 1e-13, doublePendulum)
  multivariate("cirle-parabola intersection", 1e-13 , circleParabola)
  multivariate("no-name quadratic", 1e-6, noNameQuadratic)
  multivariate("turbine", 1e-8, turbine)
  multivariate("no-name quadratic 2", 1e-10, noNameQuadratic2)

  def multivariate(name: String, tol: Double,
    func: Double => (Array[Double], Array[Interval])) = {
    println(name)
    val r = func(tol)
    println("roots: " + r._1.toList)
    println("errors: " + r._2.toList)
  }


  //-----------------------------------------
  //         Roots in QuadDouble
  //-----------------------------------------
  val tolQD = QD(1e-30)
  println(rodSystemQD(tolQD))
  println(verhulstModelQD(tolQD))
  println(predatorPreyModelQD(tolQD))
  println(carbonGasStateQD(tolQD))
  println(butlerVolmerQD(tolQD))
  println(sinXSquaredQD(tolQD))
  println(exponentsQD(tolQD))
  println(polynomial1QD(tolQD))
  println(polynomial2QD(tolQD))
  println(stressDistributionQD(tolQD).toList)
  println(sinCosSystemQD(tolQD).toList)
  println(doublePendulumQD(tolQD).toList)
  println(noNameQuadraticQD(tolQD).toList)
  println(turbineQD(tolQD).toList)
  println(noNameQuadratic2QD(tolQD).toList)



}
