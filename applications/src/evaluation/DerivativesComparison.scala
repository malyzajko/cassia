/*
 * Cerres (c) 2012-2014 EPFL, Lausanne
 */

package cerres.evaluation

import cerres.algorithms._
import cerres.macros._

import ceres.affine._
import ceres.common._

import cerres.examples._

object DerivativesComparison extends App with UnaryApplications with MultivariateApplications {
  import Macros._
  println("Automatic:")
  /*multivariate("stress distribution", 1e-10, stressDistribution)
  multivariate("sin-cos system", 1e-7, sinCosSystem)
  multivariate("double pendulum", 1e-13, doublePendulum)
  multivariate("cirle-parabola intersection", 1e-13 , circleParabola)
  multivariate("no-name quadratic", 1e-6, noNameQuadratic)
  multivariate("turbine", 1e-12, turbine)
  multivariate("no-name quadratic 2", 1e-10, noNameQuadratic2)

  def multivariate(name: String, tol: Double,
    func: Double => (Array[Double], Array[Interval])) = {
    print(name)
    val roots = func(tol)
    print(", errors: ")
    for (e <- roots._2)
      print(math.max(math.abs(e.xlo), math.abs(e.xhi)) + ":  ")
    println()
  }
  */


  //unary("polynomial 1", 1e-12, polynomial1)
  //unary("polynomial 2", 1e-5, polynomial2)
  //unary("polynomial 2_2", 1e-5, polynomial2_2)

  unary("rods", 1e-10, rodSystem)
  unary("verhulst", 1e-9, verhulstModel)
  unary("predator-prey", 1e-9, predatorPreyModel)
  unary("carbon gas state equation", 1e-12, carbonGasState)
  unary("Butler-Volmer equation", 1e-10, butlerVolmer)
  unary("x^2 - sin", 1e-10, sinXSquared)
  unary("exponents", 1e-8, exponents)


  //println("\n\n Manual")
  //println("polynomial 1: " + polynomial1Manual(1e-12))
  //println("polynomial 2: " + polynomial2Manual(1e-5))
  //println("polynomial 2_2: " + polynomial2_2Manual(1e-5))



  println("rods: " + rodSystemManual(1e-10))
  println("verhulst: " + verhulstModelManual(1e-9))
  println("predator-prey: " + predatorPreyModelManual(1e-9))
  println("carbon gas state equation: " + carbonGasStateManual(1e-12))
  println("Butler-Volmer equation: " + butlerVolmerManual(1e-10))
  println("x^2 - sin: " + sinXSquaredManual(1e-10))
  println("exponents: " + exponentsManual(1e-8))

  def unary(name: String, tol: Double, func: Double => (Double, Interval)) = {
    print(name)
    val r = func(tol)
    println(", error: " + math.max(math.abs(r._2.xlo), math.abs(r._2.xhi)))
  }


  /* ===========================================
     VERSION WITH MANUAL DERIVATIVES
   ============================================= */
 import math._
 def verhulstModelManual(tol: Double): Double = {
    val r = 4.0; val K = 1.11; val x0 = 0.1
    val x = FixedPointIteration.iterate((x: Double) => (r*x) / (1 + (x/K)), x0, tol)
    val f = (x: Double) => r / (1 + (x/K)) - 1
    val fPrime = (x: Double) => -r/ (K * (1 + (x/K))*(1 + (x/K)))
    maxAbs(errorManually(f, fPrime, x, tol))
  }

  def predatorPreyModelManual(tol: Double): Double = {
    val r = 4.0; val K = 1.11; val x0 = 0.7
    val x = FixedPointIteration.iterate((x: Double) => (r*x*x) / (1 + pow(x/K, 2)), x0, tol)
    val f = (x: Double) => (r*x) / (1 + (x/K)*(x/K)) - 1
    val fPrime = (x: Double) => r/(1 + (x/K)*(x/K)) - (2*r*x)/ (K*K*(1 + (x/K)*(x/K))*(1 + (x/K)*(x/K)))
    maxAbs(errorManually(f, fPrime, x, tol))
  }

  def rodSystemManual(tol: Double): Double = {
    val a1 = 10.0; val a2 = 13.0; val a3 = 8.0; val a4 = 10.0
    val b = 0.4; val x0 = -0.1
    val f = (x: Double) =>
      a1/a2 * cos(b) - a1/a4 * cos(x) - cos(b - x) + (a1*a1 + a2*a2 - a3*a3 + a4*a4)/(2*a2*a4)
    val fPrime = (x: Double) => (a1/a4) * sin(x) - sin(b - x)
    val alpha = NewtonsMethod.computeRoot(f, fPrime, x0, tol)
    maxAbs(errorManually(f, fPrime, alpha, tol))
  }

  def carbonGasStateManual(tol: Double): Double = {
    val T = 300; val a = 0.401; val b = 42.7e-6; val N = 1000
    val p = 3.5e7; val k = 1.3806503e-23; val x0 = 0.1
    val f = (V: Double) => (p + a * (N / V) * (N / V)) * (V - N * b) - k * N * T
    val fPrime = (V: Double) => -2*a * (N*N)/(V*V*V) * (V - N*b) + p + a * (N/V) * (N/V)
    val V = NewtonsMethod.computeRoot(f, fPrime, x0, tol)
    maxAbs(errorManually(f, fPrime, V, tol))
  }

  def butlerVolmerManual(tol: Double): Double = {
    val alpha = 0.2; val beta = 2.0; val x0 = 2.0
    val f = (x: Double) => exp(alpha*x) - exp(-(1 - alpha)*x) - beta
    val fPrime = (x: Double) => alpha * exp(alpha*x) + (1 - alpha) * exp(-(1 - alpha)*x)
    val root = NewtonsMethod.computeRoot(f, fPrime, x0, tol)
    maxAbs(errorManually(f, fPrime, root, tol))
  }

  def sinXSquaredManual(tol: Double): Double = {
    val x0 = 1.6
    val f = (x: Double) => (x/2)*(x/2) - sin(x)
    val fPrime = (x: Double) => x/2 - cos(x)
    val root = NewtonsMethod.computeRoot(f, fPrime, x0, tol)
    maxAbs(errorManually(f, fPrime, root, tol))
  }

  def exponentsManual(tol: Double): Double = {
    val x0 = 1.7
    val f = (x: Double) => exp(x)*(x - 1) - exp(-x)*(x + 1)
    val fPrime = (x: Double) => exp(x)*x - exp(-x)*x
    val root = NewtonsMethod.computeRoot(f, fPrime, x0, tol)
    maxAbs(errorManually(f, fPrime, root, tol))
  }

  def polynomial1Manual(tol: Double): Double = {
    val x0 = 1.2
    val f = (x: Double) => x*x*x/3.0 - 2 * x*x + 4.5
    val fPrime = (x: Double) => x * (x - 4)
    val root = NewtonsMethod.computeRoot(f, fPrime, x0, tol)
    maxAbs(errorManually(f, fPrime, root, tol))
  }

  def polynomial2Manual(tol: Double): Double = {
    val x0 = 6.5
    val f = (x: Double) => pow(x,6) + 4*pow(x,5) -72*pow(x,4) -214*pow(x, 3) + 1127*pow(x, 2) + 1602*x - 5040
    val fPrime = (x: Double) => x*( x*( x*( x* (6*x + 20) - 288) -642) + 2254) +1602
    val root = NewtonsMethod.computeRoot(f, fPrime, x0, tol)
    maxAbs(errorManually(f, fPrime, root, tol))
  }

  def polynomial2_2Manual(tol: Double): Double = {
    val x0 = 6.5
    val f = (x: Double) => (x-3)*(x+3)*(x+5)*(x+8)*(x-2)*(x-7)
    val fPrime = (x: Double) => x*( x*( x*( x* (6*x + 20) - 288) -642) + 2254) +1602
    val root = NewtonsMethod.computeRoot(f, fPrime, x0, tol)
    maxAbs(errorManually(f, fPrime, root, tol))
  }

 
}
