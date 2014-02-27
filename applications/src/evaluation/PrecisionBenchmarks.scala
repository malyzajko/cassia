/*
 * Cerres (c) 2012-2014 EPFL, Lausanne
 */

package cerres.evaluation

import ceres.common.{QuadDouble => QD}

import cerres.algorithms.NewtonsMethod._
import cerres.algorithms.FixedPointIteration._

object PrecisionBenchmarks extends App with ErrorComparison {
  /* ============================
      ERROR ESTIMATES
   ==============================*/
  //println("verhulst " + verhulstModelEST(1e-9))
  //println("predator-prey " + predatorPreyModelEST(1e-9))
  //println("rods " + rodSystemEST(1e-9))




  // Derivatives comparison
  //val rods = rodSystem(10.0, 13.0, 8.0, 10.0, 0.4, 0.0)
  //Interval.doubleString = "%1.16e"
  //println("Rods. Manually: " + rods._1.interval + ", auto: " + rods._2.interval)


  /*val (roots, errors) = Examples.sinCosSystem(1e-7)
  val trueErrors = getTrueErrors(sinCosSystemQD, roots)
  println("true errors: " + trueErrors.toList)
  */

  /* ============================
  //      TRUE ERRORS COMPARISON
  ===============================*/
  import cerres.examples.Examples._
  val tolQD = QD(1e-30)
  
  getTrueErrors("rods", rodSystemQD, rodSystem(1e-10)._1)
  getTrueErrors("verhulst", verhulstModelQD, verhulstModel(1e-9)._1)
  getTrueErrors("predator", predatorPreyModelQD, predatorPreyModel(1e-9)._1)
  getTrueErrors("carbon gas", carbonGasStateQD, carbonGasState(1e-12)._1)
  getTrueErrors("butler volmer", butlerVolmerQD, butlerVolmer(1e-10)._1)
  getTrueErrors("sin x", sinXSquaredQD, sinXSquared(1e-10)._1)
  getTrueErrors("exponents", exponentsQD, exponents(1e-8)._1)
  getTrueErrors("degree 3 poly", polynomial1QD, polynomial1(1e-7)._1)
  getTrueErrors("degree 6 poly", polynomial2QD, polynomial2(1e-5)._1)



  getTrueErrors("sin cos", sinCosSystemQD, sinCosSystem(1e-7)._1)
  getTrueErrors("stress distribution", stressDistributionQD, stressDistribution(1e-10)._1)
  getTrueErrors("double pendulum", doublePendulumQD, doublePendulum(1e-13)._1)
  getTrueErrors("circle parabola", circleParabolaQD, circleParabola(1e-13)._1)
  getTrueErrors("no name quadratic", noNameQuadraticQD, noNameQuadratic(1e-6)._1)
  getTrueErrors("turbine", turbineQD, turbine(1e-12)._1)
  getTrueErrors("no name quadratic 2", noNameQuadratic2QD, noNameQuadratic2(1e-10)._1)

  def verhulstModelQD: QD = {
    import QD._
    val r = 4.0; val K = 1.11; val x0 = QD(0.1)
    iterateQD((x: QD) => (r*x) / (1 + (x/K)), x0, tolQD)
  }

  def predatorPreyModelQD: QD = {
    import QD._
    val r = 4.0; val K = 1.11; val x0 = QD(0.7)
    iterateQD((x: QD) => (r*x*x) / (1 + pow(x/K, 2)), x0, tolQD)
  }

  def rodSystemQD: QD = {
    import QD._
    val a1 = 10.0; val a2 = 13.0; val a3 = 8.0; val a4 = 10.0
    val b = 0.4; val x0 = QD(-0.1)
    val f = (x: QD) =>
      a1/a2 * cos(b) - a1/a4 * cos(x) - cos(b - x) + (a1*a1 + a2*a2 - a3*a3 + a4*a4)/(2*a2*a4)
    val fPrime = (x: QD) => (a1/a4) * sin(x) - sin(b - x)
    computeRootQD(f, fPrime, x0, tolQD)
  }

  def carbonGasStateQD: QD = {
    import QD._
    val T = 300; val a = 0.401; val b = 42.7e-6; val N = 1000
    val p = 3.5e7; val k = 1.3806503e-23; val x0 = QD(0.1)
    val f = (V: QD) => (p + a * (N / V) * (N / V)) * (V - N * b) - k * N * T
    val fPrime = (V: QD) => -2*a * (N*N)/(V*V*V) * (V - N*b) + p + a * (N/V) * (N/V)
    computeRootQD(f, fPrime, x0, tolQD)
  }

  def butlerVolmerQD: QD = {
    import QD._
    val alpha = 0.2; val beta = 2.0; val x0 = QD(2.0)
    val f = (x: QD) => exp(alpha*x) - exp(-(1 - alpha)*x) - beta
    val fPrime = (x: QD) => alpha * exp(alpha*x) + (1 - alpha) * exp(-(1 - alpha)*x)
    computeRootQD(f, fPrime, x0, tolQD)
  }

  def sinXSquaredQD: QD = {
    import QD._
    val x0 = QD(1.6)
    val f = (x: QD) => (x/2)*(x/2) - sin(x)
    val fPrime = (x: QD) => x/2 - cos(x)
    computeRootQD(f, fPrime, x0, tolQD)
  }

  def exponentsQD:QD = {
    import QD._
    val x0 = QD(1.7)
    val f = (x: QD) => exp(x)*(x - 1) - exp(-x)*(x + 1)
    val fPrime = (x: QD) => exp(x)*x - exp(-x)*x
    computeRootQD(f, fPrime, x0, tolQD)
  }

  def polynomial1QD: QD = {
    import QD._
    val x0 = QD(1.2)
    val f = (x: QD) => x*x*x/3.0 - 2 * x*x + 4.5
    val fPrime = (x: QD) => x*x - 4*x
    computeRootQD(f, fPrime, x0, tolQD)
  }

  def polynomial2QD: QD = {
    import QD._
    val x0 = QD(6.5)
    val f = (x: QD) => pow(x,6) + 4.2*pow(x,5) -72.3*pow(x,4) -214.4*pow(x, 3) + 1127.1*pow(x, 2) + 1602.9*x - 5040.5
    val fPrime = (x: QD) => 6*pow(x,5) + 21.0*pow(x,4) -289.2*pow(x,3) -643.2*pow(x, 2) + 2254.2*x + 1602.9
    computeRootQD(f, fPrime, x0, tolQD)
  }

  def sinCosSystemQD: Array[QD]= {
    import QD._
    val f1 = (x1: QD, x2: QD) => x1 * cos(x2) + x2 * cos(x1) - 0.9
    val f2 = (x1: QD, x2: QD) => x1 * sin(x2) + x2 * sin(x1) - 0.1
    val f = Array(f1, f2)
    val x0 = Array(QD(0.0), QD(1.0))
    val f11 = (x1: QD, x2: QD) => cos(x2) - x2 * sin(x1)
    val f12 = (x1: QD, x2: QD) => - x1 * sin(x2) + cos(x1)
    val f21 = (x1: QD, x2: QD) => sin(x2) + x2 * cos(x1)
    val f22 = (x1: QD, x2: QD) => x1 * cos(x2) + sin(x1)
    val J = Array(Array(f11, f12), Array(f21, f22))
    computeRootQD(f, J, x0, tolQD)
  }

  def stressDistributionQD: Array[QD] = {
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
    computeRootQD(f, J, x0, tolQD)
  }

  def doublePendulumQD: Array[QD] = {
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
    computeRootQD(f, J, x0, tolQD)
  }

  def circleParabolaQD: Array[QD] = {
    import QD._
    val f1 = (x: QD, y: QD) => x*x - 2*x - y + 1.0
    val f2 = (x: QD, y: QD) => x*x + y*y - 1.0
    val f = Array(f1, f2)
    val f11 = (x: QD, y: QD) => 2*x - 2.0
    val f12 = (x: QD, y: QD) => -QD(1.0)
    val f21 = (x: QD, y: QD) => 2*x
    val f22 = (x: QD, y: QD) => 2*y
    val J = Array(Array(f11, f12), Array(f21, f22))
    val x0 = Array(QD(0.9), QD(0.2))
    computeRootQD(f, J, x0, tolQD)
  }

  def noNameQuadraticQD: Array[QD] = {
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
    computeRootQD(f, J, x0, tolQD)
  }

  def turbineQD: Array[QD] = {
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
    computeRoot3QD(f, J, x0, tolQD)
  }

  def noNameQuadratic2QD: Array[QD] = {
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
    computeRoot3QD(f, J, x0, tolQD)
  }


 private def getTrueErrors(name: String, trueRoot: QD, root: Double): QD = {
    val trueError = QD(root)- trueRoot
    println(name + ", true error: " + trueError)
    trueError
  }



  private def getTrueErrors(name: String, trueRoots: Array[QD], roots: Array[Double]): Array[QD] = {
    val trueErrors = new Array[QD](trueRoots.length)
    for (i <- 0 until trueRoots.length)
      trueErrors(i) = QD(roots(i)) - trueRoots(i)
    println(name + ", true errors: " + trueErrors.toList)
    trueErrors
  }



}
