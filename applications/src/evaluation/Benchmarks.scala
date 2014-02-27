/*
 * Cerres (c) 2012-2014 EPFL, Lausanne
 */

package cerres.benchmarks

import com.google.caliper.Param

import ceres.common.{QuadDouble => QD}
import ceres.common._

object NoErrorsBenchmark extends SimpleScalaBenchmark {
  import ExamplesNoError._
  
  //@Param(Array("10", "100", "1000", "10000"))
  //val repetitions: Int = 0

  def turbine(reps: Int) = repeat(reps) {
    turbine_noError(1e-12)(0)    
  }


/*
  //override val name = "NoErrorsBenchmark"
  //var results = 0.0
  //innerLoops = 1000
  //repetitions = 20

    def run() = {
    //results += verhulstModel_noError(1e-9)
    //results += predatorPreyModel_noError(1e-9)
    //results += rodSystem_noError(1e-10)
    //results += carbonGasState_noError(1e-12)
    //results += butlerVolmer_noError(1e-10)
    //results += sinXSquared_noError(1e-10)
    //results += exponents_noError(1e-8)
    //results += polynomial1_noError(1e-7)
    //results += polynomial2_noError(1e-8)

    //results += stressDistribution_noError(1e-10)(0)
    //results += sinCosSystem_noError(1e-7)(0)
    //results += doublePendulum_noError(1e-13)(0)
    //results += circleParabola_noError(1e-13)(0)
    //results += noNameQuadratic_noError(1e-6)(0)

    results += turbine_noError(1e-12)(0)
    results += noNameQuadratic2_noError(1e-10)(0)
  }

  def done() = println("Final accumulator: " + results)
  
  runBenchmark
  val avrgPerRepetition =
    runtimes.toList.foldLeft(0l){ (acc, x) => acc + x } / runtimes.length.toDouble
  println("average per "+innerLoops+" runs: " + avrgPerRepetition)
  */
}

/*object WithErrorsBenchmark extends Benchmark {
  import cerres.examples.Examples._
  override val name = "WithErrorsBenchmark"
  var results = 0.0
  innerLoops = 1000
  repetitions = 20

  def run() = {
    //results += verhulstModel(1e-9)._1
    //results += predatorPreyModel(1e-9)._1
    //results += rodSystem(1e-10)._1
    //results += carbonGasState(1e-12)._1
    //results += butlerVolmer(1e-10)._1
    //results += sinXSquared(1e-10)._1
    //results += exponents(1e-8)._1
    //results += polynomial1(1e-7)._1
    //results += polynomial2(1e-8)._1

    //results += stressDistribution(1e-10)._1(0)
    //results += sinCosSystem(1e-7)._1(0)
    //results += doublePendulum(1e-13)._1(0)
    //results += circleParabola(1e-13)._1(0)
    //results += noNameQuadratic(1e-6)._1(0)

    results += turbine(1e-12)._1(0)
    results += noNameQuadratic2(1e-10)._1(0)

  }

  def done() = println("Final accumulator: " + results)

  runBenchmark
  val avrgPerRepetition =
    runtimes.toList.foldLeft(0l){ (acc, x) => acc + x } / runtimes.length.toDouble
  println("average per "+innerLoops+" runs: " + avrgPerRepetition)
}

object QDBenchmark extends Benchmark {
  import cerres.examples.Examples._
  override val name = "QDBenchmark"
  var results = QD(0.0)
  innerLoops = 1000
  repetitions = 10

  def run() = {
    results += verhulstModelQD(1e-9)
    results += predatorPreyModelQD(1e-9)
    results += rodSystemQD(1e-10)
    results += carbonGasStateQD(1e-12)
    results += butlerVolmerQD(1e-10)
    results += sinXSquaredQD(1e-10)
    results += exponentsQD(1e-8)
    results += polynomial1QD(1e-7)
    results += polynomial2QD(1e-5)

    /*results += stressDistributionQD(1e-10)(0)
    results += sinCosSystemQD(1e-7)(0)
    results += doublePendulumQD(1e-13)(0)
    results += circleParabolaQD(1e-13)(0)
    results += noNameQuadraticQD(1e-6)(0)
    */
    //results += turbineQD(1e-12)(0)
    //results += noNameQuadratic2QD(1e-10)(0)
  }

  def done() = println("Final accumulator: " + results)

  runBenchmark
  val avrgPerRepetition =
    runtimes.toList.foldLeft(0l){ (acc, x) => acc + x } / runtimes.length.toDouble
  println("average per "+innerLoops+" runs: " + avrgPerRepetition)
}*/

object ExamplesNoError {
  import math._
  import cerres.algorithms.NewtonsMethod._
  import cerres.algorithms._

  def verhulstModel_noError(tol: Double): Double = {
    val r = 4.0; val K = 1.11; val x0 = 0.1
    FixedPointIteration.iterate((x: Double) => (r*x) / (1 + (x/K)), x0, tol)
  }

  def predatorPreyModel_noError(tol: Double): Double = {
    val r = 4.0; val K = 1.11; val x0 = 0.7
    FixedPointIteration.iterate((x: Double) => (r*x*x) / (1 + pow(x/K, 2)), x0, tol)
  }

  def rodSystem_noError(tol: Double): Double = {
    val a1 = 10.0; val a2 = 13.0; val a3 = 8.0; val a4 = 10.0
    val b = 0.4; val x0 = -0.1
    val f = (x: Double) =>
      a1/a2 * cos(b) - a1/a4 * cos(x) - cos(b - x) + (a1*a1 + a2*a2 - a3*a3 + a4*a4)/(2*a2*a4)
    val fPrime = (x: Double) => (a1/a4) * sin(x) - sin(b - x)
    NewtonsMethod.computeRoot(f, fPrime, x0, tol)
  }

  def carbonGasState_noError(tol: Double): Double = {
    val T = 300; val a = 0.401; val b = 42.7e-6; val N = 1000
    val p = 3.5e7; val k = 1.3806503e-23; val x0 = 0.1
    val f = (V: Double) => (p + a * (N / V) * (N / V)) * (V - N * b) - k * N * T
    val fPrime = (V: Double) => -2*a * (N*N)/(V*V*V) * (V - N*b) + p + a * (N/V) * (N/V)
    NewtonsMethod.computeRoot(f, fPrime, x0, tol)
  }

  def butlerVolmer_noError(tol: Double): Double = {
    val alpha = 0.2; val beta = 2.0; val x0 = 2.0
    val f = (x: Double) => exp(alpha*x) - exp(-(1 - alpha)*x) - beta
    val fPrime = (x: Double) => alpha * exp(alpha*x) + (1 - alpha) * exp(-(1 - alpha)*x)
    NewtonsMethod.computeRoot(f, fPrime, x0, tol)
  }

  def sinXSquared_noError(tol: Double): Double = {
    val x0 = 1.6
    val f = (x: Double) => (x/2)*(x/2) - sin(x)
    val fPrime = (x: Double) => x/2 - cos(x)
    NewtonsMethod.computeRoot(f, fPrime, x0, tol)
  }

  def exponents_noError(tol: Double): Double = {
    val x0 = 1.7
    val f = (x: Double) => exp(x)*(x - 1) - exp(-x)*(x + 1)
    val fPrime = (x: Double) => exp(x)*x - exp(-x)*x
    NewtonsMethod.computeRoot(f, fPrime, x0, tol)
  }

  def polynomial1_noError(tol: Double): Double = {
    val x0 = 1.2
    val f = (x: Double) => x*x*x/3.0 - 2 * x*x + 4.5
    val fPrime = (x: Double) => x*x - 4*x
    NewtonsMethod.computeRoot(f, fPrime, x0, tol)
  }

  def polynomial2_noError(tol: Double): Double = {
    val x0 = 6.5
    val f = (x: Double) => pow(x,6) + 4.2*pow(x,5) -72.3*pow(x,4) -214.4*pow(x, 3) + 1127.1*pow(x, 2) + 1602.9*x - 5040.5
    val fPrime = (x: Double) => 6*pow(x,5) + 21.0*pow(x,4) -289.2*pow(x,3) -643.2*pow(x, 2) + 2254.2*x + 1602.9
    NewtonsMethod.computeRoot(f, fPrime, x0, tol)
  }

  /* ====================================
      NO ERROR COMPUTATION VERSIONs
   ====================================== */
  def stressDistribution_noError(tol: Double): Array[Double] = {
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
    computeRoot(f, J, x0, tol)
  }

  // Numerical Methods with Worked Examples, C. Woodford, C. Phillips
  // @return (computed roots, certified errors)
  def sinCosSystem_noError(tol: Double): Array[Double] = {
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
    computeRoot(f, J, x0, tol)
  }

  // Numerical Methods in Scientific Computing, Dahlquist, BJoerk
  def doublePendulum_noError(tol: Double): Array[Double] = {
    import math._
    val k = 0.3
    val f1 = (x1: Double, x2: Double) => tan(x1) - k * (2*sin(x1) + sin(x2))
    val f2 = (x1: Double, x2: Double) => tan(x2) - 2*k * (sin(x1) + sin(x2))
    val f = Array(f1, f2)
    val f11 = (x1: Double, x2: Double) => 1.0 + tan(x1)*tan(x1) -2*k*cos(x1)
    val f12 = (x1: Double, x2: Double) => -k*cos(x2)
    val f21 = (x1: Double, x2: Double) => -2*k*cos(x1)
    val f22 = (x1: Double, x2: Double) => 1.0 + tan(x2)*tan(x2) -2*k*cos(x2)
    val J = Array(Array(f11, f12), Array(f21, f22))
    val x0 = Array(0.18, 0.25)
    computeRoot(f, J, x0, tol)
  }


  def circleParabola_noError(tol: Double): Array[Double] = {
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
    computeRoot(f, J, x0, tol)
  }

  def noNameQuadratic_noError(tol: Double): Array[Double] = {
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
    computeRoot(f, J, x0, tol)
  }

  def turbine_noError(tol: Double): Array[Double] = {
    import math._
    val f1 = (v: Double, w: Double, r: Double) =>
      3 + 2/(r*r) - 0.125*(3-2*v)*(w*w*r*r)/(1-v) - 4.5
    val f2 = (v: Double, w: Double, r: Double) =>
      6*v - 0.5 * v * (w*w*r*r) / (1-v) - 2.5
    val f3 = (v: Double, w: Double, r: Double) =>
      3 - 2/(r*r) - 0.125 * (1+2*v) * (w*w*r*r) / (1-v) - 0.5
    val f = Array(f1, f2, f3)

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
    val x0 = Array(0.75, 0.5, 0.5)
    computeRoot3(f, J, x0, tol)
  }

  def noNameQuadratic2_noError(tol: Double): Array[Double] = {
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
    computeRoot3(f, J, x0, tol)
  }




}
