/*
 * Cerres (c) 2012-2014 EPFL, Lausanne
 */

import ceres.common._
import ceres.affine._
import cerres.examples._
import FunctionLibrary._
import cerres.algorithms._
import ceres.common.{QuadDouble => QD}
import QD._

object CerresMain extends App with MultivariateApplications {

  val res = doublePendulumQD(1.0e-16)
  println(res(0))
  println(res(1))

  //val roots = NewtonsExamples.testP1_d(1.0e-5)
/*
  val qdTolerance = QD.pow(QD(10.0), QD(-30.0))
  val dblTolerance5 = math.pow(10.0, -5)
  val dblTolerance10 = math.pow(10.0, -10)
  val dblTolerance15 = math.pow(10.0, -15)

  val doubleString = "%1.3e"
  Interval.doubleString = doubleString

  println("(x*x*x*x*x + 3 * x*x*x*x - 11 * x*x*x - 27*x*x +10*x + 24)/20.0")
  evaluate(List(-4.5, -2.5, -0.7, 0.7, 3.8), p1_d, p1_d_1, p1_qd, p1_qd_1, p1_aa, p1_aa_1,
    dblTolerance10, dblTolerance5)

  println("\nx*x*x - x + 1")
  evaluate(List(-2.0), p2_d, p2_d_1, p2_qd, p2_qd_1, p2_aa, p2_aa_1,
    dblTolerance10, dblTolerance5)

  println("\nmath.cos(x) - 2 * x")
  evaluate(List(0.1), p3_d, p3_d_1, p3_qd, p3_qd_1, p3_aa, p3_aa_1,
    dblTolerance10, dblTolerance5)

  println("\n1/x - 0.3")
  evaluate(List(2.6), p4_d, p4_d_1, p4_qd, p4_qd_1, p4_aa, p4_aa_1,
    dblTolerance10, dblTolerance5)

  //println("\nmath.exp(x) - 2*x - 0.1")
  //evaluate(List(), p5_d, p5_d_1, p5_qd, p5_qd_1, p5_aa, p5_aa_1,
  //  dblTolerance10, dblTolerance5)
  //println("\nx * math.exp(-(x*x))")
  //evaluate(List(), p6_d, p6_d_1, p6_qd, p6_qd_1, p6_aa, p6_aa_1,
  //  dblTolerance10, dblTolerance5)

  println("\n4*x*x*x - 15*x*x + 17*x - 6")
  evaluate(List(0.5, 1.2, 2.5), p7_d, p7_d_1, p7_qd, p7_qd_1, p7_aa, p7_aa_1,
    dblTolerance10, dblTolerance5)

  println("\n3 * math.exp(x) - 4 * math.cos(x)")
  evaluate(List(-1.5, 0.5), p8_d, p8_d_1, p8_qd, p8_qd_1, p8_aa, p8_aa_1,
    dblTolerance10, dblTolerance5)

  println("\n25*x*x -10*x + 1")
  evaluate(List(0.0), p9_d, p9_d_1, p9_qd, p9_qd_1, p9_aa, p9_aa_1,
    dblTolerance10, dblTolerance5)

  println("\nmath.atan(x)")
  evaluate(List(1.0), p10_d, p10_d_1, p10_qd, p10_qd_1, p10_aa, p10_aa_1,
    dblTolerance10, dblTolerance5)

  println("\nmath.cos(x) + 2*math.sin(x) + x*x")
  evaluate(List(-1.6, -0.1), p11_d, p11_d_1, p11_qd, p11_qd_1, p11_aa, p11_aa_1,
    dblTolerance10, dblTolerance5)

  println("\n4 * math.log(x) - x")
  evaluate(List(1.9, 8.0), p12_d, p12_d_1, p12_qd, p12_qd_1, p12_aa, p12_aa_1,
    dblTolerance10, dblTolerance5)
  */
  /**
    Computes an upper bound on the true error between the computed solution xn
    and a root of the function f.
    @param f function whose roots we want
    @param xn computed solution, i.e. one possible root
    @param f' derivative of f
    @param xbar interval to consider for evaluating the slope
  */
 /* def getTrueError(f: AffineForm => AffineForm, xn: Double, fp: AffineForm => AffineForm,
    xbar: Interval): AffineForm = {
    import AffineForm._

    val xn_aa = new AForm(xn)

    val yDelta: AffineForm = abs(f(xn_aa))
    //println("yDelta: " + yDelta.interval)

    // xbar has to be an interval and f' is evaluated in affine floats
    val xbar_aa = AffineForm(xbar)
    //println("xbar: " + xbar_aa)

    val slope = fp(xbar_aa)

    val minSlope: Double = slope.interval match {
      case EmptyInterval =>
        println("we're screwed")
        return 0.0
      case NormalInterval(xlo, xhi) =>
        if (xlo <= 0.0 && xhi >= 0.0) return 0.0
        else math.min(math.abs(xlo), math.abs(xhi))
    }
    //println("minSlope: " + minSlope)

    return yDelta/minSlope
  }

  def evaluate(startPoints: List[Double], fnc: Double => Double, drv: Double => Double,
    fncQD: QD => QD, drvQD: QD => QD, fncAA:AffineForm => AffineForm, drvAA: AffineForm => AffineForm,
    dblTolerance: Double, plusMinus: Double) = {
    require(plusMinus >= 0.0, "plusMinus has to be positive")
    //val startPoints = List(-4.5, -2.5, -0.7, 0.7, 3.8)

    for (x0 <- startPoints) {
      val trueRoot = NewtonsMethod.computeRootQD(fncQD, drvQD, x0, qdTolerance)
      print("x0: " + x0 + ", true root: " + trueRoot)

      val doubleRoot = NewtonsMethod.computeRoot(fnc, drv, x0, dblTolerance)
      print(", double root: " + doubleRoot)

      val trueError = abs(trueRoot - doubleRoot)
      print(", true error: " + doubleString.format(trueError.toDouble))
      val xBar = Interval(doubleRoot - plusMinus, doubleRoot + plusMinus)
      val error = getTrueError(fncAA, doubleRoot, drvAA, xBar)
      println(", error: " + error.interval)
    }

  }
*/
}
