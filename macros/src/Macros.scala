/*
 * Cerres (c) 2012-2014 EPFL, Lausanne
 */

package cerres.macros

import scala.language.experimental.macros
import scala.reflect.macros.Context

import collection.mutable.ListBuffer
import collection.mutable.Stack

import scala.Double.{PositiveInfinity => PlusInf}
import scala.Double.{NegativeInfinity => MinusInf}

import ceres.affine._
import ceres.common._
import ceres.smartfloat.SmartFloat


object Macros extends TrueErrors {
  type D = Double
  type I = Interval
  type A = AffineForm

  //-----------------------------------
  //           DERIVATIVES
  //-----------------------------------
  def derivative(f: D => D): (D => D) = macro derivative_impl
  def jacobian(f1: Function2[D, D, D], f2: Function2[D, D, D]):
    (Array[Array[(D, D) => D]]) = macro derivative2_impl
  def jacobian(f1: Function3[D, D, D, D], f2: Function3[D, D, D, D],
    f3: Function3[D, D, D, D]): (Array[Array[(D, D, D) => D]]) = macro derivative3_impl

  /* Computes the jacobian of n equations in n variables.
   * If one want the partial derivative from just one function, one can leave the
   * rest equal to zero, for instance.*/
  def jacobian(f: List[Function4[D, D, D, D, D]]):
    Array[Array[Function4[D, D, D, D, D]]] =
    macro derivativeMult_impl[Function4[D, D, D, D, D]]

  def jacobian(f: List[Function5[D, D, D, D, D, D]]):
    Array[Array[Function5[D, D, D, D, D, D]]] =
    macro derivativeMult_impl[Function5[D, D, D, D, D, D]]

  def jacobian(f: List[Function6[D, D, D, D, D, D, D]]):
    Array[Array[Function6[D, D, D, D, D, D, D]]] =
    macro derivativeMult_impl[Function6[D, D, D, D, D, D, D]]

  def jacobian(f: List[Function7[D, D, D, D, D, D, D, D]]):
    Array[Array[Function7[D, D, D, D, D, D, D, D]]] =
    macro derivativeMult_impl[Function7[D, D, D, D, D, D, D, D]]

  def jacobian(f: List[Function8[D, D, D, D, D, D, D, D, D]]):
    Array[Array[Function8[D, D, D, D, D, D, D, D, D]]] =
    macro derivativeMult_impl[Function8[D, D, D, D, D, D, D, D, D]]

  def jacobian(f: List[Function9[D, D, D, D, D, D, D, D, D, D]]):
    Array[Array[Function9[D, D, D, D, D, D, D, D, D, D]]] =
    macro derivativeMult_impl[Function9[D, D, D, D, D, D, D, D, D, D]]


  //--------------------------------------
  //         Function Translation
  //--------------------------------------
  def double2Interval(f: Function1[D, D]): Function1[I, I] =
    macro double2Interval_impl[Function1[I, I]]
  def double2Interval(f: Function2[D, D, D]): Function2[I, I, I] =
    macro double2Interval_impl[Function2[I, I, I]]
  def double2Interval(f: Function3[D, D, D, D]): Function3[I, I, I, I] =
    macro double2Interval_impl[Function3[I, I, I, I]]
  def double2Interval(f: Function4[D, D, D, D, D]): Function4[I, I, I, I, I] =
    macro double2Interval_impl[Function4[I, I, I, I, I]]
  def double2Interval(f: Function5[D, D, D, D, D, D]):
    Function5[I, I, I, I, I, I] =
    macro double2Interval_impl[Function5[I, I, I, I, I, I]]
  def double2Interval(f: Function6[D, D, D, D, D, D, D]):
    Function6[I, I, I, I, I, I, I] =
    macro double2Interval_impl[Function6[I, I, I, I, I, I, I]]
  def double2Interval(f: Function7[D, D, D, D, D, D, D, D]):
    Function7[I, I, I, I, I, I, I, I] =
    macro double2Interval_impl[Function7[I, I, I, I, I, I, I, I]]
  def double2Interval(f: Function8[D, D, D, D, D, D, D, D, D]):
    Function8[I, I, I, I, I, I, I, I, I] =
    macro double2Interval_impl[Function8[I, I, I, I, I, I, I, I, I]]
  def double2Interval(f: Function9[D, D, D, D, D, D, D, D, D, D]):
    Function9[I, I, I, I, I, I, I, I, I, I] =
    macro double2Interval_impl[Function9[I, I, I, I, I, I, I, I, I, I]]
  def double2Interval(f: Function10[D, D, D, D, D, D, D, D, D, D, D]):
    Function10[I, I, I, I, I, I, I, I, I, I, I] =
    macro double2Interval_impl[Function10[I, I, I, I, I, I, I, I, I, I, I]]


  def double2Affine(f: Function1[D, D]): Function1[A, A] =
    macro double2Affine_impl[Function1[A, A]]
  def double2Affine(f: Function2[D, D, D]): Function2[A, A, A] =
    macro double2Affine_impl[Function2[A, A, A]]
  def double2Affine(f: Function3[D, D, D, D]): Function3[A, A, A, A] =
    macro double2Affine_impl[Function3[A, A, A, A]]
  def double2Affine(f: Function4[D, D, D, D, D]): Function4[A, A, A, A, A] =
    macro double2Affine_impl[Function4[A, A, A, A, A]]
  def double2Affine(f: Function5[D, D, D, D, D, D]):
    Function5[A, A, A, A, A, A] =
    macro double2Affine_impl[Function5[A, A, A, A, A, A]]
  def double2Affine(f: Function6[D, D, D, D, D, D, D]):
    Function6[A, A, A, A, A, A, A] =
    macro double2Affine_impl[Function6[A, A, A, A, A, A, A]]
  def double2Affine(f: Function7[D, D, D, D, D, D, D, D]):
    Function7[A, A, A, A, A, A, A, A] =
    macro double2Affine_impl[Function7[A, A, A, A, A, A, A, A]]
  def double2Affine(f: Function8[D, D, D, D, D, D, D, D, D]):
    Function8[A, A, A, A, A, A, A, A, A] =
    macro double2Affine_impl[Function8[A, A, A, A, A, A, A, A, A]]
  def double2Affine(f: Function9[D, D, D, D, D, D, D, D, D, D]):
    Function9[A, A, A, A, A, A, A, A, A, A] =
    macro double2Affine_impl[Function9[A, A, A, A, A, A, A, A, A, A]]
  def double2Affine(f: Function10[D, D, D, D, D, D, D, D, D, D, D]):
    Function10[A, A, A, A, A, A, A, A, A, A, A] =
    macro double2Affine_impl[Function10[A, A, A, A, A, A, A, A, A, A, A]]


  /* ======================================
           Unary frontends
   ======================================== */
  def errorBound(f: (Double => Double), x: Double, tol: Double): Interval = macro errorBound_impl
  def assertBound(f: (Double => Double), x: Double, tol: Double): Interval = macro assertBound_impl
  def certify(root: Double, error: Interval): SmartFloat = {
    new SmartFloat(root, maxAbs(error))
  }

  def certify(f: (Double => Double), x: Double, tol: Double): SmartFloat = macro certify_impl

  // Has manually given derivative (will be used for comparison purposes)
  def errorManually(f: (Double => Double), fPrime: (Double => Double), x: Double,
    tol: Double): Interval = macro errorManually_impl


  /* ======================================
           2D frontends
   ======================================== */
  def errorBound(f1: ((Double, Double) => Double), f2:((Double, Double) => Double), x1: Double,
    x2: Double, tol: Double): Array[Interval] = macro errorBound2_impl
  def assertBound(f1: ((Double, Double) => Double), f2:((Double, Double) => Double), x1: Double,
    x2: Double, tol: Double): Array[Interval] = macro assertBound2_impl
  //def certify(f1: ((Double, Double) => Double), f2:((Double, Double) => Double), x1: Double,
  //  x2: Double, tol: Double): Array[SmartFloat] = macro certify2_impl
  def certify(roots: Array[Double], errors: Array[Interval]): Array[SmartFloat] = {
    val result = new Array[SmartFloat](roots.length)
    for (i <- 0 until roots.length)
      result(i) = new SmartFloat(roots(i), maxAbs(errors(i)))
    result
  }
  /* ======================================
           3D frontends
   ======================================== */
  def errorBound(f1: ((Double, Double, Double) => Double), f2:((Double, Double, Double) => Double),
    f3:((Double, Double, Double) => Double), x1: Double, x2: Double, x3: Double, tol: Double):
    Array[Interval] = macro errorBound3_impl
  def assertBound(f1: ((Double, Double, Double) => Double), f2:((Double, Double, Double) => Double),
    f3:((Double, Double, Double) => Double), x1: Double, x2: Double, x3: Double, tol: Double):
    Array[Interval] = macro assertBound3_impl
  //def certify(f1: ((Double, Double, Double) => Double), f2:((Double, Double, Double) => Double),
  //  f3:((Double, Double, Double) => Double), x1: Double, x2: Double, x3: Double, tol: Double):
  //  Array[SmartFloat] = macro certify3_impl

  /* ======================================
           Any-D frontends
   ======================================== */
  def errorBoundX(f: List[Any], x: Array[Double], tol: Double): Array[Interval] = macro errorBoundX_impl
  def assertBoundX(f: List[Any], x: Array[Double], tol: Double): Array[Interval] = macro assertBoundX_impl


  /* ===========================
          DERIVATIVES
  ==============================*/
  def derivative_impl(c: Context)(f: c.Expr[Double => Double]): c.Expr[Double => Double] = {
    import c.universe._
    val (fncTree, derTree) = unaryDerivatives[Double => Double](c)(f, "math", "Double")
    //if (debug) println(derTree)
    c.Expr[Double => Double](reify(derTree.splice).tree)
  }

  def derivative2_impl(c: Context)(f1: c.Expr[(Double, Double) => Double], f2: c.Expr[(Double, Double) => Double]):
    c.Expr[Array[Array[(Double, Double) => Double]]] = {
    import c.universe._
    val (fncTree, derTree) =
      multivariate2Derivatives[(Double, Double) => Double](c)(List(f1.tree, f2.tree), "math", "Double")
    c.Expr[Array[Array[(Double, Double) => Double]]](reify(derTree.splice).tree)
  }


  def derivative3_impl(c: Context)(f1: c.Expr[(Double, Double, Double) => Double],
    f2: c.Expr[(Double, Double, Double) => Double], f3: c.Expr[(Double, Double, Double) => Double]):
    c.Expr[Array[Array[(Double, Double, Double) => Double]]] = {
    import c.universe._
    val (fncTree, derTree) =
      multivariate2Derivatives[(Double, Double, Double) => Double](c)(
      List(f1.tree, f2.tree, f3.tree), "math", "Double")
    c.Expr[Array[Array[(Double, Double, Double) => Double]]](reify(derTree.splice).tree)
  }

  def derivativeX_impl(c: Context)(f: c.Expr[List[Any]]): c.Expr[Array[Array[Any]]] = {
    import c.universe._
    val fncList = extractFunctionsFromList(c)(f.tree)
    val (dim, fncTree, derTree) = multivariate2DerivativesXD[Double](c)(fncList, "math", "Double")
    c.Expr[Array[Array[Any]]](reify(derTree.splice).tree)
  }

  def derivativeMult_impl[T](c: Context)(f: c.Expr[List[T]]):
    c.Expr[Array[Array[T]]] = {
    import c.universe._
    val fncList = extractFunctionsFromList(c)(f.tree)
    val (fncTree, derTree) =
      multivariate2Derivatives[T](c)(fncList, "math", "Double")
    c.Expr[Array[Array[T]]](reify(derTree.splice).tree)
  }

  /* ===========================
           Function translation
  ============================== */
  def double2Interval_impl[T](c: Context)(f: c.Expr[Any]): c.Expr[T] = {
    import c.universe._
    val tmp = scala2Funktion(c)(f.tree)
    val fncTree = funktion2ScalaFnc[T](c)(Funktion(tmp.params, tmp.body), "common", "Interval")
    c.Expr[T](reify(fncTree.splice).tree)
  }

  def double2Affine_impl[T](c: Context)(f: c.Expr[Any]): c.Expr[T] = {
    import c.universe._
    val tmp = scala2Funktion(c)(f.tree)
    val fncTree = funktion2ScalaFnc[T](c)(Funktion(tmp.params, tmp.body), "affine", "AffineForm")
    c.Expr[T](reify(fncTree.splice).tree)
  }


  /* ===========================
          UNARY MACRO
  ==============================*/
  def errorBound_impl(c: Context)(f: c.Expr[Double => Double], x: c.Expr[Double],
    tol: c.Expr[Double]): c.Expr[Interval] = {
    import c.universe._
    if (useIntervalsForUnary) {
      val (fncTree, derTree) = unaryDerivatives[Interval => Interval](c)(f, "common", "Interval")
      c.Expr[Interval](reify(computeTrueErrorInterval(fncTree.splice, derTree.splice, x.splice, tol.splice)).tree)
    }
    else {
      val (fncTree, derTree) = unaryDerivatives[AffineForm => AffineForm](c)(f, "affine", "AffineForm")
      c.Expr[Interval](reify(computeTrueErrorAffine(fncTree.splice, derTree.splice, x.splice, tol.splice)).tree)
    }
  }

  def assertBound_impl(c: Context)(f: c.Expr[Double => Double], x: c.Expr[Double],
    tol: c.Expr[Double]): c.Expr[Interval] = {
    import c.universe._
    if (useIntervalsForUnary) {
      val (fncTree, derTree) = unaryDerivatives[Interval => Interval](c)(f, "common", "Interval")
      c.Expr[Interval](reify(
        assertTrueError(computeTrueErrorInterval(fncTree.splice, derTree.splice, x.splice, tol.splice), tol.splice)).tree)
    }
    else {
      val (fncTree, derTree) = unaryDerivatives[AffineForm => AffineForm](c)(f, "affine", "AffineForm")
      c.Expr[Interval](reify(
        assertTrueError(computeTrueErrorAffine(fncTree.splice, derTree.splice, x.splice, tol.splice), tol.splice)).tree)
    }
  }

  def certify_impl(c: Context)(f: c.Expr[Double => Double], x: c.Expr[Double],
    tol: c.Expr[Double]): c.Expr[SmartFloat] = {
    import c.universe._
    if (useIntervalsForUnary) {
      val (fncTree, derTree) = unaryDerivatives[Interval => Interval](c)(f, "common", "Interval")
      c.Expr[SmartFloat](reify(new SmartFloat(x.splice,
        maxAbs(computeTrueErrorInterval(fncTree.splice, derTree.splice, x.splice, tol.splice)))).tree)
    }
    else {
      val (fncTree, derTree) = unaryDerivatives[AffineForm => AffineForm](c)(f, "affine", "AffineForm")
      c.Expr[SmartFloat](reify(new SmartFloat(x.splice,
        maxAbs(computeTrueErrorAffine(fncTree.splice, derTree.splice, x.splice, tol.splice)))).tree)
    }
  }

  def errorManually_impl(c: Context)(f: c.Expr[Double => Double], fPrime: c.Expr[Double => Double],
    x: c.Expr[Double], tol: c.Expr[Double]): c.Expr[Interval] = {
    import c.universe._
    val myFnc = scala2Funktion(c)(f.tree)
    val myDer = scala2Funktion(c)(fPrime.tree)

    if (useIntervalsForUnary) {
      val fncTree = funktion2ScalaFnc[Interval => Interval](c)(myFnc, "common", "Interval")
      val derTree = funktion2ScalaFnc[Interval => Interval](c)(myDer, "common", "Interval")
      c.Expr[Interval](reify(computeTrueErrorInterval(fncTree.splice, derTree.splice, x.splice, tol.splice)).tree)
    }
    else {
      val fncTree = funktion2ScalaFnc[AffineForm => AffineForm](c)(myFnc, "affine", "AffineForm")
      val derTree = funktion2ScalaFnc[AffineForm => AffineForm](c)(myDer, "affine", "AffineForm")
      c.Expr[Interval](reify(computeTrueErrorAffine(fncTree.splice, derTree.splice, x.splice, tol.splice)).tree)
    }
  }


  private def unaryDerivatives[T](c: Context)(f: c.Expr[Double => Double], pckgName: String,
    typeName: String): (c.Expr[T], c.Expr[T]) = {
    import c.universe._
    // function itself
    val myFnc = scala2Funktion(c)(f.tree)
    val fncTree = funktion2ScalaFnc[T](c)(myFnc, pckgName, typeName)

    // derivative
    val constants = pullConst(myFnc.body)
    val simplified = simplify(constants)
    val derivative = Funktion(myFnc.params, d(myFnc.params.head)(simplified))
    //val derivative = Funktion(myFnc.params, d2(myFnc.params.head)(myFnc.body))
    // TODO: preprocess trees for.spliceuation
    val derTree = funktion2ScalaFnc[T](c)(derivative, pckgName, typeName)
    (fncTree, derTree)
  }

  def assertTrueError(error: Interval, tol: Double): Interval = {
    if ((-tol > error.xhi) || (error.xlo > tol))
      throw SolutionNotIncluded("Error " + error + " is greater than tolerance " + tol)
    if (error == EmptyInterval || error.xlo == MinusInf || error.xhi == PlusInf ||
       math.max(math.abs(error.xlo), math.abs(error.xhi)) >= tol)
       throw SolutionCannotBeCertified("Error " + error + " is possibly greater than tolerance" + tol)
    error
  }

  def maxAbs(i: Interval): Double = math.max(math.abs(i.xlo), math.abs(i.xhi))

  /* ===========================
          2D MACROS
  ==============================*/
  def errorBound2_impl(c: Context)(f1: c.Expr[(Double, Double) => Double],
    f2: c.Expr[(Double, Double) => Double], x1: c.Expr[Double], x2: c.Expr[Double],
    tol: c.Expr[Double]): c.Expr[Array[Interval]] = {
    import c.universe._

    val vectorTree = doubleVector2Scala(c)(List(x1.tree, x2.tree))
    if (useIntervalsForMulti) {
      val (fncTree, derTree) = multivariate2Derivatives[(Interval, Interval) => Interval](c)(
        List(f1.tree, f2.tree), "common", "Interval")
      c.Expr[Array[Interval]](reify(front2DInterval(fncTree.splice, derTree.splice, vectorTree.splice, tol.splice)).tree)
    }
    else {
      val (fncTree, derTree) = multivariate2Derivatives[(AffineForm, AffineForm) => AffineForm](c)(
        List(f1.tree, f2.tree), "affine", "AffineForm")
      c.Expr[Array[Interval]](reify(front2DAffine(fncTree.splice, derTree.splice, vectorTree.splice, tol.splice)).tree)
    }
  }

  def assertBound2_impl(c: Context)(f1: c.Expr[(Double, Double) => Double],
    f2: c.Expr[(Double, Double) => Double], x1: c.Expr[Double], x2: c.Expr[Double],
    tol: c.Expr[Double]): c.Expr[Array[Interval]] = {
    import c.universe._

    val vectorTree = doubleVector2Scala(c)(List(x1.tree, x2.tree))
    if (useIntervalsForMulti) {
      val (fncTree, derTree) = multivariate2Derivatives[(Interval, Interval) => Interval](c)(
        List(f1.tree, f2.tree), "common", "Interval")
      c.Expr[Array[Interval]](reify(
        assertTrueErrorX(front2DInterval(fncTree.splice, derTree.splice, vectorTree.splice, tol.splice), tol.splice)).tree)
    }
    else {
      val (fncTree, derTree) = multivariate2Derivatives[(AffineForm, AffineForm) => AffineForm](c)(
        List(f1.tree, f2.tree), "affine", "AffineForm")
      c.Expr[Array[Interval]](reify(
        assertTrueErrorX(front2DAffine(fncTree.splice, derTree.splice, vectorTree.splice, tol.splice), tol.splice)).tree)
    }
  }

  // Out of order until further notice...
  def certify2_impl(c: Context)(f1: c.Expr[(Double, Double) => Double],
    f2: c.Expr[(Double, Double) => Double], x1: c.Expr[Double], x2: c.Expr[Double],
    tol: c.Expr[Double]): c.Expr[Array[SmartFloat]] = {
    import c.universe._

    val vectorTree = doubleVector2Scala(c)(List(x1.tree, x2.tree))
    if (useIntervalsForMulti) {
      val (fncTree, derTree) = multivariate2Derivatives[(Interval, Interval) => Interval](c)(
        List(f1.tree, f2.tree), "common", "Interval")
      c.Expr[Array[SmartFloat]](reify(affine2SmartFloats(List(x1.splice, x2.splice),
        front2DInterval(fncTree.splice, derTree.splice, vectorTree.splice, tol.splice))).tree)
    }
    else {
      val (fncTree, derTree) = multivariate2Derivatives[(AffineForm, AffineForm) => AffineForm](c)(
        List(f1.tree, f2.tree), "affine", "AffineForm")
      c.Expr[Array[SmartFloat]](reify(affine2SmartFloats(List(x1.splice, x2.splice),
        front2DAffine(fncTree.splice, derTree.splice, vectorTree.splice, tol.splice))).tree)
    }
  }

  def affine2SmartFloats(x: List[Double], ints: Array[Interval]): Array[SmartFloat] = {
    val result = new Array[SmartFloat](x.length)
    for (i <- 0 until x.length)
      result(i) = new SmartFloat(x(i), maxAbs(ints(i)))
    result
  }

  def assertTrueErrorX(z: Array[Interval], tol: Double): Array[Interval] = {
    for (i <- 0 until z.length) {
      val zi = z(i)
      if ((-tol > zi.xhi) || (zi.xlo > tol))
        throw SolutionNotIncluded("Error in position " + i + " from " + z.toList +
          " is greater than tolerance " + tol)
      if (zi == EmptyInterval || zi.xlo == MinusInf || zi.xhi == PlusInf ||
       math.max(math.abs(zi.xlo), math.abs(zi.xhi)) >= tol)
        throw SolutionCannotBeCertified("Error in position " + i +" from " + z.toList +
          " is possibly greater than tolerance" + tol)
    }
    return z
  }


  /* ===========================
          3D MACROS
  ==============================*/
  def errorBound3_impl(c: Context)(f1: c.Expr[(Double, Double, Double) => Double],
    f2: c.Expr[(Double, Double, Double) => Double], f3: c.Expr[(Double, Double, Double) => Double],
    x1: c.Expr[Double], x2: c.Expr[Double], x3: c.Expr[Double], tol: c.Expr[Double]): c.Expr[Array[Interval]] = {
    import c.universe._

    val vectorTree = doubleVector2Scala(c)(List(x1.tree, x2.tree, x3.tree))
    if (useIntervalsForMulti) {
      val (fncTree, derTree) = multivariate2Derivatives[(Interval, Interval, Interval) => Interval](c)(
        List(f1.tree, f2.tree, f3.tree), "common", "Interval")
      c.Expr[Array[Interval]](reify(front3DInterval(fncTree.splice, derTree.splice, vectorTree.splice, tol.splice)).tree)
    }
    else {
      val (fncTree, derTree) =
        multivariate2Derivatives[(AffineForm, AffineForm, AffineForm) => AffineForm](c)(
        List(f1.tree, f2.tree, f3.tree), "affine", "AffineForm")
      c.Expr[Array[Interval]](reify(front3DAffine(fncTree.splice, derTree.splice, vectorTree.splice, tol.splice)).tree)
    }
  }

  def assertBound3_impl(c: Context)(f1: c.Expr[(Double, Double, Double) => Double],
    f2: c.Expr[(Double, Double, Double) => Double], f3: c.Expr[(Double, Double, Double) => Double],
    x1: c.Expr[Double], x2: c.Expr[Double], x3: c.Expr[Double], tol: c.Expr[Double]): c.Expr[Array[Interval]] = {
    import c.universe._

    val vectorTree = doubleVector2Scala(c)(List(x1.tree, x2.tree, x3.tree))
    if (useIntervalsForMulti) {
      val (fncTree, derTree) =
        multivariate2Derivatives[(Interval, Interval, Interval) => Interval](c)(
        List(f1.tree, f2.tree, f3.tree), "common", "Interval")
      c.Expr[Array[Interval]](reify(
        assertTrueErrorX(front3DInterval(fncTree.splice, derTree.splice, vectorTree.splice, tol.splice), tol.splice)).tree)
    }
    else {
      val (fncTree, derTree) =
        multivariate2Derivatives[(AffineForm, AffineForm, AffineForm) => AffineForm](c)(
        List(f1.tree, f2.tree, f3.tree), "affine", "AffineForm")
      c.Expr[Array[Interval]](reify(
        assertTrueErrorX(front3DAffine(fncTree.splice, derTree.splice, vectorTree.splice, tol.splice), tol.splice)).tree)
    }
  }

  // Out of order until further notice...
  /*
  def certify3_impl(c: Context)(f1: c.Expr[(Double, Double, Double) => Double],
    f2: c.Expr[(Double, Double, Double) => Double], f3: c.Expr[(Double, Double, Double) => Double],
    x1: c.Expr[Double], x2: c.Expr[Double], x3: c.Expr[Double], tol: c.Expr[Double]): c.Expr[Array[SmartFloat]] = {
    import c.universe._

    val vectorTree = doubleVector2Scala(c)(List(x1.tree, x2.tree, x3.tree))
    if (useIntervalsForMulti) {
      val (fncTree, derTree) = multivariate2Derivatives[(Double, Double, Double) => Double,
        (Interval, Interval, Interval) => Interval](c)(List(f1, f2, f3), "interval", "Interval")
      Expr[Array[SmartFloat]](reify(affine2SmartFloats(List(x1.splice, x2.splice, x3.splice),
        front3DInterval(fncTree.splice, derTree.splice, vectorTree.splice, tol.splice))).tree)
    }
    else {
      val (fncTree, derTree) = multivariate2Derivatives[(Double, Double, Double) => Double,
        (AffineForm, AffineForm, AffineForm) => AffineForm](c)(List(f1, f2, f3), "affine", "AffineForm")
      Expr[Array[SmartFloat]](reify(affine2SmartFloats(List(x1.splice, x2.splice, x3.splice),
        front3DAffine(fncTree.splice, derTree.splice, vectorTree.splice, tol.splice))).tree)
    }
  }*/

  /* ===========================
          Any-D MACROS
  ==============================*/
  def errorBoundX_impl(c: Context)(f: c.Expr[List[Any]],
    x: c.Expr[Array[Double]], tol: c.Expr[Double]): c.Expr[Array[Interval]] = {
    import c.universe._
    val fncList = extractFunctionsFromList(c)(f.tree)
    if (useIntervalsForMulti) {
      val (dim, fncTree, derTree) = multivariate2DerivativesXD[Interval](c)(fncList, "common", "Interval")
      val vectorTree = doubleVectorXD2Scala(c)(dim, x)
      c.Expr[Array[Interval]](reify(frontXDInterval(fncTree.splice, derTree.splice, vectorTree.splice, tol.splice)).tree)
    }
    else {
      val (dim, fncTree, derTree) = multivariate2DerivativesXD[AffineForm](c)(fncList, "affine", "AffineForm")
      val vectorTree = doubleVectorXD2Scala(c)(dim, x)
      c.Expr[Array[Interval]](reify(frontXDAffine(fncTree.splice, derTree.splice, vectorTree.splice, tol.splice)).tree)
    }
  }

  def assertBoundX_impl(c: Context)(f: c.Expr[List[Any]],
    x: c.Expr[Array[Double]], tol: c.Expr[Double]): c.Expr[Array[Interval]] = {
    import c.universe._
    val fncList = extractFunctionsFromList(c)(f.tree)
    if (useIntervalsForMulti) {
      val (dim, fncTree, derTree) = multivariate2DerivativesXD[Interval](c)(fncList, "common", "Interval")
      val vectorTree = doubleVectorXD2Scala(c)(dim, x)
      c.Expr[Array[Interval]](reify(
        assertTrueErrorX(frontXDInterval(fncTree.splice, derTree.splice, vectorTree.splice, tol.splice), tol.splice)).tree)
    }
    else {
      val (dim, fncTree, derTree) = multivariate2DerivativesXD[AffineForm](c)(fncList, "affine", "AffineForm")
      val vectorTree = doubleVectorXD2Scala(c)(dim, x)
      c.Expr[Array[Interval]](reify(
        assertTrueErrorX(frontXDAffine(fncTree.splice, derTree.splice, vectorTree.splice, tol.splice), tol.splice)).tree)
    }
  }




  /* ===========================
          Helper Functions
  ==============================*/
  private def double2Range[T](c: Context)(fnc: c.Tree,
    pckgName: String, typeName: String): c.Expr[T] = {
    import c.universe._
    val tmp = scala2Funktion(c)(fnc)
    //println(tmp)
    val fncTree = funktion2ScalaFnc[T](c)(Funktion(tmp.params, tmp.body), pckgName, typeName)
    //println(fncTree)
    fncTree
  }


  private def multivariate2Derivatives[T](c: Context)(fncs: List[c.Tree],
    pckgName: String, typeName: String): (c.Expr[Array[T]], c.Expr[Array[Array[T]]]) = {
    import c.universe._

    // TODO: preprocess trees for.spliceuation
    val dim = fncs.length
    val f = new Array[Expression](dim)
    var params: List[String] = null

    var tmp = scala2Funktion(c)(fncs(0))
    params = tmp.params
    f(0) = tmp.body

    for (i <- 1 until dim) {
      tmp = scala2Funktion(c)(fncs(i))
      if (params != tmp.params) {
        c.error(c.enclosingPosition, "The parameters of the functions do not match!")
        return null
      }
      f(i) = tmp.body
    }
    val jacobian = new Array[Array[Expression]](dim)
    for (i <- 0 until dim) {
      jacobian(i) = new Array[Expression](dim)
      val f_i = simplify(pullConst(f(i)))
      //val f_i = f(i)
      for (j <- 0 until dim) {
        jacobian(i)(j) = d(params(j))(f_i)
      }
    }

    val functions = Functions(params, f)
    val derivatives = Jacobian(params, jacobian)
    val fncTree = functions2ScalaFncs[T](c)(functions, pckgName, typeName)
    val derTree = jacobian2ScalaFncs[T](c)(derivatives, pckgName, typeName)
    (fncTree, derTree)
  }

  private def multivariate2DerivativesXD[S](c: Context)(fncs: List[c.Tree], pckgName: String,
    typeName: String): (Int, c.Expr[Array[Array[S] => S]], c.Expr[Array[Array[Array[S] => S]]]) = {
    import c.universe._

    // TODO: preprocess trees for.spliceuation
    val dim = fncs.length
    val f = new Array[Expression](dim)
    var params: List[String] = null

    var tmp = scala2Funktion(c)(fncs(0))
    params = tmp.params
    f(0) = tmp.body

    for (i <- 1 until dim) {
      tmp = scala2Funktion(c)(fncs(i))
      if (params != tmp.params) {
        c.error(c.enclosingPosition, "The parameters of the functions do not match!")
        return null
      }
      f(i) = tmp.body
    }
    val jacobian = new Array[Array[Expression]](dim)
    for (i <- 0 until dim) {
      jacobian(i) = new Array[Expression](dim)
      val f_i = simplify(pullConst(f(i)))
      //val f_i = f(i)
      for (j <- 0 until dim) {
        jacobian(i)(j) = d(params(j))(f_i)
      }
    }

    val functions = Functions(params, f)
    val derivatives = Jacobian(params, jacobian)
    val fncTree = functions2ScalaXD[S](c)(functions, pckgName, typeName)
    val derTree = jacobian2ScalaXD[S](c)(derivatives, pckgName, typeName)
    (dim, fncTree, derTree)
    //return (dim, null, null)
  }

}
