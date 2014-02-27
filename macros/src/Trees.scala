/*
 * Cerres (c) 2012-2014 EPFL, Lausanne
 */

package cerres.macros

import collection.immutable.HashMap
import ceres.affine._
import ceres.common._


trait Trees extends LinearAlgebra {

  case class Jacobian(params: List[String], exprs: Array[Array[Expression]]) {
    override def toString: String = {
      var str = params.toString + " => "
      for (row <- exprs) {
        for (e <- row) str = str + "\n\t" + pretty(e)
        str = str + "\n-------------------"
      }
      str
    }
  }

  case class Functions(params: List[String], exprs: Array[Expression]) {
    override def toString: String = {
      var str = params.toString + " => "
      for (e <- exprs)
        str = str + "\n\t" + pretty(e)
      str
    }
  }

  // params are sorted in the order the parameters are given
  case class Funktion(params: List[String], body: Expression)

  object Tpe extends Enumeration {
    type Tpe = Value
    val DoubleType, IntervalType, AffineType = Value
  }
  import Tpe._

  object Expression {
    def sqrt(e: Expression) = Sqrt(e)
    def log(e: Expression) = Log(e)
    def exp(e: Expression) = Exp(e)
    def sin(e: Expression) = Sin(e)
    def cos(e: Expression) = Cos(e)
    def tan(e: Expression) = Tan(e)
    def asin(e: Expression) = ASin(e)
    def acos(e: Expression) = ACos(e)
    def atan(e: Expression) = ATan(e)
    def pow(b: Expression, e: Expression) = Pow(b, e)
    val zero = Const(0.0)
    val one = Const(1.0)
  }

  sealed abstract class Expression {
    def +(e: Expression) = Add(this, e)
    def -(e: Expression) = Sub(this, e)
    def *(e: Expression) = Mult(this, e)
    def /(e: Expression) = Div(this, e)
    def unary_- = Neg(this)
  }

  case class Const(d: Double) extends Expression
  // bound variables
  case class Var(name: String) extends Expression
  // parameters of the function defined elsewhere
  // Need name at compile time, value at runtime...
  case class Param(name: String, tpe: Tpe) extends Expression
  case class Add(l: Expression, r: Expression) extends Expression
  case class Neg(e: Expression) extends Expression
  case class Sub(l: Expression, r: Expression) extends Expression
  case class Mult(l: Expression, r: Expression) extends Expression
  case class Div(l: Expression, r: Expression) extends Expression
  case class Sqrt(e: Expression) extends Expression
  case class Pow(base: Expression, exp: Expression) extends Expression
  case class Log(e: Expression) extends Expression
  case class Exp(e: Expression) extends Expression
  case class Sin(e: Expression) extends Expression
  case class Cos(e: Expression) extends Expression
  case class Tan(e: Expression) extends Expression
  case class ASin(e: Expression) extends Expression
  case class ACos(e: Expression) extends Expression
  case class ATan(e: Expression) extends Expression

  def d(v: String)(e: Expression): Expression = {
    import Expression._
    val newExpr = e match {
      case Mult(c: Const, x) => c * d(v)(x)
      case Mult(p: Param, x) => p * d(v)(x)
      case Mult(x, p: Param) => p * d(v)(x)
      case Mult(x, c: Const) => assert(false, "const multiply"); return null
      case Const(v) => zero
      case Var(id) => if(id == v) one else zero
      case Param(n, tpe) => zero
      case Add(l, r) => d(v)(l) + d(v)(r)
      case Sub(l, r) => d(v)(l) - d(v)(r)
      case Mult(l, r) => cleanUp(d(v)(l) * r) + cleanUp(l * d(v)(r))
      case Div(f, c: Const) => Div(cleanUp(d(v)(f)), c)
      case Div(f, g) => cleanUp(cleanUp(d(v)(f) * g) - cleanUp(d(v)(g) * f)) / (g*g)
      case Sqrt(e) => d(v)(e) / (Const(2.0) * Sqrt(e))
      case Pow(base, Const(dbl)) => Const(dbl) * cleanUp(d(v)(base) * cleanUp(Pow(base, Const(dbl - 1.0))))
      case Pow(base, p: Param) => p * cleanUp(d(v)(base) * cleanUp(Pow(base, Sub(p, one))))
      case Pow(b, e) => assert(false, "non-constant powers not allowed"); return null
      case Exp(x) => Exp(x) * d(v)(x)
      case Log(x) => d(v)(x) / x
      case Sin(x) => cos(x) * d(v)(x)
      case Cos(x) => (- sin(x)) * d(v)(x)
      case Tan(x) => (one + tan(x) * tan(x)) * d(v)(x)
      case ASin(x) => d(v)(x) / sqrt(one - x*x)
      case ACos(x) => cleanUp(- d(v)(x)) / sqrt(one - x*x)
      case ATan(x) => d(v)(x) / (one + x*x)
      case Neg(x) => Neg(cleanUp(d(v)(x)))
    }
    //println("\nd("+pretty(e)+") --> " + pretty(newExpr))
    val tmp = cleanUp(newExpr)
    //println("clean up: " + pretty(tmp))
    tmp
  }

  def d2(v: String)(e: Expression): Expression = {
    import Expression._
    val newExpr = e match {
      case Mult(c: Const, x) => c * d2(v)(x)
      case Mult(p: Param, x) => p * d2(v)(x)
      case Mult(x, p: Param) => p * d2(v)(x)
      case Mult(x, c: Const) => assert(false, "const multiply"); return null
      case Const(v) => zero
      case Var(id) => if(id == v) one else zero
      case Param(n, tpe) => zero
      case Add(l, r) => d2(v)(l) + d2(v)(r)
      case Sub(l, r) => d2(v)(l) - d2(v)(r)
      case Mult(l, r) => d2(v)(l) * r + l * d2(v)(r)
      case Div(f, c: Const) => Div(d2(v)(f), c)
      case Div(f, g) => (d2(v)(f) * g - d2(v)(g) * f) / (g*g)
      case Sqrt(e) => d2(v)(e) / (Const(2.0) * Sqrt(e))
      case Pow(base, Const(dbl)) => Const(dbl) * d2(v)(base) * Pow(base, Const(dbl - 1.0))
      case Pow(base, p: Param) => p * d2(v)(base) * Pow(base, Sub(p, one))
      case Pow(b, e) => assert(false, "non-constant powers not allowed"); return null
      case Exp(x) => Exp(x) * d2(v)(x)
      case Log(x) => d2(v)(x) / x
      case Sin(x) => cos(x) * d2(v)(x)
      case Cos(x) => (- sin(x)) * d2(v)(x)
      case Tan(x) => (one + tan(x) * tan(x)) * d2(v)(x)
      case ASin(x) => d2(v)(x) / sqrt(one - x*x)
      case ACos(x) => - d2(v)(x) / sqrt(one - x*x)
      case ATan(x) => d2(v)(x) / (one + x*x)
      case Neg(x) => Neg(d2(v)(x))
    }
    newExpr
  }


  // Try to preprocess the expression so that the derivative we get does not explode
  // One run is most likely not exhaustive
  def simplify(e: Expression): Expression = e match {
    case Add(l, r) => Add(simplify(l), simplify(r))
    case Neg(e) => Neg(simplify(e))
    case Sub(l, r) => Sub(simplify(l), simplify(r))

    case Mult(l, r) => val _l = simplify(l); val _r = simplify(r)
      (_l, _r) match {
        case (x, Pow(y, Const(dbl))) => if (x == y) return Pow(x, Const(dbl + 1.0))
        case (Pow(y, Const(dbl)), x) => if (x == y) return Pow(x, Const(dbl + 1.0))
        case (Pow(x, Const(d1)), Pow(y, Const(d2))) => if (x == y) return Pow(x, Const(d1 + d2))
        case (x, y) => if (x == y) return Pow(x, Const(2.0))
      }
      //else
      Mult(_l, _r)
    case Div(l, r) => Div(simplify(l), simplify(r))
    case Sqrt(e) => Sqrt(simplify(e))
    case Pow(base, exp) => Pow(simplify(base), simplify(exp))
    case Log(e) => Log(simplify(e))
    case Exp(e) => Exp(simplify(e))
    case Sin(e) => Sin(simplify(e))
    case Cos(e) => Cos(simplify(e))
    case Tan(e) => Tan(simplify(e))
    case ASin(e) => ASin(simplify(e))
    case ACos(e) => ACos(simplify(e))
    case ATan(e) => ATan(simplify(e))
    case _ => e
  }

  // Pulls out constants to the outest level possible.
  // Since this can be done recursively in just one go, we do it separately
  def pullConst(e: Expression): Expression = e match {
    case Mult(l, r) => val left = pullConst(l); val right = pullConst(r)
      (left, right) match {
        case (Mult(Const(d1), x), Mult(Const(d2), y)) => Mult(Const(d1*d2), Mult(x, y))
        case (x, Mult(c: Const, y)) => Mult(c, Mult(x, y))
        case (Mult(c: Const, x), y) => Mult(c, Mult(x, y))
        case (Mult(p1: Param, x), Mult(p2: Param, y)) => Mult(Mult(p1, p2), Mult(x, y))
        case (x, Mult(p: Param, y)) => Mult(p, Mult(x, y))
        case (Mult(p: Param, x), y) => Mult(p, Mult(x, y))
        case _ => Mult(left, right)
      }

    case Div(l, r) => val top = pullConst(l); val bottom = pullConst(r)
      (top, bottom) match {
        case (Mult(Const(a), x), Mult(Const(b), y)) => Mult( Const(a / b), Div(x, y) )
        case ( x, Mult(Const(b), y)) => Mult( Const(1.0 / b), Div(x, y) )
        case (Mult(a: Const, x), y) => Mult( a, Div(x, y) )
        case (Mult(a: Param, x), Mult(b: Param, y)) => Mult( Div(a, b), Div(x, y) )
        case ( x, Mult(b: Param, y)) => Mult( Div(Const(1.0), b), Div(x, y) )
        case (Mult(a: Param, x), y) => Mult( a, Div(x, y) )
        case _ => Div(top, bottom)
       }

    case Add(l, r) => Add(pullConst(l), pullConst(r))
    case Neg(e) => Neg(pullConst(e))
    case Sub(l, r) => Sub(pullConst(l), pullConst(r))
    case Sqrt(e) => Sqrt(pullConst(e))
    case Pow(base, exp) => Pow(pullConst(base), pullConst(exp))
    case Log(e) => Log(pullConst(e))
    case Exp(e) => Exp(pullConst(e))
    case Sin(e) => Sin(pullConst(e))
    case Cos(e) => Cos(pullConst(e))
    case Tan(e) => Tan(pullConst(e))
    case ASin(e) => ASin(pullConst(e))
    case ACos(e) => ACos(pullConst(e))
    case ATan(e) => ATan(pullConst(e))
    case _ => e
  }


  def cleanUp(e: Expression): Expression = e match {
    // Clean up zeroes
    case Add(Const(0.0), x) => x
    case Add(x, Const(0.0)) => x
    case Sub(x, Const(0.0)) => x
    case Sub(Const(0.0), x) => Neg(x)
    case Mult(Const(0.0), x) => Const(0.0)
    case Mult(x, Const(0.0)) => Const(0.0)
    case Div(Const(0.0), x) => Const(0.0)

    // Clean up trivial expressions
    case Mult(Const(1.0), x) => x
    case Mult(x, Const(1.0)) => x
    case Neg(Const(d)) => Const(-d)
    case Neg(Neg(x)) => x
    case Pow(b, Const(1.0)) => b
    case Mult(Neg(x), Neg(y)) => Mult(x, y)
    case Neg(Mult(x, Neg(y))) => Mult(x, y)
    case Neg(Mult(Neg(x), y)) => Mult(x, y)

    // Fixme: this is probably not great if we don't have integers
    case Add(c1: Const, c2: Const) => Const(c1.d + c2.d)
    case Mult(c1: Const, c2: Const) => Const(c1.d * c2.d)
    case Mult(Const(d1), Mult(Const(d2), x)) => Mult(Const(d1 * d2), x)
    case Mult(Const(d1), Div(Const(d2), x)) => Div(Const(d1 * d2), x)
    case _ => e
  }

  def pretty(expr: Expression): String = expr match {
    case Const(d) => d.toString
    case Var(name) => name
    case Param(name, tpe) => name
    case Add(l, r) => "(%s + %s)".format(pretty(l), pretty(r))
    case Neg(e) => "-%s".format(pretty(e))
    //case Sum(terms) => terms.foldLeft(0.0)((sum: Double, y: Expression) => sum + getValue(y, input))
    case Sub(l, r) => "(%s - %s)".format(pretty(l), pretty(r))
    case Mult(l, r) => "(%s * %s)".format(pretty(l), pretty(r))
    //case Product(terms) => terms.foldLeft(1.0)((prod: Double, y: Expression) => prod * getValue(y, input)
    case Div(l, r) => "(%s / %s)".format(pretty(l), pretty(r))
    case Sqrt(e) => "sqrt(%s)".format(pretty(e))
    case Pow(b, e) => "%s^(%s)".format(pretty(b), pretty(e))
    case Log(e) => "log(%s)".format(pretty(e))
    case Exp(e) => "exp(%s)".format(pretty(e))
    case Sin(e) => "sin(%s)".format(pretty(e))
    case Cos(e) => "cos(%s)".format(pretty(e))
    case Tan(e) => "tan(%s)".format(pretty(e))
    case ASin(e) => "asin(%s)".format(pretty(e))
    case ACos(e) => "acos(%s)".format(pretty(e))
    case ATan(e) => "atan(%s)".format(pretty(e))
  }
}
