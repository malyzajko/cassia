/*
 * Cerres (c) 2012-2014 EPFL, Lausanne
 */

package cerres.macros

import scala.reflect.macros.Context
import collection.mutable.ListBuffer

import ceres.common._
import ceres.affine._


trait ASTTranslation extends Trees {
  import Tpe._

  // We're looking for definitions, i.e. val x = ... or var fnc = ...
  /*def findTree(c: Context)(tree: c.Tree, name: String): c.Tree = {
    import c.universe._
    tree match {
      case Template(parents, self, body) =>
        println("self: " + self)
        for(b <- body) {
          println("--> " + b.getClass + "  " + b)
          findTree(c)(b, name)
        }
        return null
      case ValDef(modifier, name, tpt, rhs) =>
        println("valdef -> " + name.encoded + "  ,  " + rhs)
        return null
      case _ =>
        println("#####  Other Tree: " + tree)
        return null
    }
  }*/

  // Return the function or expression tree, assumes that it can be found in this tree though...
  def findInTemplate(c: Context)(tree: c.Tree, name: String): c.Tree = {
    import c.universe._
    def isThisMyFunction(n: String)(t: Tree): Boolean = {
      //println("Looking at : " + t)
      t match {
        case ValDef(modifier, termname, tpt, rhs) =>
          return termname.encoded == n
        case _ => return false
    }}
    tree match {
      case tmpl: Template =>
        tmpl.find(isThisMyFunction(name)) match {
          case Some(valdef: ValDef) => valdef.rhs
          case None => null
        }
      case df: DefDef =>
        df.find(isThisMyFunction(name)) match {
          case Some(valdef: ValDef) => valdef.rhs
          case None => null
        }
      case md: ModuleDef =>
        md.find(isThisMyFunction(name)) match {
          case Some(valdef: ValDef) => valdef.rhs
          case None => null
        }
      case _ =>
        println("Unknown enclosing tree type: " + tree.getClass)
        return null
    }
  }

  /**
   *  Converts an AST List into a List of AST Trees.
   *  It does not actually care whether the AST in the list are Functions.
   */
  def extractFunctionsFromList(c: Context)(f: c.Tree): List[c.Tree] = {
    import c.universe._
    f match {
      case Apply(TypeApply(Select(Select(This(immutable), list), apply), tpe), args) =>
        if (immutable.encoded == "immutable" && list.encoded == "List" && apply.encoded == "apply") {
          return args
        }
        else {
          c.error(c.enclosingPosition, "Expected a list of functions, found: " + c.universe.showRaw(f))
            return null
        }
      case _ =>
        c.error(c.enclosingPosition, "Expected a list of functions, found: " + c.universe.showRaw(f))
          return null
    }
  }

  def scala2Funktion(c: Context)(f: c.Tree): Funktion = {
    import c.universe._
    //println(">>>>>>>>> enclosing class")
    //println(c.enclosingClass)

    //println("\n>>>>>>>>> enclosing method")
    //println(c.enclosingClass)

    //println("\n and we're looking for: " + f)
    f match {
      case Function(parameter, body) =>
        var params = ListBuffer[String]()
        for (p <- parameter) params += p.name.encoded
        Funktion(params.toList, scala2Expression(c)(body, params.toList))
      // Function declaration in enclosing class
      case Select(qualifier, name) =>
        // compiler crashes if we inline this
        val missingFunction = findInTemplate(c)(c.enclosingClass, name.encoded)
        if (missingFunction == null) {
          c.error(c.enclosingPosition, "Function definition not found: " + name.encoded)
          return null
        }
        val tmp = scala2Funktion(c)(missingFunction)
        return tmp
      // Function declaration within the current method
      case Ident(name) =>
        val missingFunction = findInTemplate(c)(c.enclosingMethod, name.encoded)
        if (missingFunction == null) {
          c.error(c.enclosingPosition, "Function definition not found: " + name.encoded)
          return null
        }
        val tmp = scala2Funktion(c)(missingFunction)
        return tmp
      case _ =>
        c.error(c.enclosingPosition, "Expected Function, found " + c.universe.showRaw(f))
        return null
    }
  }

  // We need to keep track which variables are given to the function, and which
  // are parameters... For those we need to find some way of getting the value
  def scala2Expression(c: Context)(e: c.Tree, vars: List[String]): Expression = {
    import c.universe._
    e match {
      case Ident(name) =>
        if (vars.contains(name.encoded)) Var(name.encoded)
        else Param(name.encoded, DoubleType)

      case Literal(cnst) => //Const(cnst.toDouble)
        cnst match {
          case Constant(b: Byte) => Const(b.toDouble)
          case Constant(s: Short) => Const(s.toDouble)
          case Constant(i: Int) => Const(i.toDouble)
          case Constant(l: Long) => Const(l.toDouble)
          case Constant(f: Float) => Const(f.toDouble)
          case Constant(d: Double) => Const(d)
          case _ => c.error(c.enclosingPosition,
            "Constant not recognized, expected Double, found: " + c.universe.showRaw(e))
            return null
        }

      case Select(This(clsName), ident) => Param(ident.encoded, DoubleType)

      case Select(lhs, op) =>
        if (op.encoded == "unary_$minus")
          Neg(scala2Expression(c)(lhs, vars))

        // TODO: Extend this to mode complex expressions
        else if (op.encoded == "mid" && lhs.isInstanceOf[Ident])
          Param(lhs.asInstanceOf[Ident].name.encoded, IntervalType)
        else {
          c.error(c.enclosingPosition, "Unknown operation: " + c.universe.showRaw(e))
          return null
        }

      case Apply(Select(Select(Select(Ident(q1), q2), q3), q4), args) =>
        if (q1.encoded == "scala" && q2.encoded == "math" && q3.encoded == "package") {
          q4.encoded match {
            case "sqrt" => Sqrt(scala2Expression(c)(args.head, vars))
            case "log" => Log(scala2Expression(c)(args.head, vars))
            case "exp" => Exp(scala2Expression(c)(args.head, vars))
            case "sin" => Sin(scala2Expression(c)(args.head, vars))
            case "cos" => Cos(scala2Expression(c)(args.head, vars))
            case "tan" => Tan(scala2Expression(c)(args.head, vars))
            case "asin" => ASin(scala2Expression(c)(args.head, vars))
            case "acos" => ACos(scala2Expression(c)(args.head, vars))
            case "atan" => ATan(scala2Expression(c)(args.head, vars))
            case "pow" => Pow(scala2Expression(c)(args.head, vars), scala2Expression(c)(args(1), vars))
            case _ =>
              c.error(c.enclosingPosition, "Unknown operation: " + c.universe.showRaw(e))
              return null
          }
        }
        else {
          c.error(c.enclosingPosition, "Unknown operation: " + c.universe.showRaw(e))
          return null
        }

      case Apply(Ident(name), args) =>
        name.encoded match {
          case "sqrt" => Sqrt(scala2Expression(c)(args.head, vars))
          case "log" => Log(scala2Expression(c)(args.head, vars))
          case "exp" => Exp(scala2Expression(c)(args.head, vars))
          case "sin" => Sin(scala2Expression(c)(args.head, vars))
          case "cos" => Cos(scala2Expression(c)(args.head, vars))
          case "tan" => Tan(scala2Expression(c)(args.head, vars))
          case "asin" => ASin(scala2Expression(c)(args.head, vars))
          case "acos" => ACos(scala2Expression(c)(args.head, vars))
          case "atan" => ATan(scala2Expression(c)(args.head, vars))
          case "pow" => Pow(scala2Expression(c)(args.head, vars), scala2Expression(c)(args(1), vars))
          case _ =>
            c.error(c.enclosingPosition, "Unknown operation: " + c.universe.showRaw(e))
            return null
        }

      case Apply(Select(lhs, name), args) =>
        name.encoded match {
          case "$plus" => Add(scala2Expression(c)(lhs, vars), scala2Expression(c)(args.head, vars))
          case "$minus" => Sub(scala2Expression(c)(lhs, vars), scala2Expression(c)(args.head, vars))
          case "$times" =>
            // Make sure constants are always left
            val right = scala2Expression(c)(args.head, vars)
            right match {
              case Const(d) => Mult(right, scala2Expression(c)(lhs, vars))
              case _ => Mult(scala2Expression(c)(lhs, vars), right)
            }
          case "$div" => Div(scala2Expression(c)(lhs, vars), scala2Expression(c)(args.head, vars))

          case _ =>
            c.error(c.enclosingPosition, "Unknown operation: " + c.universe.showRaw(e))
            return null
        }
      case _ =>
        c.error(c.enclosingPosition, "Tree not recognized: " + c.universe.showRaw(e))
        return null
    }
  }

  def doubleVector2Scala(c: Context)(v: List[c.Tree]): c.Expr[DoubleVector] = {
    import c.universe._
    val constructor = Select(New(Ident(newTypeName("DoubleVector"))), newTermName("<init>"))
    val body = array2Scala(c)(v)
    c.Expr[DoubleVector](Apply(constructor, List(Literal(Constant(v.size)), body)))
  }

  def doubleVectorXD2Scala(c: Context)(dim: Int, v: c.Expr[Array[Double]]): c.Expr[DoubleVector] = {
    import c.universe._
    val constructor = Select(New(Ident(newTypeName("DoubleVector"))), newTermName("<init>"))
    c.Expr[DoubleVector](Apply(constructor, List(Literal(Constant(dim)), v.tree)))
  }

  def functions2ScalaFncs[T](c: Context)(fncs: Functions, packageName: String,
    typeName: String): c.Expr[Array[T]] = {
    import c.universe._
    import Flag._
    val path = Select(Select(Ident(newTermName("ceres")), newTermName(packageName)),
      newTermName(typeName))

    val functionBuffer = ListBuffer[Tree]()
    for (e <- fncs.exprs) {
      val paramBuffer = ListBuffer[ValDef]()
      for (paramName <- fncs.params)
        paramBuffer += ValDef(Modifiers(PARAM), newTermName(paramName),
          Ident(newTypeName(typeName)), EmptyTree)
      val params = paramBuffer.toList
      val expr =
          if (typeName == "Interval") expr2ScalaFnc(c)(e, IntervalType, path)
          else if (typeName == "AffineForm") expr2ScalaFnc(c)(e, AffineType, path)
          else expr2ScalaFnc(c)(e, DoubleType, path)
      functionBuffer += Function(params, expr)
    }
    val body = array2Scala(c)(functionBuffer.toList)
    c.Expr[Array[T]](body)
  }

  def functions2ScalaXD[T](c: Context)(fncs: Functions, packageName: String,
    typeName: String): c.Expr[Array[Array[T] => T]] = {
    import c.universe._
    import Flag._
    val path = Select(Select(Ident(newTermName("ceres")), newTermName(packageName)),
      newTermName(typeName))

    var paramMap = new scala.collection.immutable.HashMap[String, Int]
    for (i <- 0 until fncs.params.length) {
      paramMap += ((fncs.params(i), i))
    }

    val functionBuffer = ListBuffer[Tree]()
    for (e <- fncs.exprs) {
      val params = List(ValDef(Modifiers(PARAM), newTermName("$x"),
        AppliedTypeTree(Ident(newTypeName("Array")), List(Ident(newTypeName(typeName)))), EmptyTree))
      val expr =
          if (typeName == "Interval") exprXD(c)(e, paramMap, IntervalType, path)
          else if (typeName == "AffineForm") exprXD(c)(e, paramMap, AffineType, path)
          else exprXD(c)(e, paramMap, DoubleType, path)
      functionBuffer += Function(params, expr)
    }
    val body = array2Scala(c)(functionBuffer.toList)
      c.Expr[Array[Array[T]=> T]](body)
  }

  def jacobian2ScalaXD[T](c: Context)(j: Jacobian, packageName: String,
    typeName: String): c.Expr[Array[Array[Array[T] => T]]] = {
    import c.universe._
    import Flag._
    val path = getPath(c)(packageName, typeName)

    var paramMap = new scala.collection.immutable.HashMap[String, Int]
    for (i <- 0 until j.params.length) {
      paramMap += ((j.params(i), i))
    }

    val rowBuffer = ListBuffer[Tree]()
    for (row <- j.exprs) {
      val fncBuffer = ListBuffer[Tree]()
      for (e <- row) {
        val params = List(ValDef(Modifiers(PARAM), newTermName("$x"),
          AppliedTypeTree(Ident(newTypeName("Array")), List(Ident(newTypeName(typeName)))), EmptyTree))
        val expr =
          if (typeName == "Interval") exprXD(c)(e, paramMap, IntervalType, path)
          else if (typeName == "AffineForm") exprXD(c)(e, paramMap, AffineType, path)
          else exprXD(c)(e, paramMap, DoubleType, path)
        fncBuffer += Function(params, expr)
      }
      rowBuffer += array2Scala(c)(fncBuffer.toList)
    }
    val body = array2Scala(c)(rowBuffer.toList)
    c.Expr[Array[Array[Array[T] => T]]](body)
  }

  def jacobian2ScalaFunction[T](c: Context)(j: Jacobian, packageName: String,
    typeName: String): c.Expr[Array[Array[Function2[Double, Double, Double]]]] = {
    import c.universe._
    import Flag._
    val path = getPath(c)(packageName, typeName)

    var paramMap = new scala.collection.immutable.HashMap[String, Int]
    for (i <- 0 until j.params.length) {
      paramMap += ((j.params(i), i))
    }

    val rowBuffer = ListBuffer[Tree]()
    for (row <- j.exprs) {
      val fncBuffer = ListBuffer[Tree]()
      for (e <- row) {
        val params = List(ValDef(Modifiers(PARAM), newTermName("$x"),
          AppliedTypeTree(Ident(newTypeName("Array")), List(Ident(newTypeName(typeName)))), EmptyTree))
        val expr =
          if (typeName == "Interval") exprXD(c)(e, paramMap, IntervalType, path)
          else if (typeName == "AffineForm") exprXD(c)(e, paramMap, AffineType, path)
          else exprXD(c)(e, paramMap, DoubleType, path)
        fncBuffer += Function(params, expr)
      }
      rowBuffer += array2Scala(c)(fncBuffer.toList)
    }
    val body = array2Scala(c)(rowBuffer.toList)
    c.Expr[Array[Array[Function2[Double, Double, Double]]]](body)
  }


  def jacobian2ScalaFncs[T](c: Context)(j: Jacobian, packageName: String,
    typeName: String): c.Expr[Array[Array[T]]] = {
    import c.universe._
    import Flag._
    val path = getPath(c)(packageName, typeName)

    val rowBuffer = ListBuffer[Tree]()
    for (row <- j.exprs) {
      val fncBuffer = ListBuffer[Tree]()
      for (e <- row) {
        val paramBuffer = ListBuffer[ValDef]()
        for (paramName <- j.params)
          paramBuffer += ValDef(Modifiers(PARAM), newTermName(paramName),
            Ident(newTypeName(typeName)), EmptyTree)
        val params = paramBuffer.toList
        val expr =
          if (typeName == "Interval") expr2ScalaFnc(c)(e, IntervalType, path)
          else if (typeName == "AffineForm") expr2ScalaFnc(c)(e, AffineType, path)
          else expr2ScalaFnc(c)(e, DoubleType, path)
        fncBuffer += Function(params, expr)
      }
      rowBuffer += array2Scala(c)(fncBuffer.toList)
    }
    val body = array2Scala(c)(rowBuffer.toList)
    c.Expr[Array[Array[T]]](body)
  }

  def getPath(c: Context)(packageName: String, typeName: String): c.Tree = {
    import c.universe._
    if (typeName == "Interval" || typeName == "AffineForm")
      Select(Select(Ident(newTermName("ceres")), newTermName(packageName)), newTermName(typeName))
    else
      Select(Select(Ident(newTermName("scala")), newTermName("math")), newTermName("package"))
  }

  def array2Scala(c: Context)(list: List[c.Tree]): c.Tree = {
    import c.universe._
    Apply(Select(Select(Ident(newTermName("scala")), newTermName("Array")), newTermName("apply")), list)
  }

  def funktion2ScalaFnc[T](c: Context)(f: Funktion, packageName: String,
    typeName: String): c.Expr[T] = {
    import c.universe._
    import Flag._
    val path = getPath(c)(packageName, typeName)

    val paramBuffer = ListBuffer[ValDef]()
    for (paramName <- f.params)
      paramBuffer += ValDef(Modifiers(PARAM), newTermName(paramName),
          Ident(newTypeName(typeName)), EmptyTree)
    val params = paramBuffer.toList
    val expr =
      if (typeName == "Interval") expr2ScalaFnc(c)(f.body, IntervalType, path)
      else if (typeName == "AffineForm") expr2ScalaFnc(c)(f.body, AffineType, path)
      else expr2ScalaFnc(c)(f.body, DoubleType, path)
    c.Expr[T](Function(params, expr))
  }

  // If we have function from Double to Double, we don't need to wrap...
  def expr2ScalaFnc(c: Context)(expr: Expression, fncType: Tpe, path: c.Tree): c.Tree = {
    import c.universe._
    expr match {
      case Const(d) => fncType match {
        case DoubleType => Literal(Constant(d))
        case _ => Apply(Select(path, newTermName("apply")), List(Literal(Constant(d))))
      }
      case Var(n) => Ident(newTermName(n))
      case Param(n, tpe) =>
        (tpe, fncType) match {
          case (DoubleType, DoubleType) => Ident(newTermName(n))
          case (DoubleType, _ ) => Apply(Select(path, newTermName("apply")), List(Ident(newTermName(n))))

          case (IntervalType, IntervalType) => Ident(newTermName(n))
          case (IntervalType, AffineType) => Apply(Select(path, newTermName("apply")), List(Ident(newTermName(n))))
          case (IntervalType, DoubleType) => Select(Ident(newTermName(n)), newTermName("mid"))
          case _ => c.error(c.enclosingPosition, "Parameter " + n + " has invalid type: " + tpe); null
        }
      case Neg(e) => Select(expr2ScalaFnc(c)(e, fncType, path), newTermName("unary_$minus"))
      case Add(l, r) =>
        Apply(Select(expr2ScalaFnc(c)(l, fncType, path), newTermName("$plus")),
          List(expr2ScalaFnc(c)(r, fncType, path)))
      case Sub(l, r) =>
        Apply(Select(expr2ScalaFnc(c)(l, fncType, path), newTermName("$minus")),
          List(expr2ScalaFnc(c)(r, fncType, path)))
      case Mult(l, r) =>
        Apply(Select(expr2ScalaFnc(c)(l, fncType, path), newTermName("$times")),
          List(expr2ScalaFnc(c)(r, fncType, path)))
      case Div(l, r) =>
        Apply(Select(expr2ScalaFnc(c)(l, fncType, path), newTermName("$div")),
          List(expr2ScalaFnc(c)(r, fncType, path)))

      case Sqrt(e) => Apply(Select(path, newTermName("sqrt")), List(expr2ScalaFnc(c)(e, fncType, path)))
      case Pow(b, e) =>
        Apply(Select(path, newTermName("pow")),
          List(expr2ScalaFnc(c)(b, fncType, path), expr2ScalaFnc(c)(e, fncType, path)))
      case Log(e) => Apply(Select(path, newTermName("log")), List(expr2ScalaFnc(c)(e, fncType, path)))
      case Exp(e) => Apply(Select(path, newTermName("exp")), List(expr2ScalaFnc(c)(e, fncType, path)))
      case Sin(e) => Apply(Select(path, newTermName("sin")), List(expr2ScalaFnc(c)(e, fncType, path)))
      case Cos(e) => Apply(Select(path, newTermName("cos")), List(expr2ScalaFnc(c)(e, fncType, path)))
      case Tan(e) => Apply(Select(path, newTermName("tan")), List(expr2ScalaFnc(c)(e, fncType, path)))
      case ASin(e) => Apply(Select(path, newTermName("asin")), List(expr2ScalaFnc(c)(e, fncType, path)))
      case ACos(e) => Apply(Select(path, newTermName("acos")), List(expr2ScalaFnc(c)(e, fncType, path)))
      case ATan(e) => Apply(Select(path, newTermName("atan")), List(expr2ScalaFnc(c)(e, fncType, path)))
    }
  }

  def exprXD(c: Context)(expr: Expression, paramMap: Map[String, Int], fncType: Tpe,
    path: c.Tree): c.Tree = {
    import c.universe._
    expr match {
      case Const(d) => fncType match {
        case DoubleType => Literal(Constant(d))
        case _ => Apply(Select(path, newTermName("apply")), List(Literal(Constant(d))))
      }
      case Var(n) => Apply(Select(Ident(newTermName("$x")), newTermName("apply")),
        List(Literal(Constant(paramMap(n)))))
      case Param(n, tpe) =>
        (tpe, fncType) match {
          case (DoubleType, DoubleType) => Ident(newTermName(n))
          case (DoubleType, _ ) => Apply(Select(path, newTermName("apply")), List(Ident(newTermName(n))))

          case (IntervalType, IntervalType) => Ident(newTermName(n))
          case (IntervalType, AffineType) => Apply(Select(path, newTermName("apply")), List(Ident(newTermName(n))))
          case (IntervalType, DoubleType) => Select(Ident(newTermName(n)), newTermName("mid"))
          case _ => c.error(c.enclosingPosition, "Parameter " + n + " has invalid type: " + tpe); null
        }
      case Neg(e) => Select(exprXD(c)(e, paramMap, fncType, path), newTermName("unary_$minus"))
      case Add(l, r) =>
        Apply(Select(exprXD(c)(l, paramMap, fncType, path), newTermName("$plus")),
          List(exprXD(c)(r, paramMap, fncType, path)))
      case Sub(l, r) =>
        Apply(Select(exprXD(c)(l, paramMap, fncType, path), newTermName("$minus")),
          List(exprXD(c)(r, paramMap, fncType, path)))
      case Mult(l, r) =>
        Apply(Select(exprXD(c)(l, paramMap, fncType, path), newTermName("$times")),
          List(exprXD(c)(r, paramMap, fncType, path)))
      case Div(l, r) =>
        Apply(Select(exprXD(c)(l, paramMap, fncType, path), newTermName("$div")),
          List(exprXD(c)(r, paramMap, fncType, path)))

      case Sqrt(e) => Apply(Select(path, newTermName("sqrt")), List(exprXD(c)(e, paramMap, fncType, path)))
      case Pow(b, e) =>
        Apply(Select(path, newTermName("pow")),
          List(exprXD(c)(b, paramMap, fncType, path), exprXD(c)(e, paramMap, fncType, path)))
      case Log(e) => Apply(Select(path, newTermName("log")), List(exprXD(c)(e, paramMap, fncType, path)))
      case Exp(e) => Apply(Select(path, newTermName("exp")), List(exprXD(c)(e, paramMap, fncType, path)))
      case Sin(e) => Apply(Select(path, newTermName("sin")), List(exprXD(c)(e, paramMap, fncType, path)))
      case Cos(e) => Apply(Select(path, newTermName("cos")), List(exprXD(c)(e, paramMap, fncType, path)))
      case Tan(e) => Apply(Select(path, newTermName("tan")), List(exprXD(c)(e, paramMap, fncType, path)))
      case ASin(e) => Apply(Select(path, newTermName("asin")), List(exprXD(c)(e, paramMap, fncType, path)))
      case ACos(e) => Apply(Select(path, newTermName("acos")), List(exprXD(c)(e, paramMap, fncType, path)))
      case ATan(e) => Apply(Select(path, newTermName("atan")), List(exprXD(c)(e, paramMap, fncType, path)))
    }
  }

}
