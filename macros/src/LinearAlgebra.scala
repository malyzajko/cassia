/*
 * Cerres (c) 2012-2014 EPFL, Lausanne
 */

package cerres.macros

import math._
import ceres.common._
import ceres.affine._
import ceres.common.{QuadDouble => QD}

trait LinearAlgebra {

  abstract class Vector[T](val dim: Int, val data: Array[T]) {
    def update(i: Int, y: T) = data(i) = y
    def apply(i: Int) = data(i)
    override def toString: String = data.toList.toString
  }

  class DoubleVector(d: Int, a: Array[Double]) extends Vector[Double](d, a) {
    def this(d: Int) = this(d, new Array[Double](d))

    def +(y: DoubleVector): DoubleVector = {
      assert(y.dim == dim)
      val result = new DoubleVector(dim)
      for (i <- 0 until dim) result(i) = this.data(i) + y.data(i)
      result
    }
    def unary_-(): DoubleVector = {
      val result = new DoubleVector(dim)
      for (i <- 0 until dim) result(i) = -this.data(i)
      result
    }
    def norm: Double = {
      var sum = 0.0
      for (i <- 0 until dim) sum += data(i) * data(i)
      math.sqrt(sum)
    }
    def +/-(t: Double): AffineVector = {
      assert(t >= 0.0)
      val result = new AffineVector(dim)
      for (i <- 0 until dim) result(i) = new AForm(data(i), t)
      result
    }

    def plusMinus(t: Double): IntervalVector = {
      assert(t >= 0.0)
      val result = new IntervalVector(dim)
        // Fixme: not properly rounded...
      for (i <- 0 until dim) result(i) = Interval(data(i) - t, data(i) + t)
      result
    }

    def toIntervalArray: Array[Interval] = {
      val result = new Array[Interval](dim)
      for (i <- 0 until dim) result(i) = Interval(data(i))
      result
    }

    def toAffineArray: Array[AffineForm] = {
      val result = new Array[AffineForm](dim)
      for (i <- 0 until dim) result(i) = new AForm(data(i))
      result
    }
  }


  class AffineVector(d: Int, a: Array[AffineForm]) extends Vector[AffineForm](d, a) {
    def this(d: Int) = this(d, new Array[AffineForm](d))

    def unary_-(): AffineVector = {
      val result = new AffineVector(dim)
      for (i <- 0 until dim) result(i) = -this.data(i)
      result
    }
    override def toString: String = {
      val array = new Array[String](dim)
      for (i <- 0 until dim) array(i) = "[%s, %s]".format(data(i).interval.xlo.toString, data(i).interval.xhi.toString)
      array.toList.toString
    }
    def +(y: AffineVector): AffineVector = {
      val result = new AffineVector(dim)
      for (i <- 0 until dim) result(i) = data(i) + y(i)
      result
    }
    def toInterval: IntervalVector = {
      val result = new IntervalVector(dim)
      for (i <- 0 until dim) result(i) = data(i).interval
      result
    }
  }

  class IntervalVector(d: Int, a: Array[Interval]) extends Vector[Interval](d, a) {
    def this(d: Int) = this(d, new Array[Interval](d))

    def unary_-(): IntervalVector = {
      val result = new IntervalVector(dim)
      for (i <- 0 until dim) result(i) = -this.data(i)
      result
    }
    override def toString: String = {
      val array = new Array[String](dim)
      for (i <- 0 until dim) array(i) = "[%s, %s]".format(data(i).xlo.toString, data(i).xhi.toString)
      array.toList.toString
    }
    def +(y: IntervalVector): IntervalVector = {
      val result = new IntervalVector(dim)
      for (i <- 0 until dim) result(i) = data(i) + y(i)
      result
    }
  }

  class Fnc2Vector(d: Int, a: Array[(Double, Double) => Double]) extends Vector[(Double, Double) => Double](d, a) {
    def this(d: Int) = this(d, new Array[(Double, Double) => Double](d))
    def eval(x: DoubleVector): DoubleVector = {
      assert(x.dim == this.dim)
      val result = new DoubleVector(dim)
      for (i <- 0 until dim) {
        result(i) = data(i)(x(0), x(1))
      }
      result
    }
  }

  class Fnc2AAVector(d: Int, a: Array[(AffineForm, AffineForm) => AffineForm]) extends Vector[(AffineForm, AffineForm) => AffineForm](d, a) {
    def this(d: Int) = this(d, new Array[(AffineForm, AffineForm) => AffineForm](d))
    def eval(x: DoubleVector): AffineVector = {
      assert(x.dim == this.dim)
      val result = new AffineVector(dim)
      for (i <- 0 until dim) {
        result(i) = data(i)(new AForm(x(0)), new AForm(x(1)))
      }
      result
    }
  }




  def identityMatrix(d: Int): AffineMatrix = {
    val result = new AffineMatrix(d)
    for (i <- 0 until d)
      result(i, i) = 1.0
    result
  }

  def intervalIdentityMatrix(d: Int): IntervalMatrix = {
    val result = new IntervalMatrix(d)
    for (i <- 0 until d)
      result(i, i) = 1.0
    result
  }

  // We only want square matrices
  // Fixme: constructor may not need to always initialize the arrays
  abstract class Matrix[T](val dim: Int) {
    var data = getNewArray
    def getNewArray: Array[Array[T]]
    def update(i: Int, j: Int, y: T) = data(i)(j) = y
    def apply(i: Int, j: Int) = data(i)(j)
    override def toString: String = {
      var str = ""
      for (i <- 0 until dim)
        str = str + data(i).toList + "\n"
      str
    }
  }

  class DoubleMatrix(d: Int) extends Matrix[Double](d) {
    def getNewArray = Array.fill(dim, dim)(0.0)

    def this(a: Array[Array[Double]]) = {
      this(a.length)
      data = a
    }
    // computes an approximate inverse
    def inverse: DoubleMatrix = {
      val lu = copyMatrix(data)
      val pivot = Array.fill(dim){0}
      factor(lu, pivot)  // computed LU factorization

      val inverse = new DoubleMatrix(dim)

      for (i <- 0 until dim) {
        val e = Array.fill(dim){0.0}; e(i) = 1.0  // unit vector
        solve(lu, pivot, e) // get the i-th column
        //copy e into inverse
        for (j <- 0 until dim) inverse(j, i) = e(j)
      }
      inverse
    }

    def *(v: AffineVector): AffineVector = {
      val result = new AffineVector(dim)
      for (i <- 0 until dim) {
        var sum: AffineForm = new AForm(0.0)
        for (j <- 0 until dim)
          sum = sum + AffineForm(data(i)(j)) * v(j)
        result(i) = sum
      }
      result
    }

    def *(v: IntervalVector): IntervalVector = {
      val result = new IntervalVector(dim)
      for (i <- 0 until dim) {
        var sum: Interval = Interval(0.0)
        for (j <- 0 until dim)
          sum = sum + Interval(data(i)(j)) * v(j)
        result(i) = sum
      }
      result
    }

    def *(m: AffineMatrix): AffineMatrix = {
      val result = new AffineMatrix(dim)
      for (i <- 0 until dim) {
        for (j <- 0 until dim) {
          var sum: AffineForm = new AForm(0.0)
          for (k <- 0 until dim)
            sum = sum + AffineForm(data(i)(k)) * m.data(k)(j)
          result(i, j) = sum
        }
      }
      result
    }

    def *(m: IntervalMatrix): IntervalMatrix = {
      val result = new IntervalMatrix(dim)
      for (i <- 0 until dim) {
        for (j <- 0 until dim) {
          var sum: Interval = Interval(0.0)
          for (k <- 0 until dim)
            sum = sum + Interval(data(i)(k)) * m.data(k)(j)
          result(i, j) = sum
        }
      }
      result
    }


    def -(m: AffineMatrix): AffineMatrix = {
      val result = new AffineMatrix(dim)
      for (i <- 0 until dim)
        for (j <- 0 until dim)
          result(i, j) = AffineForm(data(i)(j)) - m(i, j)
      result
    }
  }

  class AffineMatrix(d: Int) extends Matrix[AffineForm](d) {
    def getNewArray = Array.fill(dim, dim)(new AForm(0.0))
    override def toString: String = {
      var str = ""
      for (i <- 0 until dim)
        str = str + (new AffineVector(d, data(i))).toString + "\n"
      str
    }
    def mid: DoubleMatrix = {
      val result = new DoubleMatrix(dim)
      for (i <- 0 until dim) {
        for (j <- 0 until dim) {
          result(i,j) = data(i)(j).interval.mid
        }
      }
      result
    }
    def *(v: AffineVector): AffineVector = {
      val result = new AffineVector(dim)
      for (i <- 0 until dim) {
        var sum: AffineForm = new AForm(0.0)
        for (j <- 0 until dim)
          sum = sum + data(i)(j) * v(j)
        result(i) = sum
      }
      result
    }
    def -(m: AffineMatrix): AffineMatrix = {
      val result = new AffineMatrix(dim)
      for (i <- 0 until dim)
        for (j <- 0 until dim)
          result(i, j) = data(i)(j) - m(i, j)
      result
    }

    def toInterval: IntervalMatrix = {
      val result = new IntervalMatrix(dim)
      for (i <- 0 until dim)
        for (j <- 0 until dim)
          result(i, j) = data(i)(j).interval
      result
    }
  }


  class IntervalMatrix(d: Int) extends Matrix[Interval](d) {
    def getNewArray = Array.fill(dim, dim)(Interval(0.0))
    override def toString: String = {
      var str = ""
      for (i <- 0 until dim)
        str = str + (new IntervalVector(d, data(i))).toString + "\n"
      str
    }
    def mid: DoubleMatrix = {
      val result = new DoubleMatrix(dim)
      for (i <- 0 until dim) {
        for (j <- 0 until dim) {
          result(i,j) = data(i)(j).mid
        }
      }
      result
    }
    def *(v: IntervalVector): IntervalVector = {
      val result = new IntervalVector(dim)
      for (i <- 0 until dim) {
        var sum: Interval = Interval(0.0)
        for (j <- 0 until dim)
          sum = sum + data(i)(j) * v(j)
        result(i) = sum
      }
      result
    }
    def -(m: IntervalMatrix): IntervalMatrix = {
      val result = new IntervalMatrix(dim)
      for (i <- 0 until dim)
        for (j <- 0 until dim)
          result(i, j) = data(i)(j) - m(i, j)
      result
    }

  }


  class Fnc2Matrix(d: Int) extends Matrix[(Double, Double) => Double](d) {
    def getNewArray = Array.fill(dim, dim)((x: Double, y: Double) => 0.0)
    def eval(x: DoubleVector): DoubleMatrix = {
      assert(dim == x.dim)
      val result = new DoubleMatrix(dim)
      for (i <- 0 until dim) {
        for (j <- 0 until dim) {
          result(i,j) = data(i)(j)(x(0), x(1))
        }
      }
      result
    }
  }

  class Fnc2AAMatrix(d: Int) extends Matrix[(AffineForm, AffineForm) => AffineForm](d) {
    def getNewArray = Array.fill(dim, dim)((x: AffineForm, y: AffineForm) => new AForm(0.0))
    def eval(x: AffineVector): AffineMatrix = {
      assert(dim == x.dim)
      val result = new AffineMatrix(dim)
      for (i <- 0 until dim) {
        for (j <- 0 until dim) {
          result(i,j) = data(i)(j)(x(0), x(1))
        }
      }
      result
    }

  }


  def invert(A: Array[Array[Double]]): Array[Array[Double]] = {
    val lu = copyMatrix(A)
    val n = A.length
    val pivot = Array.fill(n){0}
    factor(lu, pivot)  // computed LU factorization

    val inverse = Array.fill(n, n){0.0}


    for (i <- 0 until n) {
      val e = Array.fill(n){0.0}; e(i) = 1.0  // unit vector

      solve(lu, pivot, e) // get the i-th column

      //copy e into inverse
      for (j <- 0 until n) inverse(j)(i) = e(j)
    }
    inverse
  }


  def solveSystem(A: Array[Array[Double]], b: Array[Double]): Array[Double] = {
    val lu = copyMatrix(A)
		val pivot: Array[Int] = Array.fill(b.length){0}
    factor(lu, pivot)
    solve(lu, pivot, b)
    return b
  }

  def copyMatrix(data: Array[Array[Double]]): Array[Array[Double]] = {
    val newData = Array.fill(data.length, data(0).length){0.0}
    for (i <- 0 until data.length)
      for (j <- 0 until data(0).length)
        newData(i)(j) = data(i)(j)
	  newData
	}

  /**
    Solve a linear system, using a prefactored matrix in LU form.
    (Code from Scimark).

    @param LU (in) the factored matrix in LU form.
    @param pivot (in) the pivot vector which lists the reordering used
      during the factorization stage.
    @param b    (in/out) On input, the right-hand side. On output, the solution vector.
    */
  def solve(LU: Array[Array[Double]], pvt: Array[Int], b: Array[Double]) = {
    val M = LU.length
    val N = LU(0).length
    var ii=0

    for (i <- 0 until M) {
      ///*
      val ip = pvt(i)
      var sum = b(ip)

      b(ip) = b(i)
     // */ var sum = b(i)
      if (ii==0) {
        for (j <- ii until i)
          sum -= LU(i)(j) * b(j)
      }
      else {
        if (sum == 0.0)
          ii = i
      }
      b(i) = sum
    }

    var i = N-1
    while (i >= 0) {
      var sum = b(i)
      for (j <- (i+1) until N) {
        sum -= LU(i)(j) * b(j)
      }
      b(i) = sum / LU(i)(i)
      i -= 1
    }
  }


  /**
    LU factorization (in place).  (from Scimark)
    @param A (in/out) On input, the matrix to be factored.
        On output, the compact LU factorization.
    @param pivit (out) The pivot vector records the
        reordering of the rows of A during factorization.
    @return 0, if OK, nozero value, othewise.
*/
  def factor(A: Array[Array[Double]],  pivot: Array[Int]): Int = {
    val N = A.length
    val M = A(0).length

    val minMN = math.min(M,N)

    for (j <- 0 until minMN) {
      // find pivot in column j and  test for singularity.
      var jp = j
      var t = math.abs(A(j)(j))
      for (i <- (j+1) until M) {
          val ab = math.abs(A(i)(j))
          if ( ab > t) {
              jp = i
              t = ab
          }
      }
      pivot(j) = jp

      // jp now has the index of maximum element
      // of column j, below the diagonal

      if ( A(jp)(j) == 0 )
          return 1       // factorization failed because of zero pivot

      if (jp != j) {
          // swap rows j and jp
          val tA = A(j)
          A(j) = A(jp)
          A(jp) = tA
      }
      if (j < M-1) {                // compute elements j+1:M of jth column
          // note A(j,j), was A(jp,p) previously which was
          // guarranteed not to be zero (Label #1)
          //
          val recp =  1.0 / A(j)(j)
          for (k <- (j+1) until M)
              A(k)(j) *= recp
      }

      if (j < minMN-1) {
        // rank-1 update to trailing submatrix:   E = E - x*y;
        //
        // E is the region A(j+1:M, j+1:N)
        // x is the column vector A(j+1:M,j)
        // y is row vector A(j,j+1:N)
        for (ii <- (j+1) until M) {
            val Aii = A(ii)
            val Aj = A(j)
            val AiiJ = Aii(j)
            for (jj <- (j+1) until N)
              Aii(jj) -= AiiJ * Aj(jj)
        }
      }
    }
    return 0
  }


  def solveSystemQD(A: Array[Array[QD]], b: Array[QD]): Array[QD] = {
    val lu = copyMatrixQD(A)
		val pivot: Array[Int] = Array.fill(b.length){0}
    factorQD(lu, pivot)
    solveQD(lu, pivot, b)
    return b
  }

  def copyMatrixQD(data: Array[Array[QD]]): Array[Array[QD]] = {
    val newData = Array.fill(data.length, data(0).length){QD(0.0)}
    for (i <- 0 until data.length)
      for (j <- 0 until data(0).length)
        newData(i)(j) = data(i)(j)
	  newData
	}

  /**
    Solve a linear system, using a prefactored matrix in LU form.
    (Code from Scimark).

    @param LU (in) the factored matrix in LU form.
    @param pivot (in) the pivot vector which lists the reordering used
      during the factorization stage.
    @param b    (in/out) On input, the right-hand side. On output, the solution vector.
    */
  def solveQD(LU: Array[Array[QD]], pvt: Array[Int], b: Array[QD]) = {
    val M = LU.length
    val N = LU(0).length
    var ii=0

    for (i <- 0 until M) {
      ///*
      val ip = pvt(i)
      var sum = b(ip)

      b(ip) = b(i)
     // */ var sum = b(i)
      if (ii==0) {
        for (j <- ii until i)
          sum -= LU(i)(j) * b(j)
      }
      else {
        if (sum == 0.0)
          ii = i
      }
      b(i) = sum
    }

    var i = N-1
    while (i >= 0) {
      var sum = b(i)
      for (j <- (i+1) until N) {
        sum -= LU(i)(j) * b(j)
      }
      b(i) = sum / LU(i)(i)
      i -= 1
    }
  }


  /**
    LU factorization (in place).  (from Scimark)
    @param A (in/out) On input, the matrix to be factored.
        On output, the compact LU factorization.
    @param pivit (out) The pivot vector records the
        reordering of the rows of A during factorization.
    @return 0, if OK, nozero value, othewise.
*/
  def factorQD(A: Array[Array[QD]],  pivot: Array[Int]): Int = {
    val N = A.length
    val M = A(0).length

    val minMN = math.min(M,N)

    for (j <- 0 until minMN) {
      // find pivot in column j and  test for singularity.
      var jp = j
      var t = QD.abs(A(j)(j))
      for (i <- (j+1) until M) {
          val ab = QD.abs(A(i)(j))
          if ( ab > t) {
              jp = i
              t = ab
          }
      }
      pivot(j) = jp

      // jp now has the index of maximum element
      // of column j, below the diagonal

      if ( A(jp)(j) == 0 )
          return 1       // factorization failed because of zero pivot

      if (jp != j) {
          // swap rows j and jp
          val tA = A(j)
          A(j) = A(jp)
          A(jp) = tA
      }
      if (j < M-1) {                // compute elements j+1:M of jth column
          // note A(j,j), was A(jp,p) previously which was
          // guarranteed not to be zero (Label #1)
          //
          val recp =  1.0 / A(j)(j)
          for (k <- (j+1) until M)
              A(k)(j) *= recp
      }

      if (j < minMN-1) {
        // rank-1 update to trailing submatrix:   E = E - x*y;
        //
        // E is the region A(j+1:M, j+1:N)
        // x is the column vector A(j+1:M,j)
        // y is row vector A(j,j+1:N)
        for (ii <- (j+1) until M) {
            val Aii = A(ii)
            val Aj = A(j)
            val AiiJ = Aii(j)
            for (jj <- (j+1) until N)
              Aii(jj) -= AiiJ * Aj(jj)
        }
      }
    }
    return 0
  }


  def evalVector(v: Array[(Double, Double, Double) => Double], x: Array[Double]): DoubleVector = {
    val result = new Array[Double](v.length)

    for (i <- 0 until v.length) {
      result(i) = v(i)(x(0), x(1), x(2))
    }
    new DoubleVector(3, result)
  }

  def evalMatrix(matrix: Array[Array[(Double, Double, Double) => Double]], x: Array[Double]): DoubleMatrix = {
    val dim = 3
    val result = new Array[Array[Double]](dim)

    for (i <- 0 until dim) {
      result(i) = new Array[Double](dim)
      for (j <- 0 until dim) {
        result(i)(j) = matrix(i)(j)(x(0), x(1), x(2))
      }
    }
    new DoubleMatrix(result)
  }

  def evaluateVector(v: Array[(Double, Double) => Double], x: Array[Double]): Array[Double] = {
    val result = new Array[Double](v.length)
    for (i <- 0 until v.length) {
      result(i) = v(i)(x(0), x(1))
    }
    result
  }

  def evaluateVector3(v: Array[(Double, Double, Double) => Double], x: Array[Double]): DoubleVector = {
    val result = new Array[Double](v.length)
    for (i <- 0 until v.length) {
      result(i) = v(i)(x(0), x(1), x(2))
    }
    new DoubleVector(3, result)
  }

  def evaluateMatrix(matrix: Array[Array[(Double, Double) => Double]], x: Array[Double]): Array[Array[Double]] = {
    val n = matrix.length; val m = matrix(0).length
    val result = new Array[Array[Double]](n)

    for (i <- 0 until n) {
      result(i) = new Array[Double](m)
      for (j <- 0 until m) {
        result(i)(j) = matrix(i)(j)(x(0), x(1))
      }
    }
    result
  }

  def evaluateMatrix(matrix: Array[Array[(Double, Double, Double) => Double]], x: Array[Double]): DoubleMatrix = {
    val dim = 3
    val result = new Array[Array[Double]](dim)

    for (i <- 0 until dim) {
      result(i) = new Array[Double](dim)
      for (j <- 0 until dim) {
        result(i)(j) = matrix(i)(j)(x(0), x(1), x(2))
      }
    }
    new DoubleMatrix(result)
  }

   def addVectors(x: Array[Double], y: Array[Double]): Array[Double] = {
    val res = new Array[Double](x.length)
    for (i <- 0 until x.length)
      res(i) = x(i) + y(i)
    res
  }

  def negateVector(x: Array[Double]): Array[Double] = {
    val res = new Array[Double](x.length)
    for (i <- 0 until x.length)
      res(i) = -x(i)
    res
  }

  // Euclidean norm
  def norm(v: Array[Double]): Double = {
    var sum = 0.0
    for (i <- 0 until v.length) {
      sum += v(i)*v(i)
    }
    sqrt(sum)
  }

  def evaluateVectorQD(v: Array[(QD, QD) => QD], x: Array[QD]): Array[QD] = {
    val result = new Array[QD](v.length)
    for (i <- 0 until v.length) {
      result(i) = v(i)(x(0), x(1))
    }
    result
  }

  def evaluateMatrixQD(matrix: Array[Array[(QD, QD) => QD]], x: Array[QD]): Array[Array[QD]] = {
    val n = matrix.length; val m = matrix(0).length
    val result = new Array[Array[QD]](n)

    for (i <- 0 until n) {
      result(i) = new Array[QD](m)
      for (j <- 0 until m) {
        result(i)(j) = matrix(i)(j)(x(0), x(1))
      }
    }
    result
  }

  def evaluateVectorQD(v: Array[(QD, QD, QD) => QD], x: Array[QD]): Array[QD] = {
    val result = new Array[QD](v.length)
    for (i <- 0 until v.length) {
      result(i) = v(i)(x(0), x(1), x(2))
    }
    result
  }

  def evaluateMatrixQD(matrix: Array[Array[(QD, QD, QD) => QD]], x: Array[QD]): Array[Array[QD]] = {
    val n = matrix.length; val m = matrix(0).length
    val result = new Array[Array[QD]](n)

    for (i <- 0 until n) {
      result(i) = new Array[QD](m)
      for (j <- 0 until m) {
        result(i)(j) = matrix(i)(j)(x(0), x(1), x(2))
      }
    }
    result
  }


  def addVectorsQD(x: Array[QD], y: Array[QD]): Array[QD] = {
    val res = new Array[QD](x.length)
    for (i <- 0 until x.length)
      res(i) = x(i) + y(i)
    res
  }

  def negateVectorQD(x: Array[QD]): Array[QD] = {
    val res = new Array[QD](x.length)
    for (i <- 0 until x.length)
      res(i) = -x(i)
    res
  }

  // Euclidean norm
  def normQD(v: Array[QD]): QD = {
    var sum = QD(0.0)
    for (i <- 0 until v.length) {
      sum += v(i)*v(i)
    }
    QD.sqrt(sum)
  }


}
