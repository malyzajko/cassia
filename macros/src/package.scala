/*
 * Cerres (c) 2012-2014 EPFL, Lausanne
 */
 
package cerres

import language.implicitConversions

import ceres.common._

package object macros {

  case class SolutionCannotBeCertified(s: String) extends Exception
  case class SolutionNotIncluded(s: String) extends Exception

  var debug = true

  class DoubleWrapper(d: Double) {
    def +/-(e: Double) = {
      if (e > 0) Interval(d - e , d + e)
      else Interval(d + e, d - e)
    }
  }
  implicit def double2DoubleWrapper(d: Double) = new DoubleWrapper(d)

  // Determines whether we use the interval or the affine version
  val useIntervalsForUnary = true
  val useIntervalsForMulti = true

}
