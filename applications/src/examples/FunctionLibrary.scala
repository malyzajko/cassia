/*
 * Cerres (c) 2012-2014 EPFL, Lausanne
 */

package cerres.examples

import ceres.common.{QuadDouble => QD}
import ceres.affine._
//import QD._

object FunctionLibrary {
  /* Naming convention: each function name looks like name_type[_derivative]
    name can be something like p1 for polynomial 1, type can be d for Double,
    dd for DoubleDouble, qd for QuadDouble, aa for affine form,
    derivative is optional and denotes if a function is the nth derivative of another one.
  */

  // Polynomial 1 has roots: -4, -2, -1, 1, 3
  // Computed with doubles, 1e-5: -4.000000692519392, -2.0000000033755176, -0.9999999989838327, 1.0000003092609135, 3.0000003221359677
  // Computed with doubles, 1e-10: -4.000000000000564, -2.000000000000001, -1.0, 1.000000000000051, 3.0000000000001132
  val p1_d = (x: Double) => (x*x*x*x*x + 3 * x*x*x*x - 11 * x*x*x - 27*x*x +10*x + 24)/20.0
  val p1_d_1 = (x: Double) => (5*x*x*x*x + 12*x*x*x - 33*x*x - 54*x + 10)/20.0

  val p1_qd = (x: QD) => (x*x*x*x*x + 3 * x*x*x*x - 11 * x*x*x - 27*x*x +10*x + 24)/20.0
  val p1_qd_1 = (x: QD) => (5*x*x*x*x + 12*x*x*x - 33*x*x - 54*x + 10)/20.0

  val p1_aa = (x: AffineForm) => (x*x*x*x*x + 3 * x*x*x*x - 11 * x*x*x - 27*x*x +10*x + 24)/20.0
  val p1_aa_1 = (x: AffineForm) => (5*x*x*x*x + 12*x*x*x - 33*x*x - 54*x + 10)/20.0

  // Polynomial 2 has one root at: -1.3247179572447460259609089
  // Computed with doubles, 1e-5 : -1.3247190494171253
  // Computed with doubles, 1e-10: -1.3247179572458576
  val p2_d = (x: Double) => x*x*x - x + 1
  val p2_d_1 = (x: Double) => 3 * x*x - 1
  val p2_qd = (x: QD) => x*x*x - x + 1
  val p2_qd_1 = (x: QD) => 3 * x*x - 1
  val p2_aa = (x: AffineForm) => x*x*x - x + 1
  val p2_aa_1 = (x: AffineForm) => 3 * x*x - 1

  // Polynomial 3 has one root at: 0.4501836112948735730365387
  // Computed with doubles, 1e-5 : 0.4501836542045203
  // Computed with doubles, 1e-10: 0.4501836112948739
  val p3_d = (x: Double) => math.cos(x) - 2 * x
  val p3_d_1 = (x: Double) => - math.sin(x) - 2

  val p3_qd = (x: QD) => QD.cos(x) - 2 * x
  val p3_qd_1 = (x: QD) => - QD.sin(x) - 2

  val p3_aa = (x: AffineForm) => AffineForm.cos(x) - 2 * x
  val p3_aa_1 = (x: AffineForm) => - AffineForm.sin(x) - 2


  // Polynomial 4 has one root at: 3.3333333333333334566914472
  // Computed with doubles, 1e-5 : 3.333332841403648
  // Computed with doubles, 1e-10: 3.333333333333261
  val p4_d = (x: Double) => 1/x - 0.3
  val p4_d_1 = (x: Double) => -1.0/(x*x)

  val p4_qd = (x: QD) => 1/x - 0.3
  val p4_qd_1 = (x: QD) => -1.0/(x*x)

  val p4_aa = (x: AffineForm) => 1/x - 0.3
  val p4_aa_1 = (x: AffineForm) => -1.0/(x*x)

  // Polynomial 5 has two roots: For some reason this diverges?! 
  // Computed with doubles, 1e-5 : 
  // Computed with doubles, 1e-10:
  val p5_d = (x: Double) => math.exp(x) - 2*x - 0.1
  val p5_d_1 = (x: Double) => math.exp(x) - 2

  val p5_qd = (x: QD) => QD.exp(x) - 2*x - 0.1
  val p5_qd_1 = (x: QD) => QD.exp(x) - 2

  val p5_aa = (x: AffineForm) => AffineForm.exp(x) - 2*x - 0.1
  val p5_aa_1 = (x: AffineForm) => AffineForm.exp(x) - 2

  // Function 6 has one root at 0.0  (diverges with Newton, converges with Secant) 
  // Computed with doubles, 1e-5 : -1.6999861111258868E-8
  // Computed with doubles, 1e-10: 9.926167350636332E-24
  // Aaaha: Newton gives us: -5.1020408163265305 (for both tolerances) as a possible root, starting from -5.0
  val p6_d = (x: Double) => x * math.exp(-(x*x))
  val p6_d_1 = (x: Double) => math.exp(-(x*x)) - 2*x*x*math.exp(-(x*x))

  val p6_qd = (x: QD) => x * QD.exp(-(x*x))
  val p6_qd_1 = (x: QD) => QD.exp(-(x*x)) - 2*x*x*QD.exp(-(x*x))

  val p6_aa = (x: AffineForm) => x * AffineForm.exp(-(x*x))
  val p6_aa_1 = (x: AffineForm) => AffineForm.exp(-(x*x)) - 2*x*x*AffineForm.exp(-(x*x))

  // Polynomial 7 has 3 roots at 0.75, 1.0 and 2.0
  // Computed with doubles, 1e-5 : 0.7499972689032873, 1.000000000940891, 2.0000000004806338
  // Computed with doubles, 1e-10: 0.7499999999641983, 0.9999999999999973, 2.0000000000000013
  val p7_d = (x: Double) => 4*x*x*x - 15*x*x + 17*x - 6
  val p7_d_1 = (x: Double) => 12*x*x - 30*x + 17

  val p7_qd = (x: QD) => 4*x*x*x - 15*x*x + 17*x - 6
  val p7_qd_1 = (x: QD) => 12*x*x - 30*x + 17

  val p7_aa = (x: AffineForm) => 4*x*x*x - 15*x*x + 17*x - 6
  val p7_aa_1 = (x: AffineForm) => 12*x*x - 30*x + 17

  // Function 8 has 2 roots in [-3.0, Infnt] at -1.3812060997502068247999627, 0.2548500590288502553077938 
  // Computed with doubles, 1e-5 : -1.3812065650950054, 0.2548500590506985
  // Computed with doubles, 1e-10: -1.3812060997502582, 0.25485005902885033
  val p8_d = (x: Double) => 3 * math.exp(x) - 4 * math.cos(x)
  val p8_d_1 = (x: Double) => 3 * math.exp(x)  + 4 * math.sin(x)

  val p8_qd = (x: QD) => 3 * QD.exp(x) - 4 * QD.cos(x)
  val p8_qd_1 = (x: QD) => 3 * QD.exp(x)  + 4 * QD.sin(x)

  val p8_aa = (x: AffineForm) => 3 * AffineForm.exp(x) - 4 * AffineForm.cos(x)
  val p8_aa_1 = (x: AffineForm) => 3 * AffineForm.exp(x)  + 4 * AffineForm.sin(x)

  // Polynomial 9 has one root at 0.1999999999999998223643161
  // Computed with doubles, 1e-5 : 0.199609375
  // Computed with doubles, 1e-10: 0.19999847412109376
  val p9_d = (x: Double) => 25*x*x -10*x + 1
  val p9_d_1 = (x: Double) => 50*x - 10

  val p9_qd = (x: QD) => 25*x*x -10*x + 1
  val p9_qd_1 = (x: QD) => 50*x - 10

  val p9_aa = (x: AffineForm) => 25*x*x -10*x + 1
  val p9_aa_1 = (x: AffineForm) => 50*x - 10

  // Polynomial 10 has one root at 0.0
  // Computed with doubles, 1e-5 : 1.9635306152277626E-10
  // Computed with doubles, 1e-10: 0.0
  val p10_d = (x: Double) => math.atan(x)
  val p10_d_1 = (x: Double) => 1.0/(1 + x*x)

  val p10_qd = (x: QD) => QD.atan(x)
  val p10_qd_1 = (x: QD) => 1.0/(1 + x*x)

  val p10_aa = (x: AffineForm) => AffineForm.atan(x)
  val p10_aa_1 = (x: AffineForm) => 1.0/(1 + x*x)

  // Function 11 has 2 roots at -1.2710268008159460640047188, -0.6592660457669460745373486
  // Computed with doubles, 1e-5 : -1.2710268031402283, -0.6592654415665319
  // Computed with doubles, 1e-10: -1.271026800815946, -0.6592660457664381
  val p11_d = (x: Double) => math.cos(x) + 2*math.sin(x) + x*x
  val p11_d_1 = (x: Double) => -math.sin(x) + 2*math.cos(x) + 2*x

  val p11_qd = (x: QD) => QD.cos(x) + 2*QD.sin(x) + x*x
  val p11_qd_1 = (x: QD) => -QD.sin(x) + 2*QD.cos(x) + 2*x

  val p11_aa = (x: AffineForm) => AffineForm.cos(x) + 2*AffineForm.sin(x) + x*x
  val p11_aa_1 = (x: AffineForm) => -AffineForm.sin(x) + 2*AffineForm.cos(x) + 2*x

  // Function 12 has 2 roots at 1.4296118247255556122752444, 8.6131694564413985966763966
  // Computed with doubles, 1e-5 : 1.4296077757568797, 8.613169478615234
  // Computed with doubles, 1e-10: 1.4296118247166327, 8.6131694564414
  val p12_d = (x: Double) => 4 * math.log(x) - x
  val p12_d_1 = (x: Double) => 4.0/x - 1

  val p12_qd = (x: QD) => 4 * QD.log(x) - x
  val p12_qd_1 = (x: QD) => 4.0/x - 1

  val p12_aa = (x: AffineForm) => 4 * AffineForm.log(x) - x
  val p12_aa_1 = (x: AffineForm) => 4.0/x - 1

}
