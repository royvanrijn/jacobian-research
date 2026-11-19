# Odd-cubic degree-seven safety and the sharp odd threshold nine

Date: 19 November 2026. The statements concern precisely regular winding-one odd-cubic-height phases; the sharp odd threshold is restricted to this class.

## 1. Statement and inputs

Work over any characteristic-zero field K with the regular phase

    p=t/x*(C(x)+t^2D(x)), C,D in K[x], C(0)+D(0)=0.

**Theorem OddCubicDegreeSevenSafety.** Every pure phase in this class whose least degree is odd and at most7 is eventually safe for every fixed polynomial marker. There is no pure unit-center phase of exact least degree7.

**Corollary SharpOddCubicOddThreshold.** The least odd degree of an unsafe phase in this entire class is9. The lower bound follows from the theorem and retained degree<=5 classification; the upper bound is the ALREADY RETAINED F9 witness ([ODD_PRIMITIVE_TWO_PAIR_IMAGE_WITNESSES.md](ODD_PRIMITIVE_TWO_PAIR_IMAGE_WITNESSES.md), Section1), whose dehomogenization has C=(1+x)^3,D=-(1+x)^7. This is not a new degree9 construction and not a global minimum for arbitrary Image phases.

The zero-center case is safe in all degrees by [ODD_START_HEIGHT_RATIONALITY_AND_QUINTIC_SAFETY.md](ODD_START_HEIGHT_RATIONALITY_AND_QUINTIC_SAFETY.md). The degree<=5 odd-cubic case is already retained. It remains to exclude a unit-center pure phase of exact degree7, normalizing C(0)=-D(0)=1.

## 2. Pole-degree control from exact/log separation

Retained rational-u descent gives

    u=x*sqrt(-D/C)/C in K(x), u=x+O(x^2),
    u^2=-D*x^2/C^3.

Over k=Kbar write dx/C=dR+Omega, with R rational, normalized R(0)=0, and Omega its logarithmic part. If C has degree c>=2 and h distinct roots of multiplicities m_i, R has precisely pole order m_i-1 at every repeated root and none elsewhere, including infinity. Thus its exact rational-map degree is

    D_R=c-h.                                          (S1)

Pure squarefree C is excluded, so D_R>0. Exact/log separation says

    R(v(u))-R(v(-u))=2u.

As a result E=k(u,R) is stable under the formal involution sigma:u->-u,R->R-2u. Put R_sigma=R-2u. Crucially, even if sigma does not extend to k(x),

    deg_x R_sigma=[k(x):E]*[E:k(R_sigma)]
                 =[k(x):E]*[E:k(R)]=D_R.              (S2)

The middle equality follows by applying the field automorphism sigma. Therefore

    deg_x u<=2D_R.                                    (S3)

If equality holds, R and R_sigma have disjoint pole supports: a shared pole of positive orders r,s contributes at most max(r,s)<r+s to their difference.

If u vanishes at infinity with order s>0, then its simple zero at x=0 gives

    deg_x u>=s+1.                                     (S4)

All these are exact map-degree and divisor statements, not finite-prefix tests.

## 3. Auxiliary root exclusions

One geometric root of C is completely classified in [ONE_POLE_ODD_CUBIC_CLASSIFICATION.md](ONE_POLE_ODD_CUBIC_CLASSIFICATION.md): if c>=2 its possible pure degrees are3c or5c-2. Thus none of c=3,5,7 can give degree7.

Exactly two geometric roots of C are impossible for ANY pure unit-center odd-cubic phase. Indeed degC>=2 makes the two finite residues sum0, and both cannot vanish because then1/C would have a rational primitive, contrary to the retained rational-primitive root criterion for two roots. Its logarithmic part is therefore a*dlog((x-r1)/(x-r2)), a!=0. Exact/log separation makes the ratio h(v(u))/h(v(-u)), h=(x-r1)/(x-r2), have derivative0 and constant term1. Thus it equals1 formally. Mobius injectivity forces v(u)=v(-u), contradicting u=-u. This proves the needed two-root exclusion without importing a separate two-root classification.

Squarefree nonconstant C is excluded by the retained SquarefreeCenterObstruction, or directly because R=0 contradicts the exact-part identity.

## 4. Complete exact-degree-seven case split

Put c=degC,e=degD. Rationality of sqrt(-D/C) forces c,e to have the same parity. The exact least sphere degree is max(c,e+2). If it equals7, the possibilities are

    c=1,3,5 with e=5, or c=7 with e=1,3,5.

The case c=1 is squarefree and excluded. For c=3, every root-count pattern is one-pole, two-root or squarefree, hence excluded by Section3.

For c=5,e=5, the infinity order is

    ord_zero_infinity u=(3c-e-2)/2=4.

The one-pole and two-root cases are excluded already. If h>=3 then D_R<=2, so(S3) gives degu<=4, while(S4) gives degu>=5, a contradiction.

For c=7 and e<=5 odd, u vanishes at infinity with order

    s=(19-e)/2>=7.

Again h=1,2 are excluded. If h>=4 then D_R<=3 and degu<=6, contradicting degu>=8. The only remaining root count is h=3. Here D_R=4. Bounds(S3)--(S4) force e=5,s=7 and degu=8=2D_R. Thus R and R_sigma have disjoint poles.

Every repeated C-root is a pole of R. The rational function R_sigma has no pole at infinity (both R and u are regular there), and no finite pole outside the roots of C (by u^2=-Dx^2/C^3). Disjointness therefore confines all poles of R_sigma to SIMPLE C-roots. At a simple root r of C, r!=0 and D is polynomial, so the pole order k of u satisfies2k<=3, hence k<=1. Since R is regular there, the pole order of R_sigma is also at most1. With three distinct roots and degree7, at least one root is repeated, so there are at most two simple roots. Consequently

    deg R_sigma<=2,

contradicting deg R_sigma=4 from(S2). This exhausts every degree-seven case and proves the theorem.

## 5. Scope and review limits

The safety conclusion for arbitrary original markers uses the repaired retained transfer: E2 is equivariant and kills strictly positive input weights; the negative safe branch requires fixed marker bidegree and its bounded output-weight interval. It does not use the false assertion that E2 kills every nonzero input weight.

The sharp odd threshold9 applies only to this full odd-cubic-height winding-one class. The larger height-five profile t/x*(C+t^2B+t^4D) and arbitrary quintic/septic phases remain outside this proof. No all-order conclusion is inferred from bounded moment calculations, and no historical-priority or global minimum-prime assertion is made.


The [source ledger](odd_cubic_degree_seven_safety_source_ledger.json) pins all five retained inputs, including the repaired original-marker safety transfer and the existing all-order F9 sharpness witness. The [checker](verify_odd_cubic_degree_seven_safety.py) and [receipt](odd_cubic_degree_seven_safety_verification.json) bind the exact finite algebraic controls to these proof and source bytes.

```bash
python3 research/foundations/extended-geometry/verify_odd_cubic_degree_seven_safety.py \
  --write-receipt research/foundations/extended-geometry/odd_cubic_degree_seven_safety_verification.json
```
