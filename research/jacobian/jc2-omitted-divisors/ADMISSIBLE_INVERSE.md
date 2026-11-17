# A source finite inverse equation records every omitted divisor

14 November 2026. Work over C. PF2OMIT3 gives an elimination route to the
complete omission ledger for a plane polynomial map. The extra hypotheses
on the inverse coordinate are proved available; an arbitrary resolvent
does not have this property. No unconstructed F2 map is supplied here.

## 1. Choose a primitive source coordinate with a finiteness certificate

Let F=(F1,F2):A2 -> A2 be dominant and quasi-finite, A=C[u,v], and
d=[C(x,y):C(F1,F2)]. Choose c in C and put T=x+c*y, beta=y. Require:

1. F1(T-c*beta,beta) has a nonzero constant leading coefficient in beta.
2. C(F1,F2,T)=C(x,y).

Both requirements exclude only finitely many choices of c. For the first,
the leading coefficient is the highest homogeneous part of F1 evaluated
at (-c,1). For the second, distinct embeddings of the finite separable
source field agree on x+c*y for at most one c, unless their y-images
coincide, in which case their distinct x-images never agree. Thus a
single generic constant satisfies both requirements.

After division by a nonzero constant, the equation F1(T-c*beta,beta)=u
is monic in beta. Therefore C[x,y] is finite over A[T]. Let

\[
 R(T;u,v)=a_d(u,v)T^d+\cdots+a_0(u,v)                     \tag{I1}
\]

be the primitive irreducible equation of T, normalized up to a nonzero
constant. The map U=A2 -> X_aff=Spec(A[T]/(R)) is finite birational;
its normalization is precisely U. In fact C[x,y] is finite and normal
over A[T], and both rings have the same fraction field.

The coefficient ideal is the unit ideal:

\[
 (a_0,\ldots,a_d)=A.                                     \tag{I2}
\]

If every coefficient vanished at a closed target point, the fibre of
X_aff there would be A1. The finite surjective normalization U -> X_aff
would then give a positive-dimensional fibre of F, contradicting
quasi-finiteness. Nullstellensatz proves (I2). This is stronger than
requiring only that the coefficients have gcd one.

The quasifinite hypothesis matters: the birational map (x,y)->(x,xy)
is finite over its image together with T=y, but its equation uT−v has
coefficient ideal (u,v). At zero the whole projective T-line appears,
and the source fibre is positive dimensional. This equation is rejected
by (I2), even though its coefficients have gcd one.

## 2. The entire omitted boundary is the normalized infinity divisor

Homogenize (I1) to degree d in [T:Z] and let X be its zero locus in
A2 times P1. By (I2) every fibre is a nonzero degree-d binary divisor.
Consequently X -> A2 is finite flat of degree d. Its finite normalization
N is the canonical normalization of the target in C(x,y). Normalization
commutes with restriction to the affine T-chart, giving

\[
 \boxed{U=\nu^{-1}(Z\ne0),\qquad
        N\setminus U=\nu^{-1}(Z=0).}                     \tag{I3}
\]

This identifies every omitted prime and every omitted special point.
There are no other deleted branches. The normal finite surface N is flat
over the smooth target, as in the earlier [normalization proof](PROOF.md).
Every omitted prime is a pole of T. Conversely a pole of T is omitted.
Thus a source-finite primitive coordinate reduces the simultaneous
reconstruction-pole test to a single projective infinity divisor.

The target images of all omitted primes are exactly

\[
 \boxed{S_F=V(a_d).}                                    \tag{I4}
\]

This assertion concerns omission. Ramification must still be computed on
the normalization and can differ from the discriminant of R.

The construction also gives a global source-coordinate fibration
T:N -> P1. Every finite T-fibre is exactly an affine source line A1,
and its infinity fibre as a Cartier divisor is sum_D p_D*D, with all the
omitted primes and their pole multiplicities. For a nonproper map this
fibration reaches infinity. It is a source-coordinate fibration, not the
target-line or radial pencil: it supplies no rationality statement for
the generic curve in those different pencils.

## 3. Exceptional target fibres are exact degree strata

For a closed target value y, put r_y=deg_T R(T;y). Equation (I2) ensures
R(T;y) is nonzero, allowing degree zero. Its projective fibre has finite
length r_y and infinity length d-r_y. These aggregate lengths are
preserved by the finite birational normalization N -> X.

For completeness, pass to the completed local target ring at y. The
finite flat algebra of X splits into the factors supported at its closed
fibre points. Pulling back the corresponding idempotents splits the
finite flat algebra of N. Each factor has the same generic rank before
and after normalization: both have the same generic function field.
Hence their closed fibre lengths agree. Add over the finite T-points
and, separately, over the infinity point. This proves

\[
 \boxed{\operatorname{length}(U_y)=r_y,\qquad
 L_{\rm omitted}(y)=d-r_y.}                             \tag{I5}
\]

For a Keller map every source fibre is reduced, so r_y is also its
cardinality. Distinct source points can nevertheless have the same
T-value: one must count the specialized polynomial's degree, not its
distinct roots. Normalization remains necessary for the individual
omitted local algebras and their distribution among primes.

There is also a direct source-algebra proof of the retained length. Let
M=deg_beta(F1), after the monic change in Section1. At a fixed target
(u0,v0), the algebra C[T,beta]/(F1−u0) is free of rank M over C[T].
Multiplication by F2−v0 has cokernel the source fibre algebra and
determinant Res_beta(F1−u0,F2−v0). Quasi-finiteness makes that cokernel
finite length. Smith form over C[T] identifies its C-length with the
degree of this determinant. Generically the resultant is a nonzero
constant times R^k, with k=[C(x,y):C(F1,F2,T)]; this is the norm of F2−v
over C(u,T). F2 is integral over C[u,T], so R can be normalized monic in v;
the norm differs from R^k by a nonzero constant. Specialization gives length(U_y)=k*deg R(T;y), including
primitive-coordinate collisions and nonreduced source points. For k=1
this independently proves the first equality of (I5).

The full target stratum of retained length j is therefore

\[
 a_d=\cdots=a_{j+1}=0,\qquad a_j\ne0.                    \tag{I6}
\]

At a singleton source zero fibre, with the source origin mapped to zero,

\[
 a_0(0)=0,\quad a_1(0)\ne0,\quad
 a_2(0)=\cdots=a_d(0)=0.                                \tag{I7}
\]

Thus the omitted length d−1 at zero is encoded in exact coefficient
vanishing. Properness there would require a_d(0)!=0, hence d=1 in this
singleton setting. Equation (I7) does not exclude d>1 by itself.

There is a useful primitivity certificate at a reduced singleton fibre.
If T is source-finite but not yet known primitive, let R be its minimal
equation and k=[C(x,y):C(F1,F2,T)]. Normalize X in the full source field
instead. The same affine-normalization and completed-rank argument gives

\[
 \operatorname{length}(U_y)=k\deg R(T;y).
\]

Thus length(U_0)=1 forces k=1. For a singleton-zero Keller candidate,
source finiteness itself therefore certifies primitivity. A one-point
nonreduced fibre would not give this conclusion.

## 4. Newton polygons give generic missing length and pole multiplicity

Let H be a target prime divisor. Normalize its valuation by v_H(H)=1.
For every prime D of N above H write e_D, f_D for its normal index and
residue degree. For an omitted D put p_D=-v_D(T)>0. The primitive equation
has minimum coefficient valuation zero. Its Newton polygon gives

\[
 \boxed{v_H(a_d)=\sum_{D\text{ omitted over }H}p_Df_D.}   \tag{I8}
\]

Indeed a geometric conjugate at D has T-valuation -p_D/e_D and occurs
e_D*f_D times. The final positive-slope part of the polygon has total
horizontal length

\[
 \boxed{\sum_{D\text{ omitted over }H}e_Df_D.}            \tag{I9}
\]

It counts the generic omitted scheme length. In contrast, (I8) records
the pushforward of the pole divisor. Its multiplicities are not missing
sheet counts. A general polygon does not separate all primes or determine
e and f individually; residual factorization and, if needed, further
Newton normalization data must be retained.

## 5. A complete plane control with an unramified omitted divisor

For any integer a>=1 consider

\[
 F_a(x,y)=(u=x,\ v=x^a y^2-y).
\]

This is dominant and quasi-finite of degree2. It is not Keller:
J(F_a)=2*x^a*y−1 is nonconstant. The inverse coordinate T=y satisfies

\[
 R_a(T)=u^aT^2-T-v.
\]

Its coefficients generate the unit ideal. The projective root cover has
the finite chart U=A2, and infinity chart, with z=1/T,

\[
 u^a=z+v z^2.                                           \tag{I10}
\]

The infinity chart is smooth: at z=0 its z-derivative is a unit, and at
z!=0 its v-derivative is a unit. Thus the whole cover is already normal.
The complete divisor ledger is:

| Prime | Ramifies | Omitted | Target image | e | f | Pole order of T |
|---|---|---|---|---:|---:|---:|
| D=(u=z=0) in the infinity chart | no | yes | u=0 | 1 | 1 | a |
| R=(2*x^a*y=1) in U | yes | no | 1+4*u^a*v=0 | 2 | 1 | 0 |

These are disjoint. At u=0 each target fibre has one reduced source point
and one reduced omitted point. At 1+4*u^a*v=0 the source has one point of
length2 and nothing is omitted. Elsewhere it has two reduced points.
The origin fibre ideal is (x,x^a*y^2−y)=(x,y). The path (s,s^(−a)) maps
to (s,0) and escapes along the unramified omitted branch as s tends to
zero. The retained path is (s,0).

The Newton polygon at u=0 has vertices (0,0),(1,0),(2,a): its negative
root valuation is -a, with horizontal length one. The leading coefficient
has multiplicity a, while exactly one sheet is omitted. This verifies
the distinction between (I8) and (I9). The retained ramification explains
why this source-A2 control is outside the Keller hypotheses of PF2OMIT2.

The same cover also demonstrates why primitivity and a coefficient unit
ideal are not sufficient by themselves. Choose the different primitive
function W=x^a*y. Its equation is

\[
 W^2-W-u^a v=0.
\]

It is monic, and its affine hypersurface is the full smooth finite
normalization N. In the omitted chart W=1+v*z, so W equals1 on D and has
no pole there. The distinguished source is N minus D, not the entire
finite-W chart. Reconstruction y=W/u^a has a pole on D; on the retained
u=0,W=0 branch it equals v/(W−1) and is regular. Along the escaping path,
the target coordinates tend to zero and W stays1. Hence the source is
not finite over A[W]. This primitive equation has leading coefficient one
and would hide the whole omitted divisor if the source-finiteness test
were skipped.

## 6. The degree eight marked-root control

On the n=4 zero-lower-coefficient slice, write

\[
 P(t)=qt^4+rt^3+t-b,\qquad a^2=P'(t),\qquad
 T=x+y=t/a+(1-a)/t.
\]

Eliminating a gives

\[
 (t^2-P'(t))^2-P'(t)(Tt-1)^2=0.
\]

The resultant of this equation and P(t) has coefficient content b*q^2.
After dividing it out, [inverse_equations.json](inverse_equations.json)
retains an exact degree-eight equation R_8 with 218 monomials. Its leading
coefficient is exactly b*Delta, where Delta is the binary quartic
discriminant. The coefficient ideal is the unit ideal. One squarefree
degree-eight specialization, together with the known degree-eight marked
cover, certifies that T is primitive. This uses the explicit marked-root
cover from [Fonollosa, Sections5.2–5.4](https://arxiv.org/html/2609.33491v1#S5),
not the plane-only monicity argument.

In the infinity root chart T=1/v+u−v, so its leading pole coefficient is
one along Ecrit. At the omitted zero marking a=−1, T=t/a+(1−a)/t has
leading pole coefficient two. Thus T extends to P1 on the full smooth
normalization and its inverse infinity is exactly both omitted divisors.
The normalized infinity ledger is complete, with no hidden indeterminacy.

| Target prime | Escaping Newton edge | Generic omitted length | e | f | Pole order |
|---|---|---:|---:|---:|---:|
| Delta=0 | width4, slope1/4 | 4 | 4 | 1 | 1 |
| b=0 | width1, slope1 | 1 | 1 | 1 | 1 |

The first index follows from its slope denominator4 and total width4;
the second follows from width1. Neither is read from the multiplicity of
b*Delta alone. Exceptional values include:

| Target polynomial | Specialized degree | Retained length | Omitted local lengths |
|---|---:|---:|---|
| t−b, b!=0 | 2 | 2 | 6 at infinity |
| t | 1 | 1 | 6 at infinity, 1 at the omitted zero marking |
| 2t^4−3t^3+t=t(t−1)^2(2t+1) | 3 | 3 | 4 at the double root, 1 at the omitted zero marking |
| (t−1)^3(t+1)/2 | 2 | 2 | 6 at the triple root |
| (t^2−t−1/2)^2 | 0 | 0 | 4 at each double root |

The complete occurring source cardinalities on this slice are
{0,1,2,3,4,7,8}. The partitions of four projective root multiplicities
give this upper list via 2*(number of simple roots) minus the zero-root
correction. The exact targets in the checker realize every value; retained
length4 is supplied by the cubic −(t−1)^2(t+2)/3, with simple infinity.
In particular cardinalities5 and6 do not occur. This is a complete fibre
count for the explicit degree-eight control, not for a hypothetical
degree-six plane map.

At q=r=0 the equation specializes to

\[
 R_8(T)=(-T+b)(bT+b^2-2).
\]

For b=1 both distinct retained source points have T=1, giving a double
finite root of R_8. Their normalization points remain etale of length1.
This is an explicit primitive-element collision, not ramification of
those source points. At b=0 the equation is 2T and the missing length is7.

## 7. Application to the remaining r7 and singleton-zero ledgers

For an actual degree-six plane Keller candidate, a source-finite primitive
inverse R_6 would give all target omission images by factoring a_6, and
all aggregate exceptional lengths by (I6). Generic Newton and residual
data would identify the primes over those factors. This is a concrete
completeness package, stronger than a ramification-only passport.

Under PF2OMIT2's entire single e2,f1 ramification hypothesis, every gap21
r7 image requires an additional nonproper target component. Consequently
a_6 cannot be a power of the ramified image's implicit polynomial. In
the immersed case, if pole orders are p_R,p_1,p_2,p_3, its divisor is

\[
 \operatorname{div}(a_6)=p_R[C]+\sum_{j=1}^3p_j[C_j],
 \qquad p_R,p_j\ge1.                                    \tag{I11}
\]

Repeated images contribute separately to the displayed multiplicities.
This gives an exact leading-coefficient factorization gate without
assuming that all omitted images are branch components. The separate
registered adjoint argument PF2R7EX1 excludes the two nonimmersed classes
under its stated geometric imports; their old ledgers remain conditional
controls. No global inverse R_6 for the remaining immersed family or for
an unrestricted singleton-zero candidate is currently supplied.

The polynomial (I1), the monicity and primitive-coordinate certificates,
the coefficient unit ideal, the normalized infinity primes, their e/f,
and their local Artin algebras are the remaining reviewable input package.
The checker replays explicit resultants and local controls; the general
normalization and length-transport arguments above are written proofs.
There is no claimed new global Keller map, full r7 exclusion or proof of
singleton-zero properness.
