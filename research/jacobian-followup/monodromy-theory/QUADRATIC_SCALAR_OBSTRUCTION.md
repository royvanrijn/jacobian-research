# Quadratic pole pairs with nonconstant scalars

**Scope update:** [KMFI54--55](QUADRATIC_CLASSIFICATION.md) now prove
all-weight injectivity and complete the quadratic primitive-pair reduction
for scalars of every degree. The finite certificates, bounded claims and
historical proof gate below are retained as replay evidence.

Work over algebraically closed characteristic-zero constants.
Use the primitive tensor notation of KMFI48:

P=A*[[-bc,b^2],[-c^2,bc]], gcd(b,c)=1,

with W(A,b,c)=(Ab^2)_ww+2(Abc)_wz+(Ac^2)_zz=0.
Here max(deg b,deg c)=2. No assertion below classifies every
higher-degree primitive pair.

**KMFI51 (all scalar degrees, constant leading direction).**
If the leading projective direction of the pair is constant,
every polynomial scalar solution forces, after affine conjugacy,

b=a*w+c0, c=w^2+d*w+e*z+f.

In particular rho=c/b belongs to the already completed k(w)[z]
chart. A z term in b cannot be repaired by a scalar of any degree.
The proof is a weighted positive-definite tridiagonal matrix argument.

**KMFI52 (diagonal curved case and exact finite certificates).**
In the other possible top direction, normalize to

b=w^2+beta*z+b1*w+b0,
c=kappa*w*z+c1*w+c2*z+c0, beta!=0.

The weighted leading pair is (w^2+z,kappa*w*z), after scaling,
at weights(w,z)=(1,2).
Its exact homogeneous scalar operator is the rectangular matrix
T_m(kappa) displayed below.
If it is injective, no lower terms of b,c and no scalar with
highest weight m can solve the wave equation.

The scalar can exist only if kappa is a root of C(r,s) from KMFI48
with integers r>=4,s>=0. That pair is unique, and then

deg A=r+s-4, weighted deg A=r+2s-4.

Thus a proposed pair has a finite scalar ansatz, not an unbounded
degree search. The exact retained Bezout certificates prove T_m
injective for every kappa and every m=0,...,20.
Injectivity for all m is **not proved** here.

**KMFI53 (complete quadratic class through scalar degree ten).**
Every primitive quadratic pair with a nonzero polynomial scalar
of total degree at most10 is affinely conjugate to the completed
k(w)[z] chart. Its full geometric group is natural alternating,
natural symmetric, or the binary SL2(F3) degree-eight action.
Neither PSL2(11) action nor natural A4 occurs through this entire
coefficient-unrestricted degree range.

At scalar degree exactly one the exact possibilities sharpen to
A8,A10 or binary SL2(F3) on eight letters. All three are field
actions. This is not a new Keller realization of the binary action.

[Matrix helper](quadratic_scalar_obstruction.py) ·
[checker](verify_quadratic_scalar.py) ·
[Bezout certificates](quadratic_scalar_certificates.json) ·
[model](quadratic_scalar_model.json) · [receipt](quadratic_scalar_receipt.json).

The degree-ten bound is explicit. The symbolic certificates are
universal in all lower pair coefficients and in kappa; they are
not a finite sample of such coefficients and do not prove the
unbounded matrix assertion.

## 1. Constant leading direction: a universal weighted obstruction

KMFI49 puts the leading pair in the form (0,w^2).
Write b=a*w+beta*z+c0,c=w^2+d*w+e*z+f.
If beta!=0, take weights(w,z)=(2,3).
The highest vector field is (beta*z,w^2).
A nonzero constant scaling and linear source scaling normalize
it to (z,w^2).

For a scalar monomial w^i*z^j, its highest wave is

\[
 i(i-1)w^{i-2}z^{j+2}
 +2(i+2)(j+1)w^{i+1}z^j
 +j(j-1)w^{i+4}z^{j-2}. \tag{1}
\]

Fix a weight m=2i+3j. Select the square set of output rows
w^(i+1)z^j for the same input pairs.
The diagonal is a_(i,j)=2(i+2)(j+1)>0.
For adjacent input pairs (i,j),(i-3,j+2), the off-diagonal
product, divided by the two diagonal entries, is

\[
 q_{i,j}={i(j+2)\over4(i+2)(j+3)}<1/4.
\]

A positive diagonal similarity symmetrizes the tridiagonal matrix.
Normalize its diagonal entries to1 by congruence.
Its LDL pivots satisfy p0=1,p_next=1-q/p_prev.
Induction gives every pivot>1/2.
Thus the square matrix is positive definite and invertible over Q,
hence over every characteristic-zero extension.
Every weight sector is injective.

For a general nonzero scalar A, its highest weighted component
would lie in this zero kernel. Lower pair/scalar terms have
strictly smaller weight and cannot cancel it.
Consequently beta=0. This proves KMFI51 for every scalar degree.
Once beta0, the ratio c/b belongs to k(w)[z] regardless of the
other coefficients, and KMFI46 applies.

## 2. The diagonal weighted sector

In the nonconstant projective top case, KMFI49 puts
b_top=w^2,c_top=kappa*w*z and the highest scalar on the two
eigenlines. If beta, the z coefficient of b, is zero, the pair
is already in k(w)[z]. Suppose beta!=0.
Scaling z normalizes the highest weighted pair to
v=(w^2+z,kappa*w*z) at weights1,2.

For A=w^i*z^j its wave is

\[
 D_{i,j}w^{i+2}z^j+E_{i,j}w^iz^{j+1}
                             +F_iw^{i-2}z^{j+2},
\]
\[
 D_{i,j}=(i+4)(i+3)+2\kappa(i+3)(j+1)
                         +\kappa^2(j+1)(j+2),
\]
\[
 E_{i,j}=2(i+1)(i+2+\kappa(j+2)),\qquad F_i=i(i-1).
 \tag{2}
\]

For fixed m=i+2j, T_m has floor(m/2)+1 columns and one more row.
Its column j has D on row j, E on row j+1 and F on row j+2
(the last F is zero). These are the complete output coefficients.

The first nonzero scalar coefficient has D_(i,j)=0.
But D_(i,j)=C(i+4,j), and KMFI34 proves the integer pair for a
fixed nonreal kappa unique. Hence i=r-4,j=s and
m=r+2s-4. Independently, KMFI48 fixes the highest ordinary component to a
monomial determined by the same unique residue pair. Thus it is
w^(r-4)z^s and deg A=r+s-4.
If no such r,s exist, no scalar exists.
Otherwise all scalar coefficients fit in the finite set
i+j<=r+s-4, i+2j<=r+2s-4; substituting the full pair gives
a finite linear system.

## 3. Exact low-weight matrix certificates

For each m=0,...,20, let P_m be the top maximal minor of T_m
(delete its last row) and Q_m its bottom maximal minor
(delete its first row), each divided by its nonzero rational content.
The retained certificate gives explicit polynomials U_m,V_m in Q[kappa]
such that

\[
 U_mP_m+V_mQ_m=1.
\]

Therefore those two minors never vanish simultaneously over any
characteristic-zero field, and T_m is injective for every kappa.
No random specialization or floating-point rank is used.

If total deg A<=10, its (1,2)-weight is at most20.
The diagonal curved case is impossible by these identities;
the constant top-direction curved case is impossible by KMFI51.
Every remaining pair has denominator depending on w alone and
falls within KMFI46. The quadratic pair restricts the possible
imprimitive action to the binary degree-eight exception:
a polynomial h has degree at most2, so no larger even signed
family can arise. Z-dependent ratios have only alternating or
symmetric groups. This proves the first part of KMFI53.

The top minor is a product of C(i+4,j), each with negative
rational discriminant, so it has no real root.
The bottom minors in the retained range have only simple real
roots, as certified by exact Sturm counts.
An all-m proof of real-rootedness of the bottom continuant would
therefore prove injectivity for every m. That statement is a
remaining proof gate, not a theorem asserted here.

For completeness, the degree-at-most-ten action list can be
parametrized without an unspecified group test:

| reduced primitive pair | parameter range | scalar degree | full degree |
|---|---|---:|---:|
| (1,h), deg h2 | s=0..3, epsilon0 or1 |3s+epsilon|KMFI44 with d2|
| (w,w^2+kappa*z) |r>=2,s>=0,r+2s<=12|r+2s-2|(s+2)(r+2s+2)|
| (w^2,kappa*w*z+1) |r>=s+4,r+s<=14|r+s-4|(s+2)(r-s-1)|

The latter two rows use C(r,s)=0; their group is S_N exactly
for even r and odd s, otherwise A_N.
The first row is A_N, except s0,epsilon1: it is binary degree8
when h is even about the root of lambda and A8 otherwise.
The coefficients of h and all allowed integration constants
can vary without changing this criterion.

Here is why these rows exhaust the bound. In the constant-top
case, b is constant or linear in w.
For constant b, h is quadratic. Since deg_z(lambda)<=10,
only the22 leading heat models s0..10,epsilon0 or1 need consideration.
Their exact leading w powers, retained in the checker, have degree
3s+epsilon and cannot be canceled by lower-weight data.
Only s0..3 survive the degree-ten bound.

For linear b, KMFI33 gives a z-dependent ratio with one simple
pole and quadratic R after removing affine terms.
The scalar is lambda/w^2; its highest ordinary degree is
r+2s-2 from the nonzero bottom coefficient of the j2 weighted model.
This gives the second row.
For a diagonal top pair with b depending only on w, primitivity
and the single-pole theorem force b=w^2 and
c=kappa*w*z+c1*w+c0,c0!=0, after translation/scaling.
Remove the affine term in R. Its R=c0/w has sharp bound r-s>=4,
and the scalar lambda/w^4 has degree r+s-4, giving the third row.
Each row has explicit polynomial wave witnesses by KMFI41 or44.
These are sweep fields; no new Keller occurrence is inferred.

The proof does not extend the finite certificate list to m>20.
It does, however, exclude every candidate with r+2s-4<=20,
even if its ordinary scalar degree exceeds10.

## 4. The exact scalar-linear action list

These reductions are also explicit and independent checks on the
first nonconstant degree. Translate the affine scalar to w or z.

In the diagonal top case:
- A=w has kappa^2+4kappa+10=0.
  The w*z wave coefficient forces beta0. Then KMFI33 makes
  the z-leading lambda=w*b(w)^2 a single-pole monomial, forcing
  b=w^2. The p coefficient must have a simple pole, so
  c=kappa*w*z+c1*w+c0, with c0!=0 by primitivity.
  This is r5,s0,R=c0/w up to affine R, giving degree8 and A8.
- A=z has kappa^2+2kappa+2=0.
  The wave coefficients force beta=c1=c0=0.
  Thus b(w) is quadratic and c=z*(kappa*w+c2).
  The KMFI33 single simple pole would force b to be a square
  of a linear factor canceled once by c, contradicting primitivity.
  This case has no solution.

For constant leading direction:
- If the scalar has a z term, an affine shear makes A=z.
  The wave coefficients give beta0,3a+2e0,
  a^2+2ae+3e^2=0. Hence a=e0, and the remaining coefficient
  forces c0=0, contradicting the primitive pair.
- Otherwise A=w. The coefficients give beta0,
  3a^2+2ae+e^2=0 and c0*(2a+e)=0.
  If a0, then e0,c0!=0. This is lambda proportional to w and
  quadratic h; its group is A8 unless the w term of h vanishes,
  when it is the binary degree-eight group.
  If a!=0, then c0=0 and e/a satisfies kappa^2+2kappa+3=0.
  This is r3,s0 and quadratic R, with degree10 and A10.

Examples are A=w,b=w^2,c=kappa*w*z+1 for A8;
A=w,b=w,c=w^2+kappa*z for A10;
and A=w,b=1,c=w^2 for the binary degree-eight field.
Their wave identities use the corresponding displayed quadratics.
The previous centered-origin obstruction remains valid.
[KMFI60](BRANCH_CENTERED_COMPLETION.md) now supplies an absolute
binary Keller realization after affine branch centering.

Universal weighted positive-definiteness, exact sector reductions
and finite-degree Bezout identities are written and replayable.
No formal verification or independent whole-proof review is claimed.
The general quadratic problem above the certified sector range,
higher primitive degrees and the full action spectrum remain open.
