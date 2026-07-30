# Obstructions for the canonical \(D_5\) two-mask blowdown

## 1. Outcome

The first nonlinear-classification search has a canonical target
modification.  Put

\[
\Delta(U,V)=V^2-4U^5
\]

and define

\[
\boxed{
\begin{aligned}
\beta\colon\mathbb A^4_{U,V,A,B}&\longrightarrow
\mathbb A^4_{U,V,S,T},\\
(U,V,A,B)&\longmapsto
\left(
U,V,\ VA+2UB,\ 2U^4A+VB
\right).
\end{aligned}}                                      \tag{1.1}
\]

Its mask matrix is

\[
M(U,V)=
\begin{pmatrix}
V&2U\\
2U^4&V
\end{pmatrix},
\qquad
\det M=\Delta.                                      \tag{1.2}
\]

Thus \(\beta\) is generically birational and

\[
\boxed{
A=\frac{VS-2UT}{\Delta},\qquad
B=\frac{-2U^4S+VT}{\Delta}.
}                                                   \tag{1.3}
\]

This is exactly the two-adjugate-numerator mechanism required by the
nonlinear obstruction classification.

The result of the first search is negative but strong.

1. If the two base coordinates of (1.1) remain a polynomial Dickson
   incidence \((u_0,P_5(a_0,u_0))\), polynomiality and constant Jacobian are
   incompatible, regardless of how nonlinearly \(a_0,u_0,A,B\) depend on
   the source coordinates.
2. No constant linear rechart of the four inputs to (1.3) can evade this.
   Constant Jacobian forces its first two rows to be only the standard cusp
   scalings \(U\mapsto cU,\ V\mapsto dV\), with \(d^2=c^5\).  The inverse
   divisibility then forces the two mask rows into their span, contradicting
   invertibility.
3. Exact division also rejects all \(72\) coordinate assignments for three
   primitive source charts: the unchanged incidence chart and the two
   normalized ramification-incidence charts \(R_+=0\) and \(R_-=0\).
4. At every polynomial degree, a rechart whose constant-Jacobian condition
   requires
   \(q^2-4p^5=\lambda(V^2-4U^5)\) cannot mix either mask into \(p,q\).
   Holding the other variables fixed would give a rational curve on a
   smooth genus-two generic fibre, which is impossible.

Therefore the canonical two-mask determinant is correct, but every
automorphic base-moving assembly is closed.  A surviving construction must
use a nonautomorphic log-crepant incidence rechart before applying (1.3),
or a different nonlinearly coupled blowdown.

No Keller map is constructed here.

## 2. Determinant and inverse

The derivative of \(\beta\) is block triangular because its first two
coordinates are \(U,V\).  Hence

\[
\det D\beta=\det M=\Delta.                          \tag{2.1}
\]

The adjugate identity

\[
\operatorname{adj}(M)
=
\begin{pmatrix}
V&-2U\\
-2U^4&V
\end{pmatrix}
\]

gives (1.3).  Polynomial factorization through \(\beta\) requires both
divisibilities

\[
\Delta\mid VS-2UT,\qquad
\Delta\mid -2U^4S+VT.                              \tag{2.2}
\]

Checking only \(\det M\) misses these two independent residue conditions.

## 3. Polynomial Dickson-base obstruction

Let \(X=\mathbb A^4\) have arbitrary coordinates, and suppose

\[
a_0,u_0,A,B\in\mathcal O(X)
\]

are polynomials.  Consider

\[
F=
\left(
u_0,\ P_5(a_0,u_0),\ A,\ B
\right).                                           \tag{3.1}
\]

Let

\[
\widetilde F=(a_0,u_0,A,B).
\]

The chain rule gives

\[
\boxed{
\det DF
=-J_5(a_0,u_0)\det D\widetilde F.
}                                                   \tag{3.2}
\]

Both factors on the right are polynomials.  Since \(J_5\) is nonconstant,
(3.2) cannot be a nonzero constant.

This argument does not assume that \(u_0\) is one of the original source
coordinates.  It may depend nonlinearly on both masks.  The obstruction is
that the old source coordinates \(a_0,u_0\) themselves are polynomially
visible before applying the Dickson incidence.

Consequently, factoring through the fixed-base blowdown (1.1) cannot work:
if (1.3) is polynomial, its output is exactly of the form (3.1).
Preservation of the \(D_5\) field must expose \(a_0\) only rationally, not
as a polynomial coordinate preceding \(P_5\).

## 4. Constant-linear base rigidity

One might assign constant linear combinations of four polynomial outputs
\((U,V,Z,W)\) to the inputs \((p,q,S,T)\) of (1.3).  Constant Jacobian
requires the new exceptional determinant to be proportional to the old
one:

\[
q^2-4p^5=\lambda(V^2-4U^5),
\qquad \lambda\ne0.                                \tag{4.1}
\]

Write

\[
\begin{aligned}
p&=p_UU+p_VV+p_ZZ+p_WW,\\
q&=q_UU+q_VV+q_ZZ+q_WW.
\end{aligned}
\]

Comparison of coefficients in (4.1) gives, in order,

\[
\begin{array}{c|c}
\text{coefficient}&\text{consequence}\\ \hline
Z^5,\ W^5&p_Z=p_W=0,\\
Z^2,\ W^2&q_Z=q_W=0,\\
V^5&p_V=0,\\
U^2&q_U=0,\\
V^2,\ U^5&q_V^2=\lambda=p_U^5.
\end{array}                                        \tag{4.2}
\]

Thus

\[
\boxed{p=cU,\qquad q=dV,\qquad d^2=c^5.}            \tag{4.3}
\]

After rescaling, it remains to ask whether linear mask forms \(S,T\) can
satisfy (2.2).  The first numerator has degree at most two, while
\(\Delta\) has degree five, so divisibility forces

\[
VS-2UT=0.                                          \tag{4.4}
\]

Since \(U,V\) are coprime,

\[
S=2\ell U,\qquad T=\ell V                         \tag{4.5}
\]

for a constant \(\ell\).  The four linear forms \(p,q,S,T\) are then
dependent.  They cannot be the rows of an invertible linear rechart.

### Theorem 4.1

No constant linear rechart of the canonical two-mask blowdown can make a
primitive \(D_5\) assembly simultaneously polynomial and Keller.

This theorem includes arbitrary constant linear combinations, not only
coordinate permutations.

## 5. All-degree mask rigidity

The coefficient comparison has a geometric all-degree extension.

### Theorem 5.1

Let

\[
p,q\in k[U,V,Z_1,\ldots,Z_r]
\]

satisfy

\[
q^2-4p^5=\lambda(V^2-4U^5),
\qquad \lambda\in k^\times.                         \tag{5.1}
\]

Then \(p\) and \(q\) are independent of every \(Z_i\).

### Proof

Fix \(i\), and regard all variables except \(Z_i\) as elements of the field

\[
K=k(U,V,Z_1,\ldots,\widehat Z_i,\ldots,Z_r).
\]

Equation (5.1) defines a \(K\)-morphism

\[
\mathbb A^1_{Z_i}
\longrightarrow
\mathcal C,\qquad
\mathcal C:\ y^2=4x^5+\lambda(V^2-4U^5).           \tag{5.2}
\]

The constant on the right is nonzero in \(K\), so the smooth projective
completion of \(\mathcal C\) is a genus-two hyperelliptic curve.  The map
extends from \(\mathbb A^1\) to a morphism
\(\mathbb P^1\to\overline{\mathcal C}\).  By Riemann--Hurwitz there is no
nonconstant morphism from \(\mathbb P^1\) to a genus-two curve.  Hence
\(p,q\) are independent of \(Z_i\).  Repeat for every mask. \(\square\)

Thus no polynomial automorphic incidence rechart satisfying the required
constant determinant ratio can move the branch coordinates with masks,
regardless of degree.  The first possible survivor must replace the scalar
identity (5.1) by a nonautomorphic log-crepant equation containing the
Jacobian of the rechart.

The elementary leading-face shadow is also useful.  If \(p,q\) did depend
on masks with degrees \(d,e>0\), their highest mask forms would have to
satisfy

\[
5d=2e,\qquad p_d=h^2,\qquad q_e=\pm2h^5.           \tag{5.3}
\]

In particular, equal mask-degree bounds cannot contain a survivor, and the
first formal cusp-shaped face would occur at degrees \((2,5)\).  Theorem
5.1 shows that even this face cannot integrate while the right side of
(5.1) stays fixed.

## 6. Three primitive source charts

The exact finite screen retains the primitive source Jacobian \(CQ\).

### 6.1 Unchanged chart

With source variables \((a,u,x,y)\), use

\[
\alpha_0=(a,u,CQx,y).
\]

Then

\[
(\pi_5\times\operatorname{id})\circ\alpha_0
=
(u,P,CQx,y)
\]

has determinant \(-5CQ^2=-5\Delta(u,P)\).

### 6.2 Ramification-incidence charts

For \(\gamma=\alpha\) or \(\beta\), where

\[
\alpha=\frac{3+\sqrt5}{2},\qquad
\beta=\frac{3-\sqrt5}{2},
\]

put

\[
a_0=s,\qquad u_0=\frac{s^2}{\gamma}+x.
\]

Thus \(x=0\) maps the source zero section into \(R_\gamma=0\).  With source
variables \((s,t,x,y)\), define

\[
\alpha_\gamma
=
\left(
a_0,u_0,t,C(a_0,u_0)Q(a_0,u_0)y
\right).                                           \tag{5.1}
\]

Its determinant is \(\pm CQ\), and

\[
(\pi_5\times\operatorname{id})\circ\alpha_\gamma
=
\left(
u_0,P(a_0,u_0),t,CQy
\right)                                            \tag{5.2}
\]

again has determinant \(\pm5\Delta(u_0,P)\).

For each of these three charts, assign its four displayed outputs in all
\(4!=24\) orders to \((U,V,S,T)\) in (1.3).  Exact multivariate division
gives

\[
\boxed{
\text{none of the \(72\) assignments satisfies both divisibilities (2.2).}
}                                                   \tag{5.3}
\]

Among those assignments, only the two orders per chart which retain
\((U,V)=(u_0,P)\) even pass the determinant-divisor test.  They fail the
adjugate test, as predicted by Section 3.

The finite screen is a regression for the broader automorphic theorem.  Its
value is that it checks both normalized ramification-center charts
explicitly before any nonlinear coefficients are introduced.

## 7. The first nonautomorphic cusp chart also fails

The simplest nonautomorphic incidence moves the old mask-zero section into
the normalized cusp and gives both transverse corrections one copy of the
old branch equation.  Put

\[
\begin{aligned}
D&=V^2-4U^5,\\
p&=V^2+DS,\\
q&=2V^5+DT.
\end{aligned}                                      \tag{7.1}
\]

Then \(q^2-4p^5\) vanishes modulo \(D\), so

\[
E=\frac{q^2-4p^5}{D}                               \tag{7.2}
\]

is polynomial.  To complete the rechart while retaining \(U\) as a third
coordinate, a fourth polynomial \(s(U,V,S,T)\) would have to solve the
log-crepant contraction

\[
\det\frac{\partial(p,q,s)}{\partial(V,S,T)}=E.      \tag{7.3}
\]

There is an immediate divisor contradiction.  The Hamiltonian contraction
vector is

\[
\nabla p\times\nabla q
=
\left(
D^2,\,
-2V(S+1)D,\,
-2V(T+5V^3)D
\right),                                           \tag{7.4}
\]

up to the harmless common choice of signs.  Hence the left side of (7.3)
is divisible by \(D\) for every polynomial \(s\).

The right side is not.  Modulo \(D\),

\[
\boxed{
E\equiv4V^5(T-5SV^3)\pmod D,
}                                                   \tag{7.5}
\]

which is nonzero.  Therefore (7.3) has no polynomial solution of any
degree.

This closes the first nonautomorphic normalized-cusp chart before a
coefficient ansatz is introduced.  Its failure identifies the next
necessary feature: the two transverse corrections cannot both be
divisible by the old branch equation.  Their reductions modulo \(D\) must
instead move tangentially together along the normalized cusp.

## 8. Every affine-normal tangential chart fails

The next repair is to move along the normalized cusp modulo \(D\).  Up to
rescaling the first mask, the general chart affine in the normal coordinate
is

\[
\begin{aligned}
h&=V+S,\\
p&=h^2+A(U,V,S)DT,\\
q&=2h^5+B(U,V,S)DT,\\
\rho&=(p,q,U,S),
\end{aligned}                                      \tag{8.1}
\]

where the normal coefficients do not depend on \(T\).  Its Jacobian has
degree at most one in \(T\).  For constant coefficients \(A=a,B=b\), it
specializes to

\[
\det D\rho
=2hD(5ah^3-b).                                     \tag{8.2}
\]

The log-crepant equation would be

\[
q^2-4p^5=\lambda D\det D\rho,
\qquad \lambda\ne0.                                \tag{8.3}
\]

For general \(A,B\), the right side of (8.3) has \(T\)-degree at most one.
On the left, the \(T^5\) coefficient is

\[
-4A^5D^5,
\]

so \(A=0\).  The right side is then independent of \(T\), while the
\(T^2\) coefficient on the left is

\[
B^2D^2,
\]

so \(B=0\).  The remaining map has
\(\det D\rho=0\) and is not a rechart.

Hence:

\[
\boxed{
\text{no tangential cusp chart affine in one normal coordinate has a
nondegenerate log-crepant solution.}
}                                                   \tag{8.4}
\]

This is an exact coefficient solution, not a bounded numerical search.

## 9. Updated nonlinear search frontier

The fixed-base and all automorphic base-moving routes are now closed.  The
first symmetric nonautomorphic cusp chart is also closed.  The next
candidate must replace (5.1) by a nonautomorphic incidence rechart

\[
(U,V,S,T)\longmapsto(p,q,r,s)
\]

such that:

1. \(q^2-4p^5\) satisfies the required log-Jacobian identity rather than
   merely being a scalar copy of \(\Delta(U,V)\);
2. the old zero section is moved into the exceptional divisor;
3. both adjugate numerators are divisible;
4. the rechart is nonautomorphic or otherwise avoids polynomial visibility
   of the old root coordinate \(a_0\);
5. elimination gives a finite-flat rank-five algebra.

The lowest useful coefficient search must solve the log-crepant equation
with the rechart Jacobian on its right side.  Solving only a polynomial
identity of the form (5.1), at any mask degree, is now proved insufficient.
Modulo the old branch, its base coordinates must have the form

\[
p\equiv h^2,\qquad q\equiv\pm2h^5\pmod D,           \tag{8.1}
\]

with \(h\) genuinely mask-dependent.  This supplies tangential rank along
the cusp while avoiding the common-\(D\) conormal factor in (7.4).
Section 8 closes arbitrary \(T\)-independent coefficients in the single
affine normal direction.  The
[all-degree continuation](DIHEDRAL_ALL_DEGREE_AFFINE_COMPLETION_OBSTRUCTIONS.md)
also closes every nonlinear \(T\)-dependence: a one-normal candidate is
either nonresonant and fails the degree bound, or resonant and degenerates
by factorization or the infinity valuation.  Do not launch further
one-normal coefficient searches.

The first unresolved choices are therefore:

1. two independent nonlinear tangential/normal directions, with neither
   auxiliary direction reduced to a retained-coordinate differential
   field \(K[T]\); or
2. a different birational blowdown whose exceptional determinant is
   \(\Delta\) but whose adjugate rows are not the matrix (1.2).

## 10. Reproduction

Run

```bash
.venv/bin/python scripts/verify_d5_two_mask_blowdown_obstructions.py
```

The checker verifies the determinant and adjugate formulas, the polynomial
Dickson-base chain-rule factor, the coefficient proof of constant-linear
base rigidity, the mask-row dependence, the smooth genus-two generic fibre
used by Theorem 5.1, all \(72\) exact coordinate assignments, and the
all-degree contraction-divisor failure of the first nonautomorphic cusp
chart.  It also excludes the full tangential class affine in one normal
coordinate.

The odd/even valuation ledgers and the uniform extension of these
obstructions to every \(n\ge3\) are in
[All-degree dihedral affine-completion obstructions](DIHEDRAL_ALL_DEGREE_AFFINE_COMPLETION_OBSTRUCTIONS.md).
