# Davenport quadratic-vector developability

Consider the full quadratic auxiliary vector

\[
\boxed{
\mathbf F(T,Y,U)
=\mathbf h(T,Y)+U\mathbf a(T,Y)+U^2\mathbf b(T,Y),
}
\tag{1}
\]

with three polynomial components.  This strictly contains the
constant-direction quadratic incidence already excluded.

The two highest Jacobian coefficients have a geometric meaning.  They
force the projectivized quadratic direction \([\mathbf b]\) to be a curve,
not a surface.  For affine-linear \(\mathbf b\), this gives a complete
normal-form reduction: the generic rank-two case is impossible, and only
two rank-one forms remain.

## 1. The six-coefficient determinant ledger

Write

\[
[x,y,z]=\det(x,y,z)
\]

for the scalar triple product.  Then

\[
\det D\mathbf F=\sum_{j=0}^5J_jU^j,
\]

where

\[
J_0=[h_T,h_Y,a],
\]

\[
J_1=[a_T,h_Y,a]+[h_T,a_Y,a]+2[h_T,h_Y,b],
\]

\[
\begin{aligned}
J_2={}&[b_T,h_Y,a]+[h_T,b_Y,a]+[a_T,a_Y,a]\\
&+2[a_T,h_Y,b]+2[h_T,a_Y,b],
\end{aligned}
\]

\[
\begin{aligned}
J_3={}&[b_T,a_Y,a]+[a_T,b_Y,a]\\
&+2[b_T,h_Y,b]+2[h_T,b_Y,b]+2[a_T,a_Y,b],
\end{aligned}
\]

\[
J_4=[b_T,b_Y,a]+2[b_T,a_Y,b]+2[a_T,b_Y,b],
\]

\[
\boxed{J_5=2[b_T,b_Y,b].}
\tag{2}
\]

A Keller map must have

\[
J_0\in K^*,\qquad J_1=\cdots=J_5=0.
\]

## 2. The developability equation

On the chart \(b_3\ne0\), put

\[
r_1=\frac{b_1}{b_3},\qquad r_2=\frac{b_2}{b_3}.
\]

Then

\[
[b_T,b_Y,b]
=b_3^3
\frac{\partial(r_1,r_2)}{\partial(T,Y)}.
\tag{3}
\]

Thus \(J_5=0\) says that the projective map

\[
(T,Y)\dashrightarrow[b_1:b_2:b_3]\in\mathbb P^2
\]

has differential rank at most one.  Equivalently, the ratios of the
quadratic coefficients generate a function field of transcendence degree
at most one.  The direction of the quadratic \(U\)-curve may vary along
only one rational parameter.

This does not force \(\mathbf b\) to be constant.  It divides the remaining
problem into constant-projective-direction, one-parameter projective, and
nonlinear developable cases.

## 3. Complete affine-linear classification at the top

Let

\[
\mathbf b=b_0+Tp+Yq
\]

with constant vectors \(b_0,p,q\).  Equation (2) becomes

\[
[b_T,b_Y,b]=[p,q,b_0].
\]

Hence \(J_5=0\) if and only if \(b_0,p,q\) are linearly dependent.  The
image of \(\mathbf b\) lies in a fixed target plane.

There are three nonconstant affine normal forms.

### 3.1 Rank-two image

If \(p,q\) are independent, affine source and linear target changes give

\[
\mathbf b=(T,Y,0).
\tag{4}
\]

Write \(a_3\) for the transverse component of \(\mathbf a\).  The
\(U^4\)-equation is

\[
\boxed{
(1-2T\partial_T-2Y\partial_Y)a_3=0.
}
\tag{5}
\]

On every homogeneous polynomial of degree \(d\), this operator acts by
\(1-2d\).  Since no nonnegative integer satisfies \(2d=1\), equation (5)
forces

\[
a_3=0.
\]

Now write \(h_3=H(T,Y)\).  The \(U^3\)-equation reduces to

\[
-2(T H_T+YH_Y)=0.
\tag{6}
\]

Thus \(H\) is constant.  All three vectors
\(\mathbf h_T,\mathbf h_Y,\mathbf a\) then lie in the same target plane, so
the whole Jacobian vanishes.  Therefore:

\[
\boxed{
\text{No affine-linear rank-two quadratic direction can be Keller.}
}
\tag{7}
\]

This is universal; it does not use the coefficients of the Davenport
polynomial.

### 3.2 Rank-one affine line not through the origin

If the derivative of \(\mathbf b\) has rank one and its affine image spans
a plane through the target origin, normalize

\[
\mathbf b=(1,T,0).
\tag{8}
\]

The \(U^4\)-equation gives

\[
\partial_Ya_3=0,\qquad a_3=A(T).
\]

Writing

\[
\mathbf h=(T,g_T(Y),H),\qquad
\mathbf a=(a_1,a_2,A(T)),
\]

the \(U^3\)-equation is the exact derivative

\[
\partial_Y\left(
2H+(2TA'-A)a_1-2A'a_2
\right)=0.
\tag{9}
\]

Consequently

\[
2H+(2TA'-A)a_1-2A'a_2=C(T).
\tag{10}
\]

Put

\[
v=a_2-Ta_1.
\]

Then (10) becomes

\[
H=\frac A2a_1+A'v+\frac C2.
\tag{11}
\]

After this substitution, the \(U^2\)-coefficient is again an exact
\(Y\)-derivative.  It vanishes precisely when

\[
\boxed{
\frac A2a_1^2+2A'a_1v-2A''v^2-2C'v-4A'g=D(T)
}
\tag{12}
\]

for another polynomial \(D(T)\).

Thus the first rank-one survivor is not an unrestricted coefficient
system.  It asks for a polynomial parametrization \((a_1,v)\) of a conic
whose constant term contains the Davenport polynomial \(g_T(Y)\).
Equations (11)--(12) solve \(J_3=J_2=0\) completely.

There is one more exact integration.  Eliminate \(g\) using (12).  Then

\[
J_1=\partial_Y\mathcal R_T(a_1,v),
\]

where

\[
\mathcal R_T(a,v)=-\frac{\mathcal N_T(a,v)}{24A'}
\tag{13}
\]

and

\[
\begin{aligned}
\mathcal N_T={}&
6(AA''+(A')^2)a^2v+3AC'a^2\\
&+8(A'A'''-3(A'')^2)v^3\\
&+12(A'C''-3A''C')v^2\\
&+12\bigl(
4T(A')^2-2AA'-DA''+A'D'-(C')^2
\bigr)v.
\end{aligned}
\tag{14}
\]

Consequently

\[
\boxed{\mathcal R_T(a_1,v)=E(T)}
\tag{15}
\]

for a third polynomial \(E(T)\).  The entire system
\(J_3=J_2=J_1=0\) has therefore been integrated.  Only
\(J_0\in K^*\) remains.

For fixed generic \(T\), equation (15) is a cubic of the special form

\[
a^2(\mu v+\nu)+P_3(v)=E.
\tag{16}
\]

At least one of \(a,v\) must depend nontrivially on \(Y\), since (12)
contains the degree-seven polynomial \(g_T(Y)\).  A nonconstant polynomial
map from \(\mathbb A^1_Y\) into a smooth projective cubic would extend to a
map \(\mathbb P^1\to E\) and hence be constant.  Therefore the cubic (16)
must be singular or reducible over the generic \(T\)-field.

When \(\mu\ne0\), its projective singular locus has three explicit
branches:

1. \(a=0\), with \(E\) a critical value of \(P_3\):
   \[
   P_3(v)=E,\qquad P_3'(v)=0;
   \]
2. \(v=-\nu/\mu\), with
   \[
   E=P_3(-\nu/\mu).
   \]
3. a singularity at infinity, which occurs when the cubic and quadratic
   coefficients of \(P_3\) both vanish:
   \[
   [v^3]P_3=[v^2]P_3=0.
   \]

For nonconstant polynomial \(A\),

\[
\mu=-\frac{AA''+(A')^2}{4A'}
\]

is not identically zero.  Thus these three branches exhaust the generic
singular-cubic escape.

### 3.2.1 The singularity-at-infinity branch is impossible

The third branch says

\[
A'A'''-3(A'')^2=0,\qquad
A'C''-3A''C'=0.
\tag{17}
\]

If \(n=\deg A'\ge1\), comparison of leading coefficients in the first
equation gives

\[
n(n-1)-3n^2=-n(2n+1)\ne0.
\]

Hence \(A'\) is constant and

\[
A=pT+q,\qquad p\ne0.
\]

The second equation then gives

\[
C=cT+c_0.
\]

The cubic (15) is now linear in \(v\).  Polynomiality leaves only one
exceptional possibility capable of carrying the \(Y\)-dependence of
\(g\): \(a=a(T)\), the coefficient of \(v\) in the cubic vanishes, and

\[
v=
\frac{4pg+D-\frac12Aa^2}{2(pa-c)}.
\tag{18}
\]

The vanishing coefficient fixes

\[
D'=-\frac p2a^2-2pT+2q+\frac{c^2}{p}.
\tag{19}
\]

Substitution of (18)--(19) into the remaining determinant gives

\[
J_0=
-\frac{(aa'-2)\,\Xi(T,g)\,g_Y}{4(pa-c)^2},
\tag{20}
\]

where

\[
[g]\Xi=-8p^3\ne0.
\]

If \(aa'-2=0\), then \(J_0=0\).  Otherwise (20) has \(Y\)-degree

\[
7+6=13
\]

and cannot be a nonzero constant.  Therefore the entire
singularity-at-infinity branch is closed.

### 3.2.2 The two finite singular branches are impossible

It remains to consider the two affine singularities of (16).

For the critical-value branch, let \(r\) be the critical point.  Then

\[
P_3(v)-E=(v-r)^2(\rho v+\kappa),
\]

and the cubic is

\[
a^2(\mu v+\nu)+(v-r)^2(\rho v+\kappa)=0.
\tag{21}
\]

If this cubic is irreducible, its normalization has at least two points at
infinity.  A nonconstant polynomial map
\(\mathbb A^1_Y\to(21)\) would extend to a surjective map of projective
normalizations, but the single point at infinity of \(\mathbb A^1\) cannot
cover two distinct omitted points.  Hence no such polynomial
parametrization exists.

There are two reducible cases:

1. \(\mu r+\nu=0\), when \(v-r\) is a component and the remaining factor
   is the conic from the vertical branch;
2. \(\rho v+\kappa\) is proportional to \(\mu v+\nu\), when (21) splits
   into three lines over the algebraic closure.

Thus it suffices to treat line and conic components.

For the vertical branch, the cubic factors as

\[
(v-v_0)Q_2(a,v)=0,\qquad v_0=-\nu/\mu.
\tag{22}
\]

On the line \(v=v_0\), the conic identity (12) represents the degree-seven
polynomial \(g\) by a quadratic polynomial in \(a\).  If
\(\deg_Ya=m>0\), this has degree \(2m\ne7\); if \(a\) is constant, it
cannot contain \(g\).

If the conic \(Q_2=0\) has two points at infinity, it admits no nonconstant
polynomial parametrization for the same puncture reason.  If it degenerates
into lines, substitute a line parameter \(w\) into (12).  The right side
is quadratic in \(w\):

- if its quadratic coefficient is nonzero, \(2\deg_Yw=7\), impossible;
- if it vanishes, \(w\) is affine in \(g\).  Then every output of the
  quadratic suspension is a polynomial in \((T,g,U)\), so the chain rule
  gives
  \[
  \det D\mathbf F
  =g_Y\det\frac{\partial\mathbf F}{\partial(T,g,U)}.
  \]
  The nonunit \(g_Y\) cannot divide a nonzero constant.

The only smooth conic with one point at infinity occurs when
\(\rho=0\).  As above, polynomiality forces \(A\) to be linear.  If the
conic is nondegenerate, it has the parabolic form

\[
v=\lambda a^2+\lambda_0,\qquad \lambda\ne0.
\]

Substitution in (12) has leading term \(2A'a v\), so

\[
\deg_Yg=3\deg_Ya,
\]

again incompatible with degree seven.

These alternatives exhaust both finite singular branches.  Consequently

\[
\boxed{
\text{the entire affine non-origin rank-one form } \mathbf b=(1,T,0)
\text{ is impossible.}
}
\tag{23}
\]

If \(A'=0\) and \(A\ne0\), the Davenport term disappears from (12) and the
third output is a polynomial source coordinate after a triangular change,
reducing to a relative plane Keller problem.  If \(A=0\), the constant
Jacobian equation retains a nonunit \(g_Y\) factor and is impossible.  The
genuinely new case therefore has

\[
A'\ne0.
\]

The conic identity is closely aligned with the existing moving-tangent
calculations: cancellation of its leading square terms forces a
proportional relation between \(a_1\) and \(v\), after which the next
coefficient is a tangent condition on \(g_T\).  This supplies a concrete
bridge from the quadratic-vector ansatz back to the two global Davenport
marking charts.

### 3.3 Rank-one line through the origin

The scalar-direction normal form is

\[
\mathbf b=(T,0,0).
\tag{24}
\]

Here \(J_5=J_4=0\) automatically.  The first equation is

\[
\boxed{
2T\,J(a_2,a_3)-a_2(a_3)_Y+a_3(a_2)_Y=0.
}
\tag{25}
\]

It constrains the projective linear-\(U\) direction
\([a_2:a_3]\), but it does not force it to be constant.

## 4. What survives

The general quadratic-vector opening has not been eliminated.  It has been
reduced sharply:

1. a genuinely two-parameter affine quadratic direction is impossible;
2. every affine survivor depends on only one source coordinate;
3. the non-origin line form (8) is impossible after its three integrated
   relations and singular-cubic component analysis;
4. the through-origin form (24) is governed first by the projective PDE
   (25); and
5. a nonlinear \(\mathbf b\) must parameterize a developable projective
   curve.

The best next calculation has shifted to the through-origin form (24).
Classify polynomial solutions of (25), beginning with coprime
\((a_2,a_3)\).  Its equation has a weighted projective interpretation and
should decide whether that last affine rank-one class reduces to the
constant-direction centralizer obstruction or supplies a genuinely new
candidate.

## 5. Verification

Run

```bash
.venv/bin/python scripts/verify_davenport_quadratic_vector_developability.py
```

The checker verifies all six universal determinant coefficients, the
projective developability identity, the affine-plane criterion, the
rank-two Euler obstruction, the two rank-one normal-form equations, and
the integrated conic and cubic identities, including closure of the
three singular-cubic branches.
