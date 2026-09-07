# Davenport independent marked-line opening

The product presentation of the Davenport incidence slope is not intrinsic.
Replacing the weighted coordinates

\[
\sigma=BC,\qquad \tau=cAC^2
\]

by independent slope and intercept coordinates removes the two-chart gluing
obstruction completely.  The overlap extends over the whole quadratic base
as a triangular polynomial automorphism.

This does not yet give an absolute Keller map.  It isolates the remaining
problem much more sharply: polynomialize a reciprocal modification for the
single derivative divisor

\[
D=H_W(s,W)-\sigma
\]

without making the marked root recoverable.  The elementary reciprocal
modification fails, and the first genuinely alternating quadratic Jung
pencil also fails an exact unit-ideal gate.

Work over

\[
K=\mathbb Q(a),\qquad a^2+a+2=0.
\]

## 1. The overlap is already polynomial

Let \(H_4,H_2\) be the reduced point-cover primitives from the two
proportional tangent charts.  Their root coordinates satisfy

\[
W_2=W_4+\delta,\qquad \delta=\frac{3a+5}{21},
\]

and

\[
H_4(W)-H_2(W+\delta)=m(s)W+n(s).
\]

For the inverse equation

\[
E_i(W_i)=H_i(W_i)-\sigma_iW_i+\tau_i,
\]

put

\[
\boxed{
\sigma_4=\sigma_2+m(s),\qquad
\tau_4=\tau_2-n(s)-\delta\sigma_2.
}
\]

Then

\[
E_4(W_4)=E_2(W_4+\delta)
\]

identically.  Both the root transition and the target transition are
triangular polynomial automorphisms with determinant one.  In particular,
the \(1/C\) pole found in the weighted atlas is entirely an artifact of
forcing \(\sigma\) to factor as \(BC\).

The two-chart problem is therefore solved at the marked-incidence level.
The bad polynomials of the tangent charts no longer appear in the overlap
map.  They belonged to the weighted endpoint normalization, not to the
degree-seven inverse equation itself.

## 2. The exact reciprocal-ring problem

Use

\[
\tau=H(s,W)-\sigma W.
\]

The marked-line core

\[
(W,\sigma)\longmapsto(\sigma,\tau)
\]

has Jacobian

\[
-\bigl(H_W(s,W)-\sigma\bigr)=-D.
\]

Equivalently, after replacing \(\sigma\) by the coordinate
\(D=H_W-\sigma\), the old marked source is simply

\[
R=K[s,W,D].
\]

An absolute algebraization of this opening requires a second polynomial
ring

\[
R'\cong K[x_1,\ldots,x_N]
\subset \operatorname{Frac}(R)[D^{-1}]
\]

such that:

1. the Davenport target functions \(s,\sigma,\tau\), together with any
   genuinely coupled auxiliary outputs, lie in \(R'\);
2. the volume form of \(R\) pulls back with factor \(D^{-1}\);
3. the resulting target is also affine space;
4. the generic target field still cuts out the degree-seven
   \(\mathrm{GL}_3(\mathbb F_2)\) point or line extension; and
5. no target combination recovers \(W\).

This is the precise independent-coordinate version of the missing absolute
gate.

## 3. Why the elementary reciprocal modification fails

The smallest reciprocal chart is

\[
x=DW,\qquad W=\frac{x}{D}.
\]

It has the desired Jacobian:

\[
\det\frac{\partial(s,W,D)}{\partial(s,x,D)}=D^{-1}.
\]

But its target pullbacks are

\[
\sigma=H_W\!\left(s,\frac{x}{D}\right)-D,
\]

\[
\tau=
H\!\left(s,\frac{x}{D}\right)
-\frac{x}{D}H_W\!\left(s,\frac{x}{D}\right)+x.
\]

Already at \(s=0\), the reduced Davenport primitive has degree four.
Consequently the first expression has a genuine \(D^{-3}\) pole and the
second a genuine \(D^{-4}\) pole.  Affine slope/intercept corrections alter
only degrees zero and one in \(W\), so they cannot remove these poles.

Thus the direct affine modification has the correct determinant ledger but
fails polynomiality before any boundary or Hessian comparison is needed.

## 4. The first alternating Jung pencil is also closed

The post-coordinate unit gate applies to an affine-in-\(U\) projective
coefficient pencil after choosing polynomial coordinates \((x,y)\) on the
Davenport \((T,Y)\)-plane.  The first coordinate outside the affine and
one-triangular classes has quadratic--quadratic polydegree.  Write

\[
y=T+q(Y),\qquad x=Y+p(y),
\]

where

\[
q(Y)=cY^2+dY+e,\qquad
p(y)=p_2y^2+p_1y+p_0,\qquad p_2\ne0.
\]

The inverse coordinate formulas are

\[
Y=x-p(y),\qquad T=y-q(x-p(y)).
\]

Put

\[
G(x,y)=g_{T(x,y)}(Y(x,y)).
\]

Here

\[
\deg_y T=4,\qquad \deg_y T_y=3,\qquad \deg_yG=14
\]

generically.  The unit gate would require

\[
A(x)G_y-B(x)T_y\in K^*.
\]

Since \(A\ne0\) in a nontrivial projective pencil, every coefficient of
\(G\) in \(y^{14},y^{13},y^{12},y^{11}\) must vanish.

Take all coefficients in \(x\) of those four terms and saturate their ideal
at \(p_2\).  Exact Gröbner reduction over \(K\) gives

\[
\boxed{
\left\langle
[x^j y^k]G:\ 11\le k\le14
\right\rangle:p_2^\infty=(1).
}
\]

Therefore no quadratic--quadratic alternating coordinate passes the unit
gate.  The conjugate line-cover statement follows by
\(a\mapsto-1-a\).

This closes the first live Jung polydegree proposed in the research audit;
it is stronger than checking a few normalized parabola coordinates because
all quadratic coefficients and translations are included.

## 5. All length-two Jung pencils are closed

Continue with

\[
y=T+q(Y),\qquad x=Y+p(y).
\]

If \(\deg p=m\ge2\) and \(\deg q=n\ge3\), then in the inverse coordinates
\(\deg_yY=m\) and \(\deg_yT=mn\).  The Davenport term

\[
-(5+3a)T^3Y
\]

has the unique largest degree

\[
m(3n+1).
\]

It cannot cancel with any other term, while
\(\deg_yT_y=mn-1\).  Hence the unit gate fails for every such \(m,n\), not
only for the adjacent quadratic--cubic case.

The only length-two orientation not killed by this dominance has
\(\deg q=2\).  This remaining case admits an all-degree Newton calculation.
Put

\[
r=y-e,\qquad u=Y,\qquad T=r-cu^2-du.
\]

Expanding the Davenport polynomial in \(u,r\), its three possible
Newton-leading terms are

\[
g_{r-cu^2-du}(u)
=C_7(c)u^7+C_6(c,d)u^6+C_{5,1}(c)ru^5+\text{lower terms},
\]

where

\[
\begin{aligned}
C_7={}&\frac17\bigl(
(35+21a)c^3+(-21+14a)c^2\\
&\hspace{3.8cm}-(7+7a)c+1
\bigr),\\
C_6={}&(1+a)c^3+(4a-2)c^2-(1+a)c\\
&+d\bigl((9a+15)c^2+(4a-6)c-(1+a)\bigr),\\
C_{5,1}={}&-(9a+15)c^2+(-4a+6)c+(1+a).
\end{aligned}
\]

These coefficients generate the unit ideal:

\[
\boxed{(C_7,C_6,C_{5,1})=(1)\quad\text{in }K[c,d].}
\]

There is also a short explicit certificate.  The monic gcd of \(C_7\) and
\(C_{5,1}\) is

\[
c+\frac{1+2a}{7}.
\]

At its only possible common root

\[
c=-\frac{1+2a}{7},
\]

the middle coefficient is independent of \(d\) and equals

\[
C_6=-\frac{10(2a+1)}{49}\ne0.
\]

Now let \(p\) be any polynomial of degree \(m\ge2\).  Since
\(u=x-p(y)\), the first nonzero one of the three displayed coefficients
forces

\[
\deg_yG\in\{7m,\ 6m,\ 5m+1\}.
\]

All three values exceed \(2m\), while

\[
\deg_yT=2m,\qquad \deg_yT_y=2m-1.
\]

Therefore \(A(x)G_y-B(x)T_y\) cannot be a unit.  This proves the unit-gate
obstruction for every nonlinear \(p\), without a degree bound.  The
quadratic, cubic, and quartic saturated computations above are finite
regressions of this theorem rather than the frontier.

For the reverse alternating orientation, the \(Y^7/7\) term is uniquely
dominant as soon as both shears are nonlinear, so it supplies no surviving
length-two class.  If either shear is affine, the coordinate reduces to an
affine or one-triangular case already excluded in the post-coordinate
audit.  Hence every affine-in-\(U\) coefficient pencil whose plane
coordinate has Jung length at most two is closed.

## 6. What remains live

The independent marked-line opening survives, but every elementary
reciprocal or length-two affine-in-\(U\) realization now fails.  The next
bounded search should be one of:

1. a Jung length-three coordinate, beginning with the
   quadratic--quadratic--quadratic polydegree;
2. genuinely nonlinear dependence on the auxiliary variable \(U\);
3. two coupled modification coordinates whose combined source is affine
   space and whose target does not expose \(W\); or
4. a non-elementary reciprocal link that mixes the Davenport base \(s\)
   with the polar root coordinate.

The cheapest remaining coefficient calculation is the length-three
quadratic Jung screen.  Structurally, however, nonlinear \(U\)-dependence
is the first route not governed by the affine coefficient-pencil
factorization, so it is the more promising change of mechanism.

## 7. Verification

Run

```bash
.venv/bin/python scripts/verify_davenport_independent_marked_line.py
```

The checker reconstructs the two reduced primitives, proves the polynomial
overlap identity and determinant-one transition, verifies the marked-line
Jacobian, measures the two elementary reciprocal poles, checks the first
saturated Jung ideal, and proves the all-degree three-coefficient Newton
gate.
