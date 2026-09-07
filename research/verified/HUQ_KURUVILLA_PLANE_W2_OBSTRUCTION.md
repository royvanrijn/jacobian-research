# The plane characteristic-two map has no Keller lift modulo four

This is a repository extension of Mondello's external plane theorem,
[*A Dimension-Two Counterexample to the Separable Jacobian Conjecture in
Characteristic Two*, arXiv:2608.02634v1](https://arxiv.org/abs/2608.02634).
The lifting obstruction, stabilization, and Witt-tower statements below are
not claims from Mondello's paper.

## 1. Statement

Let

\[
\begin{aligned}
 P&=x+x^2y+x^4+x^6y^2,\\
 Q&=y+x^5+x^6y+x^7y^2+x^8y^3
\end{aligned}
\]

over \(\mathbb F_2\).  This is the determinant-one, noninjective plane map
from the
[characteristic-two audit](HUQ_KURUVILLA_CHARACTERISTIC_TWO_AUDIT.md).

**Theorem (HKM2W1).**

1. For every determinant-one plane map over \(\mathbb F_2\), the cokernel of
   its first Jacobian-variation operator is canonically the top algebraic de
   Rham cohomology of the source:
   \[
   \operatorname{coker}\mathcal D_F
   \simeq H^2_{\mathrm{dR}}(\mathbb F_2[x,y])
   \simeq xy\,\mathbb F_2[x^2,y^2].
   \]
2. For the displayed map, the integral Jacobian error has class
   \[
   [K]=xy\bigl(1+x^6+x^8y^2\bigr)=xyu^2\ne0,
   \qquad u=1+x^3(1+xy).
   \]
   Hence no polynomials
   \(\widetilde P,\widetilde Q\in(\mathbb Z/4)[x,y]\) reducing to
   \(P,Q\) have constant Jacobian determinant.
3. Nonvanishing of this class is invariant under arbitrary polynomial
   left--right equivalence of plane maps over \(\mathbb F_2\).  Thus no
   polynomially left--right equivalent representative has a Keller lift
   through \(W_2(\mathbb F_2)\).
4. Every first-order Jacobian obstruction becomes exact after adjoining one
   identity coordinate.  For this map the exact stabilization extends to an
   explicit compatible Keller lift through every finite Witt ring
   \(W_n(\mathbb F_2)=\mathbb Z/2^n\).

Consequently the plane map has no compatible plane lift through higher Witt
vectors or through \(\mathbb Z_2\).  These are unrestricted polynomial
statements: no correction-degree bound is imposed.

**Theorem (HKM2W2: sharp stabilized \(W_2\) degree).** Among all polynomial
maps over \(\mathbb Z/4\) that reduce to \((P,Q,z)\) and have constant
Jacobian determinant, the minimum possible maximum total coordinate degree
is exactly \(18\).  Thus the elementary stabilized lift in Section 7 is
degree-minimal at the first Witt level, even when all three coordinates may
be corrected by arbitrary polynomials involving \(x,y,z\).

**Theorem (HKM2W3: exact stabilized \(W_3\) degree).** Among all polynomial
maps over \(\mathbb Z/8\) that reduce to \((P,Q,z)\) and have constant
Jacobian determinant, the minimum possible maximum total coordinate degree
is exactly

\[
 d_3=19.
\]

Both halves are exact and computer-assisted.  A complete necessary Boolean
coefficient system for degree 18 has 1,083 variables and 1,639 equations,
and Z3 4.15.3 certifies that it is unsatisfiable.  At degree 19, a complete
system in a \(z\)-linear first-correction ansatz has 2,685 variables and 4,513
equations and is satisfiable; a separate sparse-polynomial replay of the
returned model gives an odd constant determinant modulo eight.  The upper
bound only needs one ansatz witness, while the lower bound is unrestricted.
Multiplying one target coordinate by that odd constant gives determinant one
without changing the reduction or degree.  The preferred frozen witness
already has determinant one.  These calculations have no independent
implementation or external human review.

**Theorem (HKM2W4: fixed-representative \(W_4\) boundary).** The preferred
degree-19 representative from `HKM2W3` has no extension through
\(W_4(\mathbb F_2)=\mathbb Z/16\), in any polynomial degree.  Its next
determinant digit has a nonzero class in
\(H^3_{\mathrm{dR}}(\mathbb F_2[x,y,z])\).  This does not obstruct every
degree-19 representative: a second pinned degree-19 \(W_3\) lift has zero
next class.  For that fixed lift, however, the exact minimum extension degree
is again 52.  Its full next error has 1,027 terms and degree 52; exact affine
elimination excludes degree 51 and a 1,086-term correction attains degree 52.

For the explicit degree-25 representative (8.12), the top de Rham class
vanishes, but the exact minimum degree of an extension of that fixed
representative is 52.  A two-coefficient functional excludes every degree at
most 51, and adding a single third-coordinate correction gives determinant
one modulo 16 in degree 52.  Consequently the unrestricted minimum satisfies

\[
 \boxed{19\le d_4\le52}.
\]

No exact unrestricted value of \(d_4\) is claimed.  The preferred
support-optimized \(W_3\) point is nonextendable, while both an 818-term
degree-19 point and the degree-25 point extend first in degree 52.  These
fixed-representative results do not exclude a different degree-19
representative extending already in degree 19.

## 2. First-order Jacobian correction

Use the displayed zero-one formulas as integral representatives \(P_0,Q_0\).
Every lift modulo four is uniquely of the form

\[
 \widetilde P=P_0+2A,\qquad \widetilde Q=Q_0+2B
\]

for polynomials \(A,B\) modulo two.  Direct differentiation gives

\[
 \det D(P_0,Q_0)=1+2K,
\]

where, modulo two,

\[
\begin{aligned}
K={}&xy+x^7y+x^{10}y+x^5y^2+x^{11}y^2\\
   &+x^9y^3+x^{12}y^3+x^{13}y^4.                 \tag{2.1}
\end{aligned}
\]

In particular

\[
 [xy]K=1.                                         \tag{2.2}
\]

Modulo four, the determinant of an arbitrary lift is

\[
 \det D(\widetilde P,\widetilde Q)
 =1+2\bigl(K+\mathcal D_F(A,B)\bigr),            \tag{2.3}
\]

with the first Jacobian-variation operator

\[
 \boxed{\mathcal D_F(A,B)
 =A_xQ_y+P_xB_y+A_yQ_x+P_yB_x}.                   \tag{2.4}
\]

over \(\mathbb F_2\).  The signs become plus signs in characteristic two.

## 3. The full cokernel is top de Rham cohomology

The following argument applies to every determinant-one polynomial plane map
\(G=(G_1,G_2)\) over \(\mathbb F_2\).  Put
\(R=\mathbb F_2[x,y]\) and \(\omega=dx\wedge dy\).  The first variation
satisfies

\[
 \mathcal D_G(A,B)\,\omega
 =dA\wedge dG_2+dB\wedge dG_1
 =d\bigl(A\,dG_2+B\,dG_1\bigr).                \tag{3.1}
\]

Because \(dG_1\wedge dG_2=\omega\), the forms \(dG_1,dG_2\) are an
\(R\)-basis of \(\Omega_R^1\).  Hence every one-form is uniquely of the
form \(A\,dG_2+B\,dG_1\), and (3.1) identifies the image of
\(\mathcal D_G\) with the coefficients of all exact two-forms.  Therefore

\[
 \boxed{\operatorname{coker}\mathcal D_G
 \simeq H^2_{\mathrm{dR}}(R/\mathbb F_2).}       \tag{3.2}
\]

This quotient is elementary.  A monomial \(x^iy^j\) is an \(x\)-derivative
when \(i\) is even and a \(y\)-derivative when \(j\) is even.  Conversely,
a nonzero \(x\)-derivative has even \(x\)-exponent, and a nonzero
\(y\)-derivative has even \(y\)-exponent.  Thus every class has a unique
representative supported on monomials odd in both variables:

\[
 \boxed{H^2_{\mathrm{dR}}(R/\mathbb F_2)
 =xy\,\mathbb F_2[x^2,y^2]\,\omega.}           \tag{3.3}
\]

In particular the cokernel is rank one over the Frobenius subring
\(\mathbb F_2[x^2,y^2]\), but infinite-dimensional over \(\mathbb F_2\).
The old \(xy\)-coefficient test detects only the constant coefficient of
this Frobenius-linear obstruction.

## 4. The exact Cartier class of the Jacobian error

Projecting (2.1) to the odd--odd monomials from (3.3) gives

\[
\begin{aligned}
 [K]
 &=xy+x^7y+x^9y^3\\
 &=xy\bigl(1+x^6+x^8y^2\bigr).                  \tag{4.1}
\end{aligned}
\]

With \(r=1+xy\) and \(u=1+x^3r\),

\[
 u^2=1+x^6r^2=1+x^6+x^8y^2,
\]

so

\[
 \boxed{[K]=xyu^2\ne0.}                          \tag{4.2}
\]

Equivalently, the Cartier isomorphism sends this class to

\[
 C\bigl([K\omega]\bigr)=u\,dx\wedge dy.        \tag{4.3}
\]

This computes the entire obstruction, not only one coefficient.  In
particular \([xy]K=1\) is the lowest-coordinate witness for the nonzero
Cartier class.

For completeness, choose any integral representatives of a determinant-one
map over \(\mathbb F_2\) and write their Jacobian as \(1+2K_G\) modulo
four.  Changing the representatives by \(2(A,B)\) changes \(K_G\) by
\(\mathcal D_G(A,B)\).  Hence

\[
 \mathfrak o(G):=[K_G]\in H^2_{\mathrm{dR}}(R/\mathbb F_2) \tag{4.4}
\]

is independent of every choice.  A constant polynomial represents zero in
top de Rham cohomology, so \(G\) has a constant-Jacobian lift through
\(\mathbb Z/4\) if and only if \(\mathfrak o(G)=0\).  Equation (4.2)
therefore proves the plane lifting obstruction.

## 5. Polynomial left--right invariance

Let \(\sigma,\tau\) be polynomial automorphisms of the source and target
plane over \(\mathbb F_2\), and set

\[
 G=\tau\circ F\circ\sigma.
\]

By the Jung--van der Kulk theorem, every plane polynomial automorphism is a
composition of affine and triangular automorphisms.  Lifting the coefficients
of those generators gives polynomial automorphisms over \(\mathbb Z/4\)
with constant unit Jacobian.  We may therefore compute the obstruction of
\(G\) using compositional lifts.  If \(\widetilde\sigma,\widetilde\tau\)
lift \(\sigma,\tau\), their determinants are odd constants
\(c_\sigma,c_\tau\).  For a lift of \(F\) with determinant \(1+2K\), the
chain rule gives

\[
 \det D(\widetilde\tau\circ\widetilde F\circ\widetilde\sigma)
 =c_\tau c_\sigma\bigl(1+2K\circ\sigma\bigr)\pmod4. \tag{5.1}
\]

Write \(c_\tau c_\sigma=1+2c\pmod4\).  The constant term \(c\omega\) is
exact, while \(\det D\sigma=1\) over \(\mathbb F_2\); hence

\[
 \mathfrak o(G)=\sigma^*\mathfrak o(F).           \tag{5.2}
\]

Pullback by \(\sigma\) is an automorphism of algebraic de Rham cohomology.
Consequently

\[
 \boxed{\mathfrak o(G)=0\iff\mathfrak o(F)=0.}   \tag{5.3}
\]

Since \(\mathfrak o(F)=xyu^2\ne0\), no polynomially left--right equivalent
plane representative over \(\mathbb F_2\) has a constant-Jacobian lift
through \(\mathbb Z/4\).  The same proof works after a perfect extension of
\(\mathbb F_2\), using the corresponding Witt-vector lifts of affine and
triangular generators.

## 6. General stabilization theorem

There is an all-dimensional form of (3.1).  Let
\(R_n=k[x_1,\ldots,x_n]\), let \(G=(G_1,\ldots,G_n)\) have determinant
one, and let \(\omega_n=dx_1\wedge\cdots\wedge dx_n\).  For a correction
\(A=(A_1,\ldots,A_n)\), variation of
\(dG_1\wedge\cdots\wedge dG_n\) gives

\[
 \mathcal D_G(A)\omega_n
 =d\left(
 \sum_{i=1}^n(-1)^{i-1}A_i\,
 dG_1\wedge\cdots\widehat{dG_i}\cdots\wedge dG_n
 \right).                                         \tag{6.1}
\]

Since the \(dG_i\) form a basis, this identifies

\[
 \operatorname{coker}\mathcal D_G
 \simeq H^n_{\mathrm{dR}}(R_n/k).                \tag{6.2}
\]

Now adjoin an identity coordinate \(z\).  If \(K\omega_n\) represents any
first obstruction for \(G\), then

\[
 K\omega_n\wedge dz
 =d\bigl((-1)^n zK\omega_n\bigr).               \tag{6.3}
\]

Thus every first Jacobian obstruction becomes exact after one identity
stabilization.  Explicitly, correcting the new coordinate by \(-zK\) kills
the first-order error (the sign is immaterial in characteristic two).  This
is a general cohomological explanation of the instability seen here; it does
not by itself construct compatible lifts at all higher Witt levels.

## 7. One stabilization gives a full finite-Witt tower

The same calculation gives a sharp contrast.  Since

\[
 \det D(P_0,Q_0)=1+2K,
\]

put \(h=2K\) and, for every \(n\ge2\), define

\[
 S_n=\sum_{j=0}^{n-1}(-h)^j,
 \qquad
 \widetilde F_n(x,y,z)=\bigl(P_0,Q_0,zS_n\bigr) \pmod {2^n}.  \tag{7.1}
\]

This reduces to \((P,Q,z)\).  Its Jacobian matrix is block lower triangular.
The finite geometric-series identity gives

\[
 \det D\widetilde F_n
  =(1+h)S_n=1-(-h)^n=1\pmod {2^n}.               \tag{7.2}
\]

Moreover \(S_{n+1}\equiv S_n\pmod {2^n}\), so these polynomial maps form a
compatible tower.  Thus the obstruction class is nonzero in the plane
correction complex but becomes exact after one identity stabilization.  This
is not merely a failure of the particular \([xy]\) detector: (7.1) cancels
the entire Jacobian error at every finite Witt level.

The degrees of \(S_n\) grow with \(n\).  Their inverse limit is the restricted
two-adic power series

\[
 (1+2K)^{-1}=\sum_{j\ge0}(-2K)^j,                              \tag{7.3}
\]

not a polynomial in \(\mathbb Z_2[x,y]\).  Accordingly, (7.1) is a compatible
formal/Witt lift; it is not a finite-degree characteristic-zero polynomial
Keller map.

## 8. Sharp degree at the first stable Witt level

Let \(\deg\widetilde G\) denote the maximum total degree of the three
coordinates.  Every lift of \((P,Q,z)\) through \(\mathbb Z/4\) has the form

\[
 \widetilde G=(P_0+2A,Q_0+2B,z+2C),
 \qquad A,B,C\in\mathbb F_2[x,y,z].              \tag{8.1}
\]

Suppose \(\deg\widetilde G\le17\).  Write \(A_0,B_0\) for the coefficients
of \(z^0\) in \(A,B\), and \(C_1\) for the coefficient of \(z^1\) in
\(C\).  Then

\[
 \deg A_0,\deg B_0\le17,
 \qquad \deg C_1\le16.                           \tag{8.2}
\]

The first variation of the stabilized map is

\[
 \mathcal D_{F\times z}(A,B,C)=\mathcal D_F(A,B)+C_z. \tag{8.3}
\]

Consequently the \(z^0\)-coefficient of its half-Jacobian error is

\[
 K+\mathcal D_F(A_0,B_0)+C_1.                   \tag{8.4}
\]

Consider the two-coefficient functional

\[
 \Lambda(R)=[x^{13}y^4]R+[x^{14}y^5]R.          \tag{8.5}
\]

Using

\[
 P_x=1,\quad P_y=x^2,\quad
 Q_x=x^4+x^6y^2,\quad Q_y=1+x^6+x^8y^2,
\]

direct coefficient bookkeeping gives the all-degree identity

\[
 \boxed{\Lambda(\mathcal D_F(A_0,B_0))
 =[x^{15}y^5]A_0.}                              \tag{8.6}
\]

Indeed, the possible \([x^9y^5]A_0\), \([x^7y^3]A_0\), and
\([x^{13}y^5]B_0\) contributions occur twice and cancel; the displayed
degree-twenty coefficient is the only survivor.  Under (8.2), both the
right side of (8.6) and \(\Lambda(C_1)\) vanish.  On the other hand, (2.1)
gives

\[
 \Lambda(K)=1.                                  \tag{8.7}
\]

A constant half-Jacobian error is killed by \(\Lambda\), so (8.4)--(8.7)
contradict the assumption that \(\widetilde G\) has constant determinant.
Thus every such lift has degree at least \(18\).

The lift

\[
 (P_0,Q_0,z(1+2K))\pmod4                       \tag{8.8}
\]

has degree \(18\), and its determinant is
\((1+2K)^2=1\pmod4\).  This proves `HKM2W2`.

More generally, if \(d_n\) is the minimum degree at Witt level \(n\),
reduction to \(W_2\) and (7.1) give

\[
 18\le d_n\le17(n-1)+1,\qquad n\ge2,            \tag{8.9}
\]

with equality at \(n=2\).  Determining \(d_n\) for \(n\ge3\), or proving
that it is unbounded, remains open.

At the next level, the geometric-series upper bound \(d_3\le35\) is not
optimal.  Put

\[
 A_3=x^{15}y^2(1+x^4y^4),\qquad
 \Delta=x^2y^2+x^{10}y^4,\qquad C_3=K+\Delta.    \tag{8.10}
\]

An exact calculation over \(\mathbb F_2\) gives

\[
 \mathcal D_F(A_3,0)=K^2+\Delta.                \tag{8.11}
\]

Therefore

\[
 \widetilde F_3=
 \bigl(P_0+4A_3,Q_0,z(1+2K+4C_3)\bigr)\pmod8   \tag{8.12}
\]

has maximum coordinate degree \(25\).  Its block-triangular determinant is

\[
\begin{aligned}
 \det D\widetilde F_3
 &=(1+2K+4\mathcal D_F(A_3,0))(1+2K+4C_3)\\
 &=1+4\bigl(K+K^2+\mathcal D_F(A_3,0)+C_3\bigr)\\
 &=1\pmod8.                                      \tag{8.13}
\end{aligned}
\]

This first gives the preliminary bounds

\[
 \boxed{18\le d_3\le25.}                        \tag{8.14}
\]

The lower endpoint here is only reduction to \(W_2\).  It will be improved
below; no claim of sharpness at \(W_3\) is made.

There is nevertheless a sharp statement in the canonical first-digit gauge.
Require the \(W_3\) lift to reduce modulo four to
\((P_0,Q_0,z(1+2K))\).  Before the second correction, its determinant is

\[
 (1+2K)^2=1+4(K+K^2)\pmod8.                    \tag{8.15}
\]

The error \(K+K^2\) contains \(x^{26}y^8\), of total degree \(34\).  If a
second correction had maximum coordinate degree at most \(24\), the
\(z^0\)-coefficient of its first variation would be
\(\mathcal D_F(A_0,B_0)+C_1\), with

\[
 \deg\mathcal D_F(A_0,B_0)\le24-1+\deg Q_y=33,
 \qquad \deg C_1\le23.                           \tag{8.16}
\]

It cannot cancel \(x^{26}y^8\).  Hence (8.12) has the exact minimum degree
among extensions of the canonical \(W_2\) lift:

\[
 \boxed{d_3^{\mathrm{can}}=25.}                 \tag{8.17}
\]

Changing the first Witt digit by a noncanonical kernel correction can alter
the quadratic error, so (8.17) does not prove unrestricted \(d_3=25\).

### 8.1. Unrestricted exclusion in degree eighteen

We now allow every first- and second-Witt correction in all three
coordinates.  Suppose, for contradiction, that such a lift has maximum total
coordinate degree at most \(18\), and write it as

\[
 (P_0,Q_0,z)+2(A,B,C)+4(U,V,W)\pmod8.           \tag{8.18}
\]

Only the first two \(z\)-layers of the first correction can affect the
high-degree coefficient equations used below.  Over \(\mathbb F_2[x,y]\),
write

\[
 A=p+za+O(z^2),\qquad B=q+zb+O(z^2),\qquad
 C=c+zr+O(z^2),                                  \tag{8.19}
\]

where \(p,q,c\) have degree at most \(18\), and \(a,b,r\) have degree at
most \(17\).  The first determinant digit gives the complete equations

\[
 K+\mathcal D_F(p,q)+r\in\mathbb F_2,
 \qquad \mathcal D_F(a,b)=0.                    \tag{8.20}
\]

For clarity, the quadratic parts of the next digit in the \(z^0\) and
\(z^1\) layers are respectively

\[
\begin{aligned}
 q_0={}&p_xq_y+p_yq_x+\mathcal D_F(p,q)r\\
      &+c_x(aQ_y+bP_y)+c_y(aQ_x+bP_x),\\
 q_1={}&p_xb_y+p_yb_x+a_xq_y+a_yq_x\\
      &+a(Q_yr_x+Q_xr_y)+b(P_yr_x+P_xr_y).
                                                        \tag{8.21}
\end{aligned}
\]

There is also a coefficientwise Bockstein term.  For \(z^0\), it is the
second binary digit of the integral polynomial

\[
 K_{\mathbb Z}+p_x(Q_0)_y+(P_0)_xq_y
 -p_y(Q_0)_x-(P_0)_yq_x+\det D(P_0,Q_0)r;
\]

where \(K_{\mathbb Z}=(\det D(P_0,Q_0)-1)/2\).  For \(z^1\), it is the
second binary digit of the analogous integral
variation in \((a,b)\).  Equations (8.20) guarantee that all relevant
nonconstant coefficients are even, so these halves are well-defined modulo
two.  The formally possible term \(r\mathcal D_F(a,b)\) in \(q_1\) vanishes
by (8.20).  A \(z^2\)-coefficient of \(C\) contributes to the second digit
of the \(z^1\)-layer, but only in plane degree at most \(16\); higher
\(z\)-layers cannot affect either layer used here.

A degree-18 second correction \((U,V,W)\) can reach plane degree at most
\(27\) in the \(z^0\)-layer and at most \(26\) in the \(z^1\)-layer; the
leading contribution uses \(\deg Q_y=10\).
Therefore every coefficient of the existing second-digit error above those
two bounds must vanish.  Expanding all monomials in (8.19) gives

\[
 3\binom{20}{2}+3\binom{19}{2}=1083
\]

Boolean coefficient variables.  Equations (8.20), together with every
forced high-degree coefficient of the Bockstein terms plus (8.21), give
1,639 exact Boolean constraints.  The checked script constructs this full
system without a sparse-support ansatz; Z3 4.15.3 returns `unsat`.  Hence no
degree-18 lift exists, and (8.12) yields the refined bounds

\[
 \boxed{19\le d_3\le25.}                        \tag{8.22}
\]

This proves the unrestricted lower bound in `HKM2W3`.  It is a
solver-certified finite theorem, not an independent human proof.

### 8.2. Degree-nineteen attainment

The obstruction disappears immediately at the next degree.  Retain the
first-correction form (8.19) with

\[
 \deg(p,q,c)\le19,\qquad \deg(a,b,r)\le18,
\]

and take the second correction in the finite form

\[
\begin{aligned}
 U&=u_0+zu_1+z^2u_2,\\
 V&=v_0+zv_1+z^2v_2,\\
 W&=zw_1+z^3w_3,                                \tag{8.23}
\end{aligned}
\]

with the degree bounds forced by total degree at most \(19\).  The complete
first- and second-digit determinant equations in this ansatz have 2,685
Boolean variables and 4,513 constraints.  The SAT-specialized Z3 pipeline
returns a model.  The checker can print every exact monomial support with
`--show-model`; more importantly, it reconstructs the three polynomials

\[
 (P_0,Q_0,z)+2(A,B,C)+4(U,V,W)\pmod8
\]

from that model and separately differentiates and expands their full
three-by-three Jacobian over \(\mathbb Z/8\).  The maximum coordinate degree
is \(19\), the reduction is exactly \((P,Q,z)\), and the replay gives

\[
 \det D\widetilde F_3\in(\mathbb Z/8)^\times.
\]

Every odd residue squares to one modulo eight, so multiplying one target
coordinate by this constant converts the map to a determinant-one lift with
the same degree and the same reduction modulo two.

The preferred frozen
[support certificate](../artifacts/generated-results/huq_kuruvilla_w3_degree19_witness_first_280_second_160.json)
has determinant one directly and 440 nonzero correction coefficients:

| Witt digit | support counts by layer | total |
|---|---|---:|
| first | \(p:41,q:54,c:52,a:27,b:61,r:45\) | 280 |
| second | \(u_0:41,v_0:12,w_1:39,u_1:28,v_1:26,u_2:2,v_2:5,w_3:7\) | 160 |

The sparse highest-\(z\) tail has the compact form

\[
\begin{aligned}
 u_2&=x+x^2y^5,\\
 v_2&=x^3+x^7+x^3y^6+x^{11}+x^3y^{14},\\
 w_3&=y^2+x^2y^4+x^6+x^2y^8+x^2y^{10}
       +x^2y^{12}+y^{14}.                      \tag{8.24}
\end{aligned}
\]

The checker pins the certificate hash and replays it without invoking the
solver.  Holding its first digit fixed turns the complete second-correction
problem into an affine system over \(\mathbb F_2\): it has 1,485 variables,
709 nontrivial equations, rank 681, and hence affine dimension 804.  This
rank audit is exact; it also makes clear that support minimization is a
minimum-weight representative problem in a large affine code.

The audit verifies more: no equation mixes the three \(z\)-blocks, and each
block splits again by the connected components of its bipartite
variable--equation incidence graph.  There are respectively 59, 14, and 50
nontrivial components.  The largest \(z^0\) component has only 90 variables
and 33 equations.  Minimum Hamming weight is therefore additive across these
small affine decoding problems.

This gives a reusable **incidence-decoding rule**: for any affine
coefficient-correction problem, build the bipartite variable--equation graph
before asking for a global cardinality bound.  Connected components are
independent affine codes, so their minimum weights and certificates add.

Their equation/rank/nullity triples are

| block | variables | equations | rank | nullity |
|---|---:|---:|---:|---:|
| \(z^0:(u_0,v_0,w_1)\) | 610 | 275 | 275 | 335 |
| \(z^1:(u_1,v_1)\) | 380 | 205 | 177 | 203 |
| \(z^2:(u_2,v_2,w_3)\) | 495 | 229 | 229 | 266 |

These figures first allow any odd constant Jacobian.  Requiring determinant
one adds one independent \(z^0\) equation, giving 710 equations, rank 682,
and affine dimension 803 overall.  For the preferred 280-term first digit,
componentwise decoding gives exact block minima \(92,54,14\), both with an
arbitrary constant and with determinant one.  Therefore its exact
second-support minimum is

\[
 s_2=92+54+14=160.                              \tag{8.25}
\]

The earlier 299-term first digit has exact block minima \(111,48,6\), hence
second support 165 and total support 464.  The sparser first digit therefore
lowers both stages and gives total support \(280+160=440\).  A full nonlinear
probe finds this first digit under the bound 280; bound 275 times out, so no
global support minimum over all first corrections is claimed.

This exhibits a second reusable rule: support optimization across Witt digits
is a Pareto problem.  A first-digit kernel change alters the Bockstein and
quadratic syndrome seen by the next affine completion, so minimizing either
digit alone need not minimize their sum.

Together with the unrestricted exclusion in Section 8.1, this proves

\[
 \boxed{d_3=19.}                                \tag{8.26}
\]

The ansatz restriction in (8.23) is harmless for this existence half: one
witness suffices.  No ansatz restriction occurs in the degree-18 exclusion.

### 8.3. The next Witt digit: obstruction versus compatibility

Let \(G_{19}\) be the preferred determinant-one degree-19 representative in
Section 8.2, using the exact pinned zero--one support certificate, and put

\[
 E_{19}=\frac{\det DG_{19}-1}{8}\pmod2.
\]

Direct sparse differentiation gives 1,250 monomials in \(E_{19}\), of
maximum total degree 54.  Its projection onto the odd--odd--odd monomials has
48 terms, all of \(z\)-degree one and maximum total degree 35; in particular,

\[
 [xyz]E_{19}=1.                               \tag{8.27}
\]

For a determinant-one polynomial map in three variables over
\(\mathbb F_2\), the same exterior-form argument as in Section 3 identifies
the cokernel of the first determinant variation with

\[
 H^3_{\mathrm{dR}}(\mathbb F_2[x,y,z])
 =xyz\,\mathbb F_2[x^2,y^2,z^2]\,dx\wedge dy\wedge dz. \tag{8.28}
\]

Every extension of this fixed \(G_{19}\) to \(\mathbb Z/16\) is
\(G_{19}+8(R,S,T)\), and its next determinant digit changes by that first
variation.  Equations (8.27)--(8.28) therefore prove that \(G_{19}\) has no
\(W_4\) extension at any degree.  This is a fixed-representative obstruction:
a different degree-19 \(W_3\) lift could have a different next-digit class.

The obstruction can be imposed during the \(W_3\) search without expanding
the full determinant modulo 16.  Write its Jacobian matrix as
\(M_0+2M_1+4M_2\).  On an odd--odd--odd monomial the weight-zero determinant
vanishes.  The weight-three terms are sums of wedges of exact one-forms:
the \(M_1\)--\(M_2\) cross terms and \(\det M_1\) are exact top forms, so
their top de Rham projections vanish identically.  It is therefore enough to
compile the weight-one and weight-two terms.  This gives 241 Boolean
coefficient equations assembled from 38,760 collected symbolic terms.  On
the preferred certificate, the compiler reproduces the direct 48-term class
and its pinned hash exactly.

Holding the preferred 280-term first Witt digit fixed makes this augmented
system unsatisfiable inside the complete degree-19 existence ansatz (5,954
active constraints).  Removing that gauge restriction changes the answer:
the full 2,685-variable ansatz is satisfiable.  The resulting pinned map,
denoted \(G_{19}^{0}\), has 818 correction coefficients, split as 327 in the
first Witt digit and 491 in the second, and has constant determinant 5
modulo 8.  By construction,

\[
 \left[\frac{\det DG_{19}^{0}-5}{8}\right]
 =0\quad\text{in }H^3_{\mathrm{dR}}(\mathbb F_2[x,y,z]). \tag{8.29}
\]

Thus the proposed unrestricted cohomological UNSAT statement is false even
within the degree-19 ansatz.  Vanishing is only an unbounded extension
criterion; it does not say that a primitive has degree 19.

We also tested the genuinely joint degree-19 problem.  For each master
(W_3) model, direct expansion modulo 16 produces its next error (E), and
the subproblem

\[
 \mathcal D_{F\times z}(R,S,T)=E,
 \qquad \deg(R,S,T)\le 19,                    \tag{8.29a}
\]

is an affine system over \(\mathbb F_2\).  An inconsistent subproblem returns
a dual row combination; compiling the corresponding parity of determinant
coefficients and adding it to the master gives an exact CEGAR loop.  The
first failures have a particularly transparent form: the dual certificate is
a single coefficient outside the monomial support of the degree-19 image of
\(\mathcal D_{F\times z}\).  For example,

\[
 [y^{19}z]\,\mathcal D_{F\times z}(R,S,T)=0
 \quad\text{when }\deg(R,S,T)\le19.           \tag{8.29b}
\]

Indeed, the constant terms of \(Q_y\) and \(P_x\) would require respectively
the source monomials \(xy^{19}z\) and \(y^{20}z\), both of degree 21; all
other plane shifts have positive \(x\)-degree, and the only possible third
coordinate source is \(z^2\), whose derivative vanishes in characteristic
two.  Thus every joint lift must have
\([y^{19}z](\det DG-c)/8=0\).

The same support calculation can be done before selecting a master model.
Inside the displayed degree-19 \(W_3\) ansatz, the determinant modulo 16 has
5,396 structurally possible monomial targets.  Of these, 4,340 have zero
incidence in the nonconstant degree-19 correction operator.  The exact
layer ledger is:

| output \(z\)-degree | correction variables | nonconstant rank | determinant targets | singleton holes |
|---:|---:|---:|---:|---:|
| 0 | 500 | 275 | 1,540 | 1,223 |
| 1 | 270 | 177 | 1,485 | 1,228 |
| 2 | 405 | 229 | 1,431 | 1,165 |
| 3 | 216 | 145 | 940 | 724 |

The sorted structural-hole set has SHA-256
`de442207ad627a8202168496c37fcd2b9af7bb8cf03cbeb96bf90a662097ab99`.
This is a support statement, not a full cokernel computation: targets with
nonzero incidence can still participate in multi-row dual relations.

This first cut is satisfiable, and so are systems with 16 and 64 certified
singleton cuts; every returned master model still has an inconsistent
degree-19 completion subproblem, again witnessed by another singleton
codomain hole.  A stronger experiment compiled all 599 singleton holes
activated by the pinned class-zero model.  Their sorted target list has
SHA-256
`8faad19cd5212c598c270233f1d4407791cc75f00a689f768f690c4267c2bcb1`,
and its 599 coefficient expressions contain 6,095,343 collected determinant
terms.  Z3 returned `unknown` at the 600-second decision bound.  Separately,
the eight-hole system bit-blasts to 465,654 variables and 2,431,415 DIMACS
clauses; MiniSat remained indeterminate after 409.63 CPU seconds and
3,364,085 conflicts.

A second compiler now evaluates the quotient without expanding each cubic
Boolean conjunction separately.  For each determinant permutation it
materializes only those coefficients of a shared two-by-two Jacobian minor
that can reach a requested hole, then contracts that projected minor with
the remaining row.  Direct substitution agrees with independent sparse
determinant expansion on both pinned \(W_3\) representatives.  On the same
599-hole target set, the factored compiler uses 13,804 shared minor
coefficients, 899,303 minor products, and 623,490 final products, reducing
the construction from about twelve minutes to about two on the recorded
machine.  The 600-second master decision still returns `unknown`.

The full model-independent singleton quotient is now practical as well.  All
4,340 holes compile through 15,564 shared minor coefficients, 926,438 minor
products, and 4,200,036 final products.  Its 600-second master decision also
returns `unknown`.  Splitting by output \(z\)-degree gives the following
bounded ledger:

| layer | target hash | factored result |
|---:|:---|:---|
| 0 | `28d89469dd611fe86f59230f70774ce1d7707ad2180f32c9759390bc200d8e79` | `unknown` |
| 1 | `5a795c44391709212e1953243ab675e0abf1b9bdd2838cfe602ed8681fbe5905` | `sat`; returned full completion still inconsistent |
| 2 | `d3c20501e0284dd62938801b798a3428193377ee5e3588d2398ad7ea5b30795f` | `unknown` |
| 3 | `8d5da33464e3fab81f627186e22eeddabe8bacea283a6a7c64524a83c32c4d3e` | `sat`; returned full completion still inconsistent |

Running `solve-eqs` before bit-blasting does not decide the hard \(z^0\)
layer at the same bound.  These are bounded experiments, not an UNSAT proof.
They leave the unrestricted value of \(d_4\) open.  The remaining bottleneck
is selection of compatible \(W_3\) digits, rather than construction of the
bounded correction quotient.

For this fixed \(G_{19}^{0}\), the full error

\[
 E_{19}^{0}=\frac{\det DG_{19}^{0}-5}{8}\pmod2
\]

has 1,027 monomials and degree 52.  A completely unrestricted next-correction
search is affine over \(\mathbb F_2\).  If the correction is \((R,S,T)\), its
linearized determinant is

\[
 dR\wedge dQ\wedge dz+dP\wedge dS\wedge dz+dP\wedge dQ\wedge dT.
\]

The coefficient equations split by \(z\)-degree.  Exact elimination at
degree 51 has 10,255 active variables and 5,791 equations and is
inconsistent.  At degree 52 the system has 10,663 variables and 5,972
equations and is consistent.  One returned correction has supports

\[
 |\operatorname{supp}R|=314,\qquad
 |\operatorname{supp}S|=98,\qquad
 |\operatorname{supp}T|=674,
\]

and direct differentiation gives constant determinant 13 modulo 16.
Consequently

\[
 \boxed{d_4(G_{19}^{0})=52}.                  \tag{8.30}
\]

Now return to the block-triangular degree-25 representative
\(G_{25}=\widetilde F_3\) in (8.12), with the exact integral \(K\) from
\(\det D(P_0,Q_0)=1+2K\), and define

\[
 L=\frac{\det DG_{25}-1}{8}\pmod2.             \tag{8.31}
\]

Here \(L\in\mathbb F_2[x,y]\) is independent of \(z\), has 35 monomials and
degree 51.  Thus its three-dimensional top de Rham class vanishes, as it
must: \(L=\partial_z(zL)\).  Vanishing of the cohomology obstruction does not
yet control the degree of a primitive.

Suppose \(G_{25}+8(R,S,T)\) has constant determinant and maximum coordinate
degree at most \(d\).  Write \(R_0,S_0\) for the \(z^0\)-coefficients and
\(t_1\) for the coefficient of \(z\) in \(T\).  The necessary \(z^0\)
equation at the next digit is, up to an irrelevant constant,

\[
 L+\mathcal D_F(R_0,S_0)+t_1=0,               \tag{8.32}
\]

where \(\deg R_0,\deg S_0\le d\) and \(\deg t_1\le d-1\).  Introduce

\[
 \Lambda_4(H)=[x^{39}y^{12}]H+[x^{40}y^{13}]H. \tag{8.33}
\]

Exact coefficient bookkeeping gives the all-degree identity

\[
 \boxed{\Lambda_4(\mathcal D_F(R,S))
 =[x^{41}y^{13}]R},                            \tag{8.34}
\]

while the exact error (8.31) satisfies

\[
 \Lambda_4(L)=1.                               \tag{8.35}
\]

If \(d\le51\), the right side of (8.34) vanishes because its source monomial
has degree 54, and \(\Lambda_4(t_1)=0\) because \(\deg t_1\le50\).  Constants
are also killed by \(\Lambda_4\).  Applying the functional to (8.32)
contradicts (8.35).  Hence no constant-Jacobian extension of this fixed
representative has degree at most 51.

At degree 52 the transparent primitive works:

\[
 \boxed{G_{25}+8(0,0,zL)\pmod {16}}.           \tag{8.36}
\]

Its full three-by-three determinant is one modulo 16.  Therefore the exact
fixed-representative minimum is

\[
 d_4(G_{25})=52.                               \tag{8.37}
\]

Reduction to \(W_3\) and `HKM2W3` gives \(d_4\ge19\), while (8.36) gives
\(d_4\le52\).  The important general lesson is that a Witt-degree search is
not only a levelwise minimization problem.  One must keep a frontier whose
state records at least current degree, support, and the cohomology class of
the next determinant digit: a level-optimal representative can be
nonextendable even when a slightly larger representative continues.

## 9. Consequences and boundary of the result

The obstruction closes the direct mixed-characteristic route for this exact
plane map.  It is stronger than the earlier observation that the displayed
integer formulas themselves are not Keller: allowing arbitrary higher-degree
corrections divisible by two does not help.

It now also closes the proposed escape through polynomial plane
left--right equivalence: the obstruction is the functorial de Rham class
\(xyu^2\), not a coordinate-specific coefficient accident.  It does not
obstruct unrelated positive-characteristic Keller maps.

Stabilization at all finite Witt levels is decided by (7.1), and `HKM2W2`
settles the exact minimum degree at level two, while `HKM2W3` proves (8.26)
at the next level.  `HKM2W4` proves the fixed-representative separation in
Section 8.3 and the unrestricted bounds \(19\le d_4\le52\), but not the exact
value of \(d_4\).  Polynomial algebraization over \(\mathbb Z_2\) is not
decided.  The remaining quantitative questions include the exact level-four
degree, later degrees, asymptotic lower bounds, compatible low-support
frontiers, and impossibility of a uniformly bounded-degree or algebraic
inverse-limit family.

## 10. Exact reproduction

```bash
.venv/bin/python scripts/verify_huq_kuruvilla_plane_w2_obstruction.py
.venv/bin/python scripts/verify_huq_kuruvilla_w3_degree19_lower.py
.venv/bin/python scripts/verify_huq_kuruvilla_w3_degree19_lower.py \
  --degree 19 --seek-lift --timeout-ms 300000
.venv/bin/python scripts/verify_huq_kuruvilla_w3_degree19_lower.py \
  --degree 19 --seek-lift --require-w4-class-zero --timeout-ms 300000
.venv/bin/python scripts/verify_huq_kuruvilla_w3_degree19_lower.py \
  --degree 19 --seek-lift --require-w4-degree 19 \
  --w4-seed-cuts 32 --w4-cut-limit 64 --timeout-ms 600000
.venv/bin/python scripts/verify_huq_kuruvilla_w4_extension_obstruction.py
.venv/bin/python scripts/verify_huq_kuruvilla_w4_extension_obstruction.py \
  --search-degree19-extension
```

The checker computes the integral Jacobian error, its full odd--odd de Rham
representative \(xyu^2\), and the old \(xy\)-coefficient witness.  It checks
the first-variation/exact-form identity, the monomial description of the full
cokernel, representative-independence on dense and monomial correction
regressions, and the explicit stabilized primitive \(zK\).  It also verifies
the all-degree dual-functional identity (8.6), the sharp degree-eighteen
stable \(W_2\) lift, the degree-twenty-five \(W_3\) construction (8.12), the
universal geometric-series induction step, and the first finite levels of
the stable Witt tower.  The all-level claim is the exact identity (7.2), not
an inference from that regression.  The checker performs no bounded
correction or left--right search.

The second command constructs the unrestricted degree-18 necessary system
from (8.18)--(8.21), checks its fixed dimensions of 1,083 variables and
1,639 constraints, and asks the pinned Z3 solver for exact Boolean
satisfiability.  Its `unsat` result proves the lower bound \(d_3\ge19\).
The third command constructs the complete degree-19 system in (8.23), obtains
an exact SAT model, and directly replays its full polynomial Jacobian to prove
the matching upper bound.  Adding `--show-model` prints all support lists.
The direct replay is computationally separate from the coefficient-system
construction, but this is still one implementation and not external review.

The preferred frozen witness has a fast solver-free replay:

```bash
.venv/bin/python scripts/verify_huq_kuruvilla_w3_degree19_lower.py \
  --replay-certificate \
  artifacts/generated-results/huq_kuruvilla_w3_degree19_witness_first_280_second_160.json
```

The initial support-reduction chain uses `--minimize-second` to reach support
180.  The affine audit then proves the three-block decomposition, and
`--second-support-layer z2` and `z1` produce intermediate certificates, while
`--component-minimize-second` proves the exact component minima for any fixed
first digit.  A full nonlinear `--first-support-bound 280` search supplies
the preferred first digit; component decoding then proves minima \(92,54,14\)
and writes the pinned support-160 completion.  The full commands and hashes
are in [REPRODUCE.md](../REPRODUCE.md).  This proves the minimum only after
fixing that first digit, not a global sparsity theorem over all first
corrections.

The fourth command compiles all 241 next-class equations.  The preferred
first digit gives UNSAT when fixed, but the unrestricted degree-19 ansatz is
SAT and produces the pinned 818-term \(G_{19}^{0}\).  The fifth command runs
the joint master--subproblem loop of (8.29a); its bounded result is the
64-singleton-cut experiment described above, not an UNSAT proof.  The sixth
command independently reconstructs the preferred degree-19 and canonical degree-25
\(W_3\) representatives modulo 16.  For the former it pins and checks the
48-term top de Rham class.  For the latter it pins the 35-term error \(L\),
extracts the two-equation dual certificate (8.33), compiles the full
4,082-variable degree-51 affine \(z^0\) system as a regression, and directly
verifies the degree-52 determinant-one lift (8.36).  The final command
reconstructs \(G_{19}^{0}\), proves its full degree-51 correction system
inconsistent, finds a degree-52 correction, and directly replays its
constant determinant modulo 16.  The de Rham and dual-functional arguments
are exact and unrestricted in the omitted correction layers; each
degree-52 minimum is restricted only by holding the chosen \(W_3\)
representative fixed.
