# Complete mixed-weight (2,-1,-3) exclusion through degree seven

**Theorem (KCUB17).** Every three-variable polynomial Keller map over a
characteristic-zero field, with source/output weights (2,-1,-3) and ordinary
degree at most seven, is a polynomial automorphism. This excludes every
support size in the complete 16-slot class, including the original raw<=12
construction query. The separate (1,-1,-3) class and global simplicity minima
remain open. No smaller construction is claimed.

The distinct weights make the linear part diagonal and invertible; output
row scaling normalizes it to identity and det JF=1 without changing degree
or nonlinear occurrences. We work in that normalization and extend to an
algebraic closure for the boundary arguments.

## A double-cover quotient, not a polynomial invariant chart

The integral invariant ring has a relation, so it cannot simply be treated
as two independent polynomial coordinates. On a nonzero-a double cover put

\[
a=s^2,\quad b=u/s,\quad c=v/s^3,\qquad
F=(s^2 f,s^{-1}g,s^{-3}h).
\]

Every allowed atom gives

\[
\begin{aligned}
f&=1+\mu u^2+\nu uv+P u^4+T v^2,\\
g&=u+L v+A u^3+B u^2v+C u^5+Duv^2+Ev^3,\\
h&=v+M u^3+N u^2v+O u^5+Quv^2+R u^4v+Sv^3.
\end{aligned}
\]

There are four first-output atoms and six in each tail. The quotient f is
even, and g,h are odd under iota:(u,v)->(-u,-v). Direct differentiation,
including the factor two from a=s^2, gives

\[
J=f(g_u h_v-g_v h_u)+\tfrac12g(f_u h_v-f_v h_u)
 +\tfrac32h(g_u f_v-g_v f_u)=1. \tag{1}
\]

The independent audit reconstructs the universal original three-variable
determinant and checks (1) against it coefficient by coefficient.

## The smooth quartic boundary has paired infinities

Any f=f_u=f_v=0 point contradicts (1); hence a Keller completion requires
f=0 to be smooth. Suppose PT!=0. Set y=2Tv+nu*u. Completing the square gives

\[
y^2=-4PTu^4+(\nu^2-4\mu T)u^2-4T. \tag{2}
\]

A multiple root of the quartic would be a singular affine point, already
impossible. The smooth compactification has two points O+,O- at infinity.
In z=1/u,Y=y/u^2 coordinates its equation is

\[
Y^2=-4PT+(\nu^2-4\mu T)z^2-4Tz^4.
\]

Both points are smooth at z=0, since Y is nonzero. The differential
omega=du/y=-dz/Y is regular and nonzero there. It is also regular and
nonzero at every finite point: at a simple quartic root use y as parameter,
and elsewhere use u. Thus omega is nowhere zero on the complete curve.

Iota sends (u,y) to(-u,-y), and consequently exchanges O+ and O-. The odd
functions g,h have no finite poles. The nonzero restriction of g is not
constant: an odd constant would be zero, contradicting (1) on the boundary.
Hence g has a pole at infinity; iota forces equal pole orders m>0 at O+,O-.
Its finite zeros have total multiplicity 2m.

The boundary form of (1) is

\[
g\,dh-3h\,dg=-2\omega,
\qquad d(h/g^3)=-2\omega/g^4. \tag{3}
\]

At a finite zero of g of order k, (3) makes h/g^3 have a pole of order
4k-1. Regular h allows at most 3k, so k=1. Thus all 2m finite zeros are
simple and rho=h/g^3 has precisely 2m triple poles, total degree 6m.
At each infinity, d rho has zero order 4m. Therefore rho is regular and
rho-rho(O) has zero order 4m+1. Since h and g are odd, rho is iota-invariant;
its two infinity values agree. Subtracting this common value gives total
zero degree at least 8m+2 but pole degree 6m, impossible on a smooth complete
curve. This proves **PT=0**.

The equality of infinity values is essential; the two-infinity case is not
a direct reuse of the earlier single-infinity argument. This is a written
curve/divisor proof, not proof-assistant formalization.

## Two exact identities describe all residual parameters

Saved coefficients mu,nu,P,T are c1,c2,c3,c4. Explicit polynomial multipliers
in the full tail-coefficient ideals give

\[
\begin{array}{c|c}
P=0 &(\mu T-\nu^2/4)^4=\sum M_i I_i\\
T=0 &\nu^6=\sum N_i I_i.
\end{array}
\]

All 56 and 13 multiplier monomials are independently multiplied. These
cleared targets cover all parameter values in their branches; no generic
exceptional factor is discarded.

Thus T=0 forces nu=0 and f=1+mu*u^2+Pu^4. If T!=0, PT=0 and the first
identity give P=0,mu=nu^2/(4T). With sigma=nu/(2T),
f=1+T(v+sigma*u)^2. We now exclude every nonconstant f in these two forms.

## Complete transverse bounds and parity

On a smooth univariate boundary phi(L)=0 the restricted tails satisfy
G H'-3G'H=nonzero constant. Both source tails have transverse degree at
most five. Leading degrees can cancel only when deg H=3deg G; within
that bound the nonconstant balance is (1,3). Otherwise their degree sum
must be one. Hence deg G<=1,deg H<=3.

For f=1+mu*u^2+Pu^4 take L=u,t=v. For
f=1+T(v+sigma*u)^2 take L=v+sigma*u,t=u. A repeated root contradicts
smoothness. Every above-threshold coefficient has degree less than two in
L and vanishes on at least two boundary lines, so it vanishes globally.
The complete linear systems have ranks two and five respectively.

Both cases therefore have

\[
g=A(L)t+B(L),\quad
h=Ct^3+D(L)t^2+E(L)t+H(L),
\]

where A is even of degree<=2, B,H odd, C constant, D an odd linear polynomial,
and E even. In the first case E(0)=1. In the second case
E(0)=-sigma,A(0)=1-sigma*L0,H'(0)=1, where L0 is the original coefficient
of v in g. All of these source facts and complete linear parameterizations
are checked independently, without inverting sigma.

Up to an overall sign, the determinant equation is

\[
\mathcal J=\tfrac12\phi'(g h_t-3h g_t)
       +\phi(g_Lh_t-g_th_L)=\epsilon,\quad \epsilon\in\{1,-1\}. \tag{4}
\]

If A=0 identically, (4) factors as
(phi B'+phi'B/2)*h_t=epsilon. Nonzero odd B and degree>=2 phi make the first
factor nonconstant, impossible. Thus A is nonzero.

If C!=0, the t^3 coefficient is 3phi*A'*C, so A is a nonzero constant.
The polynomial target shear h->h-(C/A^3)g^3 removes C and preserves (4),
parity and the linear part. It may increase degree, but no degree bound on
the modified h is used. The t^2 equation then forces (phi D^2)'=0, hence
D=0. The t equation gives (phi E)'=0. If E(0)!=0 this is impossible.
Otherwise it forces E=0; in the parallel case this means sigma=0 and
H'(0)=1. The remaining equation has the factor
phi H'+3phi'H/2, whose leading degree is positive, so it cannot give a
nonzero constant. This excludes C!=0.

Now C=0 in the original source bounds. If D!=0, the t^2 equation gives
phi D^2=k A^4 for a nonzero constant k. At L=0 its left side has order two,
while the right side has order zero or eight, since nonzero even A has
degree at most two. Contradiction. Hence D=0.

The t equation next gives phi E=kA. If k=0 then E=0, ruled out by E(0)=1
in the first case. In the parallel case it again leaves the nonconstant
factor phi H'+3phi'H/2. If k!=0, degree phi=4 is impossible since deg A<=2.
For degree phi=2 it forces E=e to be a nonzero constant and A=e*phi/k.
Putting Q=B-phi H/k, the constant equation becomes

\[
e\{\phi Q'+\tfrac12\phi'Q\}=\epsilon.
\]

For any nonzero polynomial Q of degree n, the braces have degree n+1 with
nonzero leading coefficient. Zero Q gives zero. Neither is a nonzero
constant. This excludes both residual nonconstant forms of f.

The parity argument is written proof; the audit checks its coefficient
identities, degree/parity hypotheses and complete source reductions.

## The linear-first-component family is polynomially invertible

We have proved mu=nu=P=T=0, so F1=a. Over K(a), the two remaining components
are a plane Keller pair of ordinary b,c-degree at most five. The prime
cases 2,3,5 and the degree-four case of
[Nakai-Baba Theorem3, printedp408](https://ir.library.osaka-u.ac.jp/repo/ouka/all/8846/ojm14_02_08.pdf),
plus the elementary linear case, give a polynomial inverse over K(a).
This is an attributed import, not a re-proved theorem; see
[the degree-five import record](weighted-213/plane-degree5-import.json).

The constant terms in b,c are zero and the linear matrix over K[a] is
[[1,c5*a],[0,1]], with determinant one. Therefore the unique formal inverse
in b,c has coefficients in K[a], obtained recursively using this unimodular
linear matrix. It equals the polynomial inverse over K(a) by uniqueness.
Every coefficient consequently lies in K[a], and all higher coefficients
vanish, so the inverse is a polynomial globally, including at a=0.
Together with the unchanged first coordinate this proves the theorem.

## Exact evidence and retained limits

Run

```sh
.venv/bin/python research/jacobian-followup/cubic-reduction/verify_weighted_213.py
```

[Validation](weighted-213/validation.json) reconstructs all 16 source slots,
the original determinant and double-cover quotient, even/odd parity, the
even completed quartic, both exact identities, complete transverse linear
reductions, source parity/degree bounds, univariate coefficient identities,
and the generic plane family's degree/unimodular linear part. The two used
certificates contain 69 multiplier monomials, maximum multiplier degree seven.
Curve theory, the parity/divisibility proof and the plane theorem are
written/imported, not formalized or independently re-proved by this audit.

The initial unrestricted coefficient basis remains timeout20 seconds and
the raw<=12 nonlinear-first real search remains UNKNOWN30 seconds. Changed
parameter-field bases and eight concrete low-support/rank controls all give
primary units; these finite controls are not the proof. A generic full-field
certificate with an extra parameter factor is retained but unused; the
paired-infinity lemma supplies PT=0 instead of assuming that factor away.
Only the two residual identities above enter the theorem.

This closes (2,-1,-3) through degree seven. The subsequent [KCUB18 proof](WEIGHTED_113_DEGREE7_EXCLUSION.md) closes
the original (1,-1,-3) class too. No construction bound changes: raw13 in the separate
(1,-1,-2) sector, ordinary cubic12/36, cubic-linear89 and normalized3D[10,48]
remain distinct current bounds. The full global simplicity objective remains
open.
