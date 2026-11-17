# Cubic and quintic scalar string operators over algebraic curves

**QSO35.** Let F be a one-variable algebraic function field over an
algebraically closed characteristic-zero constant field C, with nonzero
derivation delta. There are no finite operators P,Q in F[D;delta] with
\[
 \operatorname{ord}P=3,\quad \operatorname{ord}Q=5,\quad
                         [P,Q]\in C^\times.              \tag{1}
\]
The exclusion descends to arbitrary characteristic-zero constants and
permits every finite algebraic coefficient cover.

Every lower constant flow is retained. The argument gives no
order-(3,7) obstruction and does not classify all cubic mates.
Together with QSO2 and QSO34, it raises the smallest possible generic
orders of a prescribed exact-R cubic Long repair to (3,7).

## Finite normalization and every lower flow

As in [QSO34's normalization](CUBIC_QUARTIC_STRING_OBSTRUCTION.md),
adjoin a cube root of the leading coefficient of P, change the
differential generator, and shift that generator by a coefficient.
This gives
\[
                 P=E^3-\tfrac32UE-\tfrac34U'+\tfrac32V.   \tag{2}
\]
Primes denote the resulting nonzero derivation epsilon. These are
finite Ore changes over an algebraic function field. The leading
coefficient of Q is constant by the highest commutator coefficient;
divide it out, retaining a nonzero constant commutator kappa.

Use the following explicitly finite operators:
\[
 \begin{split}
 B_1={}&E,\qquad B_2=E^2-U,\\
 B_4={}&E^4-2UE^2+2(V-U')E+V'+\tfrac12U^2-\tfrac56U'',\\
 B_5={}&E^5-\tfrac52UE^3+\tfrac52(V-\tfrac32U')E^2\\
       &+\tfrac54(U^2-\tfrac73U''+2V')E
       +\tfrac54(\tfrac43V''+UU'-\tfrac23U'''-2UV).
 \end{split}                                             \tag{3}
\]
They are the positive fractional powers of P, but no infinite
operator is needed below. Direct multiplication makes each
[P,B_j] have order at most one.

Subtracting B_5 from Q leaves order at most four. Its leading
coefficient must be constant, since a nonconstant leading
coefficient f at any order n>=0 contributes 3f' at order n+2
in its commutator with the monic P. Subtract a constant multiple
of B_4 and repeat with P, B_2, B_1 and 1. The order strictly
decreases at each step. Subtracting the constant multiple of P
and the scalar term preserves the string equation. Thus the full
normal form is
\[
                    Q=B_5+c_4B_4+c_2B_2+c_1B_1,
                        \qquad c_4,c_2,c_1\in C.          \tag{4}
\]
No lower flow has been silently set to zero.

The formula for B_5 agrees with
[Hayford, arXiv:2405.03260v2, Appendix A.2](https://arxiv.org/html/2405.03260v2#A2).
This is a comparison of operator formulas. The algebraic
nonexistence theorem here is not imported from that paper.

## The string equation creates its own time coordinate

Put
\[
 \begin{split}
 A={}&\tfrac12V''-\tfrac32UV,\\
 B={}&\tfrac1{12}U^{(4)}-\tfrac34UU''-\tfrac38(U')^2+
                              \tfrac32V^2+\tfrac12U^3,\\
 H_1={}&-\tfrac16U^{(4)}+\tfrac54UU''+\tfrac{15}{16}(U')^2
                              -\tfrac58U^3-\tfrac{15}{4}V^2,\\
 H_2={}&\tfrac16V^{(4)}-\tfrac54UV''-\tfrac58U'V'
                              -\tfrac58U''V+\tfrac{15}{8}U^2V,\\
 J_1={}&H_1-2c_4A-3c_2V+\tfrac32c_1U,\\
 J_2={}&H_2-c_4B+c_2(\tfrac34U^2-\tfrac14U'')-\tfrac32c_1V.
 \end{split}                                             \tag{5}
\]
Finite normal-order Ore multiplication proves the universal identity
\[
                      [P,Q]=J_1'E+\tfrac12J_1''+J_2'.    \tag{6}
\]
Hence a scalar commutator kappa!=0 forces J_1=gamma in C and
J_2'=kappa. In particular
\[
                      x=J_2/\kappa,\qquad \epsilon x=1. \tag{7}
\]
The function field is finite over C(x). Thus U,V are algebraic
functions of the constructed x and satisfy
\[
                        J_1=\gamma,\qquad J_2=\kappa x. \tag{8}
\]
No primitive was assumed before the string relation supplied it.

## Exhaustion of all finite singular balances

At a place above a finite x_0, let t=x-x_0 and use Puiseux Laurent
series. Set a=max(0,-ord_t U), b=max(0,-ord_t V), with zero for
an identically zero or bounded component.

First put R=max(1,a/2,b/3). If R>1, the only possible leading
terms of J_1 have pole order 6R and are its nonderivative U^3
and V^2 terms. Derivatives and every lower-flow term have
strictly smaller pole order. Both must occur, giving a=2R,
b=3R and u^3+6v^2=0 for nonzero principal coefficients.
But J_2 then has a unique leading U^2V term of pole order 7R:
its other terms have order at most max(3R+4,5R+2,6R)<7R.
This is impossible. Thus a<=2 and b<=3.

Suppose some pole occurs but a<2. If b=3, V^2 is uniquely
leading in J_1. If b<3 and b>a, V^{(4)} is uniquely leading
in J_2: a+4<b+4 and all its other pole orders are smaller.
If b<=a and a>0, U^{(4)} is uniquely leading in J_1.
These alternatives exhaust all poles and are impossible.
Consequently
\[
                             a=2,\qquad b\le3.            \tag{9}
\]
The derivative coefficients used here are nonzero at negative
Puiseux exponents. When a component is bounded its k-th derivative
has pole order strictly below k, which preserves all the strict
comparisons above.

Write U=u*t^{-2}+... and V=v*t^{-3}+..., allowing v=0.
The weight-six and weight-seven leading equations are
\[
 u(u-2)(u-16)+6v^2=0,\qquad v(u-4)(u-8)=0.
\]
Since u is nonzero, the complete list is
\[
                (2,0),(16,0),(4,4),(4,-4),(8,8),(8,-8).   \tag{10}
\]

## Integer resonance roots exclude all finite branching

Consider the first noninteger index r>0 in the expansions
U=sum u_r t^{r-2}, V=sum v_r t^{r-3}. Linearizing the principal
parts H_1,H_2 at a pair (u,v) in (10) gives the coefficient matrix
\[
 M_{u,v}(r)=
 \begin{pmatrix}
 A_u(r)&-15v/2\\
 -5v(r^2-8r+36-6u)/8&D_u(r)
 \end{pmatrix},                                         \tag{11}
\]
where
\[
 \begin{split}
 A_u(r)&=-\tfrac16(r-2)(r-3)(r-4)(r-5)
             +\tfrac54u(r^2-8r+18)-\tfrac{15}{8}u^2,\\
 D_u(r)&=\tfrac16(r-3)(r-4)(r-5)(r-6)
             -\tfrac54u(r^2-8r+18)+\tfrac{15}{8}u^2 .
 \end{split}
\]
The coefficients at t^{r-6},t^{r-7} in (8) say
M_{u,v}(r)*(u_r,v_r)^t=0. The c_4,c_2,c_1 terms have lower
weights, so their fractional contributions use earlier indices
r-1,r-3,r-4 and are zero. The two forcing terms only use integer
indices. Nonlinear contributions at the first fractional index
are exactly the displayed linearization.

The determinants factor as follows; each entry displays -36*det M.

| Principal pair | Factorization |
|---|---|
| (2,0) | (r-9)(r-7)(r-5)^2(r-3)^2(r-1)(r+1) |
| (16,0) | (r-15)(r-12)(r-10)(r-9)(r+1)(r+2)(r+4)(r+7) |
| (4,+/-4) | (r-10)(r-9)(r-6)(r-5)(r-3)(r-2)(r+1)(r+2) |
| (8,+/-8) | (r-12)(r-10)(r-9)(r-5)(r-3)(r+1)(r+2)(r+4) |

All roots are integers, including the repeated roots. Thus M is
invertible at the noninteger r. This excludes its first occurrence.
No resonance-compatibility assertion or full analytic Painleve
property is required for this argument about an already assumed
algebraic solution.

At a finite place where U,V are both bounded, consider the first
fractional exponent rho>0 occurring in either. If rho<4, the
coefficient of t^{rho-4} comes from the fourth derivatives alone.
Their constant coefficient matrix in (8) is
\[
              \begin{pmatrix}-1/6&0\\-c_4/12&1/6\end{pmatrix},
\]
which is invertible. All lower-derivative nonlinear terms are
less singular, so the corresponding fractional coefficients
must vanish. Hence the jets through third order are bounded.
Solving (8) for U^{(4)},V^{(4)} gives a regular polynomial
first-order system for eight jet variables. Its formal solution
at any finite initial jet has integer powers, as its first
fractional coefficient would appear one exponent earlier after
differentiation than on the polynomial right side. There is
therefore no finite branching at bounded points either.

## Rationality and the infinity contradiction

Let J=C(x,U,V). All its places over finite x are unramified by
the preceding expansions. The curve argument of
[QSO34](CUBIC_QUARTIC_STRING_OBSTRUCTION.md) applies:
Riemann--Hurwitz gives 2g-2=-d-b for d=[J:C(x)] and b>=1
points above infinity, forcing d=b=1. Thus U,V are rational.

If either has positive degree at infinity, set
R=max(deg(U)/2,deg(V)/3)>0, allowing a bounded component degree
zero in this maximum. The leading terms of J_1 again force both
U^3 and V^2 to have degree 6R; all derivatives and lower flows
have smaller degrees. Thus deg(U)=2R and deg(V)=3R.
The term U^2V of J_2 has degree 7R, above all other terms of J_2.
Its equation J_2=kappa*x would force 7R=1. This is incompatible
with the integral rational-function degrees 2R and 3R.

If both U,V are bounded, all of J_2 is bounded and cannot equal
kappa*x. Identically zero components are included. This proves (1).

## The next possible orders for a prescribed Long repair

Work in the exact-R Laurent chart with h-regular finite S,T,
[S,T]=h, and the prescribed classical pair of a-orders two
and three. Its leading coefficients are -3r and r/2 over K(r).
If generic S has a-order three, then
\[
                        \operatorname{ord}_a T\ge7.       \tag{12}
\]

Suppose instead its order were at most six. Whenever the current
mate has order divisible by three, leading commutator comparison
makes its leading coefficient ratio with the appropriate power
of S a differential constant. Subtract that constant power,
strictly lowering the order. Order two is excluded by QSO2 with
S as its odd-order mate; orders four and five by QSO34 and (1);
order zero cannot have a nonzero scalar commutator with a cubic.
The remaining possibility is therefore
\[
                   T=A(S)+T_1,\qquad \operatorname{ord}T_1=1,
                          \quad \deg A\le2.              \tag{13}
\]
Absorb A's constant term into T_1. Every nonconstant coefficient
of A is h-regular: otherwise its most negative h coefficient
would reduce to a nonconstant polynomial in the classical S,
of a-order at least two, which cannot cancel an order-one T_1.
Thus T_1 is also h-regular. The reduction of (13) has a-order
at most one or positive even order, since the classical S has
order two. It cannot have the prescribed order three.
This proves (12). It uses actual regular coefficients and
prescribed classical orders; it does not constrain an arbitrary
polar generic mate without that data.

Generic order seven and higher, genuinely deformed first images,
and nonscalar/projective carriers remain open. There is no
bounded finite Long construction or DC1/DC2 resolution here.

## The all-order cubic reconstruction obligation

Further progress should seek a uniform theorem for all cubic partner
orders, rather than extend the resonance table one order at a time.
There is a precise source boundary for an isomonodromic route.

[Miwa's 1981 Theorem 1](https://www.jstage.jst.go.jp/article/kyotoms1969/17/2/17_2_703/_pdf/-char/en)
gives a meromorphic inverse-monodromy map, but Section 2, equation
(2.22), assumes regular exponents distinct modulo integers.
A cubic spectral ramification can introduce a resonant apparent
regular pole, so that hypothesis needs attention.

[Malgrange's 1982 Theorem 3.1 and its conclusion](https://numdam.org/item/CIF_1982__17__A3_0.pdf)
allow logarithmic regular poles without that extra restriction
and produce meromorphic universal connection forms under the
stated irregular-leading hypotheses. Yet his question (3.2),
printed page 11, explicitly distinguishes this from meromorphic
generalized Schlesinger coordinates in the resonant case.

The outstanding application is: for every partner order m coprime
to three and every lower constant flow, admit the normalized
rank-three spectral connection uniformly over finite time, and
reconstruct its scalar coefficients U,V single-valuedly and
meromorphically from the universal family. This must cover the
resonant apparent pole, degenerate monodromy data, all lower flows,
and the gauge normalization; absence of finite fixed singularities
must also be proved. Meromorphic connection forms alone do not
establish this reconstruction.

This is a named application obligation, not a claim that the
general theorem is unknown in modern literature. The direct proofs
of QSO34 and QSO35 do not depend on it. The
[compact primary-source audit](cubic_hierarchy_source_audit.json)
pins both PDFs and the exact inspected pages.

## Evidence

The standalone checker directly verifies the finite normal forms
and the full all-lower-flow commutator identity. It rederives the
six principal balances and the resonance matrices by symbolic
linearization, factors the determinants, verifies the bounded-jet
leading matrix, rejects a corrupted B_5 coefficient and validates
the primary PDF pin. These exact controls support the written
all-place argument; they do not prove it by a finite scan.

Run from the repository root:

```sh
.venv/bin/python research/dixmier/dixmier-formal-lifting/check_cubic_quintic_string.py
```
