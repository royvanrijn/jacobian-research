# An easy degree-five test witness for the Hessian-root invariant

This note gives a deliberately simple pair of degree-five weighted Keller maps
for testing the stable Hessian-root invariant. Both maps are defined over
`Q`, both have Jacobian determinant `-1`, and one of them has an explicit
rational collision. Their stable inequivalence is visible before any
cross-ratio calculation: the intrinsic Hessian-root support has three points
for one map and two points for the other.

## 1. The two primitives

Put

\[
 H_{\mathrm{red}}(W)=\frac12W^2(W-1)(W^2+1)
\]

and

\[
 H_{\mathrm{dbl}}(W)=W^2(W-1)(2W^2-2W+1).
\]

For both primitives,

\[
 H(0)=H'(0)=H(1)=0,\qquad -H'(1)=-1.
\]

Thus the weighted normalization constant is `c=-1`. The remaining
admissibility data are

\[
 \frac{H_{\mathrm{red}}''(1)}c=-6,
 \qquad a_{\mathrm{red}}=-\frac54,
\]

and

\[
 \frac{H_{\mathrm{dbl}}''(1)}c=-8,
 \qquad a_{\mathrm{dbl}}=-\frac76.
\]

For either primitive define

\[
 u=1+xy,\qquad
 \gamma=1+a_Hxy+x^2z,\qquad
 W=u\gamma,
\]

\[
 C=x\gamma,
\]

\[
 B=\frac{-1+H'(W)/\gamma}{x},
\]

\[
 A=\frac{u-\bigl(WH'(W)-H(W)\bigr)/\gamma^2}{x^2}.
\]

The weighted-seed polynomiality calculation shows that these quotients are
polynomials. Denote the resulting maps by

\[
 F_{\mathrm{red}},F_{\mathrm{dbl}}:\mathbb A^3\longrightarrow\mathbb A^3.
\]

Both satisfy

\[
 \det DF=-1
\]

and have marked inverse equation

\[
 H(W)-BCW-AC^2=0.
\]

Consequently both have generic fiber degree five.

## 2. A rational collision for the double-Hessian witness

For `H=H_dbl`, take the target

\[
 (A,B,C)=(-20,20,1).
\]

The inverse polynomial is

\[
 H_{\mathrm{dbl}}(W)-20W+20,
\]

which has the two distinct simple roots `W=1` and `W=2`. Reconstruction gives

\[
 P_1=\left(-\frac1{19},20,-\frac{22990}{3}\right),
\]

\[
 P_2=\left(\frac1{44},-42,81092\right).
\]

Direct substitution gives

\[
 F_{\mathrm{dbl}}(P_1)=F_{\mathrm{dbl}}(P_2)=(-20,20,1).
\]

Hence `F_dbl` is an entirely rational explicit Keller counterexample.

## 3. Boundary-clean audit

The residual primitive factors are

\[
 R_{\mathrm{red}}(W)=W^2+1,
 \qquad
 R_{\mathrm{dbl}}(W)=2W^2-2W+1.
\]

Both have discriminant `-4`, so their two additional primitive roots are
distinct. Moreover,

\[
 \operatorname{Res}_W
 \left(R_{\mathrm{red}},-W+H_{\mathrm{red}}'(W)\right)=1,
\]

and

\[
 \operatorname{Res}_W
 \left(R_{\mathrm{dbl}},-W+H_{\mathrm{dbl}}'(W)\right)=16.
\]

Thus all additional primitive-root branches satisfy the boundary-clean
condition. For each map, the complete canonical target boundary is therefore
the intrinsically ordered pair consisting of the ramified discriminant vertex
and `C=0`, exactly as required by the stable Hessian-root theorem.

## 4. The one-line stable distinction

The Hessian polynomials are

\[
 H_{\mathrm{red}}''(W)=10W^3-6W^2+3W-1,
\]

\[
 H_{\mathrm{dbl}}''(W)=2(2W-1)^2(5W-1).
\]

The first has discriminant

\[
 \operatorname{Disc}(H_{\mathrm{red}}'')=-1080\ne0,
\]

so its Hessian-root support consists of three distinct points. The second has
support

\[
 \{1/2,1/5\}
\]

and scheme-theoretic divisor

\[
 \operatorname{div}(H_{\mathrm{dbl}}'')=2[1/2]+[1/5].
\]

On the normalization of the intrinsic discriminant open, the zeroth Fitting
ideal of the relative differential module is

\[
 \operatorname{Fitt}_0
 \Omega_{\widetilde D_H/D_H}=(H''(r)).
\]

Stable left--right equivalence preserves this module and its support after
stabilization. Hence it preserves the number of irreducible Hessian-root
components. Here that number is three versus two, which is already impossible.
The full Fitting divisor additionally records the doubled point. Therefore

\[
 \boxed{F_{\mathrm{red}}\not\sim_{\mathrm{stable}}F_{\mathrm{dbl}}.}
\]

This is the smallest practical regression for the new theorem. The essential
test is simply

```text
degree(squarefree_part(H_red'')) == 3
degree(squarefree_part(H_dbl'')) == 2
```

with the discriminants providing an equivalent quick check.

## 5. Exact reproduction

Run

```bash
.venv/bin/python scripts/verify_easy_hessian_witness.py
```

The checker verifies admissibility, polynomiality, both constant Jacobians,
generic degree five, the rational collision, both boundary-clean resultants,
and the three-component-versus-two-component Hessian-root obstruction.
