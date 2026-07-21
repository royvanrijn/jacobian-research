# An easy degree-five test witness for the Hessian-root invariant

This note gives a deliberately simple pair of degree-five weighted Keller maps
for testing the stable Hessian-root invariant.  Both maps are defined over
`Q`, both have Jacobian determinant `-1`, and one of them has an explicit
rational collision.  Their stable inequivalence is visible from the single
calculation that one Hessian polynomial is squarefree and the other is not.

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

Thus the weighted normalization constant is `c=-1`.  The remaining
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
polynomials.  Denote the resulting maps by

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

which has the two distinct simple roots `W=1` and `W=2`.  Reconstruction gives

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

The extra factor `2W^2-2W+1` has discriminant `-4`, and

\[
 \operatorname{Res}_W
 \left(2W^2-2W+1,-W+H_{\mathrm{dbl}}'(W)\right)=16.
\]

Thus its two additional primitive-root branches are distinct and satisfy the
boundary-clean condition used in the weighted boundary theorem.

## 3. The one-line stable distinction

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

while the second has the nonreduced divisor

\[
 \operatorname{div}(H_{\mathrm{dbl}}'')
 =2[1/2]+[1/5].
\]

On the normalization of the intrinsic discriminant open, the zeroth Fitting
ideal of the relative differential module is

\[
 \operatorname{Fitt}_0
 \Omega_{\widetilde D_H/D_H}=(H''(r)).
\]

Stable left--right equivalence preserves this finite divisor scheme up to an
affine change of the normalization coordinate.  In particular it preserves
reducedness.  Therefore

\[
 \boxed{F_{\mathrm{red}}\not\sim_{\mathrm{stable}}F_{\mathrm{dbl}}.}
\]

This is the smallest practical regression for the new theorem: no cross-ratio
calculation is needed.  The test is simply

```text
discriminant(H_red'') != 0
discriminant(H_dbl'') == 0
```

## 4. Exact reproduction

Run

```bash
.venv/bin/python scripts/verify_easy_hessian_witness.py
```

The checker verifies admissibility, polynomiality, both constant Jacobians,
generic degree five, the rational collision, the boundary-clean resultant,
and the reduced-versus-nonreduced Hessian divisor obstruction.
