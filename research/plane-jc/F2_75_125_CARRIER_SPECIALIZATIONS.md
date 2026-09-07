# Exact carrier specializations of the F2 lower-Laurent maps

## Result and claim boundary

The checker
[`cas/verify_f2_75_125_carrier_specializations.py`](cas/verify_f2_75_125_carrier_specializations.py)
connects the finite carrier classification to the carried lower-Laurent
presentation.  It makes an essential component distinction before doing any
coefficient substitution.

The two double-root carrier points belong to the descent-eight component and
specialize over

\[
K=\mathbf Q[\rho]/(\rho^2-3\rho+1). \tag{1}
\]

Every currently exposed zero-row, target, and layer-zero linear map is now
compiled over `K`, with exact ranks and canonical kernel/cokernel basis
digests.

The rational squarefree carrier does **not** belong to that component.  Its
cofactor

\[
R_{\rm sf}(w)=\frac{w^2-3w+3}{25} \tag{2}
\]

has discriminant `-3/625`, whereas the descent-eight presentation assumes

\[
R_{\rm dbl}(w)=\frac{(w-w_0)^2}{25(1-w_0)^2}. \tag{3}
\]

Thus substituting (2) into a parameter `w0` would be mathematically invalid.
The squarefree point must branch earlier and be compiled from the later
first-defect ledger.

For the double component, this note pins the linear maps used by the nonlinear
continuation.  That continuation is now compiled in
[`F2_75_125_NONLINEAR_FORCING.md`](F2_75_125_NONLINEAR_FORCING.md).  Neither
calculation claims that the resulting circuit ideal is the unit ideal.

The pinned result is
[`../artifacts/generated-results/jc2_f2_75_125_carrier_specializations.json`](../artifacts/generated-results/jc2_f2_75_125_carrier_specializations.json).

## 1. Exact branch routing

The carrier Wronskian leaves three geometric carrier points:

1. the squarefree rational point (2);
2. two conjugate double-root points satisfying (1).

The carried endpoint reduction in
[`F2_75_125_DERIVATION.md`](F2_75_125_DERIVATION.md) was constructed on the
earliest surviving first-defect component.  Its support prime is a movable
double root of `R`, so its base is the discriminant-zero locus (3).  Equations
(2) and (3) are disjoint strata.

For (1), both localized parameters used by the endpoint compiler are units:

\[
\rho^{-1}=3-\rho,
\qquad
(\rho-1)^{-1}=\rho-2. \tag{4}
\]

Hence the parametric unit minors remain units after specialization.

The descent-eight local ratio supplies a second quadratic equation

\[
27y^2-9y+1=0. \tag{5}
\]

Its discriminant square class is `-3`, while (1) has square class `5`.
The two quadratic fields are linearly disjoint.  With `theta=rho+y`, the
compositum has the irreducible polynomial

\[
729\theta^4-4860\theta^3+10341\theta^2-7470\theta+1756. \tag{6}
\]

Thus the double carrier row produces four lower-Laurent coefficient branches:
two embeddings of `rho` times two embeddings of `y`.  The linear maps below
are already defined over `K`; `y` first enters their nonlinear forcing and
incidence rows.

## 2. Specialized zero-row Schur maps

Let

\[
C_0=t^5(w-1)^2R_{\rm dbl}(w). \tag{7}
\]

At descent `delta`, arbitrary new `P` data are tangent directions with the
follower

\[
q_{\rm follow}=-3C_0^2p. \tag{8}
\]

Modulo these followers, the new `Q` variable is controlled by

\[
T_\delta(q)=5C_0q'-(25-\delta)C_0'q. \tag{9}
\]

In the exact factored band basis, after removing the common nonzero factor,
the column `w^j` maps to three consecutive coefficient rows:

\[
T_\delta(w^j)=A_jw^j+B_jw^{j+1}+C_jw^{j+2}, \tag{10}
\]

where, with `ell=25-delta` and exact band data `(u,nu,d)`,

\[
\begin{aligned}
A_j&=\rho(5j+u),\\
B_j&=(\rho+1)(2\ell-5j-u)-5\nu\rho,\\
C_j&=5j+u+5\nu-4\ell.
\end{aligned} \tag{11}
\]

The checker constructs these matrices over the abstract field (1), rather
than substituting a floating approximation to either root.

For descents `1..11`, their total rank is `134` and the only two `Q` kernel
directions occur at descents `5,10`.  This exactly recovers the earlier
parametric upper elimination after specialization.

For the coupled zero rows

\[
\delta=12,13,\ldots,35,37, \tag{12}
\]

every nonresonant map has a two-dimensional forcing cokernel.  Descents
`15,20,25` retain one source centralizer and have a three-dimensional
cokernel.  The direct sum of these normalized quotient targets has dimension

\[
24\cdot2+3+2=53. \tag{13}
\]

This is not the full Laurent-row cokernel.  Before common-factor division, the
coupled rows have cokernel dimension `347`; the additional `294` coordinates
are divisibility and local-jet compatibility conditions.  Thus formula (13)
is the exact quotient sequence in which the propagated forcing is reduced,
while the full nonlinear compiler retains the split

\[
294+53=347. \tag{13a}
\]

## 3. Specialized target cokernel

After the endpoint substitution, put

\[
D=w^2(w-1)^5(w-\rho)^5,
\qquad E=(w-1)(w-\rho). \tag{14}
\]

The quotient operator is

\[
\mathcal N(S)
=5wE S'+(11E+22wE')S,
\qquad \deg S\le19. \tag{15}
\]

Its matrix

\[
\mathcal N:K[w]_{\le19}\longrightarrow K[w]_{\le21} \tag{16}
\]

has shape `22 x 20`, rank `20`, and cokernel dimension two.  The complete
target image

\[
D\mathcal N(K[w]_{\le19})\subset K[w]_{\le33} \tag{17}
\]

also has rank `20`, hence a fourteen-dimensional cokernel.  Twelve of its
coordinates are exactly the local jets

\[
\begin{array}{c|c}
w=0&0,1\\
w=1&0,1,2,3,4\\
w=\rho&0,1,2,3,4.
\end{array} \tag{18}
\]

The checker verifies directly that these twelve rows annihilate (17) and,
together with the two quotient residues, span its whole cokernel.

The earlier `w=0` control and fixed `w=1` endpoint elimination remove seven
of the local rows.  Therefore the specialized target remainder is

\[
5\text{ jets at }\rho+2\text{ quotient residues}=7. \tag{19}
\]

## 4. Specialized layer-zero Hermite cokernel

The layer-zero image is

\[
K\cdot1+
w^3(w-1)^6(w-\rho)^6K[w]_{\le18}
\subset K[w]_{\le33}. \tag{20}
\]

It has rank `20` and cokernel dimension `14`.  An exact Hermite presentation
of that cokernel is

\[
\begin{aligned}
&H'(0),\ H''(0),\\
&H(1)-H(0),\ H^{(j)}(1),\quad1\le j\le5,\\
&H(\rho)-H(0),\ H^{(j)}(\rho),\quad1\le j\le5.
\end{aligned} \tag{21}
\]

The checker constructs both the `34 x 20` image matrix and the `14 x 34`
Hermite matrix over `K`, verifies ranks `20` and `14`, verifies that their
product is zero, and proves that the Hermite rows equal the full left
cokernel.

After the prior `w=0` and `w=1` eliminations, the remaining layer-zero map is
the six-row block at `rho`.  Combining (19) and (21) recovers the exact
specialized residual count

\[
\boxed{7+6=13}. \tag{22}
\]

## 5. Nonlinear continuation

The downstream compiler now performs the operations that were left open in
the first version of this note:

1. represent the ten endpoint pivot solutions as source-band circuits;
2. propagate all earlier tangent-kernel parameters through descent `12`;
3. on each successive zero row, reduce the known forcing in the exact
   cokernel basis of (10);
4. work over the quartic field (6) when the defect parameter `y` enters;
5. append the seven target and six Hermite coordinates in (22);
6. retain the five descent-eight incidence equations and the quadratic field
   relation.

It produces an exact `366`-equation arithmetic-circuit presentation:
`294+53` coupled Laurent coordinates, `7+6` final functionals, five incidence
rows, and one relative quadratic relation.  Testing that localized component
ideal remains open.

For the squarefree carrier, the downstream compiler records all `82` later
spacings `9..90` and their exact pre-target multiple census.  It routes the
carrier but does not yet eliminate any of those target-and-tail branches.

Reproduce the present result with

```bash
.venv/bin/python plane-jc/cas/verify_f2_75_125_carrier_specializations.py
.venv/bin/python plane-jc/cas/test_f2_75_125_carrier_specializations.py
.venv/bin/python plane-jc/cas/compile_f2_75_125_nonlinear_forcing.py
.venv/bin/python plane-jc/cas/test_f2_75_125_nonlinear_forcing.py
```

Intentional regeneration uses `--refresh` on the first command.
<!-- status-consumer: PF2NF1 cfd1da5136c0b6d0 -->
