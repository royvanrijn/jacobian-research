# All-degree diagonal Schur and Meng--Yang frontend theorems

## Status

This note compares three all-degree frontends.  The discovery sweeps have now
been separated from four exact degree-free statements:

- `HC4FSD1` classifies the diagonal ternary reverse-Schur equation for every
  \(D\ge4\);
- `HC4FSD2` excludes constant rank-two determinant-preserving directions on
  every normalized minimal diagonal tower of degree \(D\ge5\);
- `HC4FSD3` propagates that obstruction through every compatible lower
  homogeneous layer and also excludes ranks three and four; and
- `HC4MYGJ2` gives the exact all-order Meng--Yang normal symbol and the unique
  formal graph above every plane trace.

The proof of `HC4FSD1` is canonical in
[`HC4_DIRECT_HOMOGENEOUS_FILTRATION.md`](HC4_DIRECT_HOMOGENEOUS_FILTRATION.md),
and `HC4MYGJ2` is canonical in
[`HC4_MENG_YANG_GRAPH_OBSTRUCTIONS.md`](HC4_MENG_YANG_GRAPH_OBSTRUCTIONS.md).
The minimal-tower theorem is retained because it records the discovery face.
The stronger propagation theorem below replaces its former lower-layer
limitation.  Neither theorem excludes the underlying direct HC4 packet: they
close constant quadratic pencil admission, not constant-Hessian completion.

<!-- status-consumer: HC4MYGJ2 f28f269f8d01aa75 -->
<!-- status-consumer: HC4FSD1 aee8dfb1e21c3eb5 -->
<!-- status-consumer: HC4FSD2 d413e3020828c4b8 -->
<!-- status-consumer: HC4FSD3 1107bc6ff58456f5 -->

The original pinned regression covers potential degrees four through eight
and graph normal orders one through twelve.  It reports:

- every tested Fermat ternary-Schur radical is exactly the three-dimensional
  pure-power locus;
- every tested minimal diagonal Schur tower in degrees five through eight has
  empty constant rank-two direction scheme on the one-, two-, and
  three-active-channel strata; and
- the Meng--Yang graph equation has nonzero normal multiplier
  \(-4LN^3k(k+1)\) at every order \(k\geq1\).

The new theorem checker verifies the degree-free algebra behind all four
statements.  The finite tables remain useful regressions, not the reason for
the universal quantifiers.

## 1. Fermat ternary-Schur divisibility

Put \(m=D-2\) and normalize

\[
 \operatorname{Hess} f=\operatorname{diag}(x^m,y^m,z^m).
\]

For a ternary form \(a\in K[x,y,z]_m\), the first Schur equation asks whether

\[
 \nabla a^{\mathsf T}\operatorname{adj}(\operatorname{Hess}f)\nabla a
\]

is divisible by \(x^my^mz^m\).  Extract all coefficients of monomials that
fail this divisibility condition.  They generate a homogeneous quadratic
ideal in the coefficients of \(a\).  Singular computes its radical and
compares it in both directions with the ideal of every coefficient except
those of \(x^m,y^m,z^m\).  This was the discovery calculation.  It is now
subsumed by `HC4FSD1`, proved in the
[direct-filtration note](HC4_DIRECT_HOMOGENEOUS_FILTRATION.md): for every
\(m\ge2\), the divisibility condition holds exactly when
> \[
> a=\alpha x^m+\beta y^m+\gamma z^m,
> \qquad
> q=m^2(\alpha^2x^{m-2}+\beta^2y^{m-2}+\gamma^2z^{m-2}).
> \]
Over an algebraic closure, its coefficient-ideal radical is the pure-power
linear ideal in every degree.

| \(D\) | coefficients of \(a\) | equations | radical dimension | result |
|---:|---:|---:|---:|---|
| 4 | 6 | 12 | 3 | pure-power locus |
| 5 | 10 | 33 | 3 | pure-power locus |
| 6 | 15 | 63 | 3 | pure-power locus |
| 7 | 21 | 102 | 3 | pure-power locus |
| 8 | 28 | 150 | 3 | pure-power locus |

Thus no mixed component occurs in any degree: set-theoretically,

\[
 a=\alpha x^m+\beta y^m+\gamma z^m.
\]

This closes the completely split diagonal Hessian stratum of the ternary
Schur problem.  The next target is replacement of the diagonal Hessian with
the remaining nonsquarefree factor strata.  `HC4FSD1` does not classify those
strata or a general ternary Hessian.

## 2. Constant rank-two directions on the diagonal tower

For \(D=m+2\geq5\), use the minimal three-layer tower

\[
 \psi^{\min}_{D}=\sum_{v=x,y,z}\frac{v^{m+2}}{(m+2)(m+1)}
 +\frac{t}{m}\sum_{v=x,y,z}a_vv^m
 +\frac{t^2}{2}\sum_{v=x,y,z}a_v^2v^{m-2}.
\]

The last two layers solve the leading diagonal Schur equation.  Let \(Q\) be
an arbitrary constant symmetric \(4\times4\) matrix.  The exact result is:

> **Theorem `HC4FSD2` -- all-degree minimal diagonal rank-two
> obstruction.**  Let \(m\ge3\) and let
> \((a_x,a_y,a_z)\in\{0,1\}^3\setminus\{0\}\).  There is no constant
> symmetric matrix \(Q\) of rank two satisfying
> \[
> \det(\operatorname{Hess}\psi_D^{\min}+sQ)
> =\det\operatorname{Hess}\psi_D^{\min}
> \]
> identically in \(s,x,y,z,t\).

For a general \(Q\), the top Hessian is
\(\operatorname{diag}(x^m,y^m,z^m,0)\).  The coefficient of \(s\) first
forces \(Q_{tt}=0\).  Three distinct coefficients of \(s^2\) are then
\(-Q_{xt}^2,-Q_{yt}^2,-Q_{zt}^2\), so the complete \(t\)-row and column
vanish.  Write the remaining ternary block as \((q_{ij})\).

On \(t=0\), the coefficient equations are independent of \(m\) after their
monomial shifts.  With one active channel they force

\[
 q_{11}=q_{11}q_{22}-q_{12}^2=q_{11}q_{33}-q_{13}^2=0, \tag{2.1}
\]

so rank two requires \(q_{12}=q_{13}=0\) and a nonzero determinant on the
inactive \(2\times2\) block.  With two active channels the same equations
also force \(q_{22}=q_{12}=0\), hence \(q_{13}=q_{23}=0\) and rank at most
one.  With three active channels they force every \(q_{ij}=0\).  Thus only
the one-channel case reaches the next face.

In that case the coefficient of \(s^2\) is the inactive-block determinant
times the complementary \((x,t)\)-Hessian minor.  Direct differentiation gives

\[
 H_{xx}H_{tt}-H_{xt}^2
 =(3-m)t x^{2m-4}
 -\frac{(m-1)(m-2)}2t^2x^{2m-6}.                   \tag{2.2}
\]

The second coefficient is nonzero for every \(m\ge3\) in characteristic
zero, contradicting rank two.  This proves the theorem.

The discovery experiment extracted every coefficient of

\[
 \det(\operatorname{Hess}\psi_D^{\min}+sQ)
 -\det\operatorname{Hess}\psi_D^{\min},
\]

adjoins \(\det Q=0\), and covers rank exactly two by localizing at the six
distinct \(2\times2\) minors of \(Q\).  All six charts are unit ideals in all
twelve cases.  It is now a regression for the degree-free proof:

| \(D\) | active-channel strata | rank-two survivors |
|---:|---|---:|
| 5 | 1, 2, 3 | 0, 0, 0 |
| 6 | 1, 2, 3 | 0, 0, 0 |
| 7 | 1, 2, 3 | 0, 0, 0 |
| 8 | 1, 2, 3 | 0, 0, 0 |

For the quintic, this initially extended the negative signal from rank-one
theorem `HC4MR4`: the most economical rank-two quadratic repair is absent on
the minimal forced tower.  The displayed tower itself is not asserted to
have constant Hessian determinant.

The apparent lower-layer limitation is visible already at \(m=4\).  Adding
the permitted next homogeneous term \(t^3/6\) replaces the complementary
minor (2.2) by \(t^3\): it cancels both minimal-tower coefficients used above
but does not make the minor zero.  Thus lower layers can move the obstruction
to a later face.  The following argument shows that they can never remove it.

### 2.1 Propagation through every lower layer

Let \(D=m+2\ge5\).  Restore all compatible terms omitted from the minimal
tower:

\[
\begin{aligned}
 \psi_D&=\sum_{v=x,y,z}\frac{v^D}{D(D-1)},\\
 \psi_{D-1}&=\frac{t}{m}\sum_v\alpha_vv^m+b(x,y,z),\\
 \psi_{D-2}&=\frac{t^2}{2}\sum_v\alpha_v^2v^{m-2}
               +t\,c(x,y,z)+d(x,y,z),
\end{aligned}                                      \tag{2.3}
\]

where \(b,c,d\) have the required homogeneous degrees, and let every layer
of degree at most \(D-3\) be arbitrary.

> **Theorem `HC4FSD3` -- all-lower-layer diagonal direction
> obstruction.**  Over a characteristic-zero field, if
> \(\alpha\in K^3\setminus\{0\}\), no constant symmetric matrix \(Q\)
> of rank at least two satisfies
> \[
>  \det(\operatorname{Hess}\psi+sQ)=\det\operatorname{Hess}\psi.
>                                                               \tag{2.4}
> \]
> For \(D=5\), theorem `HC4MR4` also excludes rank one.  Hence the complete
> nonaligned diagonal quintic Schur packet admits no nonzero constant
> symmetric determinant-preserving direction of any rank.

The highest homogeneous coefficients used above are unchanged by \(b,c,d\)
or the later layers.  They first kill the \(t\)-row and column of \(Q\).  A
rank-four \(Q\) is already impossible from \([s^4]=\det Q\).  If the remaining
ternary block has rank three, the exact coefficient

\[
 [s^3]\det(H+sQ)=\det(Q_{xyz})H_{tt}               \tag{2.5}
\]

is nonzero because the leading part of \(H_{tt}\) is
\(\sum_v\alpha_v^2v^{m-2}\).  Thus only rank two remains.  With two or
three active channels, the unchanged top face makes the ternary block have
rank at most one.  With one active channel, say \(x\), it forces \(Q\) to be
supported on the inactive \(y,z\) block.  If that block is \(B\), then

\[
 [s^2]\det(H+sQ)=\det(B)(H_{xx}H_{tt}-H_{xt}^2)    \tag{2.6}
\]

exactly, including every lower layer and every cross entry of \(H\).  Rank
two makes \(\det B\ne0\), so (2.4) would force the binary Hessian determinant
in (2.6) to vanish.

Restrict to \(y=z=0\), put \(\phi(x,t)=\psi(x,0,0,t)\), and set
\(P=\phi_x,\ R=\phi_t\).  The two prescribed leading layers in (2.3) give

\[
 \deg P=D-1,\quad \operatorname{in}(P)=\frac{x^{D-1}}{D-1},
 \qquad
 \deg R=D-2,\quad \operatorname{in}(R)=
 \frac{\alpha_x\,x^{D-2}}{D-2}.                                \tag{2.7}
\]

Vanishing of the binary Hessian is \(J(P,R)=0\).  In characteristic zero
this makes \(P,R\) algebraically dependent.  The polynomial common-generator
theorem gives a polynomial \(h\in K[x,t]\) and univariate \(A,C\) with
\(P=A(h)\), \(R=C(h)\).  Degree multiplicativity and
\(\gcd(D-1,D-2)=1\) force \(\deg h=1\).  Write
\(h=ax+bt+c_0\).  The pure leading form of \(P\) in (2.7) forces \(b=0\) and
\(a\ne0\).  Hence \(P_t=0\), whereas \(R_x\ne0\), contradicting the mixed
partial identity \(P_t=R_x\).  This proves the theorem.

The common-generator result used here is Lemma 5 of I. Arzhantsev and
A. Petravchuk,
[*Closed and Irreducible Polynomials in Several Variables*](https://arxiv.org/abs/math/0608157).
The proof above needs only that lemma, the characteristic-zero Jacobian
criterion, and elementary degree arithmetic.

The practical conclusion is now definitive for this method: no further
constant-matrix search is needed on the completed diagonal quintic packet.
Only nonlinear or polynomially moving directions, collision-preserving
recharts, and non-diagonal Schur strata remain relevant there.  Constant
rank-two recognition can remain useful on non-diagonal survivors.

## 3. The Meng--Yang all-normal symbol

Let \(D_R\) denote the Hessian determinant of the Meng--Yang potential pulled
back to \(r=R(x,y,p,q)\).  For \(U\in K[y,p,q]\), theorem `HC4MYGJ2` gives

\[
 \boxed{
 D_{R+x^kU}-D_R
 \equiv -4LN^3k(k+1)x^{k-1}U\pmod{x^k}
 }
 \qquad(k\geq1).                                    \tag{3.1}
\]

The pinned run records the multipliers through \(k=12\), but (3.1) is proved
uniformly.  Since \(LN\ne0\), the graph equation is unit-triangular in the
complete normal jet.  A trace \(T=R|_{x=0}\) and target constant determine the
first normal coefficient, and every later coefficient is determined
recursively.  Hence there is a unique formal solution in
\(K[y,p,q][[x]]\) above every trace.

This changes the interpretation of the graph route.  Higher normal order is
not where the equation becomes inconsistent: it is formally unobstructed.
The remaining all-degree question is whether the unique formal branch can
terminate as a polynomial while satisfying the global marked collision and
the top-cone constraints.  A useful obstruction must therefore control
degree growth, denominators, or termination; another bounded trace family by
itself does not address the new core issue.

## 4. Route ranking after the experiment

The `route_assessment` array in the committed bounded-regression artifact is
a discovery-time record.  Its instruction to control arbitrary lower layers
was subsequently completed by `HC4FSD3`; it is preserved for provenance and
is not a live handoff.  The current ranking is the one below.

1. **Meng--Yang normal recursion: best immediate computational frontend.**
   It has an exact all-order unit symbol and replaces a growing nonlinear
   coefficient search by deterministic recursion.  The next experiment
   should compute the first several recursively forced \(U_k\) for the
   surviving quintic trace charts and look for a degree-growth or
   nontermination invariant.

2. **Ternary Schur dichotomy: best direct-theorem frontend.**  The diagonal
   radical is now classified in every degree by a support/valuation proof.
   The next step is extension to the nonsquarefree rank-three Hessian strata
   still present in `OP-HC4-D5`.

3. **Constant quadratic recognition: closed on the completed diagonal
   quintic packet.**  `HC4FSD3` excludes ranks at least two through all lower
   layers, and `HC4MR4` excludes rank one.  Nonlinear pencil directions or
   collision-preserving recharts remain possible, but no larger constant
   matrix search is needed on this packet.

The two strongest next projects are therefore complementary: extend the
uniform ternary theorem beyond the diagonal Hessian stratum on the direct
side, and derive a polynomial-termination criterion for the unique
Meng--Yang formal graph on the inherited-collision side.  Lower-layer
propagation on the diagonal constant-matrix branch is no longer an open task.

## 5. Reproduction

Run the degree-free theorem checker

```bash
.venv/bin/python scripts/verify_hc4_all_degree_frontends.py
```

For cleanup-only validation of the bounded artifact, its exact support grid,
scope warning, and generating-source hash, without SymPy or Singular replay,
run

```bash
.venv/bin/python scripts/verify_hc4_all_degree_frontends.py \
  --audit-existing-only
```

Its SHA-256 is
`ff7f2ea1f11e8aee42930e9ab90c5b061711a9d2c753e5633dfba9ee6cf24fb6`.
Replay the bounded discovery regressions with

```bash
.venv/bin/python scripts/research_hc4_all_degree_frontends.py \
  --minimum-degree 4 \
  --maximum-degree 8 \
  --maximum-normal-order 12 \
  --output artifacts/generated-results/hc4_all_degree_frontend_experiments.json
```

The generated artifact is
[`hc4_all_degree_frontend_experiments.json`](artifacts/generated-results/hc4_all_degree_frontend_experiments.json),
with SHA-256
`7176670a76f152f23c6f9b56264c39e2345d91641fdfe275130b9b32e99ffedb`.
The checker/source SHA-256 recorded inside it is
`0b6e0a1272d1201645dce96d8a039fde3c9d4c883476f40d1fed57dac446333f`.
