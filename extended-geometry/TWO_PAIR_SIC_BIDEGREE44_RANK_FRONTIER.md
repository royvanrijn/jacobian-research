# The coefficient-rank frontier in bidegree \((4,4)\)

## 1. Certified answer

Write
\[
 F=\sum_{i,j=0}^4c_{ij}
 \xi _1^i\xi _2^{4-i}z_1^jz_2^{4-j},
 \qquad C=(c_{ij}).
\]
The principal complexity invariant here is the ordinary matrix rank of
\(C\), not expanded support.

The present exact answer is the interval
\[
 \boxed{2\leq \operatorname{rank}C_{\min}\leq5.} \tag{1.1}
\]
The lower endpoint is a theorem: no rank-one bidegree-\((4,4)\) SIC
counterexample exists. The upper endpoint is attained by the known
sixteen-term witness, whose coefficient determinant is \(48\). Ranks two,
three, and four remain open. Thus (1.1) is a rigorous lower bound, not a
claim that the exact minimum is known.

No argument here uses local analysis or deformation around the existing
full-rank point.

## 2. Rank one is impossible

If \(\operatorname{rank}C=1\), then
\[
 F=A(\xi _1,\xi _2)P(z_1,z_2) \tag{2.1}
\]
for binary quartics \(A,P\). Put
\(\Lambda=A(\partial_{z_1},\partial_{z_2})\). Then
\[
 \mathcal E_2(F^m)=\Lambda^m(P^m). \tag{2.2}
\]
Over an algebraic closure every binary quartic \(A\) splits into linear
factors. The
[split-symbol theorem](SPLIT_SYMBOL_GVC_THEOREM.md) says that vanishing of
(2.2) for every \(m\geq1\) implies
\[
 \Lambda^m(QP^m)=0\qquad(m\gg0) \tag{2.3}
\]
for every fixed coordinate polynomial \(Q\).

This gives the full SIC multiplier statement. Write an arbitrary
multiplier as \(g=\sum_\nu B_\nu(\xi)Q_\nu(z)\). Then
\[
 \mathcal E_2(gF^m)
 =\sum_\nu B_\nu(\partial_z)
       \left(\Lambda^m(Q_\nu P^m)\right). \tag{2.4}
\]
Every term is eventually zero by (2.3).

> **Theorem 2.1.** Every bidegree-\((4,4)\) SIC counterexample has
> coefficient-matrix rank at least two.

This is an all-order characteristic-zero theorem, not a bounded moment
exclusion, and it uses no assumption on the root partition of \(A\).

## 3. Determinantal rank charts

For \(1\leq r\leq4\), the rank-at-most-\(r\) locus is parametrized by
\[
 C=UV^{\mathsf T},\qquad
 U,V\in\operatorname{Mat}_{5\times r}. \tag{3.1}
\]
The generic \(\mathrm{GL}_r\) gauge
\[
 (U,V)\longmapsto(UG,VG^{-\mathsf T}) \tag{3.2}
\]
has dimension \(r^2\), so
\[
 \dim\{C:\operatorname{rank}C\leq r\}=r(10-r). \tag{3.3}
\]

| rank bound | factor parameters | gauge | dimension |
|---:|---:|---:|---:|
| \(1\) | \(10\) | \(1\) | \(9\) |
| \(2\) | \(20\) | \(4\) | \(16\) |
| \(3\) | \(30\) | \(9\) | \(21\) |
| \(4\) | \(40\) | \(16\) | \(24\) |

These are global determinantal charts and do not privilege the full-rank
counterexample.

## 4. Pure and mixed moment equations

Introduce
\[
 \Phi_C(x,y)=\sum_{i,j=0}^4c_{ij}x^iy^j. \tag{4.1}
\]
Direct contraction gives
\[
 \boxed{\mu_m(C)=
 \sum_{I=0}^{4m}(4m-I)!\,I!\,
 [x^Iy^I]\Phi_C(x,y)^m.} \tag{4.2}
\]
After (3.1), this is homogeneous of degree \(2m\) in \(U,V\) and invariant
under the gauge (3.2).

For
\[
 Q_{ab}=\xi _1^a\xi _2^{1-a}z_1^bz_2^{1-b},
 \qquad a,b\in\{0,1\},
\]
the bilinear mixed tests are
\[
 \boxed{\nu_{m,a,b}(C)=
 \sum_I (I+a)!\,(4m+1-I-a)!\,
 [x^Iy^{I+a-b}]\Phi_C(x,y)^m,} \tag{4.3}
\]
omitting coefficients outside \(0,\ldots,4m\).

The apolar pairing has diagonal matrix
\[
 D=\operatorname{diag}(24,6,4,6,24). \tag{4.4}
\]
Thus \(L_C=CD\) is, up to transpose convention, the endomorphism attached
to \(C\). On a fixed-flag one-sided chart \(C\) is strictly triangular, so
\[
 \operatorname{tr}(L_C^k)=0\qquad(1\leq k\leq5). \tag{4.5}
\]
A nonzero trace is a useful semistability screen. It is not itself an SIC
defect: one must still exhibit a single fixed \(Q\) whose mixed
contractions fail to vanish eventually.

## 5. Stratum ledger

| exact rank | result | status |
|---:|---|---|
| \(1\) | excluded by Theorem 2.1 | proved |
| \(2\) | chart (3.1), tests (4.2)--(4.3) | open |
| \(3\) | chart (3.1), tests (4.2)--(4.3) | open |
| \(4\) | chart (3.1), tests (4.2)--(4.3) | open |
| \(5\) | explicit witness with \(\det C=48\) | exact counterexample |

The strata are nested: a rank-two counterexample would settle every upper
stratum at once. Conversely, failure to solve finitely many equations on
(3.1) is not a rank lower bound.

An exploratory real rank-two least-squares screen, normalized away from
zero and constrained by nonzero \(\operatorname{tr}(L_C^2)\), produced
points with very small residuals for the first twenty-four pure equations.
They are not counterexamples: there is no exact reconstruction, all-order
identity, or certified nonvanishing mixed tail. This makes exact
reconstruction on the rank-two quotient the first computational target,
not a reason to exclude rank two.

The conditioning audit explains why that screen produced deceptively
small residuals. On a rank-two factor chart, the singular values of the
first-twelve-moment Jacobian at the numerical point range from about
\(3.1\) to \(2.2\cdot10^{-15}\). The seventh normalized residual is the
active one; later raw moments are small largely because of scale. The
numerical point is therefore not suitable for exact reconstruction as it
stands.

## 6. Exact rank-two invariant-quotient obstruction

There is nevertheless an exact reason for truncated semistable rank-two
leads. Let
\[
 X_2=\{C\in\operatorname{Mat}_5:\operatorname{rank}C\leq2\}.
\]
It has dimension \(16\), and its generic diagonal
\(\mathrm{SL}_2\)-quotient has dimension \(13\). On the global chart with
pivot rows of \(U\) equal to the identity, take
\[
 U=\begin{pmatrix}
 1&0\\0&1\\14&17\\18&4\\6&13
 \end{pmatrix},\qquad
 W=\begin{pmatrix}
 8&10&1&8&4\\19&1&4&6&17
 \end{pmatrix},\qquad C=UW.
 \tag{6.1}
\]
Exact evaluation modulo \(1000003\) gives
\[
 \operatorname{rank}
 d(\mu_1,\ldots,\mu_{13})|_C=13. \tag{6.2}
\]
Thus the first thirteen moments are algebraically independent on \(X_2\)
and attain the invariant-quotient dimension.

They cannot be a homogeneous system of parameters. The degree-\(n\)
coordinate ring of the rank-two determinantal variety has the Cauchy
decomposition
\[
 \mathbb Q[X_2]_n
 =\bigoplus_{\substack{\lambda\vdash n\\\ell(\lambda)\leq2}}
 S_\lambda(\operatorname{Sym}^4)^*
 \otimes S_\lambda(\operatorname{Sym}^4). \tag{6.3}
\]
For an \(\mathrm{SL}_2\)-module with even weights, invariant multiplicity
is weight-zero multiplicity minus weight-two multiplicity. Applying this
to (6.3) gives the invariant Hilbert coefficients
\[
 1,1,5,13,53,149,483,1274,3370,7994,18398,39472,81962,
 161896,\ldots. \tag{6.4}
\]
The invariant ring is Cohen--Macaulay: \(X_2\) is the
\(\mathrm{GL}_2\)-quotient of the polynomial factor space \((U,W)\), and
the diagonal \(\mathrm{SL}_2\) action commutes with that gauge, so the
ring is an invariant ring of the reductive product
\(\mathrm{GL}_2\times\mathrm{SL}_2\) acting on a polynomial ring.

If degrees \(1,\ldots,13\) were parameters, the Hilbert numerator would
have nonnegative coefficients. Exact expansion instead gives
\[
 \boxed{
 [t^{69}]H_{\mathbb Q[X_2]^{\mathrm{SL}_2}}(t)
 \prod_{m=1}^{13}(1-t^m)=-5266.} \tag{6.5}
\]
Consequently:

> **Proposition 6.1.** The common zero fiber of
> \(\mu_1,\ldots,\mu_{13}\) on \(X_2\) contains a semistable point.

By itself this proposition does not determine whether the point has rank
one or two.  The finite rank-one analysis below closes that ambiguity.
A zero of thirteen moments still does not satisfy the all-order SIC
premise.

The boundary problem is itself tightly calibrated. On the rank-one Segre
cone, \(\mu_1,\ldots,\mu_6\) have exact Jacobian rank six, the dimension
of its generic invariant quotient. Through degree \(100\), their required
Hilbert numerator is
\[
\begin{aligned}
1+t^2+t^3+4t^4+2t^5+7t^6+5t^7+8t^8+5t^9
 +7t^{10}+2t^{11}\\
{}+4t^{12}+t^{13}+t^{14}+t^{16},
\end{aligned} \tag{6.6}
\]
with zero coefficients after degree \(16\) through the checked cutoff. It
is nonnegative and has coefficient sum \(50\). Thus the first six moments
are a Hilbert-compatible parameter candidate on rank one.

The collided-root strata have short exact cutoffs.  If the symbol has at
most two roots, normalize it to \(u^rv^{4-r}\).  For \(r=0,4\), the first
moment cuts out the one-sided hyperplane.  For \(r=1,2,3\), exact
Gröbner reduction proves that the first four moments have radical equal
to the expected union of the two one-sided linear loci; eighth powers of
the radical generators lie in the moment ideal.  The three-root partition
\((2,1,1)\) is already cut out by the first five moments by the existing
three-root theorem.

The expected annihilator sections are in the diagonal
\(\mathrm{SL}_2\)-nullcone.  Indeed, write \(P=\ell^4\) and choose the
one-parameter subgroup adapted to \(\ell\).  The equation \(A(\ell)=0\)
removes the unique extremal coefficient of \(A\); every remaining weight
of \(A\otimes\ell^4\) is then strictly positive.  The at-most-two-root
and \((2,1,1)\) root-partition strata are therefore closed by finite
prefixes.  It remains here to treat every squarefree cross-ratio.

There is now a smaller exact squarefree chart. Normalize
\[
A=uv(u-v)(u-\lambda v),\qquad e=1,
\]
and eliminate \(b\) with \(\mu_1\). The second moment is linear in \(a\);
after removing a harmless scalar its pivot is
\[
p=5c-5(\lambda+1)d+2\lambda^2+8\lambda+2. \tag{6.7}
\]
On \(p=0\), exact Gröbner reduction over \(\mathbb Q(\lambda)\) shows that
\(\mu_2,\ldots,\mu_6\) generate the unit ideal. Thus the generic
squarefree fiber lies entirely on the principal chart \(p\ne0\), where
\(\mu_2\) eliminates \(a\).

In the remaining \((c,d)\)-plane the three expected annihilator sections
have the remarkably small radical ideal
\[
\boxed{8c-3d^2=0,\qquad d(d-4)(d-4\lambda)=0.} \tag{6.8}
\]
After imposing the first relation, the substituted moments
\(\mu_3,\mu_4,\mu_5,\mu_6\) contain the cubic factor in (6.8) with
multiplicities \(1,1,2,2\), respectively. The unresolved generic step is
settled by a smaller resultant decomposition.  The gcd of all six
pairwise \(c\)-resultants is
\[
d(d-4)(d-4\lambda)q_\lambda(d)^6, \tag{6.9}
\]
where \(q_\lambda\) is quadratic in \(d\).  On the component cut out by
\(q_\lambda\), exact reduction gives \(p^3=0\), so this component is
absent from the principal chart.  On the other three branches the exact
Gröbner bases force, respectively,
\[
(d,c)=(0,0),\quad(4,6),\quad(4\lambda,6\lambda^2). \tag{6.10}
\]
Thus over \(\mathbb Q(\lambda)\) the squarefree six-moment fiber consists
exactly of the annihilator sections.  Two full generic Gröbner variants,
a direct
saturation, and a large fraction-free Bézout solve were slower failed
routes; they are not mathematical evidence.

Three distinguished special squarefree orbits are also closed.  The
harmonic orbit has the existing rational representative \(\lambda=2\).
Over each of the quadratic fields
\[
\lambda^2+4\lambda+1=0,\qquad
\lambda^2-\lambda+1=0, \tag{6.11}
\]
the projective moment fiber again has degree four and radical equal to
the four annihilator sections; every generator of the six-element
radical basis has its eighth power in the moment ideal.  The first orbit
is where the pivot meets an expected section, and the second is the
equianharmonic orbit.  The remaining squarefree task is now precise:
the pivot-boundary resultant gcd is exactly
\[
\lambda^4(\lambda-1)^4
(\lambda^2+4\lambda+1)(\lambda^2-6\lambda+6)
(6\lambda^2-6\lambda+1), \tag{6.12}
\]
where the three quadratic factors form the already-closed
pivot-annihilator orbit.  On the three expected \(d\)-branches, shifting
the expected \(c\)-value and removing the common powers
\((1,1,2,2)\) gives complete pairwise-resultant gcds.  Besides
\(\lambda=0,1\) and the same pivot orbit, their only new factors form
the \(S_3\)-orbit represented by
\[
22\lambda^4-54\lambda^3+\lambda^2-54\lambda+22=0. \tag{6.13}
\]
An exact computation over this quartic field again gives a degree-four
fiber with the expected radical and eighth-power certificates.

The specialization problem admits a final, simpler compatibility test.
Put \(h=8c-3d^2\) and
\(g=d(d-4)(d-4\lambda)\).  Verified modular reconstruction produced
candidate rational standard bases implying
\[
\begin{aligned}
(f_3,f_4,f_5,f_6,\;zph-1)&\ni
  \lambda^4(\lambda-1)^4,\\
(f_3,f_4,f_5,f_6,\;h,\;zpg-1)&=(1).
\end{aligned}
\]
The second membership is now exact without the large four-variable
calculation.  Substitute \(c=3d^2/8\), divide the four moment numerators
by the invertible factors \(g,g,g^2,g^2\), and retain the Rabinowitsch
equation for \(pg\).  The resulting ideal in
\(\mathbb Q[z,d,\lambda]\) has exact reduced standard basis \((1)\).

For the first membership, the same output occurs modulo \(101,103,107\).
Singular's verified
modular standard-basis routine checks that the proposed output is a
standard basis containing the input, but in this nonhomogeneous global
order it does not by itself certify the reverse ideal containment.
Attempts to produce an exact rational lift or syzygy identity exceeded
the available resource envelope.  Therefore the first membership remains
the precise characteristic-zero gate; its modular output is reproducible
evidence, not a theorem.

If the remaining membership is certified, then every squarefree six-moment zero
is one of the four annihilator sections. Combined with the exact
collided-root cutoffs, this would imply that the semistable point of
Proposition 6.1 can be chosen with exact coefficient-matrix rank two.
At present Proposition 6.1 only guarantees rank at most two.

The least degree repair is
\[
 \boxed{\mu_1,\ldots,\mu_{12},\mu_{14}.} \tag{6.14}
\]
At the same point (6.1), these corrected moments also have exact Jacobian
rank thirteen. Their Hilbert numerator is nonnegative through degree
\(100\), with last observed nonzero coefficient in degree \(82\).
Moreover, its coefficients through degree \(82\) are exactly
palindromic.  This is the predicted top degree: the square
rank-at-most-two determinantal ring has \(a\)-invariant \(-5\cdot2=-10\),
unchanged after taking \(\mathrm{SL}_2\)-invariants, while the degrees in
(6.14) sum to \(92\).  Thus (6.14) is a Gorenstein-perfect parameter
candidate.  This remains a necessary-test result, not a proof that
(6.14) is a system of parameters.

## 7. Next exact gates

The efficient order is now:

1. finish the remaining squarefree Rabinowitsch membership above, preferably by
   sparse rational reconstruction of modular target-only lift identities;
2. conditional on that boundary closure, start from the exact-rank-two
   semistable thirteen-moment survivor;
3. on rank two, whose determinantal quotient has expected affine
   dimension \(16-3=13\), normalize a non-null lowest
   Clebsch--Gordan component and remove the
   residual torus before elimination;
4. use the Gorenstein-perfect corrected degree set
   \(\mu_1,\ldots,\mu_{12},\mu_{14}\), rather than the impossible
   consecutive set, to organize the remaining finite-moment fiber;
5. seek an exact component or an effective holonomic recurrence upgrading
   finite pure vanishing to all orders;
6. apply (4.3) with fixed low-bidegree multipliers on every all-order
   survivor; and
7. only after a global rank-two exclusion, repeat for ranks three and
   four.

Promotion to a counterexample requires an exact coefficient field,
all-order pure vanishing, and one fixed mixed defect for infinitely many
\(m\). Promotion of the lower bound beyond two requires a global proof on
the entire determinantal variety.

## 8. Reproduction

Run

```bash
python3 scripts/verify_two_pair_sic_bidegree44_rank_frontier.py
python3 scripts/verify_two_pair_sic_bidegree44_rank_two_invariants.py
.venv/bin/python scripts/verify_two_variable_quartic_squarefree_pivot.py
.venv/bin/python scripts/verify_two_variable_quartic_two_root_finite.py
```

The dependency-free checker verifies the chart dimensions, exact
representatives of ranks one through four, the pure and all four bilinear
mixed formulas through order four, the nilpotent trace screen, and the
known rank-five determinant and moment formulas. It does not turn the
exploratory rank-two residuals into a theorem.

The second checker constructs the rank-two moment Jacobian directly,
proves both displayed rank-thirteen statements modulo a good prime, and
computes the invariant Hilbert coefficients from the Cauchy and
Jacobi--Trudi formulas. It verifies the exact obstruction (6.5), the
corrected necessary test through degree \(100\), and palindromy through
the predicted Gorenstein top degree \(82\).

The third checker proves the generic squarefree pivot-boundary unit certificate,
the two-generator expected affine ideal (6.8), the four exact cubic-factor
multiplicities, and the generic resultant decomposition (6.9)--(6.10).
It also closes the pivot-annihilator and equianharmonic quadratic orbits;
the branch-quartic orbit is closed over its degree-four field, and the
harmonic orbit is supplied by the separate \(\lambda=2\) anchor.  It
extracts the complete pivot-boundary and expected-branch exceptional
gcds and proves the \(h=0,\ g\ne0\) Rabinowitsch chart empty by a smaller
exact unit certificate. It does not certify the remaining \(h\ne0\)
uniform Rabinowitsch membership.
Collided-root quartics are supplied separately by the two-root and
three-root theorems.
