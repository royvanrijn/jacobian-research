# Binary GVC adelic factorial tomography and the six-step packet family

## 1. Status and scope

This note records a finite experiment on the projected two-colour return
semigroup that occurs in the binary Hall--carry quotient.  The calculation is
exact in that model.  It does **not** prove unrestricted
\(\operatorname{GVC}(2)\), and the collisions found below are not GVC
counterexamples.

The useful outcome is sharper.  A strengthened census through radial span
seven uses the unrelated primes \(5,7,11,13\), three prime-power layers, and
units modulo \(p^3\).  Every primitive relation with unequal scalar factorial
partitions is valuation-separated; there is no accidental finite-window
adelic collision.  The only fully decorated \(C_2,C_3\)-blind relations are
three symmetry orbits of one exact all-span family.  A \(C_4\) character
separates every member of that family.  Consecutive-residue rigidity now
closes arbitrary nonfree factorizations after fixed marked promotion.  The
remaining question is whether an actual Hall--jet shell promotes to one such
fixed marked packet at all.

## 2. Finite packet model

Fix a radial span `r`.  Use the marked columns

\[
 R_i=(1,0,i),\qquad B_j=(0,1,j),
 \qquad 0\leq i,j\leq r.
\tag{2.1}
\]

A state `x=(x^R,x^B)` has marked counts

\[
 c_R=|x^R|,\qquad c_B=|x^B|
\]

and total level

\[
 w(x)=\sum_i i x^R_i+\sum_j j x^B_j.
\tag{2.2}
\]

Thus a fibre fixes `(c_R,c_B,w)`.  Its radial vector is represented by

\[
 \rho(x)=\bigl(w,\ r(c_R+c_B)-w\bigr).
\tag{2.3}
\]

The positive and negative parts of an integer kernel element of (2.1) are
two states in one fibre.  A conformally indecomposable kernel element is a
Graver relation.  Since the rank is three, every circuit has support at most
four.  The experimental packet envelope therefore keeps mixed-colour Graver
relations of support at least five.

This deliberately over-approximates the unresolved Hall packets.  It retains
every projected primitive obstruction in the declared span, including some
which an eventual Hall lift may remove by support loss, an additional
character, or a separator.

## 3. Recorded packet and signature

For each side of a Graver relation, the script records:

1. the radial vector (2.3);
2. the positive multiplicity partition on the `R` side;
3. the positive multiplicity partition on the `B` side;
4. the exact marked level--multiplicity lists;
5. the base-`p` digits of the scaled multiplicities;
6. every outgoing Kummer carry, with its digit position;
7. the Legendre layers of both radial factorials;
8. the `p`-free units of the two multinomials, the radial factorial, and
   their product modulo `p^k`; and
9. the group-ring elements in \(\mathbb Z[C_t]\) whose Fourier evaluations
   give all marked-side \(C_t\)-character traces, together with the two
   monomial character classes.

At scale \(N=q p^e\), the scalar weight attached to a state is

\[
 \mathcal W_x(N)=
 \binom{Nc_R}{Nx^R}
 \binom{Nc_B}{Nx^B}
 (N\rho_x)!(N\rho_y)!.
\tag{3.1}
\]

The radial factorial is common inside a fibre.  Kummer's theorem gives the
two multinomial valuations as the sums of the recorded carry amounts.
Legendre's formula gives the remaining valuation.  Removing those powers of
`p` gives the stored unit modulo `p^k`.

Two checked-in runs are complementary.  The full span-six pilot stores every
probe row and uses

\[
 p\in\{2,3,5,7,11\},\qquad 1\leq e\leq2,qquad 1\leq q\leq3,
 \qquad k=2,qquad t\in\{2,3\}.
\tag{3.2}
\]

The exact marked levels are retained as metadata but are not inserted into
the comparison key: doing that would identify a state by definition.  The
comparison instead asks how much of that marking is recovered by the stated
torsion characters.

The strengthened primitive census uses

\[
 p\in\{5,7,11,13\},\qquad 1\leq e\leq3,\qquad 1\leq q\leq3,
 \qquad k=3,
\tag{3.3}
\]

again with \(C_2,C_3\) in the configured decorated signature.  It separates
the scalar signature

\[
 \Sigma_W(x)=
 \left(v_p\mathcal W_x(qp^e),
 p^{-v_p\mathcal W_x(qp^e)}\mathcal W_x(qp^e)\bmod p^3
 \right)_{p,e,q}
\tag{3.4}
\]

from the marked digit/carry/unit and torsion decorations.  Thus a packet can
be reported as valuation-separated or unit-separated without first inserting
its exact partition into the comparison key.  Exact Stirling data are stored
as an independent all-scale separator.  Primitive-only mode omits global
packet-sum grouping, but retains a first-separator certificate for every
Graver relation and every probe row for each survivor.

## 4. Exact Graver computation

Let `A_r` be the matrix with columns (2.1).  The Graver basis is computed as
the nontrivial part of the Hilbert basis of the Lawrence lifting

\[
 \{(u,v)\in\mathbb N^{2n}:A_ru-A_rv=0\}.
\tag{4.1}
\]

The diagonal generators are discarded.  Every remaining Hilbert generator
has disjoint positive and negative supports and is a Graver relation.

One full-support calculation suffices for all projected supports.  A Graver
relation of a coordinate subconfiguration remains Graver after extension by
zero, and every full relation supported in that subconfiguration restricts
back.  The artifact therefore records the universal basis once and obtains
the basis of each nonfree projected support by exact support filtering.

Normaliz 3.10.2 performs the Hilbert-basis calculation.  Translation and a
common dilation of all active levels are removed, after which side exchange,
colour exchange, and level reversal are canonicalized.

## 5. Separator funnels through span seven

The strengthened run gives the following exact comparison.  The span-six
column is recomputed with (3.3), so the two columns use the same adelic window.

| stage | span six | span seven |
|---|---:|---:|
| full two-colour Graver basis | 8,559 | 34,890 |
| after translation/dilation and finite symmetries | 1,584 | 6,601 |
| mixed primitive packet envelope, support at least five | 1,490 | 6,401 |
| exact projected support semigroups represented | 868 | 3,107 |
| entropy- and valuation-separated | 1,037 | 4,750 |
| finite adelic scalar collisions | 453 | 1,651 |
| accidental finite-window scalar collisions | 0 | 0 |
| separated by marked-side digits or carries | 60 | 221 |
| separated next by configured torsion characters | 391 | 1,427 |
| equal full configured \(C_2,C_3\) signature | 2 | 3 |

In both spans, the finite adelic scalar collision set is exactly the set of
relations whose combined multiplicity partitions agree.  Every unequal pair
is already separated by a total valuation in the configured window; no unit
residue is needed after that step.  This is a bounded statement, but it has no
false positive against the exact scaled-factorial test.

The first-prime distribution shows that this is genuinely adelic rather than
a disguised one-prime test:

| first separating prime | unequal relations |
|---:|---:|
| \(5\) | 4,557 |
| \(7\) | 148 |
| \(11\) | 45 |
| \(13\) | 0 |

All first separators occur already at \(e=1\).  Thus unrelated primes are
essential in this range, while higher prime powers and unit residues add no
new separation through span seven.  They remain recorded because that
redundancy is not proved uniformly in the span or in actual Hall lifts.

For the 1,430 span-seven relations whose operator and polynomial partitions
also agree separately, the first separating torsion orders are

| first order | relations |
|---:|---:|
| \(2\) | 1,244 |
| \(3\) | 183 |
| \(4\) | 3 |

Consequently \(C_2,C_3,C_4\) separate every primitive relation in this
span-seven projected census.

## 6. The universal six-step collision family

The three \(C_2,C_3\)-blind orbits are instances of the all-span relation

\[
 R_{s+6}B_aB_{a+1}=R_sB_{a+3}B_{a+4},
 \qquad s,a\geq0.
\tag{6.1}
\]

Whenever the displayed levels lie in the declared span, (6.1) is a Graver
relation.  Indeed a conformal subrelation must use both operator marks.  One
polynomial mark on each side cannot compensate their level difference six;
using both polynomial marks gives the whole relation.  Its exact projected
fibre has only two states: with `R_(s+6)`, the unique polynomial level sum is
`2a+1`, while with `R_s` it is `2a+7`.

After the configured affine-level and finite symmetries, span six has two
orbits and span seven has three.  Stored span-seven representatives are

\[
\begin{aligned}
 R_6B_2B_3&=R_0B_5B_6,\\
 R_6B_1B_2&=R_0B_4B_5,\\
 R_7B_0B_1&=R_1B_3B_4.
\end{aligned}
\tag{6.2}
\]

Writing `r` for the declared span, both sides have marked partitions
\(\lambda_R=(1)\), \(\lambda_B=(1,1)\), and radial vector

\[
 \rho_{r,s,a}=(s+2a+7,\ 3r-s-2a-7).
\]

Consequently, for every \(N\geq1\), both scalar weights equal

\[
 \binom{2N}{N}
 \bigl(N(s+2a+7)\bigr)!
 \bigl(N(3r-s-2a-7)\bigr)!.
\tag{6.3}
\]

Thus every prime has identical Legendre/Kummer valuations, multiplicity
digits, and factorial units.  Modulo two, a level shift by six fixes the
operator residue and a shifted adjacent pair has the same two-element
histogram.  Modulo three, shifts by six and three are invisible.  Hence all
\(C_2,C_3\) group-ring traces agree at every scale.

A \(C_4\) character always distinguishes the operator marks because their
difference is six.  The artifact records `4` as the first separating torsion
order for all three span-seven orbits.  Level reversal sends the relative
offset \(\delta=s-a\) to \(-\delta-2\), explaining the orbit growth as the
span increases.

## 7. Fixed-character Hall termination of the six-step family

The projected collision itself can be closed more strongly than the census
suggests.  Subtracting `s` times the operator count and `a` times the
polynomial count from the level row identifies every support in (6.1) with

\[
 R_6,R_0,B_0,B_1,B_3,B_4.
\tag{7.1}
\]

At scale \(N\), write a state in this order as
\((r_6,r_0,b_0,b_1,b_3,b_4)\).  Fix the color counts \((N,2N)\), total
level \(7N\), and the \(C_2,C_3\) marked histograms of either endpoint.
The polynomial-side histogram equations are

\[
 b_0+b_4=b_1+b_3=N,
 \qquad
 b_0+b_3=b_1+b_4=N.
\tag{7.2}
\]

They give \(b_3=b_4=t\) and \(b_0=b_1=N-t\).  The level equation then
gives \(r_6=N-t\), \(r_0=t\).  Thus the complete blind fibre is exactly

\[
 z_t=(N-t,t,N-t,N-t,t,t),
 \qquad 0\leq t\leq N.
\tag{7.3}
\]

Let \(A_+,A_-\) be the coefficients of \(R_6,R_0\), let
\(B_0,B_1,B_3,B_4\) be the polynomial-channel coefficients, and put

\[
 U=A_+B_0B_1,
 \qquad
 V=A_-B_3B_4.
\tag{7.4}
\]

The operator multinomial of \(z_t\) is \(\binom Nt\), while its
polynomial multinomial is

\[
 \frac{(2N)!}{(N-t)!^2t!^2}
 =\binom{2N}{N}\binom Nt^2.
\]

After removing the common radial factorial, the complete \(C_2,C_3\)-blind
row is therefore

\[
 \mathcal F_N(U,V)=
 \binom{2N}{N}
 \sum_{t=0}^N\binom Nt^3U^{N-t}V^t.
\tag{7.5}
\]

The sum is a two-variable Franel polynomial.  Equivalently,

\[
 \binom{2N}{N}^{-1}\mathcal F_N(U,V)
 =\operatorname {CT}_{y,z}
 \left((1+y)(1+z)(U+Vy^{-1}z^{-1})\right)^N.
\tag{7.6}
\]

This gives the promised terminal theorem.

> **Theorem 7.1 (fixed-character six-step packet termination).**  Work over
> a characteristic-zero field.  Let a fixed homomorphism to a finite abelian
> marking group refine the \(C_2,C_3\)-blind fibre (7.3), and let \(h\) be
> the order of the relative endpoint class.  If the normalized marking-fibre
> rows containing both endpoints vanish at scales \(h\) and \(2h\), then
> \(U=V=0\); hence the six-step support strictly drops.  If the relative
> class is nonzero, the same marking already separates the primitive
> endpoints at scale one.

Indeed, if the marking is \(\pi\), then along (7.3)
\(\pi(z_t)=\pi(z_0)+t\delta\), where
\(\delta=\pi(z_1-z_0)\).  Let \(h\) be the order of \(\delta\), with
\(h=1\) when the marking is blind.  At scale \(N=h\), the
endpoint-containing class consists of \(t=0,h\); at scale \(N=2h\), it
consists of \(t=0,h,2h\).  Removing the nonzero common radial and central
binomial factors gives

\[
 U^h+V^h=0,
\tag{7.7}
\]

and

\[
 U^{2h}+\binom{2h}{h}^{\!3}U^hV^h+V^{2h}=0.
\tag{7.8}
\]

Set \(X=U^h\), \(Y=V^h\).  Equation (7.7) gives \(Y=-X\), so (7.8)
becomes

\[
 \left(2-\binom{2h}{h}^{\!3}\right)X^2=0.
\tag{7.9}
\]

Since \(\binom{2h}{h}\geq2\), its coefficient is nonzero in
characteristic zero.  Hence \(X=Y=0\), and therefore \(U=V=0\).  Each
product in (7.4) then loses at least one active channel.  This proves the
theorem.

The canonical \(C_4\) operator marking has
\(6-0=2\pmod4\), so its relative character has order two.  It separates
the primitive endpoints at odd scales; if an even-scale class groups them,
Theorem 7.1 closes that class using the rows at scales two and four.  With no
further character at all, \(h=1\), and the first two Franel rows are

\[
 U+V,
 \qquad
 U^2+8UV+V^2.
\]

Their common zero is again only \(U=V=0\).

> **Corollary 7.2.**  The six-step family is not a nonzero obstruction after
> promotion to any fixed finite-character Hall packet.  Such a packet is
> character-separated, Hall-terminal in two character-periods, or loses
> support.

There is also an exact affine-carry refinement which does not assume that the
borrow state is independent of the prime.  Remove the common central binomial
factor from (7.5) and write

\[
 f_N(U,V)=\sum_{t=0}^N\binom Nt^3U^{N-t}V^t.
\tag{7.10}
\]

> **Proposition 7.3 (two-digit Franel carry transform).**  Let \(p\) be an
> odd prime, \(0\leq r<p\), \(1\leq d<p\), and \(N=pd+r\).  The
> valuation-zero part of (7.10), reduced modulo \(p\), is
> \[
>  f_d(U^p,V^p)f_r(U,V).
> \tag{7.11}
> \]
> Every remaining coefficient has valuation exactly three.  After division
> by \(p^3\), that carry shell is
> \[
>  d^3 f_{d-1}(U^p,V^p)
>  \sum_{k=1}^{p-1-r}
>  \frac{(-1)^{k+1}U^{p-k}V^{r+k}}
>       {\left(k\binom{r+k}{k}\right)^3}
>  \pmod p.
> \tag{7.12}
> \]

Write \(t=pj+s\).  If \(s\leq r\), Lucas gives valuation zero and
\[
 \binom{pd+r}{pj+s}^3
 \equiv\binom dj^3\binom rs^3\pmod p,
\]
which sums to (7.11).  If \(s=r+k>r\), subtraction borrows once, so
Kummer gives valuation one for the binomial and valuation three for its
cube.  The two-digit factorial-unit formula gives

\[
 p^{-1}\binom{pd+r}{pj+r+k}
 \equiv
 \frac{(-1)^{k+1}d\binom{d-1}{j}}
      {k\binom{r+k}{k}}
 \pmod p.
\tag{7.13}
\]

Cubing and summing proves (7.12).

The first affine digit already detects the blind Hall branch.  Take
\(d=r=1\), so \(N=p+1\), and suppose \(V=-U\).  The four no-borrow
parameters are \(t=0,1,p,p+1\).  Their exact sum is

\[
 2\left(1-(p+1)^3\right)U^{p+1},
\]

and hence its first normalized unit is

\[
 p^{-1}\sum_{t\in\{0,1,p,p+1\}}
 \binom{p+1}{t}^3U^{p+1-t}(-U)^t
 \equiv-6U^{p+1}\pmod p.                       \tag{7.14}
\]

This is nonzero for every good \(p\geq5\) while \(U\) is active.  For
completeness, the normalized borrow shell at every \(p\geq7\) is

\[
 U^{p+1}\sum_{k=1}^{p-2}\frac1{k^3(k+1)^3}
 =20U^{p+1}\pmod p.                             \tag{7.15}
\]

Indeed

\[
 \frac1{k^3(k+1)^3}
 =\frac6k-\frac3{k^2}+\frac1{k^3}
  -\frac6{k+1}-\frac3{(k+1)^2}-\frac1{(k+1)^3},
\]

and summing from \(1\) to \(p-2\), using the vanishing of the complete
second and third inverse-power sums in \(\mathbb F_p\), gives \(12+6+2\).
Thus neither the no-borrow correction nor the borrow shell supplies an
internal affine-carry escape for the complete six-step fibre.

The calculation is not special to the levels in the span-seven census.  Its
actual hypothesis is a rank-one marked lattice.

> **Theorem 7.4 (squarefree rank-one Hall packets).**  Let a marked return
> matrix on ordered columns
> \(R_+,R_-,B_0,B_1,B_2,B_3\) have integer kernel
> \[
>  \ker_{\mathbb Z}A
>  =\mathbb Z(-1,1,-1,-1,1,1).                  \tag{7.16}
> \]
> Suppose its endpoint fibre has operator count \(N\), polynomial count
> \(2N\), and contains
> \((N,0,N,N,0,0)\).  Then its complete nonnegative fibre is (7.3), its
> normalized row is (7.5), and all conclusions of Theorem 7.1 and
> Proposition 7.3 hold.  Thus no squarefree rank-one marked Hall packet can
> remain nonzero after its two endpoint-containing rows are inherited.

Every integral point in the fibre differs from the displayed endpoint by an
integer multiple of the primitive generator (7.16).  Nonnegativity restricts
that multiple to (0\leq t\leq N\), proving (7.3).  The two color
multinomials depend only on the six multiplicities, not on the level or
marking rows, so their product is again
(\binom{2N}{N}\binom Nt^3\).  The fixed-character and affine-carry proofs
therefore apply without change.

The odd-shift family gives a uniform all-span application.

> **Corollary 7.5 (odd-shift Franel family).**  For every odd \(g\geq3\),
> the projected primitive relation
> \[
>  R_{s+2g}B_aB_{a+1}=R_sB_{a+g}B_{a+g+1}        \tag{7.17}
> \]
> has, after level normalization, a complete \(C_2,C_g\)-blind scaled
> fibre equal to (7.3).  Consequently Theorem 7.1 and Proposition 7.3
> apply verbatim: every fixed finite-character promotion is separated or
> loses support, and the prime-dependent two-digit carry split has no
> internal Hall survivor.  A \(C_{g+1}\) marking separates the primitive
> endpoints.

Normalize (7.16) to

\[
 R_{2g},R_0,B_0,B_1,B_g,B_{g+1}.
\]

The \(C_g\) polynomial histogram gives

\[
 b_0+b_g=N,\qquad b_1+b_{g+1}=N.
\]

Write \(t=b_g\), \(u=b_{g+1}\).  Since \(g\) is odd, the \(C_2\)
histogram gives

\[
 b_0+b_{g+1}=N,\qquad b_1+b_g=N,
\]

and hence \(t=u\).  The total level \((2g+1)N\) then forces the
\(R_{2g}\)-count to be \(N-t\).  This is exactly (7.3), and its
multinomial coefficient is still
\(\binom{2N}{N}\binom Nt^3\).  All subsequent arguments depend only on
that line and its weights.  Finally \(2g\not\equiv0\pmod{g+1}\), so
\(C_{g+1}\) distinguishes the operator endpoints.  Primitivity on the
displayed support follows because a proper subrelation using both operator
marks cannot bridge their level gap \(2g\) with only one polynomial mark
from each side.

## 8. Nonfree factorizations after marked promotion

The primitive two-state census does not by itself control a nonfree scaled
fibre: one state can have several factorizations into Hilbert atoms.  The
[nonfree-factorization tomography theorem](BINARY_GVC_NONFREE_FACTORIZATION_TOMOGRAPHY.md)
now closes that ambiguity in every fixed marked span.

For radial span \(s\), put \(q=\lceil s/2\rceil\).  The complete residue
histograms modulo \(q\) and \(q+1\) are the incidence matrix of a bipartite
path on levels \(0,\ldots,2q-1\), hence are injective there.  Adding level
\(2q\) creates one cycle, with kernel

\[
 (\underbrace{1,\ldots,1}_{q},0,
  \underbrace{-1,\ldots,-1}_{q}).
\tag{8.1}
\]

For two colours, equality of total radial level forces the two kernel
coefficients to be opposite.  The resulting return is a conformal sum of

\[
 R_iB_{i+q+1}=R_{i+q+1}B_i,\qquad0\leq i<q,
\tag{8.2}
\]

which are the already-safe four-point beta circuits.  Thus the
\(C_q,C_{q+1}\) atom signature is injective modulo safe circuits in every
span, and any factorization relation with equal signature multisets reduces
atomwise to those circuits.  Together with factorial-trace separation, this
removes the general same-vector/nonfree-semigroup obstruction after a fixed
marked packet has been inherited.

The exact Hilbert-basis experiment independently checks the smaller
\(C_2,C_3,C_4\) labels through span six.  It finds the first genuine
factorial-only square, but \(C_3\) separates it, and no all-signature
factorization collision remains in 65 span-four, 400 span-five, or 1,469
span-six profiles.  These bounded counts are evidence; the all-span
conclusion is the incidence-forest proof above.

## 9. What remains open

Corollary 7.2 proves the fixed-character version of the proposed Hall-lift
theorem.  It also shows that neither distinct-ray factorial independence nor
a higher torsion order is needed once this packet has actually been exposed.
Proposition 7.3 further proves that allowing the carry subset itself to vary
with \(p\) does not create an internal cancellation on the complete six-step
line: the first affine digit has the nonzero correction (7.14).
Theorem 7.4 closes every squarefree rank-one marked lift, and Corollary 7.5
closes the same internal mechanism for every odd shift.  Thus
increasing the radial span merely to find a higher torsion-blind adjacent-pair
relation will not move the frontier.

It does **not** prove that the prime-dependent Hall--jet initial shell exposes
one fixed packet.  Before scale-compatible promotion, an affine carry state
may depend on the exposing prime and may select different proper subsets of
(7.3) at different scales.  Pure vanishing is known only for the sum of all
tied initial shells, so (7.7)--(7.8) cannot be assigned to one six-step block
without a factorial-compatible inheritance argument.  In particular, another
tied packet could cancel (7.14) before the six-step block is separately
exposed.  The subsequent
[translation-tangent theorem](BINARY_GVC_TRANSLATION_TANGENT_RIGIDITY.md)
shows that a coefficient-blind Hilbert-module inheritance principle is false,
but proves that the primitive one-direction tangent kernel and every
prime-power-character collision with sufficiently large underlying prime are
flat.  Ruling out the remaining exceptional-small-prime, mixed-prime, and
two-dimensional cross-packet cancellation over one common high quotient is
precisely the corrected
Cartesian affine-carry promotion problem `(SC)` in the uniform binary note.

Thus a binary counterexample route through (6.1) would now require all of:

1. a prime-dependent affine carry selection which never promotes to a fixed
   finite-character refinement of (7.3);
2. cancellation with other tied shells which prevents the two Franel rows
   from being inherited by the six-step block;
3. an actual pure-zero binary operator/polynomial pair realizing that shell;
   and
4. a polynomial multiplier with a persistent nonzero mixed tail.

No such lift is supplied by the census.  The surviving global theorem is
therefore not six-step terminality, which is now proved, but fixed-packet
exposure from the prime-dependent affine Hall filtration.

## 10. Reproduction

Run from the repository root:

```bash
.venv/bin/python scripts/research_binary_gvc_prime_power_tomography.py \
  --output artifacts/generated-results/binary_gvc_prime_power_tomography.json.gz \
  --summary-output artifacts/generated-results/binary_gvc_prime_power_tomography_summary.json
```

The compressed artifact contains every Graver relation, each projected
support basis, all packet metadata, and every probe row.  The summary contains
the census and the complete records for the two collision orbits.  Their
logical result hash and whole-file hashes are recorded in
[`REPRODUCE.md`](../REPRODUCE.md) and the generated-results index.

Run the strengthened primitive census with:

```bash
.venv/bin/python scripts/research_binary_gvc_prime_power_tomography.py \
  --radial-degree 7 --primitive-only \
  --primes 5,7,11,13 --max-exponent 3 --max-quotient 3 \
  --unit-power 3 --torsion-orders 2,3 \
  --output artifacts/generated-results/binary_gvc_adelic_tomography_span7.json.gz \
  --summary-output artifacts/generated-results/binary_gvc_adelic_tomography_span7_summary.json
```

Primitive-only mode records the universal Graver basis, the exact basis of
each represented support semigroup, every primitive classification and
first-separator certificate, and full tomography rows for the three
survivors.  It deliberately omits global repeated-state collision groups,
which are not needed to separate primitive semigroup relations.

Replay the exact six-step fibre and Franel-weight identities with:

```bash
python3 scripts/verify_binary_gvc_six_step_packet_termination.py
```

The script is a dependency-free bounded regression for the general proof in
Theorem 7.1; it is not substituted for the characteristic-zero argument
(7.2)--(7.9).

Replay the consecutive-residue incidence theorem without invoking Normaliz:

```bash
.venv/bin/python scripts/research_binary_gvc_nonfree_factorization.py \
  --verify-consecutive-residues
```

Run the complete span-four factorial/Graver census and the larger
Hilbert-atom signature censuses with:

```bash
.venv/bin/python scripts/research_binary_gvc_nonfree_factorization.py

.venv/bin/python scripts/research_binary_gvc_nonfree_factorization.py \
  --radial-degree 5 --signature-only --torsion-orders 2,3,4 \
  --output artifacts/generated-results/binary_gvc_nonfree_factorization_span5_signature.json

.venv/bin/python scripts/research_binary_gvc_nonfree_factorization.py \
  --radial-degree 6 --signature-only --torsion-orders 2,3,4 \
  --output artifacts/generated-results/binary_gvc_nonfree_factorization_span6_signature.json
```

The first command writes the span-four artifact at its default generated
results path.  The exact logical and whole-file hashes are recorded in the
generated-results index and in REPRODUCE.md.
