# Structural closure targets for Strong Factorial

Internal research synthesis, 21 November 2026. Inspected repository base and
latest GitHub main: `19eb9912fb60142f8b620dac24896931d4ff3e46`.

The research target is now a theorem uniform in the support parameters and
moment order. Previously completed individual supports remain useful checks.
They are recorded in
[the two-support synthesis](STEPANOV_SFC_PROGRESS_AND_LOWER_THREE_CLOSURE.md).
This note records the structural results and the exact implications that
remain open. **Neither unrestricted SFC nor all two-support SFC is closed.**

The fixed-window continuation proves a finite window at every starting order
for all monomial-times-linear frames, an exact all-support initial resultant,
and an actual support-size SFC theorem for rational block averages under an
explicit size condition. The latest audit adds an exact moving-caustic
certificate for common-factor descent, proves obstructions to two proposed
universal scalar-operator shortcuts, and verifies shifted-resultant positivity
for all binary exponent parameters at starting orders one through five.
The sharp general monomial-linear window and unrestricted SFC remain open.
The ranked next proof targets and their scopes are in
[SFC_GENERAL_CLOSURE_RESEARCH_PROGRAMME.md](SFC_GENERAL_CLOSURE_RESEARCH_PROGRAMME.md).

## 1. Results established in this investigation

The ordinary factorial functional is
`L(X^alpha)=prod_i alpha_i!`. All ordinary coefficient-field conclusions below
hold over every field of characteristic zero. Restrictions on auxiliary
rational Gamma parameters or rational linear-form coefficients are part of
their respective statements.

| Structural result | Quantifiers and conclusion | Remaining limitation |
|---|---|---|
| Balanced private amplification | For arbitrary support size `s`, shared monomials with identical exponent profiles, arbitrary coefficients and every `D>=6s`, `sum c_j X_j^D M_j` satisfies SFC at every starting order. Every finite nonnegative integer support admits an explicit private enlargement of this form. | The private factors and profile balancing change the moments. Their removal is unproved. |
| Separating shared rows | Balanced base monomials with partition layers admitting a rational separating-set packing of total capacity at least `6s` satisfy all-order SFC. Every intermediate within-layer merger certified by the final packing is permitted. | The prescribed layers remain distinct from each other and the base variables; loss of separation is not covered. |
| Intrinsic zero-pattern amplification | Balanced columns with distinct zero/nonzero patterns, separated by `ell` rows, satisfy all-order SFC whenever every positive exponent is at least `6s ell`. One may take `ell<=s-1`. No new variables are needed. | Identical zero patterns cannot satisfy any positive uniform near-pure suppression exponent; the general case requires a different argument. |
| Ramified polynomial pullback | A normalized polynomial cover whose critical values lie in the prescribed seed fibers preserves every normalized factorial moment of a monomial-linear seed, including repeated coefficients and zero exponents. A degree-`d` cover gives `dS` actual terms with the original seed window. | Target coefficients are prescribed reciprocals of preimages. A nontrivial cover degree divides `gcd(S_target,1+sum r_i)`, so primitive cases remain. |
| Monomial times linear form: fixed windows | For every nonnegative exponent profile and all active coefficients, put `d=1+sum r_i`, `e=s-ceil(s/d)`. Every window of length `de+s d(d+1)/2` contains a nonzero moment, at every starting order. When the grouped primitive is polynomial, the stronger fixed length `s d(d-1)/2` works. | These general bounds exceed the required SFC length `s`. |
| Reduced-divisor contact | With `nu=deg rad(D)` and `T0=deg(J/gcd(J,rad(D)))`, every window of length `(nu+T0-1)n-T0+1` works. All forced neutral, repeated-root, and shared-divisor factors are removed exactly. | The final numerator has exactly the computed degree; sharp closure requires control of its order at zero. |
| Exact initial resultant | For all support sizes and formal exponent parameters, the resultant of `g1,...,gs` is the explicit product of all nonempty subset sums `1+sum_I r_i`, with multiplicity `(|I|-1)!(s-|I|)!`. Uniform ternary certificates also cover orders 2–4 for `(a,a,b)` and orders 3–5 for `(1,1,b)`. | The initial formula does not extend unchanged to shifted windows. The ternary certificates have fixed starting orders. |
| Shifted binary resultant positivity | For each starting order 1 through 5, the exact resultant is coefficientwise positive in both common exponents; its displayed residual is Schur-positive. Each identity covers all nonnegative exponent parameters and all coefficient vectors. | Starting order and support size remain bounded. Universal Newton/Schur positivity is conjectural. |
| Pointed moving-caustic detector | For every characteristic-zero polynomial pair with `a'=H'D`, `deg D=S>=1`, the moving resultant is a nonzero scalar times an absolutely irreducible polynomial to power `e=[K(z):K(a,H)]`. This detects the maximal common right factor and, for actual frames, the existing support reduction to `S/e`. | The proposed implication from a forbidden sharp gap to `e>1` is unproved; polynomial `H` is required. |
| Scalar-operator branch obstructions | Two exact positive-exponent binary frames force late Euler resonances in every support-width homogeneous annihilator. The second also forces a nonconstant polynomial into its kernel. | These refute universal operator hypotheses, not SFC. The actual marked branch must remain distinguished from conjugate combinations. |
| Rational block-average SFC | For group sizes `p_j=b v_j` and exponent totals `R_j=a v_j`, actual support-size SFC holds at every starting order whenever `2(b-floor(b/d))>=a h d`, where `d=1+a sum v_j` and `h` is the number of active groups. | Coefficients agree within each group; the rational-average and size conditions are retained. |
| Cubic determinant stability | For every order `n`, all monic perturbations of the stated consecutive cubic models with coefficient errors of three-adic valuation greater than one have resultant valuation exactly `deg gcd(R_n,R_(n+1))`. | General support moment polynomials need not have this normalization. |
| Prime-power scalar transfer | All actual nonresonant leaves and their seeds are controlled on the closed root disk, uniformly in every prime power `q`. This gives explicit mixed-offset criteria and unconditional offset-one and offset-two coprimality. | Higher upper-response degeneracies and gaps with several prime factors remain. |
| Prime-power divisibility closure | For every prime-power gap `q`, every positive divisor `a` of `q`, and every starting order, `X^a(A+B X^q)` has a nonzero moment in each consecutive pair. The same holds for every integer Gamma shape `1<=beta<=a`. | This does not cover every `a<q`, or gaps with several prime factors. Parts of the ordinary family overlap earlier theorems; the proof is one uniform mechanism. |
| Actual upper-response degree bounds | For every prime gap `p` and offset `1<=h<=p`, the mixed-pair gcd has degree at most `2 floor((h-1)/2)`. Every prime power also has a denominator-depth degree bound and the stronger criterion `v_p(kappa)+2s+D>0`. | A possible nonreal common factor can remain; the prime translation does not extend unchanged to proper prime powers. |
| Identical disjoint monomial pairs | For every nonconstant monomial `M`, arbitrary coefficients, and every starting order, `A M(X)+B M(Y)` on disjoint copies has a nonzero moment in each consecutive pair. | The two monomials have the same exponent template and disjoint variables. |
| Overlapping swap pairs | Mirrored exponent rows `(b_i+d_i,b_i)` and `(b_i,b_i+d_i)` give all-order SFC(2), with arbitrary coefficients, whenever `sum_(b_i=0) d_i>=2`. Equal rows and specified symmetric mergers are allowed. | The private-degree threshold is sharp for this normalized-reciprocal certificate; smaller private degree is not an SFC counterexample. |
| Cyclic monomial blocks | Arbitrary degrees and support sizes, with disjoint cyclic copies: the moment zero set is exactly the complement of the semigroup generated by the cycle lengths. Powers and disjoint products also satisfy SFC. | The root-of-unity coefficient structure and disjoint-variable assumptions are retained. |
| Cyclic polynomial blocks with a common monomial | The exact semigroup zero set extends to nonnegative rational polynomial templates and an overlapping common monomial whose exponents agree across each cycle. Symmetric within-copy identifications, powers and disjoint products are covered. | The cyclic coefficient pattern and identical-copy symmetry remain essential; arbitrary cross-copy merging is open. |

Full proofs, exact scope and executable evidence are in:

- [SFC_BALANCED_PRIVATE_AMPLIFICATION.md](SFC_BALANCED_PRIVATE_AMPLIFICATION.md).
- [SFC_SEPARATING_ROW_AMPLIFICATION.md](SFC_SEPARATING_ROW_AMPLIFICATION.md).
- [SFC_DEAMPLIFICATION_BRIDGES.md](SFC_DEAMPLIFICATION_BRIDGES.md).
- [SFC_MONOMIAL_LINEAR_CONTACT_BOUND.md](SFC_MONOMIAL_LINEAR_CONTACT_BOUND.md).
- [SFC_ALGEBRAIC_UNIFORM_WINDOWS.md](SFC_ALGEBRAIC_UNIFORM_WINDOWS.md).
- [SFC_REDUCED_DIVISOR_CONTACT_BOUND.md](SFC_REDUCED_DIVISOR_CONTACT_BOUND.md).
- [SFC_INITIAL_RESULTANT_AND_SHIFTED_WINDOWS.md](SFC_INITIAL_RESULTANT_AND_SHIFTED_WINDOWS.md).
- [SFC_SHIFTED_RESULTANT_POSITIVITY_AUDIT.md](SFC_SHIFTED_RESULTANT_POSITIVITY_AUDIT.md).
- [SFC_POINTED_PENCIL_AND_CAUSTIC_DESCENT.md](SFC_POINTED_PENCIL_AND_CAUSTIC_DESCENT.md).
- [SFC_DIFFERENTIAL_BRANCH_OBSTRUCTIONS.md](SFC_DIFFERENTIAL_BRANCH_OBSTRUCTIONS.md).
- [SFC_INVERSE_POWER_CONTACT_AUDIT.md](SFC_INVERSE_POWER_CONTACT_AUDIT.md).
- [SFC_RAMIFIED_PULLBACK_CLOSURE.md](SFC_RAMIFIED_PULLBACK_CLOSURE.md).
- [SFC_CUBIC_UNIT_DEFORMATION_TRANSFER.md](SFC_CUBIC_UNIT_DEFORMATION_TRANSFER.md).
- [SFC_PRIME_POWER_SCALAR_TRANSFER.md](SFC_PRIME_POWER_SCALAR_TRANSFER.md).
- [SFC_PRIME_POWER_UPPER_RESPONSE.md](SFC_PRIME_POWER_UPPER_RESPONSE.md).
- [SFC_PRIME_POWER_DIVIDED_DIFFERENCE_CLOSURE.md](SFC_PRIME_POWER_DIVIDED_DIFFERENCE_CLOSURE.md).
- [SFC_IDENTICAL_MONOMIAL_PAIR.md](SFC_IDENTICAL_MONOMIAL_PAIR.md).
- [SFC_OVERLAPPING_SWAP_PAIRS.md](SFC_OVERLAPPING_SWAP_PAIRS.md).
- [SFC_CYCLIC_MONOMIAL_BLOCKS.md](SFC_CYCLIC_MONOMIAL_BLOCKS.md).
- [SFC_CYCLIC_POLYNOMIAL_TILTS.md](SFC_CYCLIC_POLYNOMIAL_TILTS.md).

These are internal written proofs with exact computational controls. The
finite controls do not establish their universal quantifiers; the written
arguments do. No literature-priority, external-refereeing or formal-verification
claim is made.

### 1.1 A proved enlargement for every finite support

The support-size direction now has a quantitative theorem. Let the base
monomials `M_j` have the same multiset of positive exponents, with arbitrary
overlaps, and let `X_j` be private variables. Put

\[
 F=\sum_{j=1}^s c_jX_j^D M_j,\qquad D\ge6s.
\]

For the common endpoint moment `H_m=Lambda(M_j^m)`, define
`P_m=Lambda(F^m)/(H_m(Dm)!)`. Its pure coefficients are one. At a
composition `k` of `m`, integer factorial convexity and a labelled-coloring
injection give

\[
 0<[c^k]P_m\le\binom{m}{k_1,\ldots,k_s}^{1-D}.
\]

Consequently `P_m` is uniformly close to the power sum `sum c_j^m` on
`max |c_j|=1`, with error at most
`binom(m+s-1,s-1)m^(1-D)` for `m>=2`; at `m=1` the error is zero.
The inspected classical Turán power-sum theorem supplies an index in
every `s`-term window at which the power sum has modulus at least
`2(s/(8e(n+s-1)))^s`. The exact rational comparison in the proof shows
that `D=6s` makes the error smaller, uniformly in `s,n`. The worst bound
used is `243/256`, with a separate smaller bound for the initial window.
This is a proof for all orders, not an eventual or sampled assertion.

An arbitrary support can be balanced explicitly: for column `j`, add private
variables carrying all exponent entries from every other column. All resulting
columns then have the same exponent multiset. Add `X_j^(6s)` and apply the
theorem. Characteristic-zero descent uses only the finitely generated field
of a hypothetical counterexample. No positivity of its coefficients is used.

This gives a known regular endpoint for every support. It does not show that
removing the new private variables preserves regularity. In particular,
identifying all private marker variables into one destroys the multinomial suppression
that proves the theorem. That is the precise outstanding transfer step for
this route.

The same proof permits nonzero rational polynomial templates `P_j` which
are nonnegative on the positive orthant and have equal pure factorial
moments. Weighted AM–GM gives the required mixed-moment bound. Coordinate
permutations of one nonnegative polynomial are a concrete sufficient
condition; negative individual monomial coefficients are allowed, for
example in sums of squares. Fresh private markers distinguish the blocks,
so the active-block window is no longer than the actual support-size
window. This positivity assumption concerns fixed auxiliary templates,
not arbitrary coefficients of the final polynomial.

### 1.2 Shared rows and actual merger preservation

The [separating-row theorem](SFC_SEPARATING_ROW_AMPLIFICATION.md) replaces
private markers by partition layers. A column chooses one symbol in each
layer, with symbol-variable exponent `D_l`. A set of layers separates if
its joint symbols distinguish all columns. Choose nonnegative rational
weights `w_B` on separating sets, with

\[
 \sum_{B\ni l}w_B\le D_l,\qquad W=\sum_Bw_B\ge6s.
\]

An injective map from colored words to their layer words yields the exact
bound `H_marker(k)/H_marker(m e_j)<=C_k^(-W)`. The existing quantitative
comparison then proves SFC at every order. If the final symbol coarsening
retains this certificate, every intermediate coarsening does as well.
This proves actual factorial row-merger preservation for the stated
class, including families in which every variable is shared by several
terms. It does not permit merging layers with one another or with base
variables.

For example the four binary words `00,01,10,11`, with two exponent-24
layers, use four shared variables. The redundant code `000,011,101,110`
works with exponent 16 in each layer and weight 8 on each separating
pair. The theorem covers arbitrary sizes and codes; these are examples.

The same note proves an intrinsic result: balanced exponent columns with
distinct binary zero patterns satisfy SFC if their minimum positive
exponent is at least `6s ell`, where `ell` selected rows distinguish all
columns. Conditional factorial convexity gives suppression exponent
`D/ell`; no additional marker variables are introduced. Always
`ell<=s-1`.

The exact boundary is also proved. For two balanced columns `a,b`, the
composition `(m-1,1)` has normalized mixed moment asymptotic to
`C_(a,b) m^(-beta_(a,b))`, where
`beta_(a,b)=sum_(a_i=0)b_i` and `C_(a,b)>0`. Any positive uniform
suppression exponent is at most every such beta. Equal zero patterns
give beta zero, so fixed exponent dilation cannot extend this particular
argument to those supports.

## 2. A precise contact target for an unbounded support family

Let

\[
 F(Y)=\left(\prod_{i=1}^sY_i^{r_i}\right)\sum_{i=1}^sc_iY_i,
 \qquad r_i\ge0,\quad c_i\ne0.
\]

Set

\[
 P(z)=\prod_i(1-c_i z),\qquad
 A(z)=\prod_i(1-c_i z)^{r_i},\qquad a(z)=zA(z),
\]
\[
 H(0)=0,\quad H'=A/P,\qquad B=H\circ a^{-1},
 \qquad D=P+zP A'/A.
\]

The primitive `H` is polynomial if every exponent is positive, and otherwise
is a well-defined formal series with rational derivative. In both cases
`a'(0)=H'(0)=1`, `deg D=s`, and

\[
 B'(a(z))=\frac1{D(z)},\qquad
 \frac{\mathcal L(F^m)}{m!\prod_i(r_i m)!}
 =(m+1)[t^{m+1}]B(t).
\tag{1}
\]

The desired all-order theorem for this entire family is exactly the bound

\[
 \boxed{\operatorname{ord}_0(H-b\circ a)\le\deg b+s}
\tag{2}
\]

for normalized polynomial truncations `b`. This is an **open contact
statement**, not an imported theorem. The retained first bound was

\[
 \operatorname{ord}_0(H-b\circ a)\le(2n-1)s+1,
 \qquad \deg b\le n.
\tag{3}
\]

It is independent of all the exponents. Its proof differentiates (1), clears
only the fixed-degree factors `P,D`, and obtains a nonzero numerator of
degree at most `(2n-1)s-n` for `B^(n+1)(a)`. A nonzero polynomial cannot
vanish to order exceeding its degree. The numerator cannot be identically
zero, because that would make `B` polynomial and force the impossible
polynomial identity `D B'(a)=1`.

For `n=1`, (3) equals the desired `s+1` contact bound. Directly, vanishing of
moments `1,...,s` would give `H-a=O(z^(s+2))`, while
`(H-a)'=H'(1-D)` has order at most `s`. Thus the initial support-size window
is proved for every choice of exponents and coefficients in this family.

The [reduced-divisor proof](SFC_REDUCED_DIVISOR_CONTACT_BOUND.md) now improves
(3). Group equal coefficients into multiplicities `p_j` and exponent totals
`R_j`, put `J=prod_(R_j!=p_j)(1-c_jz)`, and define

\[
 \nu=\deg\operatorname{rad}(D),\qquad
 T_0=\deg\left(J/\gcd(J,\operatorname{rad}(D))\right).
\]

The sharper bound is

\[
 \operatorname{ord}_0(H-b\circ a)\le n\nu+(n-1)T_0+1.
\tag{3a}
\]

There is also a bound independent of the starting order. Put
`d_a=1+sum r_i` and `e=s-ceil(s/d_a)`. The
[uniform algebraic-window proof](SFC_ALGEBRAIC_UNIFORM_WINDOWS.md) gives

\[
 \operatorname{ord}_0(H-b\circ a)\le n+C_{s,R},\qquad
 C_{s,R}=d_a e+\frac{s d_a(d_a+1)}2.
\tag{3b}
\]

It uses the actual algebraic observable `rho=B'=1/D(a^{-1})`. A resultant
clears its poles, and an exact connection matrix gives a minimal scalar
differential operator with bounded coefficient degrees. Every conjugate
of `rho` tends to zero at infinity, so the operator has no nonzero
polynomial solution. Applying it to the truncation proves (3b), including
mixed zero/positive exponents. No generic-rank or positivity inference is used.

If every group satisfies `R_j>=p_j`, then `H` is polynomial. Subtracting
its constant trace gives a differential operator of order at most `d_a-1`
with no nonzero rational solution, improving (3b) to the fixed window
`s d_a(d_a-1)/2`. Each of these proved bounds may be combined by taking
the smallest. None supplies the sharp length `s` in all cases.

The initial window also has a complete determinant formula. With the
standard normalization `Res(c1,c2^2,...,cs^s)=1`,

\[
 \operatorname{Res}(g_1,\ldots,g_s)=
 \prod_{\varnothing\ne I\subseteq[s]}
 \left(1+\sum_{i\in I}r_i\right)^{(|I|-1)!(s-|I|)!}.
\tag{3c}
\]

The [proof](SFC_INITIAL_RESULTANT_AND_SHIFTED_WINDOWS.md) uses a triangular
change to weighted elementary symmetric polynomials and multiplicativity
of the resultant, proving all multiplicities and the positive normalization.
Nonnegative exponents make every factor nonzero. This certifies the initial
window, not an arbitrary starting order. A negative formal-exponent control
shows that the same subset-resonance test cannot simply be shifted.

There is useful structure beyond these estimates:

\[
 \frac{dH}{a}=\frac{dz}{z\prod_i(1-c_i z)},
\tag{4}
\]

whose poles and residues are independent of the exponents. The exponent
information remains in the powers of `a`; it has not disappeared from the
moments. The sharper bound must use the prescribed residues or the special
zero-fiber profile, rather than treating `H` as an arbitrary observable.

There is an explicit algebraic form of the remaining target. With
`E=P A'/A`, define the actual derivative numerators

\[
 M_2=-D',\qquad
 M_{k+1}=PD M_k'
 +\bigl(P'D-(k-1)ED-(2k-1)PD'\bigr)M_k.
\]

The sharp statement for this family is equivalent to
`z^s` not dividing `M_(n+1)` for every `n>=1`. Indeed the denominator in
the exact formula for `B^(n+1)(a)` is a unit at zero. Thus the next proof
can target the first `s` coefficients of a specified polynomial recurrence,
with the prescribed integer residues in `E/P`. Merely renaming this
nondivisibility condition does not prove it. The reduced-divisor proof
bounds the entire numerator degree, which it proves is exactly the stated
value. A further universal leading-degree cancellation cannot close that
argument. The fixed differential-operator bound is also larger than the
sharp target; it excludes long gaps without proving this local nondivisibility.

An exact rank-`s` residue module with invertible shifts is also available.
Its consecutive residue vectors can nevertheless become dependent at
admissible complex coefficients and recover independence at the next order.
The proof file gives the exact specialization `t=(1+3i)/5`. Consequently,
rank alone does not justify scalar propagation. A successful argument must
use the distinguished residue at zero, which represents the actual
factorial moment.

### 2.1 A valid redistribution corollary

Some unequal exponents do transfer. Partition variables into blocks of
sizes `p_j` carrying the same coefficient `d_j`, put `S_j=sum_i X_(j,i)`
and `R_j=sum_i r_(j,i)`, and let
`F=(prod X_(j,i)^r_(j,i)) sum_j d_j S_j`. Rising-factorial Vandermonde
within each block gives the exact identity

\[
 \frac{\Lambda(F^m)}{m!\prod_{j,i}(r_{j,i}m)!}
       =[z^m]\prod_j(1-d_jz)^{-R_jm-p_j}.
\]

Indeed the contribution of block `j` receiving linear degree `k` is
`prod_i(r_(j,i)m)! (R_jm+p_j)_k`; the outer multinomial expansion
then gives the displayed coefficient. Thus redistribution of the
individual exponents within a block preserves the entire zero pattern
when its total `R_j` is fixed.

If `R_j=q p_j` for one integer `q>=1` on every active block, comparison
with the Gamma-shape-`p_j` phase
`A(Y)=prod_j Y_j^(q p_j) sum_j d_jY_j` gives

\[
 \Lambda(F^m)=
 \frac{\prod_{j,i}(r_{j,i}m)!}{\prod_j(p_j)_{q p_jm}}
                       \Lambda_{\mathbf p}(A^m).
\]

All displayed factors are nonzero. The
[weighted radial theorem](WEIGHTED_RADIAL_LINEAR_SFC_TRANSFER.md)
therefore supplies an all-order window of length equal to the number
of distinct nonzero `d_j`, no greater than the actual support size
`sum_(d_j!=0) p_j`. Inactive blocks are independent monomial factors
and require no balance. This is an explicit corollary of that theorem,
not a separate closure mechanism. Different coefficient values cannot
be merged: `Lambda(Y^2(3X-Y))=0`, whereas `Lambda(XY(3X-Y))=4`.

Taking powers does not propagate the initial-window theorem: the
support of `F^k` is `binom(k+s-1,s-1)`, and its moments sample the
orders `km` instead of a consecutive block starting at `k`. Both
changes must be accounted for in any proposed amplification.

### 2.2 Exact ramified pullback and common-factor descent

The [ramified pullback theorem](SFC_RAMIFIED_PULLBACK_CLOSURE.md) now
preserves the entire distinguished moment sequence, not just its first
window. For grouped seed data

\[
 Q_0=z\prod_j(1-c_jz)^{p_j},\qquad
 a_0=z\prod_j(1-c_jz)^{R_j},
\]

take a polynomial `h` of degree `d` with
`h(0)=0,h'(0)=1,h'|Q_0(h)`. Equivalently all its finite critical values
are roots of `Q_0`. Put

\[
 Q=Q_0(h)/h',\qquad a=a_0(h).
\]

At a preimage `beta` over `1/c_j`, with local degree `e`, the target
coefficient `1/beta` has multiplicity `e(p_j-1)+1` and aggregate
common exponent `eR_j`. A nonzero preimage of zero has coefficient
multiplicity one and exponent `e`. Arbitrarily allocating each aggregate
exponent among its equal-coefficient variables produces an ordinary
monomial-linear phase of actual support `d sum_j p_j`.

For every `m>=0`, the exact identity is

\[
 \frac{\Lambda(F_h^m)}{m!\prod_i(mr_i^{\rm target})!}
 =\operatorname{Res}_0a_0(h)^{-m}\frac{h'\,dz}{Q_0(h)}
 =\frac{\Lambda(F_0^m)}{m!\prod_i(mr_i^{\rm seed})!}.
\]

Thus every zero and every seed window transfer, over a finite splitting
extension when necessary. The mechanism is substitution at the actual
residue at zero; it makes no inference from the rank of a larger residue
module. Successive admissible covers compose, with multiplicative support
size. The second cover's admissibility is checked against the updated
root set.

For example, for arbitrary integers `q,d>=2,b>=1`, use
`P_0=1-z^q,A_0=P_0^b` and `h=1-(1-z/d)^d`. The target has `dq` terms,
common exponents `db`, `d(q-1)` copies of `b`, and `d-1` copies of one.
Its complete normalized zero pattern is

\[
 g_m=0\ (q\nmid m),\qquad g_{qk}=\binom{(bq+1)k}{k}\ne0.
\]

Its optimal window is `q`, and every variable occurs in every monomial.
This supplies unequal-exponent families outside the identical-zero-pattern
limitation of the preceding suppression argument. The coefficients are
constrained by the cover; they cannot be made arbitrary in this statement.
The note also gives a conditional plane-tree construction via specified
Shabat polynomials and an elementary one-term seed.

The reduction has an exact limit. Writing `S` for support and `R` for
the total common exponent, every degree-`d` pullback satisfies

\[
 S=dS_0,\qquad R+1=d(R_0+1),\qquad d\mid\gcd(S,R+1).
\]

When `H` is polynomial, every existing common polynomial right factor
of `a,H` descends to another canonical grouped frame; a rational common
factor can be normalized to a polynomial one. This classifies the
composition reduction without asserting a factor exists. In particular,
`r=(1,1,2)` has degrees `deg a=5,deg H=2` for any nonzero coefficients,
and admits no nontrivial descent. The current exact certificates give its
three-moment window at starting orders 1,2,3. For all starting orders, the
reduced-divisor and trace-centered arguments give `min(3n,30)`; the desired
length 3 at arbitrary later orders remains open here. General formal
observables also retain the sharp target. The degree test does not imply SFC failure; known equal-exponent
families also pass this primitiveness test.

An exact counterexample further shows that the weaker condition
`a'=H'D`, `deg D=s`, by itself does not imply sharp contact: the note
exhibits `n=s=2` and contact five, whose shared critical point maps outside
the required zero fiber. This rules out dropping the root alignment from
the next contact argument; it is not an SFC counterexample.

### 2.3 Rational block averages now close under a size condition

Let the active coefficient groups have `p_j=b v_j` variables and total
common exponents `R_j=a v_j`, for positive integers `a,b,v_j`. Write
`h` for the number of distinct active coefficients, `V=sum v_j`,
`d=1+aV`, and `s=bV`. Within each group the exponents may be distributed
arbitrarily, including zeros. Part C of the
[uniform-window proof](SFC_ALGEBRAIC_UNIFORM_WINDOWS.md) proves actual
support-size SFC at every starting order whenever

\[
 \boxed{2\left(b-\left\lfloor\frac bd\right\rfloor\right)\ge ahd.}
\tag{5}
\]

The exact bridge is a power of a polynomial inverse. For
`R(z)=prod_j(1-c_jz)^(v_j)` and `V_*=(zR(z^a))^{-1}`,

\[
 [t^{am+b}]V_*^b=\frac b{am+b}\,g_m.
\]

Trace centering bounds contacts of `V_*^b` by the approximant degree plus
`(d-1)floor(b/d)+ahd(d-1)/2`. The inverse power has only one coefficient
class modulo `a`; inequality (5) makes that constant no larger than `as`.
Thus a zero block of `s` moments would exceed the proved contact bound.
This is an all-start proof, not a finite search.

For example, for arbitrary `c,d` not both zero,

\[
 X_1Y_1\bigl(c(X_1+X_2+X_3+X_4)
             +d(Y_1+Y_2+Y_3+Y_4)\bigr)
\]

satisfies SFC. With two distinct nonzero group coefficients the proof gives
a seven-moment window for its eight actual terms; the one-group boundaries
have every moment nonzero. More generally, take any `h>=2`, `a=1`,
`v_j=1`, and `b=(h+1)^2+1`. These phases satisfy the theorem, including all coefficient faces. On the
full active stratum, `gcd(hb,h+1)=1`, so they cannot be nontrivial
support-multiplying polynomial pullbacks of smaller canonical frames.
This is an infinite primitive family, with coefficient equality required
inside each group; it does not close independent coefficients or arbitrary
unequal positive exponents.

The [inverse-power audit](SFC_INVERSE_POWER_CONTACT_AUDIT.md) separately
records why the sharper target
`ord(V^b-Q)<=deg Q+b #CV^times(A)` remains open. For a genuine truncation
of `(z-z^2)^{-1}` squared, the proposed divided-polynomial construction has
five nonzero critical values when its claimed budget is four. The full
critical-value resultant certifies this failure. It refutes that shortcut,
not the sharp target or SFC.

## 3. Determinant amplification that is uniform in the order

Define

\[
 Q_n(x)=\sum_{k=0}^n(-1)^k\binom{n+2k}{3k}x^{n-k},
 \qquad R_n=Q_n\pmod3.
\]

For monic `P,T` of degrees `n,n+1`, respectively, with integral coefficients
in a complete algebraically closed valued extension of `Q_3`, assume
every coefficient error in `P-Q_n,T-Q_(n+1)` has valuation greater than one.
Then

\[
 \boxed{v_3\operatorname{Res}(P,T)
       =\deg\gcd(R_n,R_{n+1}).}
\tag{5}
\]

Thus the same precision works at every order even when the residual common
factor grows. The proof does not bound its degree. At each common residue
root, the lower model derivative is a unit. Lift that root to the integer
model in an unramified extension. An exact cubic trace identity gives
valuation exactly one for the next model value there. A perturbed root
stays closer than that scale. Each common residue disk therefore contributes
exactly one to the resultant valuation; other disks contribute zero.

The perturbed roots themselves are allowed to be ramified. Only the
integer-model lift is used in the unramified trace step. This distinction
is indispensable.

For rational units `u_L/u_0=zeta^L mod9`, the complete moment kernel
`Lambda_u(X^L)=L!u_L` admits a common normalization to this model. It follows
that every `X(A+BX^2)` has a nonzero weighted moment in each positive pair.
This transfers to ordinary multivariate families. For example, put

\[
 S=X+3Y-3Z.
\]

For every `r>=0`, every nonzero phase

\[
 S^{3^r}\bigl(A+B S^{2\cdot3^r}\bigr)
\tag{6}
\]

has a nonzero ordinary factorial moment in every positive pair. The proof
allows arbitrarily many linear-form variables when one coefficient is a
three-adic unit, the others are divisible by three, and the sum of those
other coefficients is divisible by nine. It uses the exact complete
homogeneous-polynomial identity for the moments of `S`, not an arbitrary
linear substitution preserving the factorial functional.

The precision is sharp for unrestricted deformations: congruence merely
modulo three fails. Moreover, the direct all-prime analogue modulo `p^2`
is false. At `p=5`, the model orders `156,157` admit integral perturbations
modulo `25` with a common root. The full rational unit-kernel realization
is retained in the certificate. These are counterexamples to proposed
transfer rules, **not ordinary SFC counterexamples**.

## 4. The prime-power gap problem is reduced to its actual upper response

Write

\[
 T_{n,\gamma,q}(t)=\sum_{k=0}^n\binom nk(\gamma)_{qk}t^k,
 \qquad q=p^s.
\]

For positive rational `gamma` integral at `p`, the full scalar-leaf argument
retains all actual seeds on the closed disk. The final normalized polynomial
is a unique restricted family

\[
 Q(X;v)=Q_0(X)+qv\sum_j\eta_j(v)G_j(X),
 \qquad \eta_j\in\mathbb Z_p\langle v\rangle.
\]

In particular, with `B_s=(q-1)/(p-1)-s`,

\[
 Q(-qn;v)=Q_0(-qn)+p^{B_s+1}vA(v),
 \qquad v_p(Q_0(-qn))=B_s.
\tag{7}
\]

This is the uniform substitute for recomputing one support's matrix. The
mixed-offset pair uses an explicit integral linear functional `T_h` on
the same polynomial family. The complete transfer gives sufficient
nondegeneracy conditions in terms of `T_h(Q_0)`, `T_h(G_j)` and
`T_h(Q_0')`. The subsequent upper-response proof also handles vanishing
or cancellation of the actual first response.

Offsets one and two close for all these parameters. Offset three has at
most two distinct common roots, and the retained squarefreeness theorem
gives gcd degree at most two. A possible irreducible quadratic remains a
real obstruction to completing that argument; a two-root bound alone does
not imply coprimality.

Put `C=gamma+q(n+1)`, `kappa=q(n+1)/(C)_h`, and
`D=min_j v_p(T_h(G_j))`. The stronger uniform exclusion is

\[
                         v_p(\kappa)+2s+D>0,
 \qquad1\le h\le q.
\tag{8}
\]

If all the functional values vanish, the pair is coprime as well.
The scalar auxiliary has a nonzero constant term of valuation `B_s`;
under (8), every coefficient above degree one has larger valuation.
Only the actual linear coefficient can compete. A common root would
therefore be unique on the whole integral root disk. Rational conjugacy
and the exact rational Gamma Gram factorization then exclude it. No
positivity of the phase over a complex field is assumed.

### 4.1 Divided differences cancel the factorial loss

The upper functional has the exact polynomial form

\[
 \mathcal T_h(G_j)=(-1)^{q+1}(X)_h[C,C-q,A_r,A_j],
 \qquad A_j=j-q+1,
\]

where brackets denote an order-three divided difference, including
coincident nodes. Its product rule gives

\[
 D\ge B_p(h):=\sum_{\ell\ge1}
                  \max\{\lfloor h/p^\ell\rfloor-3,0\}.
 \tag{9}
\]

Each product-rule term has four consecutive factorial blocks. Counting
multiples of every `p^ell` in those blocks proves (9) for arbitrary
integral nodes. This supplies the divisibility that was lost by treating
the upper functional as merely integral.

For `a=p^w<=q=p^s`, `h=a`, and `gamma=an+beta` with integer
`1<=beta<=a`, exact denominator cancellation gives

\[
 v_p(\kappa)=s-v_p(a!)-\delta,
 \quad B_p(a)\ge v_p(a!)-3w+2,
\]

where `delta=1` only when `p=2,w=s`, and is zero otherwise. Thus

\[
 v_p(\kappa)+2s+D\ge3(s-w)+2-\delta>0.
 \tag{10}
\]

The offset-one case handles `w=0`. This proves the full divisor-of-a-
prime-power family in the table, for all orders and coefficient fields.
The proof note identifies the overlap with earlier relative-gap and
dyadic results; it does not claim every case in the family is new.

For every ordinary `1<=a<q`, write `a=p^w A`, `p` not dividing `A`,
and `N=n+1`. The same mechanism gives

\[
 v_p\binom{(A+p^{s-w})N-1}{A-1}<K,
 \quad
 K=3s-\sum_{\ell\ge1}\min\{\lfloor a/p^\ell\rfloor,3\}\ge3
 \quad\Longrightarrow\quad\gcd=1.
 \tag{11}
\]

For `A=1` this holds at every order. For `A>1`, the surviving congruence
classes remain unresolved by this criterion. They are not counterexamples.

### 4.2 What higher responses still leave open

The actual leaf coefficients satisfy a stronger slope-one filtration.
It bounds the mixed-pair gcd by an even integer depending only on
`v_p(kappa)+2s+D`. At prime gaps, a formal Artin--Schreier translation
identifies a guaranteed unit response at degree `h-2`, giving

\[
 \deg\gcd(T_{n,\gamma,p},T_{n+1,\gamma+h,p})
                         \le2\lfloor(h-1)/2\rfloor.
 \tag{12}
\]

This is independent of the moment order and the denominator depth.
The actual `p=5,gamma=11,h=5` degeneracy is retained as a control;
the proof does not assume a nonzero first variation. For proper prime
powers the naive translation is false already in its first coefficient.
The proof records that failure explicitly. An imaginary quadratic can
split at a chosen prime, so a two-root bound or a split local picture
does not by itself complete rational descent. Explicitly, for every
odd prime `p` and `r>=1`, `X^2-X+p^r` has negative discriminant over
the reals and is irreducible over the rationals, but Hensel lifting at
`0,1 mod p` gives two roots in `Q_p`. This refutes the proposed
split-to-rational shortcut, not SFC.

## 5. Two further structural closures with explicit hypotheses

For disjoint identical copies of a nonconstant monomial `M`, put

\[
 a_k=\Lambda(M^k)/k!,\qquad
 P_n(z)=\sum_{k=0}^n a_k a_{n-k}z^k.
\]

These are the normalized moment polynomials for `A M(X)+B M(Y)`.
For total monomial degree at least two, the sequence
`delta_k=k!/a_k` is log-concave. Consequently the coefficients of
the coefficientwise reciprocal of `P_n`, after division by `binom(n,k)`,
are log-concave and symmetric. The inspected coefficient version of
Mao--Wang--Zheng's reciprocal interlacing theorem separates the zeros of
`P_n` and `P_(n+1)`. Degree one is the elementary geometric-sum case.
This is a new application in this investigation of a credited external
theorem; its precise version, proof audit and hypotheses are in the
[monomial-pair note](SFC_IDENTICAL_MONOMIAL_PAIR.md).

For arbitrarily many terms, the separate cyclic-block theorem uses only
an elementary formal Kaluza coefficient identity. A nonconstant monomial
has positive log-convex coefficients `a_k`, so all coefficients of
`log(sum a_k t^k)` are strictly positive rational numbers. Root-of-unity
filtering therefore gives an exact noncancellation certificate.
For disjoint cycles of lengths `q_1,...,q_r`, with a common scalar and
positive rational orbit scales,

\[
                \Lambda(F^n)=0
          \quad\Longleftrightarrow\quad
                n\notin\langle q_1,\ldots,q_r\rangle.
 \tag{13}
\]

This identifies the whole zero set and proves SFC for those sums, their
powers and their disjoint products. It does not require convergence of
the factorial series. The assumptions cannot be discarded: positive
formal coefficients alone fail for `exp(t)`, whose nontrivial cyclic
product is identically one. Arbitrary independent coefficients remain
outside this cyclic argument.

### 5.1 Two actual merging extensions

The [overlapping-pair theorem](SFC_OVERLAPPING_SWAP_PAIRS.md) permits
mirrored rows `(b_i+d_i,b_i)` and `(b_i,b_i+d_i)` with arbitrary shared
exponents, provided the total private degree is at least two. Its
binomial-normalized reciprocal coefficients stay log-concave because the
private factors alone control the derivative of their ratio. Equal rows
factor out, and specified symmetric row mergers retain the theorem. The
imported reciprocal-interlacing theorem is the same inspected version used
for identical disjoint pairs; no new analytic input is assumed.

The [polynomial cyclic theorem](SFC_CYCLIC_POLYNOMIAL_TILTS.md) replaces
monomial templates by nonnegative rational polynomials and permits a common
monomial overlapping all copies, with equal exponents on corresponding
coordinates of a cycle. At each fixed target order `n`, that monomial
changes the Gamma shapes to `na_i+1`. The connected-permutation expansion
gives strictly positive formal cumulants for these shapes. The cyclic
filter then gives exactly the same semigroup zero set as before. The
proof fixes `n` before changing the shapes; it does not invent one global
generating function from order-dependent kernels. Exposed faces control
support collisions in powers. Symmetric variable identifications within
copies preserve the complete zero pattern.

There is also a new obstruction on an actual integer support. For the
features `X_jY_j`, the reciprocal moment form is
`(sum c_j)^m/(m!)^2`, which is real stable. Merging the `X_1,X_2` rows
keeps all monomials distinct but changes the quadratic reciprocal to
`((sum c_j)^2-c_1c_2)/4`, which is not real stable or Lorentzian. Its
two-variable restriction fails normalized log-concavity. Thus those
auxiliary properties cannot supply universal row-merger preservation,
even though the two-term endpoints themselves satisfy SFC.

## 6. What would reach arbitrary supports

For a distinct integer support matrix `A=(a_1,...,a_s)`, let

\[
 M_m(c)=\mathcal L\left(\left(\sum_jc_jX^{a_j}\right)^m\right),
 \qquad I_{A,n}=(M_n,\ldots,M_{n+s-1}).
\]

Unrestricted SFC is the assertion
`sqrt(I_(A,n))=(c_1,...,c_s)` for every such matrix and every `n>=1`.
This regular-system formulation is already in Edo--van den Essen,
[*The Strong Factorial Conjecture*, Remark 2.9(b)](https://arxiv.org/pdf/1304.3956v2).
Calling the ideal regular is a restatement until a proof of regularity is
provided.

The retained
[unrestricted mechanism audit](UNRESTRICTED_SFC_MECHANISMS.tex) supplies a
concrete sufficient induction scheme with two open preservation lemmas:

1. Adjoining a unit row to a nonconstant support column preserves the
   projective regularity of each consecutive moment block. Its coefficient
   multiplier is `k_j!`.
2. Merging two integer exponent rows preserves that regularity whenever the
   resulting columns remain distinct. Its coefficient multiplier is
   `binom((u+v) dot k, u dot k)`.

Starting with independent coordinate monomials, these operations construct
every integer support. The retained audit proves that construction and its
base case; it does not prove either preservation lemma. Positivity of a
Gram matrix, arbitrary moment-measure interpolation and generic regularity
do not supply them. Explicit countercontrols to those broader arguments
are already recorded there.

The new amplification theorem adds a second sufficient route: every
support has an explicit enlarged regular family, so it remains to justify
removing those particular private decorations. This implication is not
formal and is not proved here. The known broader positive-measure
obstructions prevent replacing the missing integer-toric argument by
positivity alone.

The new [exact merger obstruction](SFC_ROW_MERGING_OBSTRUCTIONS.md) makes
this limitation explicit. Exchangeable positive features `U,V` can have
`E((1+3i)U+(1-3i)V)^2=E((1+3i)U+(1-3i)V)^3=0` while retaining a
strict polynomial Gram form. The privately amplified features
`X^D U,Y^D V` satisfy every consecutive two-moment window for `D>=12`,
by the same quantitative theorem. Merging the private variables into one
`Z` restores those exact two zeros. The endpoint columns remain distinct.
The other features do not have the ordinary independent-exponential law,
so this is not an ordinary factorial counterexample. It rules out a
de-amplification proof using only equal pure moments, Hölder bounds,
strict Gram positivity and integer decorations. The proof also records
the exact inverse Beta operator and correlated-Gaussian/Jacobi formulas.

The [de-amplification bridge audit](SFC_DEAMPLIFICATION_BRIDGES.md) now
gives an exact finite-range arithmetic approximation. Fix a maximum
moment order `M`, choose `p>M`, and put
`D_r=sum_(i=1)^phi(p^r) p^(ri)`. Then
`(D_r k)!/(D_r!)^k=1 mod p^r` for every `k<=M`.
Actual decorated forms, divided by `(D_r!)^m`, therefore approach the
undecorated forms. The prime and finite window are fixed before taking
this limit. No fixed prime works at all orders.

Two separate gaps remain. Full-rank minors can tend to zero; an eventual
bound on their valuations is equivalent to full rank of the limiting
matrix, so it is not supplied by regularity at positive exponents.
Moreover, exactly balanced lifts satisfy a closed pure-coefficient
identity even after diagonal coefficient rescaling and degreewise
normalization. The unbalanced support `X,Y^2` violates that identity by
`-16` at order two. The arithmetic construction removes equal new
private powers from a fixed base; it cannot silently remove the fixed
profile-balancing factors from an arbitrary base.

The same audit proves exact algebraic countercontrols to ordinary ideal
transport under inverse factorial maps, including a Laurent resonance.
It also supplies an unbounded family already satisfying the amplification
theorem whose reciprocal quadratic has two positive Hessian directions.
Neither invertible coefficient operators nor universal reciprocal stability
closes the remaining ordinary-ideal implication. General private-power
removal for distinct integer supports is equivalent to full SFC.

The current work therefore has three general obligations: sharpen the
observable contact bound, exclude the surviving upper-response factors
beyond the proved prime-power divisibility cases, and establish a valid
support-transfer preservation theorem. The divided-difference identity
and factor-degree bounds sharpen the second obligation. The new uniform
amplification, separating shared rows, overlapping pairs, symmetric cyclic
mergers and exact ramified pullbacks sharpen the third, with each exact
restriction retained. Common-right-factor descent isolates part of the
first obligation. The fixed algebraic windows and rational block-average
closure now improve it further, but the general sharp primitive contact
problem still requires proof. None proves unrestricted SFC. Each remaining obligation is stated independently of an individual
support.

## 7. Verification boundary

The proof files contain their exact replay commands and source pins.
The cubic determinant has a separately implemented FLINT resultant audit;
the primary checker uses standard-library exact arithmetic. Observable
identities, residue obstructions, stage normalizations and upper responses
have exact controls with their finite domains recorded in the receipts.
The ramified-pullback checker adds 45 exact polynomial covers, 51 direct
factorial-moment expansions including order-zero controls, and explicit
composition, degree and failed-contact certificates. The current continuation
adds six standard-library checkers: an all-support initial-resultant proof
with 960 complete-interpolation determinants, reduced-divisor cancellations,
an exact failed critical-value budget, rational-observable connections,
trace-centered minimal operators, and rational block-average identities.
Their receipts distinguish finite identity controls from written universal
proofs. Independent team replays and invalid-source rejection checks are
recorded in the manifest.

The existing manifest is extended with the current bytes and the new
scopes. The status registry and generated navigation have not been
promoted in this research pass. The pre-existing navigation/source-pin
failures recorded in the manifest remain separate from these replayed
results. No stale retained proof hash has been replaced merely to make an
old checker pass. No publication or remote repository update is made.
