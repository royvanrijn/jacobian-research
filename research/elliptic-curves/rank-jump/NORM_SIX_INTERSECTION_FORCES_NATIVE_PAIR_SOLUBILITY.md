# A norm-six trace class forces simultaneous rational solubility

**Triple follow-up:** the [minimal marked carrier calculation](MINIMAL_CARRIER_AND_RATIONAL_SPLITTING_OF_A_TWO_DIRECTION_BLOCK.md)
finds genus one for the +7 two-direction basis, but the three-cover relation
has intersection degree twelve. Its rational component recovers the known
parameter; the remaining degree-eleven component is irreducible. The
degree-one hypothesis for that triple is therefore excluded.

There is now a generic geometric explanation for one of the simultaneous
square events on the +8 fibre. The native trace classes of `0911e` and
`1795d` admit a norm-six difference after translation. Their corresponding
rational bisections consequently intersect once. A unique geometric
intersection of two curves defined over Q is rational, so both cover
conditions become soluble together.

The exact equations recover the previously observed parameter
\(t=-4112/1937\), without inserting it into the elimination. The construction
also forces a group relation: the two resulting points sum to a generic
section. It therefore explains **one** quotient direction, not two and
not the full +8 jump. This is a positive global-solubility mechanism with
an intrinsic dimension limit.

The pair and translate were first identified retrospectively. The general
norm-six criterion below uses only the generic marked lattice and bisection
classes; its discovery and this three-pair comparison are not a prospective
validation or an authorized new search.

## The generic criterion

Work on the smooth rootless elliptic K3 of the published R17 family, with
Mordell–Weil lattice M, all marked sections defined over Q, and no reducible
geometric fibres. Write \(h(v)=\langle v,v\rangle\). The existing
[bisection lattice proof](../../elkies-k3/BISECTION_COLLISION_SEARCH.md)
gives the class of a rational bisection with trace w as

\[
B_w=\left(\frac{h(w)-2}{4},\,2,\,w\right)
\quad\text{in }U\oplus M(-1).
\]

Thus

\[
B_v\cdot B_w=\frac{h(v-w)}2-2.
\]

Inversion and translation by a rational section S send a bisection of
trace w to one of trace \(2S-w\). These operations preserve the base and
extend to automorphisms of the smooth minimal K3. For two distinct native
quadratic fields, the two resulting bisections are distinct curves.

**Sufficient solubility criterion.** Suppose the two native trace classes
\(\tau_i,\tau_j\) have a representative R satisfying

\[
R\equiv\tau_i+\tau_j\pmod {2M},\qquad h(R)=6.
\]

Then \(S=(R+\tau_i+\tau_j)/2\) is a rational generic section, and

\[
B_{\tau_i}\cdot(S-B_{\tau_j})
=\frac{h(2S-\tau_i-\tau_j)}2-2=1.
\]

The two smooth rational curves have no common component. Their intersection
is a zero-dimensional Q-scheme of length one, hence a reduced Q-point.
Its two preimages give rational points on both native covers over the same
base parameter. If the branch divisors of the two covers are disjoint,
their fibre product is smooth, so this is a rational point on the smooth
simultaneous genus-one carrier itself. This proves global solubility, not
merely local solubility or a rational divisor class.

The criterion does not guarantee that the intersection parameter has a
smooth elliptic fibre. When it does, the constructed points satisfy

\[
P_i+P_j=S(t_0),\qquad [P_i]=-[P_j]\text{ in }E_{t_0}(\mathbf Q)\otimes\mathbf Q/G.
\]

Thus this construction provides at most one quotient direction. Its
nonvanishing requires a separate certificate. Large generic character rank
does not remove this specialization relation.

## The observed +8 pair, derived from generic equations

Use P1,...,P17 for the published generic basis. The previous exact
specialization relation transports to

\[
S=P_{12}-P_{14}+P_{15}-P_{16}.
\]

For the native traces \(\tau_G\) of `0911e` and \(\tau_F\) of `1795d`,

\[
2S-\tau_G-\tau_F=P_1+P_3+P_6+P_{12}=R,\qquad h(R)=6.
\]

This equality and norm use only integer vectors in the exact generic
height lattice. Both traces have norm ten. The norm-six coset has exactly
the two signed shortest representatives \(\pm R\); no other vector of
norm at most six lies in it. The opposite representative gives the
conjugate relation with translate \(\tau_G+\tau_F-S\).

The [equation protocol](NATIVE_PAIR_COLLAPSE_LOCUS_PROTOCOL.json) fixes this
pair and S. It reconstructs all required generic sections from the published
equations, forms \(Q=S-P_G(t,u)\), and imposes the residual chord for F.
Over \(u^2=q_G(t)\), its residual is

\[
a(t)+u b(t).
\]

On the explicitly certified open set where \(b\ne0\), a common point forces
\(u=-a/b\). Taking the norm gives \(a^2-q_Gb^2=0\). Imposing F's residual
x-quadratic removes the third intersection of its chord. Their polynomial
gcd, after excluding singular, branch and denominator-zero parameters, is

\[
\boxed{t+4112/1937}.
\]

The chord norm initially has degree 21; the residual-quadratic gcd already
has degree one, and saturation removes nothing further. The exact rational
functions for both cover roots and both point maps are retained. They
satisfy both square equations and the elliptic group relation in the finite
parameter algebra. The open-set calculation alone makes no assertion about
excluded places; the independent lattice intersection number one, together
with this exhibited point, proves that there are no additional geometric
intersections elsewhere for these two translated curves.

Only the 32 already frozen R17 parameters were evaluated. The sole hit is
`08234-009`. Its parameter is obtained from the generic intersection
polynomial, while its label is attached afterwards. The
[previous independent witness calculation](PAIRED_SOLUBILITY_AND_SPECIALIZATION_COLLAPSE.md)
proves that this common quotient line is nonzero. No new specialization,
rank search, or exceptional-point input was used by the elimination.

This distinction matters: an equation computable from generic data can still
have been chosen using retrospective evidence. This particular S was. The
norm-six existence criterion, however, can be tested directly on a supplied
pair of generic trace classes without exceptional points.

## Three fixed carrier comparisons

The [lattice protocol](NORM_SIX_CARRIER_SOLUBILITY_PROTOCOL.json) tests only
the three previously compared pairs. For their trace sum w, the lattice
\(2M+\mathbf Z w\) is the union of \(2M\) and \(w+2M\). Since M has minimum
four, every nonzero vector of norm at most six in this lattice lies in the
latter coset. This reduces each test to a small exact lattice calculation;
the full atlas or global short-vector population is not regenerated.

| Pair | Norm-six representative? | Criterion proves | Independently known carrier solubility |
|---|---|---|---|
| F,G = 1795d,0911e | Yes: exactly one pair ±R | Globally soluble | Yes; carrier Jacobian rank 2 |
| F,D = 1795d,11278 | No short representative | UNKNOWN | Yes; carrier Jacobian rank 3 |
| A,D = 030cb,11278 | No short representative | UNKNOWN | No: labelled nonzero Sha class |

The [global carrier results](GLOBAL_CARRIER_SOLUBILITY_AND_SPECIALIZATION.md)
and [A,D obstruction](NATIVE_PAIR_CARRIER_HAS_A_SHA_OBSTRUCTION.md) are
inherited from their prior certificates. F,D proves that the criterion is
not necessary. A,D is consistent with the sufficient criterion, but failure
of the norm-six test does not prove its insolubility. The norm residues of
the trace-sum cosets are respectively 2, 2, and 0 modulo four; A,D therefore
cannot have norm six in its coset at all.

Each worker is bounded at 60 seconds. A second implementation enumerates
the three cosets with exact rational LDL bounds and integer endpoint
rounding, independent of PARI's short-vector enumeration. It visits only
164, 204 and 208 recursion nodes, and recovers the same two, zero and zero
signed vectors. The same verifier checks the generic translate and rational
intersection with a separate explicit group law.

## Mechanisms, limits, and the next falsifiable test

1. **Solubility — strongest positive mechanism here:** a norm-six trace
   congruence forces a degree-one intersection of rational bisections,
   which forces simultaneous rational cover points. This is a genuine
   implication before specialization, not a square-value correlation.
2. **Incidence — its dimension limit:** that intersection imposes a generic
   sum relation, so it cannot explain two independent quotient directions.
   A graph of such relations at one parameter leaves at most one quotient
   direction per connected component, before further dependencies. More
   soluble cover labels do not by themselves constitute a larger block.
3. **Weak or excluded explanations:** norm-six absence is not a global
   obstruction; F,D is a counterexample to necessity. Carrier Jacobian
   rank and generic capacity remain insufficient explanations of the full
   jump. Five +8 witness directions remain beyond the exact rank-three
   quartet, as before.
4. **Next experiment:** the +7 fibre has a relation among three soluble
   covers while every pair remains independent. Form the curve traced on
   the K3 by a translated sum/difference of two of those cover points,
   and intersect it with the third bisection. Compute its generic degree,
   divisor class, and intersection number first. A degree-one intersection
   could force three rational lifts with only one relation, hence leave
   room for a two-dimensional soluble block. A higher intersection degree,
   a non-birational map hiding a lifting obstruction, or failure to recover
   the fixed parameter would falsify that proposed explanation. This is a
   bounded test on the existing triple, not a new parameter campaign.
5. **Information for Agent 1:** the norm-six condition is interpretable as
   a sufficient **carrier solubility** feature. Any future use must retain
   the sum-relation certificate and distinguish the intersection parameter
   from an arbitrary rational point on the same genus-one carrier. It is
   not a high-rank predictor. No visibility feature or search-policy change
   is proposed.

The desired large-jump chain remains incomplete. The new result explains
why a particular pair becomes soluble together, but its geometry itself
limits the gain to one direction. Extending degree-one intersection
constructions to larger, independent blocks is the next mathematical issue.

## Reproduction

From the repository root, with Sage 10.9 and PARI 2.17.3:

```sh
sage -python elliptic-curves/rank-jump/native_pair_collapse_locus.py check
sage -python elliptic-curves/rank-jump/norm_six_carrier_solubility.py check
sage -python elliptic-curves/rank-jump/verify_native_intersection_solubility.py check
```

Immutable evidence: [generic equation input](../../artifacts/generated-results/elliptic-curves/rank_jump_native_pair_collapse_locus_inputs_v1.json),
[intersection polynomial and maps](../../artifacts/generated-results/elliptic-curves/rank_jump_native_pair_collapse_locus_v1.json),
[three-pair lattice input](../../artifacts/generated-results/elliptic-curves/rank_jump_norm_six_carrier_solubility_inputs_v1.json),
[exact short cosets](../../artifacts/generated-results/elliptic-curves/rank_jump_norm_six_carrier_solubility_v1.json),
and [independent verification](../../artifacts/generated-results/elliptic-curves/rank_jump_native_intersection_solubility_verification_v1.json).
The scripts and inputs are hash-bound. All changes are confined to rank-jump
analysis; the active high-rank search and mathematical status registry are
untouched.

The concurrent commit `84aecacc` picked up these analysis scripts and
certificates before the rank-jump note was committed. Their exact hashes
remain the certified ones; no shared history was rewritten to separate them.
