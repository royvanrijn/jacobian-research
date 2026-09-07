# Exceptional branch-divisor code

For a finite set of exact square conditions over \(\mathbf Q(T)\),

\[
 z_i^2=g_i(T),
\]

the branch-divisor code records `div(g_i) mod 2` as a binary row. A closed
\(\mathbf Q\)-branch place represented by an irreducible factor of degree
`d` is one Galois-orbit column of weight `d`; this is equivalent to keeping
its `d` geometric conjugate places separately for every branch count below.
The place at infinity is an additional weight-one column whenever the finite
divisor has odd degree.

Every nonzero codeword gives an elementary quadratic quotient

\[
 z^2=\prod_i g_i(T)^{e_i},\qquad e_i\in\mathbf F_2.
\]

If its geometric branch count is `B>0`, its genus is `(B-2)/2`. A zero word
is reported separately as complete squareclass cancellation: geometrically it
is an unramified split condition, not a genus-zero double cover. For the
connected multiquadratic cover with code dimension `r` and branch-union size
`B`, tame Riemann--Hurwitz gives

\[
 g=1+2^{r-2}(B-4).
\]

The implementation also reports pairs of directions with common branch
support and code automorphisms. It keeps a compressed column for each closed
branch place, but expands its degree formally over \(\overline{\mathbf Q}\):
every geometric place in one irreducible factor has the same incidence
column. It reports the product of the symmetric groups on those equal columns.
This is the complete coordinate-permutation symmetry visible to the incidence
code; it does not assert that such permutations are induced by a rational
\(PGL_2\) automorphism of the base line or by an automorphism of the cover.

## Pinned calculation

Run:

```sh
PYTHONPATH=elliptic-curves:elliptic-curves/cas .venv/bin/python \
  elliptic-curves/cas/analyze_exceptional_branch_code.py
```

The generated artifact is
[`elliptic_exceptional_branch_code.json`](../../artifacts/generated-results/elliptic-curves/elliptic_exceptional_branch_code.json).
It uses only exact factorization over \(\mathbf Q\) and no rational-point
search.

The currently pinned inputs are the ten independent E22 directions `P13`
through `P22`, the eight Fermigier rank-20 directions, their combined
eighteen-direction code over the common Fermigier `T`-line, and the eight
certified exceptional directions each on ICARM 245 and ICARM 275. In every
one of these data sets, all direction rows are disjoint irreducible sextics.
Hence every code is the weighted direct sum of six-place rows:

| Input | code dimension | branch-union size | smallest quotient genus | low-genus words / cancellations |
| --- | ---: | ---: | ---: | ---: |
| Fermigier E22 | 10 | 60 | 2 | 0 / 0 |
| Fermigier rank-20 anchor | 8 | 48 | 2 | 0 / 0 |
| Fermigier E22 plus rank-20 anchor | 18 | 108 | 2 | 0 / 0 |
| ICARM 245 | 8 | 48 | 2 | 0 / 0 |
| ICARM 275 | 8 | 48 | 2 | 0 / 0 |

Thus this computation exhausts all `2^18-1=262143` elementary quotients of
the combined Fermigier code, in addition to the separate `2^10-1` E22 and
`2^8-1` anchor codes. Their genus histograms are exactly the binomial
distributions on branch counts `6k`; the combined code has `C(18,k)`
quotients of genus `3k-1`. No combination of three or more of these declared
directions produces a new genus-zero or genus-one quadratic quotient.

The artifact also contains two deliberately segregated supplemental inputs:
the 16 non-generic preimages observed in the declared height-200,000 search
at ICARM 243 and the 19 observed at ICARM 226. The analyzer reconstructs the
two six-root families and verifies each listed value is a rational square at
its stated anchor before factoring its branch divisor. Their individual codes
have dimensions 16 and 19, branch-union sizes 96 and 114, respectively, and
again have minimum quotient genus two with no cancellation or shared branch
support. This is an exact analysis of those finite lists, **not** a claim that
either bounded list exhausts all exceptional directions or proves a new rank
statement.

ICARM 262, the remaining sub-cutoff rank-20 fibres, and bisections of the
eventual rootless K3 are deliberately not given a negative conclusion here:
this checkout does not yet contain a pinned, independently certified list of
their square conditions. The reusable
[`branch_divisor_code.py`](../cas/branch_divisor_code.py) accepts arbitrary
exact numerator/denominator polynomials once those records are available. In
particular, `analyze_exceptional_branch_code.py --conditions-file FILE` reads
one or more JSON files with

```json
{"families":[{"name":"rootless_bisections","parameter":"u","conditions":[
  {"label":"B1","numerator_coefficients_ascending":[...],
   "denominator_coefficients_ascending":[...]}
]}]}
```

and appends their complete branch-code reports to a fresh artifact. This
interface is exact; no numerical roots or point searches are used.
