# Strong Factorial after balanced private amplification

Internal research proof, 21 November 2026. Repository base:
`19eb9912fb60142f8b620dac24896931d4ff3e46`.

## Result and scope

Let `M_1,...,M_s` be monomials in any shared set of variables. They may
overlap arbitrarily and need not be distinct. Assume their exponent columns
have the same multiset of positive entries. Introduce fresh independent
variables `X_1,...,X_s`. For every integer `D>=6s` and every characteristic-zero
field, every nonzero member of the coefficient family

\[
 F=\sum_{j=1}^s c_j X_j^D M_j
\tag{1}
\]

satisfies Strong Factorial at every starting order, with the actual number
of nonzero coefficients as its window length. The coefficients are arbitrary;
they are not restricted to real, positive or root-of-unity values.

Every finite nonnegative integer support has an explicit enlargement satisfying these
hypotheses: privately complete each column's exponent profile to the union
of all the original column profiles, then add the fresh powers `X_j^D`.
The exponent `6s` depends only on the support size, not on the original
degrees, overlaps or moment order. No assertion that the decorations can be
removed is made. **Unrestricted SFC and general row-merging preservation
remain open.** No novelty or external-refereeing claim is made.

The one external analytic input is a classical quantitative power-sum
inequality. All factorial estimates below are exact inequalities between
positive rational structure constants, established before complex
coefficients are substituted. They do not assume positivity of a
characteristic-zero coefficient field or of a complex factorial moment.

## 1. The precise power-sum input

For complex numbers `z_1,...,z_s` with `max |z_j|=1`, and every `n>=1`,
Turán's second main theorem gives

\[
 \max_{n\le m<n+s}\left|\sum_{j=1}^s z_j^m\right|
 \ge T_{s,n},\qquad
 T_{s,n}=2\left(\frac{s}{8e(n+s-1)}\right)^s.
\tag{2}
\]

The inspected source is Hugh L. Montgomery and Robert C. Vaughan,
*Multiplicative Number Theory III*, Appendix K, Theorem K.3, printed
pages 334–337, in the
[author-hosted PDF](https://personal.science.psu.edu/rcv4/571s25/Vol3.pdf).
Its general weighted statement orders the roots by modulus and uses the
least absolute value of the partial sums of the weights. Here every
weight is `1`, that minimum is `1`, and its parameters are `M=n-1`, `N=s`.
The theorem explicitly permits `M=0`; repeated or zero roots are allowed.
No hypothesis about zeros of the Riemann zeta function is used.

The complete source bytes and the interpolation/Chebyshev proof of K.3
were inspected. Source PDF SHA-256:
`ef098f0a0e5e27fd16721e44fb2197b163357de45b16111e7b82d179f36e7855`.
The local audit copy is
`research/artifacts/local/stepanov_merging/root_montgomery_vaughan_turan.pdf`.
The [source ledger](sfc_balanced_private_amplification_source_ledger.json)
records the exact input and scope. This is a source pin, not a formal
verification of Turán's theorem.

## 2. Exact suppression of the mixed coefficients

Write `M_j=prod_i Y_i^(a_ij)`. The balanced exponent profiles give a common
pure moment

\[
 H_m=\prod_i (m a_{ij})!\quad\text{independent of }j.
\]

For a composition `k=(k_1,...,k_s)` of `m`, put

\[
 C_k=\binom{m}{k_1,\ldots,k_s},\qquad
 H(k)=\prod_i\left(\sum_j a_{ij}k_j\right)!.
\]

**Lemma 1.** For `m>=1`, `H(k)<=H_m`; the case `m=0` is trivially equal to one.

**Proof.** Fix a row `i`. Take a list of `m` integers containing `k_j`
copies of `m a_ij`. Their average is the integer `sum_j a_ij k_j`.
Among nonnegative integer lists of fixed length and sum, the product of
factorials is minimized when the entries differ by at most one: replacing
`u,v` with `u+1,v-1` when `v>=u+2` multiplies that product by
`(u+1)/v<1`. Here the average is integral, so

\[
 \left(\sum_j a_{ij}k_j\right)!^m
 \le \prod_j (m a_{ij})!^{k_j}.
\]

Multiply over the rows and use the common pure moments to obtain
`H(k)^m<=H_m^m`. All quantities are positive integers. QED.

**Lemma 2.** For every positive integer `D`,

\[
 \frac{\prod_j(Dk_j)!}{(Dm)!}\le C_k^{-D}.
\tag{3}
\]

**Proof.** Partition `Dm` labelled positions into `D` fixed labelled blocks
of size `m`. There are `C_k^D` colorings with `k_j` positions of color `j`
inside each block. These inject into all colorings with `Dk_j` positions
of each color, whose number is `(Dm)!/prod_j(Dk_j)!`. QED.

For `m>=1`, normalize the actual moment of (1) by

\[
 P_m(c)=\frac{\Lambda(F^m)}{H_m(Dm)!},\qquad
 p_m(c)=\sum_jc_j^m.
\]

Every pure coefficient of `P_m` is exactly one. Lemmas 1 and 2 give the
coefficient bound

\[
 0< [c^k]P_m
 =C_k\,\frac{H(k)}{H_m}\,
       \frac{\prod_j(Dk_j)!}{(Dm)!}
 \le C_k^{1-D}.
\tag{4}
\]

For a nonpure composition and `m>=2`, `C_k>=m`: choose one nonempty color
of size `a` with `1<=a<=m-1` and use `C_k>=binom(m,a)>=m`.
Consequently, on `max_j |c_j|=1`,

\[
 |P_m(c)-p_m(c)|\le B_{s,D}(m),\qquad
 B_{s,D}(m)=\binom{m+s-1}{s-1}m^{1-D}\quad(m\ge2).
\tag{5}
\]

The composition count intentionally includes the pure terms, so (5) is
a convenient upper bound. At `m=1` the error is exactly zero.

## 3. The uniform exponent `D=6s`

Assume `s>=2` and `D>=6s`.

First, `B_(s,D)(m)` decreases for real `m>=2`. Indeed it is

\[
 \frac{m^{s-D}}{(s-1)!}
       \prod_{j=1}^{s-1}(1+j/m),
\]

a product of positive nonincreasing factors when `D>=s`.
Likewise `B_(s,D)(n)(n+s-1)^s` is decreasing for `n>=2`, because it is

\[
 \frac{n^{2s-D}}{(s-1)!}
   \prod_{j=1}^{s-1}(1+j/n)(1+(s-1)/n)^s.
\tag{6}
\]

For `n>=2`, divide the maximum error in the window by (2), use (6),
and then use `e<3`:

\[
 \frac{\max_{n\le m<n+s}B_{s,D}(m)}{T_{s,n}}
 \le \frac{s(s+1)}2
       \left(\frac{8e(s+1)}s\right)^s 2^{-D}
 < R_s,
\]

\[
 R_s=\frac{s(s+1)}2
       \left(\frac{3(s+1)}{8s}\right)^s.
\tag{7}
\]

The comparison is rational after the elementary bound `e<3`:

\[
 R_2=\frac{243}{256}<1,
\]

\[
 \frac{R_{s+1}}{R_s}
 =\frac{3(s+2)^2}{8s(s+1)}
   \left(\frac{s(s+2)}{(s+1)^2}\right)^s<1\quad(s\ge2).
\tag{8}
\]

For the first factor in (8), the required inequality is exactly
`8s(s+1)-3(s+2)^2=(s-2)(5s+6)>=0`; the second factor is strictly less
than one. Thus the maximum error is less than `(243/256)T_(s,n)`.

For `n=1`, `P_1=p_1`. All remaining errors in the first window are at
most `B_(s,D)(2)=s(s+1)2^(-D)`. Hence

\[
 \frac{\max_{1\le m\le s}|P_m-p_m|}{T_{s,1}}
 < Q_s=\frac{s(s+1)}2(3/8)^s\le Q_2=27/64.
\tag{9}
\]

Here `Q_(s+1)/Q_s=3(s+2)/(8s)<=3/4` for `s>=2`.

Choose the index supplied by (2). Equations (5)–(9) make its actual
normalized moment nonzero. In fact, for every `max |c_j|=1`,

\[
 \max_{n\le m<n+s}|P_m(c)|
 >\frac{13}{256}T_{s,n}
 >\frac{13}{128}
       \left(\frac{s}{24(n+s-1)}\right)^s.
\tag{10}
\]

The one-term case is immediate. If some coefficients vanish, delete them
and apply the proof with the active count `r`; the remaining columns are
still balanced and `D>=6s>=6r`. This proves the actual support-size window.

For an arbitrary characteristic-zero field, a hypothetical zero block
uses only finitely many coefficients. Their finitely generated field over
`Q` embeds into `C`. The embedding preserves nonzero coefficients and all
factorial-moment equalities. Normalize the embedded coefficients by their
maximum absolute value and apply the complex proof. This descent does not
assert that an arbitrary entire coefficient field embeds into `C`.

## 4. Every support has an explicit enlargement

Let `A=(a_ij)` be any finite nonnegative integer exponent matrix, with
`s` columns and arbitrary overlaps. For each column `j`, introduce fresh
variables `Z_(j,h,i)` for every `h!=j` and every positive entry `a_ih`, and
set

\[
 \widetilde M_j(Y,Z)=Y^{a_j}
       \prod_{h\ne j}\prod_{i:a_{ih}>0}Z_{j,h,i}^{a_{ih}}.
\tag{11}
\]

The exponent profile of every `Mtilde_j` is the multiset union of the
profiles of all columns of `A`. Each `Z` variable appears in just one
column. The theorem therefore applies to

\[
 \widetilde F=\sum_j c_j X_j^{6s}\widetilde M_j.
\tag{12}
\]

This constructs an all-order SFC family from every support without any
bound on its initial degree or dimension. It is achieved by actual
nonnegative integer exponent rows. A private power can also be built
from unit rows and mergers involving only that column. The fresh markers
make all final monomials distinct, even if the input columns coincided.

More generally, any variable identifications among the base monomials
which leave their exponent profiles equal retain the theorem, provided
the fresh marker variables stay private. An arbitrary base identification
can be followed by another explicit private balancing step. This is a
proved class of compensated identifications. It is not preservation for
an unmodified arbitrary row merger.

## 5. Positive polynomial templates with equal pure moments

The private-amplification argument also applies to polynomial blocks.
Let `P_1,...,P_s in Q[Y_1,...,Y_d]` be nonzero polynomials that are
nonnegative on the positive real orthant. Assume that for every integer
`m>=1` their ordinary factorial moments have one common value

\[
                  \Lambda(P_j^m)=H_m>0\qquad(1\le j\le s).
 \tag{P1}
\]

Introduce fresh variables `X_1,...,X_s` that occur in none of the
templates. For every integer `D>=6s`, every characteristic-zero field
`K`, and every nonzero coefficient vector `(c_1,...,c_s) in K^s`, put

\[
                       F=\sum_{j=1}^s c_jX_j^D P_j(Y).
 \tag{P2}
\]

**Corollary.** If `r` coefficients are nonzero, every block of `r`
consecutive positive moment orders contains a nonzero moment of `F`.
In particular `F` satisfies ordinary SFC. Its actual support size is

\[
       \#\operatorname{supp}(F)
           =\sum_{j:c_j\ne0}\#\operatorname{supp}(P_j)\ge r.
 \tag{P3}
\]

The templates may overlap arbitrarily, may have nonzero constant
terms, and need not be distinct. Their individual monomial coefficients
may have either sign; pointwise nonnegativity is the hypothesis. The
final coefficients `c_j` are unrestricted elements of `K`.

**Proof.** The ordinary factorial functional on rational polynomials
has the integral representation

\[
 \Lambda(H)=\int_{(0,\infty)^d}H(y)
                              e^{-y_1-\cdots-y_d}\,dy.
 \tag{P4}
\]

This follows termwise from the elementary factorial integral; every
polynomial is absolutely integrable against this density. Fix a
composition `k_1+...+k_s=m>=1`. Weighted AM--GM, applied pointwise to
the nonnegative numbers `P_j(y)^m` with weights `k_j/m`, gives

\[
 0\le\prod_jP_j(y)^{k_j}
       \le\sum_{j:k_j>0}\frac{k_j}{m}P_j(y)^m.
\]

Indices with weight zero are omitted. Integrating and using (P1) yields

\[
                 0\le H(k):=\Lambda\!\left(\prod_jP_j^{k_j}\right)
                       \le H_m.
 \tag{P5}
\]

After expanding (P2) by its `s` blocks, the coefficient of `c^k` in
its normalized moment is exactly

\[
 \binom{m}{k_1,\ldots,k_s}
       \frac{H(k)}{H_m}
       \frac{\prod_j(Dk_j)!}{(Dm)!}.
 \tag{P6}
\]

Its pure coefficients are one. Inequality (P5) replaces the balanced
monomial-profile lemma, and the private-marker coloring inequality
bounds (P6) by `binom(m;k_1,...,k_s)^(1-D)`. The same mixed-mass
estimate, Turan power-sum bound, and `D>=6s` comparison therefore
prove the corollary over `C`. Removing zero coefficients replaces `s`
by `r` while keeping `D>=6r`.

Distinct private markers distinguish all monomials belonging to different
active blocks. Within each block multiplication by `c_jX_j^D` preserves
the nonzero coefficients of `P_j`. This proves (P3), even when a
template has negative coefficients. The block-count window is no
greater than the actual support-size window.

For general `K`, the finitely generated field over `Q` generated by
the finitely many `c_j` embeds into `C`. The rational templates and
all factorial-moment identities are preserved by that embedding.
The complex conclusion therefore excludes any hypothetical zero block
over `K`. Positivity is used only for the fixed rational templates in
(P4)--(P5), before arbitrary field coefficients are substituted. QED.

A simple sufficient condition for (P1) is to choose one nonzero rational
polynomial `P` that is nonnegative on the positive orthant and let each
`P_j` be obtained from `P` by a permutation of the common coordinate
variables. The product exponential measure in (P4) is permutation
invariant, so all pure moments agree. They are strictly positive:
a nonzero nonnegative polynomial is positive on a nonempty open subset
of the orthant. Polynomials with nonnegative rational coefficients are
an immediate concrete subfamily, as are nonzero sums of squares.

Equal exponent profiles for monomials are the monomial instance of
this coordinate-permutation construction. The extension enlarges the
base blocks, but it keeps the fresh private powers and does not supply
a way to remove them or to infer SFC for an undecorated polynomial.

## 6. Why this does not yet imply unrestricted SFC

The decoration map changes mixed coefficients, not just each moment by a
nonzero scalar. Removing the private factors would require a separate
regularity-preservation statement. No such statement follows from (4),
because decreasing `D` increases the permitted mixed mass and can exceed
the quantitative power-sum lower bound. A continuous deformation or
nonzero endpoint resultant is insufficient to bridge that gap.

Balancing is also essential to this particular comparison. Without it,
the pure coefficients cannot in general all be normalized to one with a
single scalar for each degree and a fixed rescaling of the variables.
For example, `M_1=1`, `M_2=Y` have pure moments `1` and `m!`; their ratio
is not a fixed scalar to the power `m`. Thus (2) cannot be applied by
silently using degree-dependent coefficient variables.

Finally, the fresh markers must remain private. Merging all `X_j` into
one variable makes `X^(Dm)` a common nonzero scalar factor in the moments,
so it removes exactly the multinomial suppression in (3). When the base
monomials are distinct, this merger keeps the columns distinct, so the
missing implication is a substantive instance of the original merging
problem. If base monomials coincide, the merger can also cause column
collisions; those are excluded from the general preservation obligation.

## 7. Computational assurance

The [exact checker](verify_sfc_balanced_private_amplification.py) checks
the coefficient identities, integer factorial inequalities, rational
constant comparisons, balancing construction and failure controls. It
also constructs selected finite projective-regularity certificates by
exhibiting full-rank homogeneous Macaulay matrices modulo a prime. A
nonzero full minor modulo that prime certifies a nonzero rational minor;
this certifies those specified finite instances in characteristic zero.
The [receipt](sfc_balanced_private_amplification_verification.json) records
the precise domains and witnesses.

The universal conclusion comes from the written proof and the stated
classical theorem, not from the finite certificates. Root, the arithmetic reviewer and the factorial reviewer independently
reconstructed all inequalities, the
source hypotheses, the initial window and coefficient-boundary cases,
the balancing construction and characteristic-zero descent. The positive
polynomial-template extension was separately reconstructed by root and
both reviewers.

Replay from the repository root:

```sh
python research/factorial/image-rigidity/verify_sfc_balanced_private_amplification.py
```

The source proof is the inspected classical theorem, not an RH assumption
or a positive-measure regularity-preservation claim. No remote publication
or change to the global status registry is made by these files.
