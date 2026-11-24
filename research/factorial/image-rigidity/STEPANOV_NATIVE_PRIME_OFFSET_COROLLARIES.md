# An extra p-adic digit and explicit prime-gap corollaries

Internal derivation, 21 November 2026, from base revision
`19eb9912fb60142f8b620dac24896931d4ff3e46`.

## Result and exact status

**Continuation update.** The residual classes in this note describe the
output of its extra-digit criterion. Later proofs have closed both lower-three
prime gaps five and seven, and the lower-four exceptional prime gap eleven,
at every order. See [the current synthesis](STEPANOV_SFC_PROGRESS_AND_LOWER_THREE_CLOSURE.md),
[the lower-three prime closures](STEPANOV_LOWER_THREE_PRIME_CLOSURE.md), and
[the lower-four gap-eleven closure](STEPANOV_LOWER_FOUR_ELEVENTH_GAP_CLOSURE.md).
The original congruence calculations remain valid intermediate reductions;
those three families are no longer open.

The proposed integral-correction argument is valid and already covered
by a stronger retained pointwise theorem. A further ordinary refinement
proved here improves its residual congruences from modulo `p^2` to modulo
`p^3`: the lower reciprocal root has valuation at least one, and the
primitive jet forces another positive valuation when the Hasse scalar
vanishes. This additional digit uses no positivity argument.

The note also extracts an explicit harmonic-number corollary of the
retained theorem. It gives, for example, the ordinary all-terminal family

\[
 X^6(A+BX^p),\qquad p>6\text{ prime},\quad p\notin\{13,137\}.
\]

The retained relative-gap theorem removes `p=137`, leaving only `p=13`
among primes greater than six. The gap-five case is separately retained
as the prime-plus-one family. No assertion at gaps two or three is made
here. For lower exponent six, the newly extracted content relative to
the relative-gap theorem is at prime gaps seven, eleven and seventeen.

For fixed lower exponent, the bare statement “all but finitely many prime
gaps” already follows from `d>=3a`. The additional content here is an
explicit algebraic exclusion rule that also applies within `a<p<3a`.
No literature-priority claim is made.

There are two distinct dependencies concerning positivity:

* The integral-correction/Krylov argument in Section 2 is algebraic and
  uses no positivity principle. The additional ordinary digit in
  Theorem 3a is also algebraic.
* The stronger harmonic corollary uses the retained native-prime theorem.
  Its nonintegral-correction branch proves uniqueness of a possible
  common root and excludes a rational real root using strict positivity
  of an actual even-order real Gamma integral. That positivity is proved
  for these actual moments, not assumed for a determinant or for arbitrary
  coefficient fields. Rational Bezout identities then transfer the
  ordinary conclusions to every characteristic-zero coefficient field.

## 1. Exact imported statements and notation

Let

\[
 T_{n,\gamma,d}(t)=\sum_{k=0}^n\binom nk(\gamma)_{dk}t^k.
\]

The actual ordinary pair at lower exponent `a`, gap `p`, and terminal
`N=n+1>=2` is

\[
 T_{n,an+1,p},\qquad T_{n+1,a(n+1)+1,p}.
 \tag{1}
\]

The following sources were read, including their relevant proof steps:

* [DIAGONAL_PRIME_GAP_SFC.tex](DIAGONAL_PRIME_GAP_SFC.tex),
  Theorem `thm:offset-rank` and its proof: a necessary and sufficient mixed
  rank criterion for `1<=h<=d`, with integral polynomial correction
  coefficients and a unit triangular change of columns.
* [LOWER_TWO_PRIME_GAP_SFC.tex](LOWER_TWO_PRIME_GAP_SFC.tex),
  proof of `thm:mixed`, especially concomitant equation (3): the algebraic
  rising-jet columns and the full-rank restriction argument.
* [HIGHER_GAP_ADJOINT_RANK.tex](HIGHER_GAP_ADJOINT_RANK.tex),
  equations (2)–(8) and proof of `thm:rank`: the lower differential equation
  is monic of order `d+1` at every nonzero point; independent adjoint
  solutions give a full-rank concomitant with kernel the actual nonzero
  lower jet. No simplicity assumption on a lower root is made.
* [PRIME_PLUS_ONE_EXPONENT_GAP_SFC.tex](PRIME_PLUS_ONE_EXPONENT_GAP_SFC.tex),
  `lem:kernel` and `thm:criterion`: for prime `p`, `1<=h<=p+1`, `n>=1`,
  and positive real algebraic `gamma` with a `p`-integral embedding, set

  \[
  C=\gamma+pn+p,\qquad
  \kappa_h=\frac{p(n+1)}{(C)_h},\qquad
  \sigma_h(x)=[Y^2](x+Y)_h.
  \]

  The lower and upper moments are coprime if
  `v_p(kappa_h)>=0` **or** `sigma_h(gamma mod p)!=0`.
  The second condition is a Hasse derivative, with no division by two
  in the residue field. The proof uses one fixed row set of the reduced
  concomitant, independent of the putative common root.

For the ordinary problem, `gamma=an+1` and `h=a`, so the embedding and
positivity hypotheses on the shape hold exactly. A missing phase
coefficient leaves a monomial, whose positive factorial moments are
nonzero. Otherwise the exact normalization is

\[
 \mathcal L_1((X^a(A+BX^p))^m)
       =A^m(am)!T_{m,am+1,p}(B/A).
 \tag{2}
\]

## 2. The proposed algebraic sieve is valid at every terminal

**Theorem 1.** Let `p` be prime, `1<=a<p`, and `N>=2`. The pair (1) is
coprime whenever

\[
 (a+p)N\not\equiv r\pmod {p^2}
           \quad\text{for every }1\le r\le a-1.
 \tag{3}
\]

The excluded congruence classes are potential exceptions to this
criterion, not demonstrated common-root cases. There are exactly `a-1`
distinct excluded classes modulo `p^2`; when `a=1` there are none.

**Proof.** Work in an algebraically closed complete extension of
`Q_p`. Put `n=N-1`, `gamma=an+1`, and

\[
 C=(a+p)n+p+1=(a+p)N-a+1.
\]

Since `C+a-1=(a+p)(n+1)`, exact cancellation gives

\[
 \kappa_a=\frac{p(n+1)}{(C)_a}
   =\frac{p}{(a+p)\prod_{j=0}^{a-2}(C+j)}.
 \tag{4}
\]

The empty product is one. The factor `a+p` is a `p`-adic unit. The
remaining factors occupy an interval of length `a-1<p`, so at most one
is divisible by `p`. Therefore

\[
 v_p(\kappa_a)\ge0
 \quad\Longleftrightarrow\quad
 p^2\nmid C+j\quad(0\le j\le a-2).
 \tag{5}
\]

This is exactly (3), with `r=a-1-j`. Distinct `r` give distinct classes
because `a+p` is invertible modulo `p^2`.

For completeness, the integral branch can be proved directly from the
mixed-rank criterion. A lower root `t` gives `u=epsilon/t`, with
`epsilon=(-1)^(p+1)`. Every coefficient at degree `k` of the lower `T`
has valuation at least `k`: `(an+1)_(pk)` contains at least `k` multiples
of `p`, and the binomial coefficient is integral. The lower monic
reciprocal is

\[
 u^nT_{n,an+1,p}(\epsilon/u).
\]

If `v_p(u)<1`, its leading term has strictly smaller valuation than
every other term, which is impossible at a root. Thus `v_p(u)>=1`.

In the offset criterion, every correction to a column is
`epsilon*kappa_a*u*v_k*B_n`, where `v_k` is an integer polynomial in
`C`. Under (5), each correction reduces to zero. The reduced cyclic
space is `k[z]/(z^p)`. On this space put

\[
 A=z-zD_z,\qquad U=A+an.
\]

The auxiliary recurrence reduces to `B_m=U^m1`, always at the same
lower shape `gamma`. The matrix of `U` in `1,z,...,z^(p-1)` has
diagonal entries `an-j` and successive subdiagonal entries one.
Its eigenvalues are the distinct elements of `F_p`, and `U^p=U`.
The seed `1` is cyclic: its first `p` iterates are monic in distinct
successive degrees. Hence its projection onto every eigenspace is
nonzero.

On the `p-1` nonzero eigenspaces, the vectors

\[
 U^{n+1}1,\ldots,U^{n+p-1}1
\]

give a Vandermonde matrix multiplied by nonzero diagonal entries.
They have rank `p-1` for every `n>=1`. The retained triangularization
from these columns to the `L_j` has integral entries and diagonal
`1,-1,...,(-1)^(p-2)`, all units. The corrected matrix therefore has
rank `p-1` modulo the maximal ideal and in characteristic zero.
This contradicts the necessary mixed-rank deficiency of a common root.

Fourier inversion is used only in characteristic zero to derive the
cyclic criterion; no Fourier matrix is reduced modulo `p`. All
denominators removed in the rank criterion are nonzero at the actual
positive integer shape. The argument requires no extra seed condition
and no lower-root simplicity hypothesis. Coprimality over `Q` supplies
a rational Bezout identity and thus transfers to every characteristic-zero
coefficient field. This proves the theorem.

## 3. The stronger retained theorem gives an exact harmonic sieve

Put `H_0=0` and `H_m=sum_(i=1)^m 1/i`. Define the explicit integers

\[
 s_{a,j}=\sigma_a(-j)
   =(-1)^j j!(a-1-j)!
                 \bigl(H_{a-1-j}-H_j\bigr),
       \quad 0\le j\le a-1.
 \tag{6}
\]

**Lemma 2.** Formula (6) holds as a rational identity and its value
is an integer. Moreover `s_(a,j)=0` if and only if `a` is odd and
`j=(a-1)/2`.

**Proof.** Factor `(x)_a=(x+j)G_j(x)`, where
`G_j(x)=product_(k!=j)(x+k)`. The coefficient of `Y^2` in
`(-j+Y)_a` is `G_j'(-j)`, an integer. Its product-rule formula is

\[
 G_j'(-j)=G_j(-j)\sum_{k\ne j}\frac1{k-j},
\]

which is exactly (6). The factorial factor is nonzero, and harmonic
numbers are strictly increasing, so the only zero is the stated
central index.

**Theorem 3 (retained baseline residual classes).** For `1<=a<p`, any failure
of ordinary consecutive coprimality must lie in a terminal class

\[
 (a+p)N\equiv a-1-j\pmod {p^2},
 \qquad 0\le j\le a-2,\qquad p\mid s_{a,j}.
 \tag{7}
\]

Here every prime divides the integer zero. This is a necessary
condition for an actual exception, not an assertion that an exception
exists in any listed class.
This `p^2` condition records the retained baseline; the current stronger
ordinary condition is (11), modulo `p^3`.

**Proof.** By Theorem 1 an unexcluded terminal has a unique `j` with
`p^2 | C+j`. Then `gamma mod p=C mod p=-j`. If `p` does not divide
`s_(a,j)`, the Hasse-derivative branch of the retained pointwise theorem
excludes the common root. This proves (7).

The imported Hasse branch is substantive. Briefly, its primitive jet
has reduction proportional to the derivative vector
`k_m=((x)_m)'|_(x=gamma mod p)`, with `k_p=-1`. The upper relation
evaluates on this kernel to `sigma_a(gamma mod p)`. If that scalar is
nonzero, one fixed determinant has affine reduction with nonzero slope,
so it has at most one integral zero. All conjugates of a common root
obey the same determinant, forcing a common factor over `Q` to be
linear. An actual even-index moment is strictly positive for real
phase ratio, by its Gamma integral, excluding that rational root.
This established property of the real moment sequence is the
positivity input inherited here. It is not a proposed determinant
positivity assumption or an ordering imposed on an arbitrary field.

### One additional digit from the actual lower-root valuation

**Lemma 3a (rational-shape integral boundary).** Let `p` be prime,
`1<=h<=p+1`, `n>=1`, and let `gamma>0` be rational and `p`-integral.
With the notation of the retained pointwise criterion, suppose

\[
 \sigma_h(\bar\gamma)=0,\qquad v_p(\kappa_h)\ge-1.
 \tag{9}
\]

Then `T_(n;gamma,p)` and `T_(n+1;gamma+h,p)` are coprime. This
implication uses no positivity argument.

**Proof.** In every block of `p` factors `gamma+j`, exactly one is
zero modulo `p`, and hence belongs to `p Z_p`. Consequently the
coefficient at degree `k` of the lower moment has valuation at least
`k`, exactly as in Section 2. Every lower reciprocal root satisfies
`v_p(u)>=1`. This argument applies to a rational `p`-integral shape,
not only an integer shape.

Suppose there were a common root. Divide its actual nonzero rising-jet
vector by a component of minimum valuation to obtain a primitive
integral vector `w=(w_1,...,w_p)`. This scaling preserves all its
homogeneous equations. The retained reduced-kernel lemma gives

\[
 \bar w=c k,\qquad c\ne0,\qquad
 k_m=((x)_m)'|_{x=\bar\gamma},\quad k_p=-1.
\]

In particular `w_p` is a unit. Let `v=(v_1,...,v_p)` be the integral
rising-basis coefficient row of `V_h`; its constant coefficient is
omitted because `R=0`. The exact upper equation and kernel evaluation
are

\[
 w_p=\epsilon\kappa_h u\,(v w),\qquad
 \overline{v w}=c\,\sigma_h(\bar\gamma)=0.
\]

Thus `v w` belongs to the maximal ideal. By (9) and `v_p(u)>=1`,
`epsilon*kappa_h*u` is integral. The right side belongs to the maximal
ideal while the left side is a unit, a contradiction. Fractional
valuations of the individual normalized jet entries cause no problem:
only strict positivity of the valuation of `v w` is needed.

Combining this lemma with the retained theorem gives the strengthened
pointwise criterion, for the stated rational shapes,

\[
 \sigma_h(\bar\gamma)\ne0
       \quad\text{or}\quad v_p(\kappa_h)\ge-1.
 \tag{10}
\]

The first branch still uses the retained theorem's proof. The new
algebraic implication is specifically the zero-Hasse branch (9).
Mere integrality of a shape in an arbitrary ramified extension of
`Q_p` does not supply the lower-root bound `v_p(u)>=1`; this note
does not assert (10) under that weaker hypothesis.

There is an exact scope countercontrol. Take `p=3`, `n=1`, and the
positive real shape `gamma=1+sqrt(3)`, embedded at a three-adic square
root of three. This shape is integral, but

\[
 (\gamma)_3=(1+\sqrt3)(2+\sqrt3)(3+\sqrt3)
             =24+14\sqrt3.
\]

The two summands have distinct three-adic valuations `1` and `1/2`.
The lower monic reciprocal is `u+(gamma)_3`, so its root has valuation
`1/2`. This disproves the stronger root bound for merely integral
ramified shapes; it is not a counterexample to Strong Factorial.

**Theorem 3a (ordinary sieve modulo the prime cube).** For `1<=a<p`, an actual
ordinary common root at terminal `N` can exist only if

\[
 (a+p)N\equiv a-1-j\pmod {p^3},\qquad
 0\le j\le a-2,\qquad p\mid s_{a,j}.
 \tag{11}
\]

**Proof.** Theorem 3 supplies `sigma_a(gamma mod p)=0` and a unique
index `j` with `p^2 | C+j`. Lemma 3a rules out `v_p(kappa_a)>=-1`.
The rational number `kappa_a` has integer valuation, so its valuation
must be at most `-2`. Formula (4), with at most one `p`-divisible
denominator factor, makes this exactly `p^3 | C+j`. Substituting
`C=(a+p)N-a+1` gives (11).

There is no iteration to an arbitrary number of digits in this proof.
At valuation `v_p(kappa_a)=-2`, the product `kappa_a*u` need not be
integral, so the same contradiction no longer follows. Further
refinement requires a new quantitative jet-valuation input.

## 4. All-terminal even offsets and the surviving odd-offset progression

For `a>=2`, define the nonzero integer

\[
 D_a=\prod_{\substack{0\le j\le a-1\\s_{a,j}\ne0}}|s_{a,j}|.
\]

**Corollary 4.** Let `a>=2` be even. For every prime `p>a` with
`p` not dividing `D_a`, every nonzero ordinary binomial
`X^a(A+BX^p)` satisfies Strong Factorial at every terminal `N>=2`,
over every characteristic-zero coefficient field.

**Proof.** Lemma 2 makes every `s_(a,j)` nonzero. The stipulated
prime divides none of them, so (7) is empty. Apply (2) and rational
Bezout transfer.

**Corollary 5.** Let `a>=3` be odd. For every prime `p>a` with
`p` not dividing `D_a`, any actual consecutive common root must occur at
a terminal satisfying

\[
 (a+p)N\equiv (a-1)/2\pmod {p^3}.
 \tag{8}
\]

**Proof.** By Lemma 2 only `j=(a-1)/2` has zero Hasse value. It lies
in `0,...,a-2`. Substitute it in (11).

Corollary 4 is an explicit consequence of the retained theorem.
Corollary 5 additionally uses the new ordinary digit refinement;
neither is an independent proof of the retained theorem. The finite
calculations below certify the displayed constants, not their
all-terminal scope.

## 5. Exact small-offset consequences and remaining classes

The full Hasse-value lists and their prime divisors above `a` are:

| `a` | `s_(a,0),...,s_(a,a-1)` | Prime divisors above `a` of nonzero values |
|---|---|---|
| 3 | `3,0,-3` | none |
| 4 | `11,-1,-1,11` | `11` |
| 5 | `50,-5,0,5,-50` | none |
| 6 | `274,-26,4,4,-26,274` | `13,137` |
| 7 | `1764,-154,28,0,-28,154,-1764` | `11` |
| 8 | `13068,-1044,188,-36,-36,188,-1044,13068` | `11,29,47` |
| 10 | `1026576,-69264,11016,-2664,576,576,-2664,11016,-69264,1026576` | `13,17,37,7129` |

The retained
[relative-gap-three theorem](RELATIVE_GAP_THREE_STRONG_FACTORIAL.md)
states all-terminal ordinary SFC for `a>=2,d>=3a`. Combining that
statement with Corollary 4 gives:

* Lower exponent six: all prime gaps `p>6` except possibly `13`.
  Gap `137` is covered by `137>=18`; gap `5` is separately covered
  by the retained `a=p+1` family. The new near-diagonal prime-gap
  consequences are `7,11,17`.
* Lower exponent eight: all prime gaps `p>8` except possibly `11`;
  `29,47` are covered by `p>=24`.
* Lower exponent ten: all prime gaps `p>10` except possibly `13,17`;
  `37,7129` are covered by `p>=30`.

The phrase “except possibly” describes this combination of criteria,
not proven counterexamples or a claim that every other repository
method fails on those supports.

The same exact sieve gives the following intermediate prime-gap domains:

| Lower exponent and gap | Retained criterion alone (baseline) | Additional ordinary digit |
|---|---|---|
| `a=4,p=11` | `N=97 (mod 121)` | `N=1065 (mod 1331)` |
| `a=6,p=13` | `N=18,89 (mod 169)` | `N=694,1272 (mod 2197)` |
| `a=8,p=11` | `N=115 (mod 121)` | `N=841 (mod 1331)` |
| `a=10,p=13` | `N=147,162 (mod 169)` | `N=669,2006 (mod 2197)` |
| `a=10,p=17` | `N=75,118 (mod 289)` | `N=364,1274 (mod 4913)` |

These are computed from (7) and (11), so they are necessary conditions
only. The first row has since been closed at every order by the linked
gap-eleven theorem; the other rows retain their stated intermediate scope.

For the requested lower-three examples, the unrefined and refined
terminal classes are exactly:

| Gap | Integral-correction baseline | Retained native-prime baseline | Additional ordinary digit |
|---|---|---|---|
| `p=5` | `N=19,22 (mod 25)` | `N=22 (mod 25)` | `N=47 (mod 125)` |
| `p=7` | `N=5,10 (mod 49)` | `N=5 (mod 49)` | `N=103 (mod 343)` |

The removed classes correspond to `j=0`, whose Hasse value is `3`.
The classes surviving this criterion have `j=1` and Hasse value zero.
Their first candidate terminals are respectively `47` and `103`, without
requiring a finite-prefix certificate. The continuation closes both entire
progressions. In particular
the former residues represented by `N=22` at gap five and `N=5` at gap
seven are excluded by the new algebraic lemma: each has
`v_p(kappa_a)=-1`. The checker includes these as negative controls
against incorrectly retaining every old `p^2` residue. It does not
execute the remaining infinite progressions.

## 6. Verification and novelty limits

Run:

```text
python -B research/factorial/image-rigidity/verify_stepanov_native_prime_corollaries.py --repo-root . --output /tmp/stepanov-prime.json
```

The checker uses exact integer/rational arithmetic and finite-field
Gaussian elimination. It verifies the symbolic cancellation in (4),
the Hasse/harmonic formulas for the recorded offsets, complete
factorizations of the displayed values, the exact residue-class lists,
finite valuation/congruence guards, and all cyclic shifts and positive
power classes for selected primes. It also compares the actual integer
auxiliary recurrence with the reduced Krylov model in finite cases.
It also verifies the primitive-kernel Hasse evaluation in finite
characteristics, the refined `p^3` residue tables, and exact boundary
examples with `v_p(kappa_a)=-1` that the new lemma excludes. The general
all-prime and all-terminal conclusions follow from the written arguments
and the declared retained inputs, not extrapolation of these guards.

Repository searches located the general native-prime pointwise theorem
and its explicit lower-four corollary. They did not locate the general
harmonic formula/corollary or the explicit lower-six and lower-eight
consequences stated here. Searches also did not locate the specific
zero-Hasse, valuation-at-least-minus-one improvement or this general
ordinary `p^3` sieve in the base revision. Other `p^3` congruences in
the repository concern different families and are not superseded.
This is an account of the inspected repository, not a claim of novelty
in the mathematical literature. No global status files have been changed.

The [retained internal replay](stepanov_native_prime_corollaries_verification.json) passed on
21 November 2026 with SymPy 1.14.0. It contains 3,654 complete
valuation/congruence guards, 336 finite-field Krylov rank checks,
52 comparisons with the actual integer auxiliary recurrence, and
273 lower-coefficient valuation guards, in addition to the symbolic
cancellation and exact Hasse-value/factorization certificates.
The completed extra-digit checks add 418 primitive-kernel/Hasse
evaluations, 110 explicitly eliminated valuation-minus-one lifts,
90 rational-shape coefficient guards, and the ramified-shape scope
countercontrol. The additional digit proof and all displayed refined
residue values were independently reconstructed by the parallel
determinant reviewer; this is internal AI review.
The all-order corollaries have the written proofs above; these finite
guards do not replace those proofs or the declared retained dependencies.
