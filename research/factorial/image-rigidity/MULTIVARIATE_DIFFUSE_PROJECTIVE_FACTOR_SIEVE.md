# A projective common-factor sieve at physical diffuse contacts

Consider ordinary two-monomial consecutive factorial moments
\(P_n,P_N\) with \(N=n+1\) and \(n\ge2\). If they have a nonconstant
common rational factor, every irreducible common factor has even
degree, because the even moment is strictly positive on the real
line. There is a scale-invariant polynomial that records the shape
of a common factor and has an exact denominator formula in terms of
its root valuations. This gives two restrictions near any fixed
core's positive diffuse contacts. Bounded valuation spread reduces
the possible common factors to finitely many projective shapes,
including at contacts of higher degree or multiplicity. At primes
that force a binomial residue on a single shared Newton slope,
the reduction of the projective polynomial must have roots of
bounded multiplicative order.

The local tests below allow different prime sets for different
profiles. The control is on the spread of the shared valuations,
or on the primes whose residues force torsion. The row heights,
counts and gaps can vary independently subject to the diffuse
hypotheses in the last section.

## Exact projective denominators

Let \(H\in\mathbb Q[z]\) be any monic common factor of even degree
\(d>0\), with roots \(\alpha_1,\ldots,\alpha_d\) counted with
multiplicity. Write \(c=\prod_i\alpha_i\) and set
\[
           B_H(Y)=\prod_{i=1}^d(Y-\alpha_i^d/c).
\]
The moment polynomials have nonzero constant terms, so \(c\ne0\).
In the factorial case all roots occur in nonreal conjugate pairs
and \(c>0\). The coefficients of \(B_H\) are rational by symmetry,
and its constant term is one. A common nonzero rescaling of the
roots leaves \(B_H\) unchanged. The projective multiset of roots
of \(H\) can be recovered from \(B_H\) up to finitely many choices:
choose a \(d\)-th root of each root of \(B_H\) to obtain the possible
roots of \(H\) divided by \(c^{1/d}\).

Let \(L_H\) be the least positive integer with \(L_H B_H\in\mathbb Z[Y]\).
For a prime \(p\), sort the valuations \(s_1\le\cdots\le s_d\) of the
roots of \(H\) in \(\overline{\mathbb Q}_p\), with \(v_p(p)=1\), and
put \(S=\sum_i s_i=v_p(c)\). Then the denominator is exact:
\[
 v_p(L_H)=
    \max_{0\le k\le d}\left(kS-d\sum_{i\le k}s_i\right)
    =\sum_{i=1}^d (S-ds_i)_+.                                   \tag{1}
\]
Indeed the valuations of \(\alpha_i^d/c\) are \(ds_i-S\) in the
same order. The coefficient at distance \(k\) from the leading
term has valuation at least the sum of the \(k\) least valuations.
When some are negative, take \(k\) to be their number. Their subset
is then the unique one of minimal valuation in the elementary
symmetric sum, giving equality. If none are negative, all
are zero and the coefficients are integral. For \(d=2\) this
reduces to the absolute difference of the two valuations, and
\[
    B_H(Y)=Y^2-\left((\alpha_1+\alpha_2)^2/c-2\right)Y+1 .
\]
Thus (1) extends the reduced-denominator calculation of the
[quadratic trace sieve](MULTIVARIATE_QUADRATIC_NEWTON_TRACE_SIEVE.md)
to arbitrary common factor degree and unequal root moduli.

For each prime \(p\), let \(W_p(d)\) be the maximum of the
right side of (1) over size-\(d\) multisets of shared root
valuations of \(P_n,P_N\). The multiplicity of a valuation
\(a/b\) in lowest terms must be a multiple of \(b\), and is
bounded by the smaller of its two Newton edge widths. The
[shared stepped-residue capacities](MULTIVARIATE_COMMON_FACTOR_FIELD_DISCRIMINANT.md#the-local-branch-and-discriminant-ledger)
may further restrict these multisets. The multiplicity condition
makes \(W_p(d)\) an integer. If the permitted list is empty at
any prime, a degree-\(d\) common factor is excluded already.
Otherwise \(L_H\mid \prod_p p^{W_p(d)}\). The product is finite,
since once the two polynomials are made integral, their
nonzero endpoint coefficients are units at all but finitely
many primes. A simple bound valid for each allowed multiset is
\[
 \max_k\left(kS-d\sum_{i\le k}s_i\right)
         \le \left\lfloor d^2/4\right\rfloor(s_d-s_1).
\]
In particular this tests spread independent of a large common
valuation change at all roots.

## A gcd bound from adjacent coefficient ratios

The denominator also has an explicit bound that needs neither
Newton face enumeration nor prime factorization. For the
degree-\(j\) coefficients \(a_{j,k}\) in the physical variable,
let
\[
 E_j=\operatorname{lcm}_{1\le k,l\le j}
   \operatorname{num}\left(\frac{a_{j,k}a_{j,l-1}}
                                {a_{j,k-1}a_{j,l}}\right),
\]
where \(\operatorname{num}\) is the positive reduced numerator
and the coefficient ratios are positive in this setting.
This integer is invariant under a rational change of argument
scale and a common scaling of the coefficients. It satisfies
\[
               L_H\mid \gcd(E_n,E_N)^{\lfloor d^2/4\rfloor}.
                                                                    \tag{2}
\]
At each prime the valuation of \(E_j\) is the maximum minus
the minimum of the first differences
\(v_p(a_{j,k})-v_p(a_{j,k-1})\). A Newton edge slope is the
average of these consecutive first differences over the
edge's index interval. Hence the Newton slopes of either
parent all lie in its interval of first differences. The
shared root valuations lie in the negative of the
intersection of the two intervals, whose diameter is at
most \(\min(v_p(E_n),v_p(E_N))\). Apply (1) and its spread
bound to obtain (2). The intersection diameter can replace
this minimum when first-difference valuations are known.
In particular every irrational projective branch in the
diffuse estimate below forces growth of the explicit
overlap integer \(\gcd(E_n,E_N)\).

## Projective residues from a uniform shared slope

Suppose the roots of \(H\) all have one valuation \(a/b\) at
a prime, in lowest terms. There \(b\mid d\), and the monic
stepped residue \(R\in\mathbb F_p[Z]\) of \(H\) is a divisor
of the two stepped residues on the shared slope. In particular
\(R\) has degree \(d/b\) and nonzero constant term. Then
\[
     \overline{B_H}(Y)=
       \prod_{R(w)=0}
          \left(Y-\frac{w^{d/b}}{R(0)}\right)^b .               \tag{3}
\]
The product counts roots of \(R\) with multiplicity. To check
this, take \(\pi^b=p\) and change variables by \(\alpha=\pi^a x\).
The monic polynomial of the rescaled roots reduces to
\(R(x^b)\). Factorization of initial forms follows from Gauss
normalization and multiplication in the residue field, which
also gives the required divisibility of the stepped residues.
Equivalently one can apply the
[Newton polygon product rule of Guardia, Montes and Nart](https://arxiv.org/pdf/0807.2620).
The product of the rescaled roots reduces to \(R(0)\), since
\(d\) is even. Taking \(d\)-th powers and dividing by this
product proves (3), including when \(p\mid b\).

If \(R\mid Z^g-A\) with \(A\in\mathbb F_p^\times\), every root of
\(\overline{B_H}\) has multiplicative order dividing \(2g\).
In fact for \(w\) a root of \(R\) one has \(w^g=A\) and
\(\left((-1)^{d/b}R(0)\right)^g=A^{d/b}\) by multiplying over
the roots of \(R\). Hence the \(g\)-th power of
\(w^{d/b}/R(0)\) is a sign. This applies, for instance, when
one of the two stepped face polynomials is a binomial of
degree \(g\). The containment
\[
              \overline{B_H}\mid (Y^{2g}-1)^d                  \tag{4}
\]
also holds for inseparable residual polynomials.

There is a quantitative consequence over primes that can
vary with or within the profile. Assume the coefficients of
\(B_H\) lie in a fixed bounded real box. Fix \(g\ge1\) and
write \(S_g(Y)=(Y^{2g}-1)^d\). For any set \(\mathcal T_g\) of
primes satisfying the uniform-slope binomial hypothesis with
this \(g\), either \(B_H\mid S_g\) in \(\mathbb Q[Y]\), or
\[
               \prod_{p\in\mathcal T_g}p
                        \le C\,L_H^{\,2gd-d+1},                \tag{5}
\]
with a computable \(C\) depending on \(g,d\) and the real box.
Indeed divide \(S_g\) by the monic \(B_H\). Its remainder has
coefficients bounded on the box, and clearing denominators
by \(L_H^{2gd-d+1}\) gives an integer polynomial by the
long-division steps. At the stated primes \(p\nmid L_H\)
by (1), and the entire remainder reduces to zero by (4).
If the remainder is nonzero, any one nonzero coefficient
after clearing denominators proves (5).

For a test against a possible degree-\(d\) factor, the set
\(\mathcal T_g\) can be selected without finding the factor.
At each included prime require every permitted slope multiset
to be uniform, and every permitted monic common stepped
residue of degree \(d/b\) there to divide some binomial of
degree \(g\) with nonzero constant term. This is a finite
local test. It can also require (3) directly for any
particular candidate shape. At primes where the two full
parent polynomials have unit endpoints, the usual normalized
reductions are the stepped residues with \(b=1,a=0\).

## Consequences at finite physical diffuse contacts

Use the notation \(\eta,\tau,\lambda,h,q,v\) of the
[physical diffuse comparison](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md).
Fix the ordinary core and a compact positive strength
interval, and let \(\eta\to0\) along actual collisions. At
the finitely many orders below the absolute cutoff, the
integer \(h\) has bounded range, and all limiting
\((q_*,v_*)\) belong to the finite algebraic reference list.

To retain a physical norm coordinate with rational
coefficients, use the lower tail first-factor scale
\(R_{n,1}\) of the
[common rational normalization](MULTIVARIATE_DIFFUSE_CLOCK_PRIME_HEIGHT.md#the-two-eliminants-in-rational-coordinates).
The polynomials in that common variable converge to
\[
 H_n(q_*,q_*^{-N(n-1)}z),\qquad
 H_N(q_*,\Lambda^h v_* q_*^{-N(n-1)}z),
\]
where \(H_j\) is the raw reference polynomial in the
physical diffuse comparison and \(\Lambda=N/n\). The
actual roots are the physical roots multiplied by
\(R_{n,1}>0\), which leaves \(B_H\) unchanged.

At each fixed reference contact, let \(\mathcal H_d\) be
the finite list of monic polynomials of degree \(d\)
obtained by selecting roots from the reference gcd with
multiplicity. It suffices to keep selections invariant
under complex conjugation. At nearby sufficiently diffuse
ordinary collisions, every monic degree-\(d\) common factor
in the common rational variable has coefficients, and
projective coefficients \(B_H\), within
\(O(\eta^\gamma)\) of those of one member \(H_*\) of
\(\mathcal H_d\). In general one may take
\(\gamma=1/[N(nN)^4]\). If the contacting roots are simple
in both polynomials, one may take \(\gamma=1/m_*\) using
the angular branch order \(m_*\) from the local bound in
the physical diffuse comparison.

To verify the claim, the coefficient and center estimates
and the algebraic confinement of \(q,v\) give an
\(O(\eta^{1/(nN)^4})\) coefficient error relative to the
fixed contact. The upper and lower leading coefficients
stay away from zero and the roots stay in a compact
annulus. The product of the distances from an actual root
to the roots of a reference parent is bounded by a
constant times the coefficient error, evaluated at the
actual root. Thus some distance is at most a constant
times the error to the reciprocal parent degree power.
For sufficiently small error, a common root must be
near the same root of both parents. Use disjoint disks
at the fixed reference roots. The root counts in those
disks agree with the reference multiplicities by Rouche's
theorem, so the common factor cannot exceed the smaller
multiplicity. In the simple case use the sharper local
parameter confinement and linear root perturbation on
the simple branch disks. Symmetric powers and division
by the product \(c\) are uniformly Lipschitz in the
matched roots on this compact annulus. These steps give
the assertion, including the norm coordinate \(c\) in
the common rational variable.

The following statements result.

1. If all projective coefficient candidates in
   \(\mathcal H_d\) have some irrational coordinate, then a
   degree-\(d\) common factor forces, on a branch whose such
   coefficient has algebraic degree \(e\),
   \[
       \sum_p W_p(d)\log p\ \ge\
             \log L_H\ \ge\ \frac{\gamma}{e}\log(1/\eta)-O(1).
                                                                    \tag{6}
   \]
   For a rational projective candidate, either \(B_H\) equals
   the candidate exactly, or \(\log L_H\ge
   \gamma\log(1/\eta)-O(1)\).
   To prove the irrational case, evaluate the fixed integral
   minimal polynomial of the candidate coefficient at the
   rational actual coefficient. The resulting nonzero rational
   number has absolute value at least \(L_H^{-e}\), and is
   \(O(\eta^\gamma)\) by the mean-value bound. In the rational
   case use separation of two distinct rationals. All constants
   can be bounded using the finite reference data. In particular
   a uniform bound on \(L_H\), allowing its prime divisors to
   vary, leaves only exact rational projective candidates.

2. The same argument applies to the physical norm \(c\) in the
   common rational normalization. Its denominator is exactly
   \(\prod_p p^{\max(0,-S)}\) with \(S\) formed from the shared
   valuations in that variable. For a given choice of slopes
   in the original physical variable, these are obtained by
   adding \(v_p(R_{n,1})\). If both the norm denominator and
   the projective denominator are bounded, then for sufficiently
   small \(\eta\) each common factor in that variable is itself
   one of the rational reference factors. Indeed both coordinates
   \((B_H,c)\) must stabilize to their reference values, and
   for fixed \((B_H,c)\) there are only finitely many root
   multisets obtained by taking \(d\)-th roots. Root matching
   in the fixed annulus identifies the reference member itself.
   In particular if none of the selected reference factors is
   rational, the two bounded-denominator conditions exclude
   that branch.

3. At a branch where \(B_{H_*}\nmid S_g\), this nondivisibility
   persists for sufficiently small \(\eta\), by continuity of
   the monic remainder coefficients. Thus its binomial primes
   obey (5) uniformly. In particular, if the prime product in
   (5) tends to infinity and the projective denominator stays
   bounded, the limiting projective roots must all be roots
   of unity of order dividing \(2g\). Here \(g\le N\) when it
   is taken from one parent binomial face. Even with unbounded
   projective denominators, (5) forces their growth relative to
   the logarithmic mass of these primes.

If \(L_H=1\) and all roots of a reference selection have
equal modulus, then for sufficiently small \(\eta\) the rational
projective shape chosen in (6) must be cyclotomic. All its roots
have modulus one, and it is monic integral. This assertion follows directly by forming
the monic polynomials of successive powers of those roots.
The polynomials are integral, with coefficients bounded
in terms of \(d\), so only finitely many occur. The set of
powers of any one root is therefore finite. This version of
the torsion restriction requires only equal moduli in the
reference and zero slope spread; the actual root moduli can
vary.

These restrictions can be combined for every possible even
common-factor degree below the absolute cutoff, with the higher
orders already covered by that cutoff. Free arithmetic growth
of the common valuation spreads can satisfy (6), and an exact
rational limiting shape can survive the shape test. Thus these
inequalities constrain the remaining mixed-tail problem without
settling all profiles.


## Exact arithmetic fixture

For the polynomial \(H(z)=(z^2+3)(z^2+6)\), the degree is four
and \(c=18\). Its projective polynomial is
\[
  B_H(Y)=(Y-\tfrac12)^2(Y-2)^2
         =Y^4-5Y^3+\tfrac{33}{4}Y^2-5Y+1.
\]
At \(p=2\) the slopes give root valuations
\(0,0,\tfrac12,\tfrac12\), so (1) gives \(v_2(L_H)=2\). At
\(p=3\) all four valuations are \(1/2\); the monic stepped
residue is \(R(Z)=(Z+1)(Z+2)=Z^2-1\) in \(\mathbb F_3[Z]\).
Formula (3) gives \(\overline{B_H}(Y)=(Y+1)^4\), and (4)
holds with \(g=2\). The rational shape itself has non-torsion
roots. This fixture concerns the local formulas for a monic
factor, with no assertion that \(H\) divides a factorial pair.
