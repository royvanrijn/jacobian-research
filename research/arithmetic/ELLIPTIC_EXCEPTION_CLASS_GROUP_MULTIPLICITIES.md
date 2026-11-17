# Class-group multiplicities recover the residual Sha dimension

The existing local two-isogeny descent proves that this elliptic
curve has trivial two-primary Tate--Shafarevich group. The proved
all-prime Selmer-parity theorem then makes `dim Sha[p]` even for
every prime, without assuming that `Sha` is finite. In the
strict-Selmer converse range, the number of actual `E[p]` quotient
copies in the division-field class group determines this dimension
exactly. Two copies with no third copy certify a deep formal image
and vanishing `Sha[p]`.

## 1. The retained descent also computes the two-isogeny Selmer groups

Use the two curves from
[EEA Theorem 2.2](ELLIPTIC_EXCEPTION_ARITHMETIC.md):

\[
 \begin{aligned}
 E_0 &: y^2=x^3-10656x^2+27091730x,\\
 E_0'&: y^2=x^3+21312x^2+5183416x,
 \end{aligned}
\]

where `q_0=2709173`, `27091730=10q_0`, and
`5183416=2^3 7^3 1889`. They are linked by dual degree-two
isogenies `phi:E_0->E_0'`, `psi:E_0'->E_0`, with
`psi phi=[2]` and `phi psi=[2]`.

For `E:y^2=x^3+ax^2+bx`, the squareclass map is

\[
 \alpha(O)=1,\quad\alpha((0,0))=b,\quad
                  \alpha((x,y))=x\bmod F^{*2}.                  \tag{1.1}
\]

Here `F` may be `Q` or any completion. The corresponding map on
`E'` is `alpha'`. These are the connecting Kummer maps for
`psi` and `phi`, respectively. The explicit isogenies and
`ker alpha=psi(E'(F))`, `ker alpha'=phi(E(F))` are Dujella,
*Number Theory*, section 15.5, printed pp.507--509:
[retained section](sources/dujella-number-theory-2021-section15-5.pdf).
The formulas are algebraic and have the same meaning over every
completion. To identify the connecting character explicitly, a
preimage of a point with nonzero `x` is defined over `F(sqrt(x))`;
the corresponding field at `(0,0)` is `F(sqrt(b))`.

Write `S_psi` and `S_phi` for the two isogeny Selmer groups,
viewed in `Q^*/Q^{*2}` through these maps.

**Theorem ESM1 (the two-primary Sha group vanishes).** One has

\[
 \begin{aligned}
 S_\psi&=\alpha(E_0(\mathbb Q))
                    =\{1,2,5q_0,10q_0\},\\
 S_\phi&=\alpha'(E_0'(\mathbb Q))
                    =\{1,14,1889,26446\}.
 \end{aligned}                                                 \tag{1.2}
\]

Consequently

\[
 \boxed{\Sha(E_0/\mathbb Q)[2^\infty]
              =\Sha(E_0'/\mathbb Q)[2^\infty]=0.}               \tag{1.3}
\]

**Why the existing exclusions are local Selmer exclusions.** A
class in `S_psi` must be represented by a signed squarefree divisor
`d` of `b`: for every prime not dividing `b`, the valuation of the
`x` coordinate of a local point is even. This is the same unique
smallest-valuation argument as EEA Theorem 2.1, applied over `Q_l`;
the identity and `(0,0)` satisfy the assertion too. Local membership
in the Kummer image is equivalent to local solubility of

\[
       N^2=dM^4+aM^2e^2+(b/d)e^4,\qquad(M,e)\ne(0,0).          \tag{1.4}
\]

When `e!=0` the map back to `E` is
`(x,y)=(dM^2/e^2,dMN/e^3)`; the special cases give the identity
and the rational two-torsion class. Over `Q_l` scale `M,e` to
`Z_l` with at least one a unit. The integral coefficients make
`N` integral, so the primitive projective residue exclusions apply
to every local solution, including those initially having
denominators. Local membership is a subgroup condition, since
(1.1) is a connecting homomorphism.

For `E_0`, the real place excludes all negative classes. Its eight
positive candidate classes are `H` and `5H`, with
`H={1,2,5q_0,10q_0}`. Every class in `H` has a rational
representative among `O,P,T,P+T`, as proved in EEA. Any locally
soluble class in `5H` would therefore give local membership of
`5`, by multiplying by one of these global classes. Its quartic
is

\[
           N^2=5M^4-10656M^2e^2+5418346e^4.
\]

At two the primitive right-side residues modulo sixteen are
`{5,10,15}`, all nonsquares. Hence `5` is absent even locally at
two. This proves the first equality in (1.2).

For `E_0'`, candidates are the sixteen signed squarefree divisors
of `2*7*1889`. The complete local exclusions in EEA Theorem 2.2 are

| Local modulus | Excluded classes |
| --- | --- |
| `16` at two | `2,7,3778,13223,-1,-14,-1889,-26446` |
| `25` at five, among the survivors | `-2,-7,-3778,-13223` |

The only remaining classes are `1,14,1889,26446`. They all have
rational representatives: use `O,T'`, the point
`Q'=(3584/25,4301248/125)` of class `14`, and `Q'+T'`.
Therefore no further local condition removes any survivor, proving
the second equality in (1.2). This uses the written valuation and
primitive-normalization arguments plus the complete residue controls
already retained in EEA, not a new rational-solubility assumption.

**From the isogeny Selmer groups to Sha.** The isogeny Kummer
sequences and their local diagrams give

\[
 \begin{aligned}
 0&\to E_0(\mathbb Q)/\psi E_0'(\mathbb Q)
         \to S_\psi\to\Sha(E_0'/\mathbb Q)[\psi]\to0,\\
 0&\to E_0'(\mathbb Q)/\phi E_0(\mathbb Q)
         \to S_\phi\to\Sha(E_0/\mathbb Q)[\phi]\to0.
 \end{aligned}                                                 \tag{1.5}
\]

These follow exactly as the multiplication-by-`p` Selmer--Sha
sequence in ECS2, replacing that multiplication with the isogeny.
The first two terms in each row coincide by (1.2), so both Sha
kernels are zero. If `xi in Sha(E_0)[2]`, then `phi(xi)` belongs
to `Sha(E_0')[psi]`, hence is zero; now `xi` belongs to
`Sha(E_0)[phi]`, hence is zero. The same argument with the curves
interchanged gives `Sha(E_0')[2]=0`. A nonzero two-primary
torsion group would contain an element of order two. This proves
(1.3). \(\square\)

The local quartic data in this argument are replayed by the existing
command

```sh
.venv/bin/python research/arithmetic/verify_elliptic_exception_arithmetic.py
```

Its 36 exact controls include the full primitive residue tables and
the global witnesses. The new Sha conclusion is the written
isogeny-cohomology argument (1.5), not a numerical Sha calculation.

## 2. Even residual Sha dimension without finiteness of Sha

Write `E` for either integral model of the original curve. For any
prime `p`, let `D_p` be the maximal divisible subgroup of
`Sha(E/Q)[p^infinity]`, and put

\[
 D_p\simeq(\mathbb Q_p/\mathbb Z_p)^{d_p},\qquad
 W_p=\Sha(E/\mathbb Q)[p^\infty]/D_p,\qquad
 t_p=\dim_{\mathbb F_p}\Sha(E/\mathbb Q)[p].
                                                               \tag{2.1}
\]

Each `d_p` is finite and `W_p` is finite. This follows from
finiteness of the bounded Selmer groups; the precise retained
interface is Milne ADT, I Remark 6.7, printed p.77/PDF p.85:
[retained source](sources/milne-arithmetic-duality-2006.pdf),
[updated locator metadata](sources/milne-arithmetic-duality-2006.json).
It asserts the cofinite structure, not finiteness of the whole
p-primary group.

**Theorem ESM2 (even p-primary corank and even p-torsion dimension).**
For this fixed curve, for every prime `p`,

\[
                          \boxed{d_p\equiv t_p\equiv0\pmod2.}   \tag{2.2}
\]

Neither `d_p=0` nor finiteness of `Sha[p^infinity]` is asserted at
an odd prime.

**The parity input and the divisible part.** Tim and Vladimir
Dokchitser, *On the Birch--Swinnerton-Dyer quotients modulo squares*,
arXiv:math/0610290v3, Theorem 1.4, printed/PDF p.2, and
Theorem 4.19, printed/PDF pp.26--27, prove the p-parity theorem for
every elliptic curve over `Q` and every prime `p`. Their notation
on p.4 is

\[
        \operatorname{rk}_p(E/\mathbb Q)=
                   \operatorname{rank}E(\mathbb Q)+d_p.
\]

The parity of this integer is independent of `p`, being the
analytic-rank parity or equivalently the root-number sign. The
all-prime theorem is unconditional; it is separate from their
preceding theorem about algebraic rank under finiteness of Sha:
[primary arXiv v3](https://arxiv.org/pdf/math/0610290v3),
[retained PDF](sources/dokchitser-dokchitser-2008-parity-v3.pdf),
[version metadata](sources/dokchitser-dokchitser-2008-parity-v3.json).
This deep parity theorem is imported, not reproved.

By ESM1, `d_2=0`. The established rational rank is two, so
`rk_2=2`. Comparing the parity at `p` with the parity at two
makes `2+d_p` even. Hence `d_p` is even, without knowing whether
it is zero. This does not compute the analytic rank; it only
fixes its parity.

**The finite reduced part.** The Cassels--Tate pairing has left
and right radicals equal to the divisible subgroup by Milne ADT,
I Theorem 6.13(a), printed p.82/PDF p.90. The canonical elliptic
polarization is induced by the rational divisor `(O)`, so the
pairing is alternating, by Remark 6.10(b), printed p.79/PDF p.87.
It therefore induces a perfect alternating pairing on the finite
group `W_p`.

Such a finite abelian p-group has paired cyclic elementary divisors.
Indeed choose an element of maximal order `p^n`; perfection gives
a second element pairing with it to exact order `p^n`. These two
generate a nondegenerate `(Z/p^n)^2` summand. Split off its
orthogonal complement and repeat. Consequently `dim W_p[p]` is
even. This uses the pairing on `W_p`, not a claimed nondegenerate
pairing on `Sha[p]` or on `W_p[p]`; those restrictions can be
degenerate.

Since `D_p` is divisible, the p-torsion sequence gives

\[
              0\to D_p[p]\to\Sha(E/\mathbb Q)[p]
                                      \to W_p[p]\to0.
\]

Thus `t_p=d_p+dim W_p[p]` is even, completing the proof of
(2.2). \(\square\)

For a general elliptic curve, Cassels--Tate pairing alone only gives
`t_p congruent d_p mod2`. The extra input here is the proved
vanishing at two and the all-prime parity theorem. We have not
replaced an unknown divisible subgroup with zero. Finiteness of
`Sha[p]` holds unconditionally; finiteness of `Sha[p^infinity]`
is a different assertion and is not established here at odd primes.

## 3. Actual quotient multiplicities recover the Sha dimension

Retain the notation of
[ECS1--ECS3](ELLIPTIC_EXCEPTION_CLASS_GROUP_SELMER_CONVERSE.md):

\[
 K_p=\mathbb Q(E[p]),\quad G_p=\operatorname{Gal}(K_p/\mathbb Q),\quad
 c_p=\dim_{\mathbb F_p}
       \operatorname{Hom}_{G_p}(\operatorname{Cl}(K_p)/p,E[p]).
                                                               \tag{3.1}
\]

By [EOI2](ELLIPTIC_EXCEPTION_ADELIC_OPEN_IMAGE.md), all ECS image
and local hypotheses hold outside one fixed finite set of primes.
Every assertion below is in that finite-exception range.

**Theorem ESM3 (residual Sha dimension from the actual class-group
quotient multiplicity).** The integer `c_p>=1` is the largest `m`
for which `Cl(K_p)/p` has a `G_p`-equivariant quotient `E[p]^m`.
It determines

\[
                  \boxed{t_p=2\left\lfloor\frac{c_p-1}{2}
                                      \right\rfloor.}          \tag{3.2}
\]

In particular:

| Actual multiplicity | Exact conclusion |
| --- | --- |
| `c_p=1` | `e_p=1` and `Sha[p]=0` |
| `c_p=2` | `e_p>1` and `Sha[p]=0` |
| `c_p>=3` | `Sha[p]!=0`; the formal depth can require further information |
| `c_p` even | `e_p>1` |

The first three rows are equivalences. Thus at a nondeep prime a
second quotient copy is accompanied by at least a third copy.

**Proof of the multiplicity interpretation.** Under the ECS image
hypothesis, `V=E[p]` is simple and `End_(G_p)(V)=F_p`. Choose a
basis `f_1,...,f_(c_p)` of the Hom space and evaluate to obtain
`Cl(K_p)/p -> V^(c_p)`. If its image were proper, its nonzero
quotient inside this semisimple isotypic module would admit a map
to `V`. Composing would give a nontrivial scalar linear relation
among the `f_i`, a contradiction. Thus evaluation is surjective.
A quotient `V^m` conversely gives `m` independent maps, so
`m<=c_p`. This concerns actual quotients; it is not an assertion
about all composition factors of `Cl(K_p)/p`.

**Proof of the recovery formula.** ECS2 gives

\[
 c_p=\begin{cases}
       1+t_p,&e_p=1,\\
       2+t_p-\rho_p,&e_p>1,
      \end{cases}
 \quad\rho_p=\operatorname{rank}\delta_p\in\{0,1\}.             \tag{3.3}
\]

Here `delta_p` is exactly the localization obstruction defined in
ECS equation (3.2); it is not set to zero. In all cases `c_p` is
one of `t_p+1,t_p+2`. By ESM2, `t_p` is even. These two adjacent
integers determine that even integer by (3.2).

If `t_p=0`, then `delta_p=0`, so (3.3) gives multiplicities one
and two exactly in the nondeep and deep cases. If `t_p!=0`,
evenness gives `t_p>=2`, hence `c_p>=3`; conversely `c_p>=3`
forces `t_p!=0` by (3.2). Finally a nondeep multiplicity is
`1+t_p`, which is odd, so an even multiplicity forces a deep
formal image. \(\square\)

More precisely an even `c_p` has `e_p>1` and `rho_p=0`.
An odd `c_p>=3` may be nondeep, or deep with `rho_p=1`.
No actual nonvanishing of `delta_p` is asserted or excluded.
For example, the abstract exact-sequence data alone allow
`M=F_p^2`, `X=F_p^2`, `L=F_p`, `S=M directsum X`, with
`loc(M)=0` and either zero localization on `X` or a nonzero
coordinate map. Both have a perfect alternating pairing on `X`.
This is only a linear-algebra illustration, not an arithmetic
realization or a substitute for a Poitou--Tate calculation.

## 4. What this changes in the density problem

Outside the finite exceptional set, the exact nonzero-Sha prime set
is `{p:c_p>=3}`. The primes with a deep formal image and zero
`Sha[p]` are exactly `{p:c_p=2}`. The two-copy set in EWC3 is
the disjoint union of these two multiplicity ranges. Its density-zero
condition still requires control of both rational simultaneous
Wieferich primes and the Sha contribution, as proved in EOI2.
No density estimate for either multiplicity range is obtained here.

The parity result does not prove that `Sha` is finite, that its
odd-primary divisible coranks vanish, or that infinitely many primes
satisfy `e_p=1`. It also does not decide the branch-offset condition
needed for fatal primes, any conditional correlation estimate, or
the all-primes Hasse measure. The distinction between fixed-prime
module structure and varying-prime density remains essential.

All new cohomological and group-theoretic steps are the written
arguments above. The p-parity theorem and Cassels--Tate duality
are primary imports with exact locators. No new artificial finite
checker, formal verification, external review or analytic-rank
computation is claimed.
