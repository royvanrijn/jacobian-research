# Deep rational formal images force unramified division extensions

A simultaneous elliptic Wieferich prime for the two rational directions
forces a specific large quotient of the ideal class group of the
corresponding division field. The implication works at every power
`p^r`, uniformly outside one fixed finite prime set. It is an
unconditional arithmetic restriction on the remaining formal-depth
problem, not a density theorem for its primes. Rank two already forces
one full torsion-module quotient at all sufficiently large primes;
additional formal depth supplies part or all of a second copy.

## 1. The curve and uniform Kummer input

Retain the existing curve

\[
 E:y^2=X^3-10758382X+6600759744,\qquad X=x-3552,
 \qquad \Gamma=E(\mathbb Q),
\]

and the two independent rational points `P,Q` of ETC3a. Their
independence and the non-CM assertion are retained in
[EKX1](ELLIPTIC_EXCEPTION_KUMMER_BRANCH_EXCLUSION.md).
Consequently `(P,Q)` generates a Zariski-dense cyclic group in
`E^2`: a proper connected subgroup would give a nonzero integer
endomorphism relation, as in EKX1's nondegeneracy proof. A full
rational Mordell--Weil basis is not assumed.

At a good prime `p>=7`, rational two-torsion and Hasse's bound give
`p` prime to `N_p=#E(F_p)`. Let `e_p>=1` be the actual rational
formal image depth of EHD3:

\[
 \overline\Gamma\cap E_1(\mathbb Q_p)=E_{e_p}(\mathbb Q_p).
                                                               \tag{1.1}
\]

Here `E_n` is the usual formal filtration with
`v_p(-X/y)>=n`. The closure is taken in `E(Q_p)`.

For `r>=1`, define

\[
 K_{p,r}=\mathbb Q(E[p^r]),\qquad
 M_{p,r}=K_{p,r}([p^r]^{-1}P,[p^r]^{-1}Q),\qquad
 G_{p,r}=\operatorname{Gal}(K_{p,r}/\mathbb Q).
                                                               \tag{1.2}
\]

The notation includes the full division fibres. Corvaja--Demeio--
Javanpeykar--Lombardo--Zannier, Theorem 6.9, applied to `E^2/Q`
and `(P,Q)`, gives a constant `c_2` such that

\[
 [E[p^r]^2:\operatorname{Gal}(M_{p,r}/K_{p,r})]\le c_2
                 \quad\text{for every }p,r.                   \tag{1.3}
\]

The identification is its Kummer translation map. In particular,
for every `p>c_2` and every `r>=1`,

\[
 \operatorname{Gal}(M_{p,r}/K_{p,r})\simeq E[p^r]^2:
                                                               \tag{1.4}
\]

a proper subgroup of this finite `p`-group would have index at
least `p>c_2`. The prime bound here does not depend on `r`.
The directly read primary interface is Theorem 6.9, printed p.2144,
and equation (13), printed pp.2142--2143, in the
[retained published PDF](sources/cdjlz-2022-ramified-covers.pdf),
with [version metadata](sources/cdjlz-2022-ramified-covers.json).
No Frobenius-at-`p` inference is made from this global degree theorem.

For each bad rational prime `ell`, let

\[
 m_\ell=|E(\mathbb Q_\ell)/E_1(\mathbb Q_\ell)|.
\]

These are finite: the local group is compact and the formal subgroup
is open. Fix a finite set `Sigma_E` containing the bad primes,
all primes dividing these finitely many `m_ell`, all primes at most
`c_2`, and the primes below seven. This set is independent of `r`.
We assert existence of this finite set, not a computed Kummer bound
or a displayed complete numerical list.

## 2. The class-group quotient at every level

**Theorem EWC1 (formal depth to a class-group module).** For every
prime `p notin Sigma_E` and every integer `r>=1`,

\[
 e_p\ge r+1
 \quad\Longrightarrow\quad
 M_{p,r}/K_{p,r}\text{ is everywhere unramified and abelian},
                                                               \tag{2.1}
\]

with Galois group `E[p^r]^2`. There is a
`G_(p,r)`-equivariant surjective homomorphism

\[
 \boxed{\operatorname{Cl}(K_{p,r})/p^r\operatorname{Cl}(K_{p,r})
                          \twoheadrightarrow E[p^r]\oplus E[p^r].}
                                                               \tag{2.2}
\]

In particular

\[
                  \boxed{p^{4r}\mid h(K_{p,r}).}                \tag{2.3}
\]

**Local divisibility at the residue prime.** Since `p` does not
divide `N_p`, the local group is the direct sum of its canonical
prime-to-`p` torsion lift and its formal group. Multiplication by
`p^r` is invertible on the first summand. On the formal part use
`t=-X/y`. Milne II, Lemma 5.9, printed p.67/PDF p.75 of
[the retained elliptic-curves book](sources/milne-elliptic-curves-2006.pdf),
gives `t(U+V)=t(U)+t(V)+O(p^(5n))` for `U,V in E_n` on
this integral short model. Adding `p` copies of a point in `E_n`
therefore gives

\[
 t([p]U)=p\,t(U)+O(p^{5n}).
\]

If `U` has exact depth `n`, its image has exact depth `n+1`.
Moreover the induced map `E_n/E_(n+1) -> E_(n+1)/E_(n+2)`
is an isomorphism: the error vanishes modulo `p^(n+2)`, and the
leading map is multiplication by `p` on the parameter. Successive
lifts in these graded quotients and completeness produce a preimage
of every point in `E_(n+1)`; exact depth also gives injectivity.
Thus `[p]E_n=E_(n+1)` for every `n>=1`. Iterating gives

\[
 [p^r]E_1(\mathbb Q_p)=E_{r+1}(\mathbb Q_p).
                                                               \tag{2.4}
\]

Equation (1.1) consequently gives the exact equivalence

\[
       e_p\ge r+1\quad\Longleftrightarrow\quad
                   \Gamma\subseteq[p^r]E(\mathbb Q_p).         \tag{2.5}
\]

To see the torsion splitting explicitly, lift any reduced point,
multiply it by `N_p`, and uniquely divide that formal point by
the `p`-unit `N_p`. Subtracting the result gives its unique
`N_p`-torsion lift. Thus taking the closure in (1.1) loses no local
divisibility information.

Under the left side of (2.5), both `P,Q` have `p^r`-division points
in `E(Q_p)`. Every other division point differs from these by
`E[p^r]`. The completion of `K_(p,r)` already contains that
torsion. Hence `M_(p,r)/K_(p,r)` is locally trivial at every prime
above `p`.

This step uses the retained formal addition estimate and the graded
lifting argument above. Dainobu, section 4.2, printed/PDF pp.16--17,
gives an alternative proof using the formal logarithm; that is
comparison context here, not an additional premise. We use `p>=7`.

**Local divisibility at the bad primes.** At a fixed bad prime
`ell`, multiplication by `p^r` is an automorphism of the formal
pro-`ell` group because `p!=ell`. It is also an automorphism of
the finite quotient of order `m_ell` because `p` does not divide
`m_ell`. Lifting in the quotient and then dividing in the formal
subgroup proves surjectivity on `E(Q_ell)`; the same two steps
prove its kernel is trivial. Thus `P,Q` have local division points
already in `E(Q_ell)`. As above, adjoining the full division fibres
over `K_(p,r)` creates no extension of its completions at these
bad primes. In particular it is unramified there.

**All other finite primes.** At a good prime `ell!=p`, the curve
extends to a proper abelian scheme and `[p^r]` is finite etale.
Properness extends `P,Q` to sections even if their affine
coordinates have denominators. Its division fibres and the torsion
field are therefore unramified at `ell`. In particular
`M_(p,r)/K_(p,r)` is unramified there. The exact primary statement
is Milne AVs Proposition 20.7 and Corollary 20.8, printed/PDF p.44
of [the retained source](sources/milne-abelian-varieties-2022.pdf).

**Infinite primes, reciprocity and equivariance.** The Weil pairing
and canonical elliptic selfduality give `mu_p subset K_(p,r)`.
Since `p` is odd, this field is totally imaginary, so there is no
infinite-prime ramification to check. The retained pairing interface
is Milne AVs, section 16 opening, printed/PDF p.29.

Together with (1.4), the preceding local arguments prove (2.1).
Global unramified class field theory now supplies a surjective Artin
map from the ideal class group to this abelian Galois group. We
use Milne, *Class Field Theory*, version 4.03, Theorem 0.4 and its
following Hilbert class field paragraph, printed p.4/PDF p.13,
and Theorem 0.8, printed p.9/PDF p.18. The latter identifies the
canonical map sending prime ideals to Frobenius elements and its
factorization through a ray class group with ramified support.
Here there is no ramification, so the modulus is one and the map
factors through the ordinary ideal class group:
[primary source](https://www.jmilne.org/math/CourseNotes/CFT.pdf),
[retained PDF](sources/milne-cft403.pdf),
[metadata](sources/milne-cft403.json).

Both fields in (1.2) are Galois over `Q`, since `P,Q` are rational.
Under Kummer translation, conjugation by a lift of
`g in G_(p,r)` is the diagonal action of `g` on `E[p^r]^2`:
the affine choices of division points cancel. Artin reciprocity
commutes with the action on prime ideals and conjugation of their
Frobenius elements. Thus its quotient map is equivariant.
The target has exponent `p^r`, so the map factors through the
left side of (2.2). Its order is `p^(4r)`, proving (2.3).
\(\square\)

## 3. Consequences for the remaining fatal-prime problem

**Corollary EWC2 (a class-group sufficient condition for safety).**
For an active prime `p notin Sigma_E`, absence of a
`G_(p,1)`-equivariant quotient `E[p]^2` from
`Cl(K_(p,1))/p Cl(K_(p,1))` implies `e_p=1`.
In particular

\[
 v_p(h(K_{p,1}))<4\quad\Longrightarrow\quad e_p=1.             \tag{3.1}
\]

Such a prime is safe in both unresolved classes by EDC2--EDC3 in
[the dependency-core note](ELLIPTIC_EXCEPTION_DEPENDENCY_CLOSED_CORE.md).
Conversely, every sufficiently large fatal prime must have the
specific quotient (2.2) at level one. More generally, absence of
the level-`r` quotient implies `e_p<=r`.

**Proof.** These are contrapositives of EWC1, followed by the
core-adjusted local safety criterion. \(\square\)

There is no converse asserting that such a class-group quotient
forces a deep formal image or a fatal numerator branch. The quotient
may have other origins; the reduction support and branch offset
condition `k_p<e_p`, `3` not dividing `k_p`, remain separate.
No frequency or summability estimate for these class-group modules
as `p` varies is established. Thus the bridge neither proves that
there are finitely many simultaneous elliptic Wieferich primes nor
settles all-finite-set positivity or the all-primes Hasse measure.

## 4. A conditional route from class-group density to weighted means

**Conditional corollary EWC3 (density propagation).** For each good
prime `p`, put `K_p=Q(E[p])`, `G_p=Gal(K_p/Q)`, and define

\[
 \mathcal B=\left\{p:\operatorname{Cl}(K_p)/p\operatorname{Cl}(K_p)
       \text{ admits a }G_p\text{-equivariant quotient }E[p]^2\right\}.
                                                               \tag{4.1}
\]

If `B` has natural density zero among primes, then the primes with
`e_p>1`, and therefore the active primes with `3|e_p`, have natural
prime density zero. Under this unproved hypothesis, every fixed
positive finite-place event `C_0={alpha=a} cap intersection_(q in S) G_q`
satisfies

\[
 \mathbb E\!\left(\sum_{\substack{p\le Y\\p\text{ active}}}h_p
                           \,\middle|\,C_0\right)
                         \gg_{S,a}\log\log Y.                  \tag{4.2}
\]

Here `h_p` is the actual predictable hazard of ETP4--ETP5, set to
zero before its fixed cutoff, as in
[EAM4](ELLIPTIC_EXCEPTION_AFFINE_COSET_POLE_MEANS.md).

**Proof.** Outside the fixed finite set `Sigma_E`, EWC1 at level one
gives `{p:e_p>1} subset B`. Density zero passes to subsets and is
unchanged by a finite set. Since `{p:3|e_p}` is contained in
`{p:e_p>1}`, the hypothesis of EAM4 holds. Its conclusion is (4.2).
\(\square\)

The class-group density hypothesis is not established here. Even
(4.2) would be only a first-moment statement: no variance estimate,
almost-sure recurrence, single-vector tail bound or all-primes Hasse
conclusion is inferred.

## 5. The baseline quotient and the extra contribution of depth

Let `Gamma_0=<P,Q,T>` and let `I=[Gamma:Gamma_0]`. The index is
finite because `P,Q` are independent in the rank-two rational group
and `T` generates its torsion. Enlarge the fixed finite exclusion set
to `Sigma'_E=Sigma_E union {p:p divides I}`. This enlargement is still
independent of `r`; no numerical value of `I` or a full displayed
Mordell--Weil basis is needed.

**Theorem EWC4 (mixed-level quotient).** For every prime
`p notin Sigma'_E` and every integer `r>=1`, put
`s=min(r,e_p-1)`. Then there is a `G_(p,r)`-equivariant surjection

\[
 \boxed{\operatorname{Cl}(K_{p,r})/p^r\operatorname{Cl}(K_{p,r})
                      \twoheadrightarrow E[p^r]\oplus E[p^s].}
                                                               \tag{5.1}
\]

Interpret `E[p^0]` as the zero group. In particular,

\[
                  \boxed{v_p h(K_{p,r})\ge 2r+2\min(r,e_p-1).}  \tag{5.2}
\]

Thus `E[p^r]` is a quotient and `p^(2r)` divides `h(K_(p,r))`
for every such `p,r`, including when `e_p=1`. At `e_p>=r+1`,
(5.1) recovers the two full copies in EWC1.

**A procyclic formal coordinate.** The graded proof of (2.4) also
constructs a group isomorphism `f:E_1(Q_p)->Z_p`, without importing
a formal logarithm. Choose `U` of exact formal depth one. For any
`V in E_1`, its class modulo `E_2` determines a unique digit
`a_0 in {0,...,p-1}` such that `V-[a_0]U in E_2`. Uniquely divide
this difference by `[p]`, and repeat. Completeness yields a unique
expansion

\[
                    V=\lim_{n\to\infty}
                           [a_0+a_1p+\cdots+a_np^n]U.
\]

The inverse of `a -> [a]U` is the required continuous group
isomorphism, and `f(E_n)=p^(n-1)Z_p`. Compose it with the canonical
formal projection in the torsion splitting of EWC1, writing this
homomorphism as `f_0:E(Q_p)->Z_p`.

Put `u=f_0(P)`, `v=f_0(Q)`. The closure of the actual rational
formal image is `E_(e_p)`. Since `[I]Gamma subset Gamma_0` and
`p` does not divide `I`, the formal closures of `Gamma_0` and
`Gamma` agree: in `Z_p`, multiplication by the unit `I` preserves
every closed subgroup. Also `f_0(T)=0`. Consequently

\[
                    u\mathbb Z_p+v\mathbb Z_p
                              =p^{e_p-1}\mathbb Z_p.            \tag{5.3}
\]

**The coefficient lattice.** Work modulo `p^r`. A Smith basis for
the row `(u,v)` gives an invertible coefficient matrix
`A in GL_2(Z/p^r Z)` for which the first new row entry is zero,
and the second has valuation `e_p-1` when `e_p-1<r`; if
`e_p-1>=r`, both entries are zero. Lift the two columns of `A`
to integers and form the corresponding rational combinations
`D_1,D_2` of `P,Q`. Their formal coordinates satisfy

\[
 f_0(D_1)\in p^r\mathbb Z_p,\qquad
 f_0(D_2)\in p^s\mathbb Z_p.
\]

The finite prime-to-`p` torsion summand is divisible by both powers,
so (2.4) gives

\[
                  D_1\in[p^r]E(\mathbb Q_p),\qquad
                  D_2\in[p^s]E(\mathbb Q_p).                  \tag{5.4}
\]

Equivalently, the coefficient lattice for local `p^r` divisibility
has Smith form `(Z/p^r Z) directsum p^(r-s)(Z/p^r Z)`.

**The global quotient and its uniform degree.** Inside the fixed
original division extension `M_(p,r)` of (1.2), take

\[
 F_{p,r}=K_{p,r}([p^r]^{-1}D_1,[p^s]^{-1}D_2),                 \tag{5.5}
\]

omitting the second fibre when `s=0`. Use the full Kummer group
`E[p^r]^2` of the original fixed points `P,Q`, supplied uniformly
by (1.4). The matrix `A` acts invertibly on this group. Projecting
its second new coordinate by multiplication by `p^(r-s)` gives a
surjection onto

\[
              \operatorname{Gal}(F_{p,r}/K_{p,r})
                            \simeq E[p^r]\oplus E[p^s].        \tag{5.6}
\]

This is a quotient of the fixed-point Kummer extension; the deep
Kummer theorem is not reapplied to the varying points `D_1,D_2`.
The bound on `p` is therefore independent of `r` and of their
integer coefficients. Both points are rational, so their fields are
Galois over `Q`, and the action in (5.6) is the usual diagonal
torsion action.

By (5.4), this extension is locally trivial above `p`. At each bad
prime the local multiplication maps by `p^r` and `p^s` are
bijective by the same finite-quotient argument as in EWC1. At all
other finite primes, the division fibres are finite etale on a
proper good model. Thus `F_(p,r)/K_(p,r)` is everywhere unramified;
infinite primes cause no issue because `K_(p,r)` is totally
imaginary. The equivariant Artin argument of EWC1 now gives (5.1).
Its target has order `p^(2r+2s)`, proving (5.2). \(\square\)

The density hypothesis in EWC3 concerns the existence of a **second**
full `E[p]` copy at level one. One copy is already unconditional
for every sufficiently large prime. In particular, vanishing of the
`p`-part of these class groups is not a possible density criterion
here. The mixed quotient does not prove a converse or any density
estimate for the additional copy.

## 6. Relation to the established class-group framework

Dainobu's *Ideal class groups of division fields of elliptic curves
and everywhere unramified rational points*, Journal of Number
Theory 264 (2024), 211--232, DOI
[10.1016/j.jnt.2024.05.007](https://doi.org/10.1016/j.jnt.2024.05.007),
introduces everywhere unramified rational divisibility in
Definition 1.3 and maps it into class-group homomorphisms in
Proposition 2.1. Theorem 4.1 gives class-number bounds under its
cohomological and irreducibility assumptions; Proposition 4.5
isolates the residue-prime local condition under its Tamagawa
hypothesis. These are directly read in
[arXiv v3](sources/dainobu-2024-unramified-rational-points.pdf),
with [version metadata](sources/dainobu-2024-unramified-rational-points.json).

EWC1 is a specialization of this established Kummer/class-field
mechanism, proved directly with the retained uniform Kummer input.
It uses the stronger local divisibility in `E(Q_p)` and the fixed
finite exclusion set at the other bad primes. We do not import
Dainobu's cohomology, Tamagawa or irreducibility hypotheses without
checking them, and make no novelty claim. The direct construction
identifies an actual module quotient, rather than relying only on
a class-number divisor or a semisimplified component.

The uniform Kummer theorem is used as a degree statement for each
`p^r`. Chebotarev is not applied at the ramified residue prime of
that division field. No artificial finite checker is added to these
written local and class-field arguments, and no numerical Kummer
threshold or division-field class number is computed.
