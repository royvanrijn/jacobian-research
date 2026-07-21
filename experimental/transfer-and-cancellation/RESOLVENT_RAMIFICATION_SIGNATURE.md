# Resolvent--ramification signature and polynomial equivalence

This note separates two objects which are easy to conflate:

1. the intrinsic finite extension and its boundary valuations; and
2. a chosen one-variable inverse resolvent, its derivative, and displayed
   reconstruction formulas.

Only the first object is automatically invariant under polynomial left--right
equivalence.  A resolvent becomes an invariant certificate only after its
critical and polar divisors have been identified with valuations in the
canonical finite model.

Throughout, `k` has characteristic zero, `X=Y=A^d_k`, and
`F:X -> Y` is a dominant generically finite Keller map.  Geometric statements
are made after base change to an algebraic closure.  Arithmetic residue
degrees are retained separately.

## 1. The intrinsic signature

Put `K=k(Y)` and `L=k(X)`, with the inclusion induced by `F`, and let

\[
 \bar X_F=\operatorname{Norm}_Y(L)\longrightarrow Y                 \tag{1}
\]

be the normalization of `Y` in `L`.  Since `Y` is excellent, (1) is finite.
The map `F` factors through (1).  The induced map
`j_F:X -> bar(X)_F` is a quasi-finite birational morphism to a normal scheme,
hence an open immersion by Zariski's Main Theorem.  Thus

\[
 \partial_F=\bar X_F\setminus X                                  \tag{2}
\]

is canonical as a pair with `X`, even though a projective closure or a
particular sequence of blow-ups is not.

### Definition

The **intrinsic resolvent--ramification signature** `Sigma(F)` consists of:

- the generic degree `nu=[L:K]`;
- the geometric monodromy permutation group of the Galois closure of `L/K`;
- for every geometric prime divisor `Z` of `Y` met by `partial_F`, the
  multiset
  \[
    \mathcal B_Z(F)=\{(e(E/Z),f(E/Z)):E\subset\partial_F, E\mapsto Z\};
                                                                    \tag{3}
  \]
- the generic sheet-loss length
  \[
    \ell_Z(F)=\sum_{E\subset\partial_F, E\mapsto Z}e(E/Z)f(E/Z);  \tag{4}
  \]
- for every boundary prime `E`, the valuation functional
  `v_E:k[X]^* -> Z`, considered up to automorphism of the `k`-algebra
  `k[X]`.

Over an algebraically closed field all geometric residue degrees in (3) are
one.  Over a nonsplit coefficient field, one geometric orbit is a single
arithmetic prime and its residue degree is the orbit size.  Formula (4) is a
length: it counts nearby generic sheets specializing to boundary.  It should
not be confused with the number of points in a reduced special fiber.
Indeed, extension of the DVR at the generic point of `Z` gives
`nu=sum_E e(E/Z)f(E/Z)` (there is no defect in characteristic zero); (4) is
the part of this sum supported outside the affine-source open.

The last entry is the coordinate-free meaning of **reconstruction pole
divisors**.  Choosing coordinates `(x_1,...,x_d)` replaces the functional by
the vector `(v_E(x_1),...,v_E(x_d))`; that vector is not invariant under a
source automorphism.  Its support, as a boundary valuation of the affine
coordinate ring, is intrinsic.

### Resolvent decoration

Suppose a primitive element `theta in L` gives a one-variable equation

\[
 \Phi(T)=0,qquad \Phi\in K[T].                                  \tag{5}
\]

For a generic target parameter, write the divisor of `Phi_T` on the
`T`-line as a partition

\[
 (a_1,\ldots,a_s),                                                \tag{6}
\]

so the corresponding local resolvent indices are `(a_i+1)`.  Also record the
polar divisors of the displayed rational reconstruction formulas.  These are
the **presented resolvent data**.

They join `Sigma(F)` only after the following checks at every generic critical
prime:

1. the resolvent hypersurface is normal there (or has been normalized);
2. the critical prime maps to a boundary prime `E` of (1);
3. a transverse uniformizer has local equation `z=u^{a_i+1}` up to a unit;
4. every displayed reconstruction pole is classified as either a genuine
   boundary valuation or a removable chart pole.

Under these hypotheses `a_i+1=e(E/Z)`, so the index is intrinsic.  Without
them, neither (6) nor a denominator list is an invariant of `F`.

### Adversarial examples

There are two useful failure modes.

- In `k(s) subset k(t)`, `t^3=s`, the primitive element `u=t+t^2` has two
  conjugates which coincide over `s=1`.  Its power-basis resolvent therefore
  acquires a discriminant factor at `s=1`, although the normalized extension
  is unramified there.  A naïve resolvent discriminant records failure of the
  chosen generator, not just field ramification.
- In C24, `s=x/A` has a pole on the affine source divisor `A=0`.  In the
  projective resolvent chart `U=Ps`, exactly one nonzero infinity root is that
  finite divisor.  Counting all poles of `s` as dicritical would therefore
  overcount the boundary by one.  The two-chart normalization is essential.

Polynomial target automorphisms give a third elementary warning: the degree
and equation of a displayed discriminant hypersurface can change.  Its prime
and the data (3), not its coordinate degree, are invariant.

## 2. Invariance theorems

### Theorem 2.1 (source automorphisms)

If `alpha:X -> X` is a polynomial automorphism, then `F` and `F alpha` have
the same intrinsic signature.  Their presented resolvent data need not agree.

**Proof.**  The automorphism `alpha^*:L -> L` identifies the two extensions
over `K`.  Functoriality and uniqueness of normalization identify their
finite models and carry the open copy of `X` to the open copy of `X`.
Therefore boundary primes, extensions of divisorial valuations, `e`, `f`,
and (4) are carried bijectively.  The Galois closures are isomorphic as
permutation extensions.  A chosen source coordinate or primitive element is
not fixed by `alpha`, which is why raw reconstruction vectors and (6) are not
asserted invariant.  QED

### Theorem 2.2 (target automorphisms)

If `beta:Y -> Y` is a polynomial automorphism, then `F` and `beta F` have
isomorphic intrinsic signatures, with every target divisor transported by
`beta`.  No equation, ordinary polynomial degree, or named coordinate
hyperplane is preserved in general.

**Proof.**  The automorphism `beta^*:K -> K` transports the inclusion
`K subset L`.  Integral closure commutes with this isomorphism.  The induced
isomorphism of finite normal models preserves extensions of DVRs and hence
`e`, `f`, and the fundamental length sums (4).  It also identifies the Galois
closures and their actions on the `nu` embeddings.  QED

### Corollary 2.3 (left--right equivalence)

If `G=beta F alpha`, then `Sigma(G)=Sigma(F)` after transporting divisors by
`beta`.  This is exactly the combination of Theorems 2.1 and 2.2.

### Theorem 2.4 (adjoining identity coordinates)

For `F^+=F times id_(A^a)`, the nontrivial part of `Sigma(F^+)` equals
`Sigma(F)`: the degree and geometric monodromy are unchanged, every boundary
prime is `E times A^a` over `Z times A^a`, its `e` and `f` are unchanged, and
its sheet-loss length is unchanged.  Consequently this core is invariant
under polynomial stable left--right equivalence.

**Proof.**  The new field extension is the purely transcendental base change
`L(t_1,...,t_a)/K(t_1,...,t_a)`.  Linear disjointness from a purely
transcendental extension preserves degree and the geometric Galois closure.
Because `bar(X)_F` is normal in characteristic zero,
`bar(X)_F times A^a` is normal and is the integral closure after this base
change.  The identity variables have valuation zero at `E times A^a`; the
extensions of the corresponding DVRs have the same value groups and residue
extensions.  This proves all assertions.  An arbitrary stable left--right
equivalence is then handled by Theorems 2.1 and 2.2.  QED

This theorem does **not** say that a particular one-variable resolvent survives
mixing old and new coordinates.  Only the intrinsic finite-cover data do.

## 3. C01 and a generic weighted seed

For C01 the projective cubic resolvent is

\[
 cU^3-2U^2V+bUV^2-2aV^3.                                \tag{7}
\]

It has degree `3`, monodromy `S_3`, and generic critical partition `1^2`.
There is one geometric boundary divisor over the irreducible discriminant,
with `(e,f)=(2,1)` and sheet loss `2`.  There is no second generic
nonproperness component.  The two affine/projective root charts show that no
root at projective infinity is lost.

Now let `H` be a generic admissible degree-`n` weighted seed.  Geometrically

\[
 H(W)=hW^2(W-1)\prod_{j=1}^{n-3}(W-\rho_j),               \tag{8}
\]

with distinct nonzero `rho_j`.  Its inverse pencil is

\[
 H(W)-BCW+cAC^2.                                         \tag{9}
\]

The universal pencil theorem gives degree `n` and monodromy `S_n`.  The
generic critical partition is `1^(n-1)`.  The boundary table is

| target prime | geometric boundary primes | `(e,f)` | sheet loss |
|---|---:|---|---:|
| discriminant `Delta_H` | 1 | `(2,1)` | 2 |
| `C=0` | `n-3` | `(1,1)` each | `n-3` |

Over a nonsplit field, the second row is grouped according to the irreducible
factors of `H/(W^2(W-1))`; their degrees are the residue degrees.  This is the
generic, not the canonical, weighted family.  For
`H=W^(n-1)(1-W)` the `C=0` row is instead one zero-cluster prime with
`(e,f)=(n-2,1)` and loss `n-2`.

C01 is the `n=3` minimal case of the generic table.

## 4. C24 for arbitrary `(m,r)`

Put

\[
 N=r(m+1)+1,
\quad D(T)=1-T(Q-PT)^m,
\quad \Psi(T)=C\int_0^T D(t)^r\,dt-R.                    \tag{10}
\]

The reconstruction formulas in the master construction prove that (10) is a
primitive inverse resolvent and that the generic degree is `N`.

### 4.1 Critical divisor, ramification, and monodromy

One has

\[
 \Psi_T=C D(T)^r.                                        \tag{11}
\]

The divisor `D=0` is geometrically integral: with `Y=Q-PT`, its coordinate
ring is

\[
 k[P,Q,T]/(1-TY^m)\cong k[Q,Y,Y^{-1}],
\quad T=Y^{-m},\quad P=(Q-Y)Y^m.                         \tag{12}
\]

For generic `(P,Q)`, `D` has `m+1` distinct roots, each of multiplicity `r`
in (11).  Hence the critical partition is

\[
 r^{m+1},                                                \tag{13}
\]

and the local ramification index is `r+1`.

The critical values are generically distinct.  It suffices to specialize
`Q=0`, `P!=0`: the critical points differ by `(m+1)`-st roots of unity and
their critical values are a common nonzero scalar times those roots.  Thus a
generic finite branch cycle is one `(r+1)`-cycle.  Irreducibility of `Psi`
makes the branch-cycle action transitive.  There are `m+1` finite branch
cycles and their defects sum to

\[
 (m+1)r=N-1.                                             \tag{14}
\]

The supports therefore form a connected uniform hypertree.  Here is the
elementary group lemma needed below.  If connected `k`-cycles have support
union of size `1+s(k-1)`, then the group they generate is `S` on that union
for even `k` and `A` for odd `k`, provided `s>=2`.  To prove it, order the
hyperedges so that each new edge meets the preceding union in one letter;
equality in the union-size bound forces such an order.  For the first two
edges, relabel their cycles
`a=(0,1,...,k-1)` and `b=(0,k,...,2k-2)`.  The commutator
`a^(-1)b^(-1)ab` is the 3-cycle `(0,k,1)`.  Conjugating this commutator by
powers of `a` and `b`, and multiplying two conjugates with one common pair,
gives every 3-cycle on the union.  Thus the generated group contains the
alternating group; adjoining each later edge repeats the same argument.
When `k=2`, the edges are a connected ordinary tree and generate the full
symmetric group.  Finally, a `k`-cycle is odd exactly when `k` is even, which
proves the lemma.

Apply the lemma with `k=r+1` and `s=m+1>=2`.  If `r+1` is even, a branch
cycle is odd and the group is `S_N`; if `r+1` is odd, all finite branch cycles
and the `N`-cycle at infinity are even and the group is `A_N`.  Therefore

\[
 \operatorname{Mon}(C24_{m,r})=
 \begin{cases}
 S_N,&r\text{ odd},\\
 A_N,&r\text{ even}.
 \end{cases}                                             \tag{15}
\]

The integral divisor (12) maps generically one-to-one to the discriminant,
because generic critical values are distinct.  Its boundary data are
`(e,f)=(r+1,1)` and sheet loss `r+1`.

### 4.2 The projective resolvent chart over `P=0`

Set `T=U/P`.  Multiplication of (10) by `P^(r+1)` gives the regular chart

\[
 C J(U,P)-RP^{r+1}=0,
\quad
 J(U,P)=\int_0^U\{P-V(Q-V)^m\}^r\,dV.                    \tag{16}
\]

At `P=0`, after `U=Qw`, its nonzero roots are the roots of

\[
 K_{m,r}(w)=w^{-(r+1)}\int_0^w v^r(1-v)^{mr}\,dv,        \tag{17}
\]

a degree-`mr` polynomial.  It is squarefree: a common zero of `K` and `K'`
would be a nonzero zero of both the integral and its integrand; the only
candidate is `w=1`, where the beta integral is nonzero in characteristic
zero.

Let `q=h(0)` be the selected root of the cancellation polynomial
`M_(m,r)`.  Exact substitution gives

\[
 K_{m,r}\left(-{q\over1-q}\right)
 =(1-q)^{-mr}\int_0^1u^r(1-q+qu)^{mr}\,du=0.             \tag{18}
\]

This root is not boundary: it is the affine divisor `A=0`, since there
`U=Ps=Bx=-qy` and `Q=(1-q)y`.  The restriction `A=0 -> P=0` is generically
degree one.  Indeed `q!=0,1`, so `Q` determines `(x,y)`; étaleness of the
three-dimensional map and `P=AB` with `B!=0` then makes the remaining
`z -> R` map linear with nonzero coefficient.  The `U=0` cluster is the
affine divisor `B=0` and contributes `r+1` sheets.  The other `mr-1` simple
roots of (17) are genuine boundary primes.  Consequently the complete
geometric table is

| target prime | geometric boundary primes | `(e,f)` | sheet loss |
|---|---:|---|---:|
| discriminant `Delta_(m,r)` | 1 | `(r+1,1)` | `r+1` |
| `P=0` | `mr-1` | `(1,1)` each | `mr-1` |

Over `K_q=Q(q)`, the arithmetic residue degrees over `P=0` are the degrees of
the irreducible factors of `K_(m,r)` after removing the distinguished linear
factor (18).  Equivalently, they are the orbit sizes of the other roots of
`M_(m,r)` under the stabilizer of `q`, after the Möbius change
`w=-q/(1-q)`.

The raw pole `s=x/A` sees all `mr` nonzero roots in (17); the intrinsic table
sees only `mr-1`.  This is the promised adversarial chart check.

## 5. Four exact small cases

All residue degrees in the geometric column are one.

| `(m,r)` | `N` | monodromy | critical partition | `Delta`: `(e;loss)` | `P=0`: geometric primes/loss | arithmetic `f` over `K_q` |
|---|---:|---|---|---|---|---|
| `(1,1)` | 3 | `S_3` | `1^2` | `2;2` | `0/0` | none |
| `(1,2)` | 5 | `A_5` | `2^2` | `3;3` | `1/1` | `1` |
| `(2,1)` | 4 | `S_4` | `1^3` | `2;2` | `1/1` | `1` |
| `(2,2)` | 7 | `A_7` | `2^3` | `3;3` | `3/3` | `1,2` |

The arithmetic `1,2` in the last row reflects the `D_4` splitting-field
action of the quartic parameter polynomial: after fixing one root, the other
three split into orbits of sizes one and two.

Run

```bash
.venv/bin/python scripts/verify_resolvent_ramification_signature.py
```

It checks (11), degrees, squarefreeness of (17), identity (18), the displayed
arithmetic factor degrees, and the four finite permutation groups.  These are
exact regressions for the four rows, not proofs of the all-`m,r` arguments.

## 6. Equivalence consequences

### Theorem 6.1

For `r>=2`, `C24_(m,r)` is not polynomially left--right equivalent, nor
stably polynomially left--right equivalent, to a generic weighted-seed map.

**Proof.**  Degrees first force the weighted inverse degree to be
`n=N`.  The C24 discriminant boundary prime has ramification index `r+1>=3`
and generic sheet loss `r+1`; the generic weighted discriminant prime has
index two and loss two.  Theorems 2.1--2.4 preserve this boundary valuation
data.  Independently, when `r` is even, (15) gives `A_N` rather than the
weighted `S_N`.  Therefore no such equivalence exists.  QED

This proves the conclusion asserted informally in
`MASTER_CANCELLATION_CONSTRUCTION.md`, but repairs its justification: raw
resolvent indices are not enough; the proof requires their identification
with the canonical boundary primes.

### The `r=1` classification problem

Let `r=1`, `m>1`, so `N=m+2`.  The presently computed intrinsic data agree
with a generic degree-`N` weighted seed:

- degree `N` and monodromy `S_N`;
- critical partition `1^(N-1)`;
- one discriminant prime with `(e,f)=(2,1)` and loss two;
- `N-3=m-1` unramified geometric boundary primes over one additional target
  divisor, with total loss `m-1`.

Thus the resolvent--ramification signature developed here gives **no
obstruction** for `C24_(m,1)`.  The case `m=1` is explicitly C01 after the
linear changes already recorded in the master construction.  For `m>1`, no
polynomial left--right transformation is known, and no stronger invariant in
this repository currently separates the maps.

A precise next condition is to decide whether there are polynomial target
coordinates identifying the two nonproperness components together with an
isomorphism of the finite normal pairs

\[
 (\bar X_{C24},X)\cong(\bar X_{G_H},A^3)                 \tag{19}
\]

whose restriction to the affine opens is induced by a polynomial source
automorphism.  Equality of the tables above is necessary but not sufficient:
it omits conductor ideals, intersections of boundary primes over
`Delta cap V(P)`, and the filtered valuation semigroups
`v_E(k[X]\setminus0)`.  Computing any one of those structures is a concrete
route to an obstruction; constructing (19) with polynomial restrictions is a
route to an equivalence.

## 7. Concise comparison

Here “generic C04” means (8), and `n=N` when it is compared with C24.

| family | degree | monodromy | critical partition | discriminant boundary | second boundary component |
|---|---:|---|---|---|---|
| C01 | 3 | `S_3` | `1^2` | one prime, `e=2`, loss 2 | none |
| generic C04, degree `n` | `n` | `S_n` | `1^(n-1)` | one prime, `e=2`, loss 2 | `n-3` geometric `e=1` primes; loss `n-3` |
| C24 `(m,r)` | `r(m+1)+1` | `S_N` if `r` odd; `A_N` if `r` even | `r^(m+1)` | one prime, `e=r+1`, loss `r+1` | `mr-1` geometric `e=1` primes; loss `mr-1` |

The table proves inequivalence for `r>=2`.  For `r=1` every displayed entry
matches generic C04, leaving the polynomial equivalence question open.

## External inputs and hypotheses

The precise external inputs are the following standard results.

1. A scheme of finite type over a field is excellent (hence Nagata), so its
   normalization in a finite function-field extension is finite.  We apply
   this only to the normal affine space `Y` and the finite separable extension
   `L/K`.
2. Zariski's Main Theorem in the form: a separated quasi-finite birational
   morphism to a normal noetherian integral scheme is an open immersion.  It
   is applied to `j_F`, which is quasi-finite because its fibers lie in the
   fibers of the Keller map.
3. The fundamental equality for the integral closure of a DVR in a finite
   separable extension, `sum e_i f_i=[L:K]`; characteristic zero rules out
   inseparability and defect in this divisorial situation.
4. A finite separable extension and its Galois closure keep their degree and
   geometric permutation group after adjoining algebraically independent
   variables.  Normality is preserved by polynomial extension.
5. For a finite cover of `P^1_C`, a small loop around a tame local equation
   `z=u^e` gives an `e`-cycle, the finite branch cycles generate the connected
   monodromy group, and their product is inverse to the cycle at infinity.
   This is the standard branch-cycle consequence of the Riemann existence
   theorem; the comparison reference used elsewhere in this repository is
   SGA 1, Expose XII, Theorem 5.1 and Corollary 5.2.

The minimal-transitive-factorization/hypertree lemma itself is proved in
Section 4.1 and is not an external classification theorem.  For an arbitrary
characteristic-zero base, the coefficients and finite cover descend to a
finitely generated subfield and the geometric group is computed after an
embedding in `C`.  No bounded symbolic calculation is used for an all-degree
assertion.
