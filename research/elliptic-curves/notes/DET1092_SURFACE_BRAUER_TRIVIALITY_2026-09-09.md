# No nonconstant Brauer classes on the determinant-1092 parent over Q

## The theorem and its arithmetic scope

**New deduction, with independently verified arithmetic inputs.** For the
completed determinant-1092 K3 parent `X/Q`, pullback identifies

\[
 \boxed{\operatorname{Br}(X)=\operatorname{Br}(\mathbf Q).}
\]

This is a statement about **all torsion orders**, not only the2-primary
part. Equivalently, normalize a surface Brauer class by subtracting its
evaluation at a rational point on the zero section: the normalized class
is zero. There is no nonconstant surface Brauer class over Q whose values
could distinguish302 from the controls.

The geometric Brauer group is not zero: the already proved geometric
Picard rank19 gives `Br(Xbar)=(Q/Z)^3`. We do not determine its subgroup of
Galois invariants here. A geometric invariant, a class descending to
`Br(X)`, a fibre Kummer class, and a fibre Selmer/Sha class must not be
identified with each other.

This strengthens the earlier
[section-tree and local blindness theorem](DET1092_BRAUER_SECTION_AND_LOCAL_GATE_2026-09-08.md).
That theorem left transcendental surface classes at other primes open.
The present theorem excludes every nonconstant class on this surface over Q,
without constructing a Brauer algebra or enlarging the local-place panel.

## 1. Existing evidence, not a new point-count campaign

**Verified applications already completed.** The
[full parent geometry theorem](CURVE302_PARENT_SEARCH_AND_GEOMETRY_2026-09-07.md)
provides:

- a full integral geometric Picard lattice `L=U+(-M17)`, determinant1092,
  with all19 generators represented by divisors over Q;
- smooth proper K3 reduction at the already fixed primes149 and151;
- full Frobenius polynomials obtained from exact counts over `Fp` and `Fp^2`.

The full characteristic polynomials are

\[
 R_p(T)=(T-p)^{19}(T+p)(T^2-a_pT+p^2),
 \qquad a_{149}=-248,\quad a_{151}=-244.
\]

| Prime | Count over `Fp` | Previously certified count over `Fp^2` |
|---|---:|---:|
|149|24636|493345524|
|151|25276|520355556|

The new constructor independently recounts the base-field points by direct
quadratic-character summation. It deliberately **reuses** the immutable
extension-field counts; those are not represented as new independent counts.
Good reduction is rechecked from the degree24 squarefree discriminant,
coprimality with the short `a4`, and smooth infinity.

## 2. The arithmetic reduction lattices are already complete

**New verified application of the existing lattice obstruction.** At either
prime the eigenvalue `p` has multiplicity19. Tate's theorem therefore gives
arithmetic Picard rank19. The geometric Picard rank is20, but that is not
the lattice to use in Artin--Tate over `Fp`.

The inherited rank19 lattice embeds with some finite index `e` in the
arithmetic reduction lattice. Since its determinant is

\[
1092=2^2\cdot3\cdot7\cdot13,
\]

integrality forces `e=1` or2. Its discriminant group is cyclic of order1092.
The unique possible order-two extension has a half-word of square `-45`.
Every divisor on a K3 has even self-intersection, so this extension cannot
be a Picard lattice. Consequently `e=1` at both primes, and

\[
 \operatorname{Pic}(X)\longrightarrow\operatorname{Pic}(X_p)
 \quad\text{is an isomorphism of integral lattices},
 \qquad |\det\operatorname{Pic}(X_p)|=1092.
\]

Over a finite field the invariant geometric divisor classes descend;
`Br(Fp)=0` removes the descent obstruction. The already rational generators
of `Pic(X)` thus account for the entire arithmetic Picard group, not only
its rational span or its rank.

## 3. Artin--Tate gives two trivial finite-field Brauer groups

**Established literature, verified application.** In the monic
characteristic-polynomial convention, Artin--Tate for this K3 gives

\[
 \left.\frac{R_p(T)}{(T-p)^{19}}\right|_{T=p}
  =p^{22-19-1}\,\#\operatorname{Br}(X_p)\,
                  |\det\operatorname{Pic}(X_p)|.
\]

Tate and Artin--Tate are established for these K3 reductions; see the
[Artin--Tate normalization in Kedlaya's notes, section18.2](https://kskedlaya.org/weil-cohom/chapter-18.html)
and the Tate references in the completed parent proof.

**New exact calculation.** Substitution reduces the order to

\[
 \#\operatorname{Br}(X_p)=\frac{4p-2a_p}{1092}=1
 \qquad(p=149,151).
\]

The numerator is1092 in both cases. The independent checker also uses
the reciprocal convention `P_p(u)=det(1-u Frob_p)` and verifies

\[
 \frac{p}{1092}
 \left.\frac{P_p(u)}{(1-pu)^{19}}\right|_{u=1/p}=1.
\]

This guards against confusing the arithmetic rank19 formula over `Fp`
with the geometric rank20 discriminant-squareclass calculation over `Fp^2`.

## 4. Why two good reductions rule out global nonconstant classes

**New deduction from established cohomology.** The following argument
also gives a reusable criterion. Let `ell` be any prime and choose
`p in {149,151}` with `p != ell`.

The finite-field Kummer sequence and the preceding calculation give

\[
 H^2(X_p,\mu_\ell)=\operatorname{Pic}(X_p)/\ell.
\tag{1}
\]

A geometric K3 has `H^1(Xbar_p,mu_ell)=0`. The finite field has
cohomological dimension one for torsion prime to `p`. Hence
Hochschild--Serre identifies

\[
 H^2(X_p,\mu_\ell)
 \simeq H^2(X_{\overline{\mathbf F}_p},\mu_\ell)^{\mathrm{Frob}_p}.
\tag{2}
\]

Smooth proper base change identifies the geometric cohomology here with
`H^2(Xbar,mu_ell)`, compatibly with specialized divisor classes. Since all
of `Pic(X_p)` comes from `Pic(X)`, equations(1)--(2) imply:

> Every Frobenius-fixed geometric Kummer lift at this prime is the
> reduction of a global divisor class modulo `ell`.

Now take `alpha in Br(X)[ell]` and choose a Kummer lift
`b in H^2(X,mu_ell)`. Its geometric restriction is Galois-invariant, hence
Frobenius-fixed at this prime. It therefore equals the class of a divisor
on `X` modulo `ell`. Subtract that divisor class from `b`. The geometric
restriction of the Brauer image is now zero, so `alpha in Br1(X)`.

The full geometric Picard group is the trivial Galois module `Z^19`.
Thus `H^1(Q,Pic(Xbar))=0`, and the algebraic Brauer exact sequence gives
`Br1(X)=Br(Q)`. This proves that every global order-`ell` class is constant.
The ingredients are the Kummer sequence, smooth proper base change and
Hochschild--Serre; see
[Colliot-Thelene--Skorobogatov, sections2.2,3.1 and4.3](https://www.imo.universite-paris-saclay.fr/~jean-louis.colliot-thelene/BGgroup_book.pdf).

Finally evaluation at a rational point splits off `Br(Q)`. If the
normalized subgroup were nonzero, a nonzero finite-order element would
have a nonzero multiple of prime order. The preceding argument excludes
every such prime. Hence the normalized subgroup is zero, proving the theorem.

One prime would leave its residue-characteristic primary part untreated.
The two existing distinct primes suffice for every prime; no additional
prime or cohomology computation is required.

## 5. Consequences for the seed problem, and what is not proved

**New obstruction.** A proposed incidence mechanism that requires a
nonconstant class in `Br(X)` over Q cannot work on this parent, irrespective
of its order, presentation, degree, or chosen elliptic fibration. In
particular there is no remaining transcendental surface Brauer class to
evaluate at extra primes in search of a302/control distinction.

**Verified distinctions, not new rank claims.** The first302 seed still has
its certified nongeneric elliptic Kummer class. Its specialized cubic-field
ideal class remains subject to the previously identified norm-equation fork.
Neither object is a nonconstant class on the original surface. Likewise,
the conic rank18 construction uses a pulled-back surface, and the nine
blinded M16 recoveries change the reference subgroup. Nothing here changes
their certified outcomes or forces specialized elliptic Sha groups to vanish.
Ramified classes on open subsets, classes on new base changes, and new
rational points on individual fibres are not excluded.

**Audited literature boundary.** The tempting identification of the
trigonal2-division curve's Jacobian2-torsion with surface Brauer2-torsion
in [van Geemen, Theorem7.6](https://sites.unimi.it/vangeemen/0408006.pdf)
assumes Picard rank two. The more general Theorem6.2 assumes that the
covering involution acts trivially on Picard. Original elliptic inversion
negates our17-dimensional MW summand, so that hypothesis also fails.
Neither theorem is applied to this rank19 parent. No genus10 Jacobian
torsion group or generic Selmer dimension has been computed or inferred.

**Open constructive goal.** This closes a broad arithmetic route, not the
seed-incidence problem itself. A prospective construction through302 still
needs an actual rational class outside generic MW17, with its original
elliptic independence checked. No such new construction is claimed here.

## 6. Checkpoints and reproduction

**Verified computation.** Both new scripts finish in under one second in
the observed runs, within25-second caps. The initial independent checker
made a `Factorization`-versus-list comparison error before any reduction
conclusion; its exact line, source hash and reconstructible correction are
retained. No mathematical input, prime, count or limit changed.

The [packet](../../artifacts/generated-results/elliptic-curves/det1092_surface_brauer_triviality_v1/)
binds the old point-count certificate, the parent, both scripts, the
arithmetic-lattice obstruction and both Artin--Tate orders. The proof is
hybrid: exact arithmetic certificates plus the stated established
cohomological theorems, not a formal proof assistant verification.

```sh
timeout 25s sage -python research/elliptic-curves/cas/verify_det1092_surface_brauer_triviality.sage
```

No exceptional point, later V3 artifact, parameter sweep, point search,
class-group calculation, production mutation or detached job was used.
