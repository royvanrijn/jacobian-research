# Arithmetic density of omitted-partition strata

This note combines the
[contact-strata theorem](COINCIDENT_ROOT_REBUILD.md), the
[exceptional-component normalizations](COMPONENT_NORMALIZATION.md), and the
[finite-field symmetric-monodromy theorem](FINITE_FIELD_CHEBOTAREV.md).
For each fixed inverse degree it turns the geometric omitted-value phase
diagram into a finite-field rarity theorem.  It also gives a joint version in
which the omitted partition of the seed and the Frobenius cycle type of a
squarefree target are simultaneously prescribed.

The two partitions in the joint statement have different meanings.  The
omitted partition records the multiple roots of the unique omitted member of
the pencil.  A Frobenius partition records the residue degrees of an ordinary
squarefree member of the same pencil.

## Merged arithmetic phase-diagram theorem

> **Theorem.**  Fix `N>=4`.  There is a nonzero integer `B_N` such that, for
> every finite field `F_q` of characteristic not dividing `B_N`, the following
> statements hold.
>
> 1. A positive-dimensional full-contact stratum of dimension `d` has
>    `c(q)q^d+O_N(q^(d-1/2))` rational seeds, where `c(q)` is the number of
>    top-dimensional geometric components fixed by Frobenius.  In particular,
>    every geometrically integral atomic component has leading constant one.
> 2. The full nonsurjective seed locus has
>    `q^(floor(N/2)-1)+O_N(q^(floor(N/2)-3/2))` rational seeds.  For `N>=5`,
>    its density in the normalized admissible seed space is
>    `q^(-floor((N-3)/2))(1+O_N(q^(-1/2)))`.
> 3. On a geometrically integral stratum normalization, prescribing a
>    squarefree Frobenius cycle type `mu` multiplies the seed--target main term
>    by `|C_mu|/N!`.
> 4. At a fixed nonzero bridge scale, the same counts hold in the reduced
>    mixed-moment coordinates.

Parts 1 and 2 merge the contact-strata, component-normalization, enumerative,
and phase-diagram theorems with Lang--Weil.  Part 3 merges the universal
`S_N` theorem with the twist argument from finite-field Chebotarev, now over
the function field of a seed stratum.  Part 4 is the arithmetic base change
of the Gaussian invariant-transport theorem.  The remaining sections prove
the four parts and state their exact scope.

## 1. Integral models and the counting convention

Fix `N>=4`.  The ambient normalized affine slice is

\[
 \mathcal S_N=
 \{H:H(0)=H'(0)=H(1)=0,\ H'(1)=-1,\ \deg H\le N\}
 \simeq\mathbb A^{N-3}.
\]

Let `A_N` be the open on which the degree is exactly `N` and the weighted
admissibility conditions hold.  For a full-contact partition

\[
 \lambda=(\lambda_1,\ldots,\lambda_\ell),\qquad
 \lambda_i\ge2,\qquad \sum_i\lambda_i=N,
\]

write `E_lambda` for the exact omitted-value stratum and put

\[
 d_\lambda=\ell(\lambda)-1.
\]

All equations involved have integer coefficients after clearing finitely
many denominators.  Generic flatness, openness of geometric integrality, and
the characteristic-zero results therefore give an integer `B_N!=0` such
that the seed slice, the retained contact models, the component
normalizations, and the classification by a unique omitted value all have
good reduction in every characteristic not dividing `B_N`.  We enlarge
`B_N` by `N!`; in particular every retained characteristic is greater than
`N` after adding the finitely many smaller primes.

The integer `B_N` is degreewise, not seedwise.  It is distinct from the
seed-specific bad ideal in the effective Chebotarev theorem.  Its existence
is enough for the asymptotic statements below.  The displayed equations,
Jacobian saturations, and projective degrees in the geometric notes make it
effective in principle.

All point counts mean counts of geometric seed parameters defined over
`F_q`.  Thus "nonsurjective" means that the base change of the weighted map
to an algebraic closure has an omitted value.  This is not the different
question whether the function on the finite set `F_q^3` is surjective.  In
fact, the [finite-field detection-limit theorem](FINITE_FIELD_DETECTION_LIMITS.md)
shows that every good weighted marked-root reduction is nonbijective on
`F_q^3`, regardless of whether its seed lies in `N_N`.

### Descent of the omitted value

There is one arithmetic point that is not supplied by dimension alone.  Let
`H in A_N(F_q)` be geometrically exceptional.  Over `bar(F_q)` it has a
unique omitted pencil value `(s,t)`.  Frobenius sends `(s,t)` to another
omitted value of the same seed, so uniqueness forces

\[
 (s^q,t^q)=(s,t).
\]

Hence `(s,t) in F_q^2`.  Moreover, for every multiplicity `m`, the set of
roots having multiplicity `m` is Frobenius-stable.  The monic polynomial
`Q_m` supported on that set therefore belongs to `F_q[W]`.  Thus an
`F_q`-rational exceptional seed has an `F_q`-rational point in its exact
quotient root model, and conversely.  This identifies rational seed points
with rational exact-stratum points; it rules out a hidden contribution from
omitted values defined only over an extension field.

The use of uniqueness here is valid in every retained characteristic, not
merely by specialization from characteristic zero.  Since `p>N`, every
nonconstant polynomial of degree at most `N` has nonzero derivative;
polynomial Mason--Stothers therefore applies to the two degree-`N` pencil
members exactly as in the characteristic-zero proof.  The remaining
all-double case uses only unique factorization and `p!=2`.  Hence the unique
omitted-value argument repeats verbatim in these characteristics.

## 2. Lang--Weil on a single stratum

Let `Z` be a geometrically irreducible component of the reduced retained
model of `E_lambda`, and let `tilde Z -> Z` be its normalization.  Suppose
the component and its normalization are defined over `F_q`.  For every good
`q`,

\[
 \boxed{
 \#Z(\mathbb F_q)
 =q^{d_\lambda}+O_{N,\lambda}(q^{d_\lambda-1/2}).}
 \tag{1}
\]

Indeed, `tilde Z` is geometrically integral of dimension `d_lambda`, so
Lang--Weil applies to it.  A finite birational normalization differs from
its image only above a closed subset of dimension at most `d_lambda-1`.
Deleting the root diagonals and the inadmissible divisors also changes the
count by only `O(q^(d_lambda-1))`.  Both changes are absorbed by the error in
(1).

For a reducible stratum the correct leading coefficient is

\[
 c_\lambda(q)=
 \#\{\text{top-dimensional geometric components fixed by Frobenius}\}.
 \tag{2}
\]

Consequently

\[
 \boxed{
 \#\mathcal E_\lambda(\mathbb F_q)
 =c_\lambda(q)q^{d_\lambda}
  +O_{N,\lambda}(q^{d_\lambda-1/2}).}
 \tag{3}
\]

Components in a nontrivial Frobenius orbit contribute rational points only
through their lower-dimensional intersections.  Thus `c_lambda(q)` can
depend on the Frobenius class; over extensions of one fixed finite base field
it is periodic in the extension degree.  On a geometrically irreducible
normalization defined over `F_q`, it is exactly one.  A labelled-root model
has instead the expected finite-cover factor, such as `a!b!`; that factor
does not occur when seeds themselves are counted.

When `d_lambda=0`, (3) must be replaced by the exact finite Frobenius fixed-
point count.  An error of order `q^(-1/2)` cannot describe a varying integer.
This qualification affects maximally collided zero-dimensional strata, but
not the top-dimensional density theorem for `N>=5`.

## 3. Counts for every exceptional component

For an atomic partition `2^a3^b`, with `2a+3b=N`, the geometric results are
stronger: the stratum is irreducible and its component has the explicit
smooth geometrically integral normalization

\[
 \widetilde{\mathcal C}_{a,b}
 =\{\Phi(Q^2R^3)=0\}\cap D(DA).
\]

It has dimension `a+b-1`, and the normalization map has generic degree one.
After enlarging `B_N` once more if necessary, (1) therefore gives

\[
 \boxed{
 \#\mathcal C_{a,b}(\mathbb F_q)
 =q^{a+b-1}+O_N(q^{a+b-3/2}).}
 \tag{4}
\]

The same estimate holds for the exact atomic stratum
`E_(2^a3^b)`, because its complement in the component has smaller
dimension.  The exact projective degrees and graph Chow classes computed in
the enumerative theorem bound the complexity of the varieties in (4); they
can be inserted into an explicit Lang--Weil bound if a numerical constant is
required.

## 4. Quantitative rarity of nonsurjective seeds

Set

\[
 D_N=\left\lfloor\frac N2\right\rfloor-1,
 \qquad
 m_N=N-3,
 \qquad
 r_N=m_N-D_N=\left\lfloor\frac{N-3}{2}\right\rfloor.
 \tag{5}
\]

There is a unique component of dimension `D_N`: it is indexed by

\[
 (a,b)=
 \begin{cases}
 (N/2,0),&N\text{ even},\\
 ((N-3)/2,1),&N\text{ odd}.
 \end{cases}
\]

Every other component has dimension at most `D_N-1`, and every intersection
is lower-dimensional.  The disjoint omitted-partition stratification and
(4) now imply

\[
 \boxed{
 \#\mathcal N_N(\mathbb F_q)
 =q^{D_N}+O_N(q^{D_N-1/2}).}
 \tag{6}
\]

Since `A_N` is a nonempty open of affine `m_N`-space,

\[
 \#\mathcal A_N(\mathbb F_q)
 =q^{m_N}+O_N(q^{m_N-1}).                              \tag{7}
\]

Dividing (6) by (7) gives the promised phase-diagram density:

\[
 \boxed{
 \frac{\#\mathcal N_N(\mathbb F_q)}
      {\#\mathcal A_N(\mathbb F_q)}
 =q^{-\lfloor(N-3)/2\rfloor}
  \bigl(1+O_N(q^{-1/2})\bigr),\qquad N\ge5.}
 \tag{8}
\]

Thus the leading constant is one in the normalized seed space.  More
general constants arise only from taking unions with several
top-dimensional arithmetic components, changing the quotient convention,
or counting labelled factorizations.

For `N=4`, the same calculation has `r_N=0`: the exceptional component is
dense, agreeing with the geometric transition rather than giving rarity.

## 5. Joint omitted-partition and Frobenius statistics

Let `Z` be as in Section 2 and let `mu` be a partition of `N`.  Write

\[
 z_\mu=\prod_i i^{m_i}m_i!,
\]

where `m_i` is the number of parts of `mu` equal to `i`.  Over
`Z x A^2_(s,t)`, remove the discriminant of

\[
 H(W)-sW+t.
\]

For `x in Z`, write `H_x` for its image in the normalized seed space.  Replace
`Z` here by its geometrically integral normalization, which changes the
eventual seed count only over a lower-dimensional locus.  The open

\[
 B=Z\times\mathbb A^2_{s,t}\setminus V(\Delta)
\]

is then normal and geometrically integral.  The universal monodromy theorem,
applied over `\overline{\mathbb F_q(Z)}(s,t)`, gives `S_N`.  Consequently the generic fiber
of the ordered-root `S_N`-torsor over `B` is geometrically connected.  A
finite etale cover of a normal integral base is the normalization in its
generic etale algebra, so the total ordered-root cover is geometrically
integral.  Every constant twist becomes isomorphic to it over `bar(F_q)` and
is geometrically integral as well.

This is the extra relative-family step beyond the seedwise Chebotarev
theorem.  After spreading out and enlarging `B_N`, Lang--Weil on each twist
in dimension `d_lambda+2`, followed by the exact twist point-count identity,
gives

\[
 \boxed{
 \begin{aligned}
 &\#\{(x,s,t):x\in Z(\mathbb F_q),\
       H_x(W)-sW+t\text{ squarefree of Frobenius type }\mu\}\\
 &\hspace{35mm}
 =\frac1{z_\mu}q^{d_\lambda+2}
  +O_{N,\lambda}(q^{d_\lambda+3/2}).
 \end{aligned}}
 \tag{9}
\]

Here `1/z_mu=|C_mu|/N!`.  The discriminant boundary has dimension at most
`d_lambda+1`, so it does not alter the main term.

Applying (9) to the unique top component and absorbing all lower-dimensional
components yields

\[
 \boxed{
 \#\{(H,s,t):H\in\mathcal N_N(\mathbb F_q),\
       \operatorname{Frob}_{H,s,t}\text{ has type }\mu\}
 =\frac1{z_\mu}q^{D_N+2}
  +O_N(q^{D_N+3/2}).}
 \tag{10}
\]

Relative to all admissible seed--target pairs, (10) has density

\[
 \boxed{
 \frac1{z_\mu}q^{-\lfloor(N-3)/2\rfloor}
 \bigl(1+O_N(q^{-1/2})\bigr).}
 \tag{11}
\]

Equation (11) is the precise sense in which exceptional-seed rarity and the
universal symmetric-group cycle law coexist.  The unique omitted target
itself lies on the discriminant and is not one of the squarefree targets
counted in (9)--(11).

For weighted Keller targets with `C!=0`, each pencil target `(s,t)` has
exactly `q-1` lifts

\[
 (A,B,C)=\left(\frac{t}{cC^2},\frac{s}{C},C\right).
\]

Multiplying (9) or (10) by `q-1` gives the corresponding three-dimensional
target count.

## 6. Gaussian moment coordinates

Fix a bridge scale `eta!=0`.  In characteristic greater than `N`, factorial
normalization is invertible and the mixed-moment coordinate map has
triangular Jacobian `eta^(N-3)`.  Away from the finitely many primes already
placed in `B_N`, it is a polynomial automorphism

\[
 \Psi_{N,\eta}:\mathcal S_N\xrightarrow{\sim}\mathbb A^{N-3}
\]

over `F_q`.  It therefore preserves `F_q`-point counts exactly.  Equations
(3), (4), (6), and (8) hold unchanged for the corresponding exceptional
partition strata in the mixed-moment coordinates
`(M_3,...,M_(N-1))`.

Here the finite-field coordinates are the reductions of the algebraic
characteristic-zero moment polynomials.  No probability measure or Gaussian
expectation over `F_q` is being introduced.

This arithmetic invariance is affine.  As in the geometric moment theorem,
it makes no assertion that the nonlinear coordinate change preserves the
degrees of chosen projective closures.

## 7. What the existing geometry supplies

The proof uses no new component classification.  Its inputs are precisely:

1. the disjoint classification by a unique omitted partition;
2. `dim E_lambda=ell(lambda)-1` and the exact closure order;
3. geometric integrality and degree-one normalization for every atomic
   component;
4. the unique top component and the lower-dimensional intersection lattice;
5. full `S_N` monodromy for the universal pencil; and
6. the polynomial moment-coordinate isomorphism.

The only arithmetic operation is spreading these fixed-degree varieties and
covers away from finitely many primes and applying Lang--Weil, with the twist
identity for the joint Chebotarev statement.  The result is therefore a
direct arithmetic consequence of the proved geometry.

## 8. Remaining gaps and effectivity upgrades

The qualitative rarity theorem (6)--(8) has no remaining geometric gap once
the descent and relative-cover arguments above are included.  The following
stronger refinements are not yet supplied by the existing theorems.

1. **An explicit degreewise bad integer.**  The proof constructs `B_N` by
   spreading out finitely many fixed-degree schemes and morphisms.  It does
   not yet compute the product of discriminants, resultants, Jacobian minors,
   and generic-degree certificates that would give a numerical `B_N`.
2. **Constants for every non-atomic stratum.**  Formula (3) is complete in
   terms of Frobenius on geometric components, but the repository proves
   geometric irreducibility uniformly only for the atomic `2^a3^b` strata.
   Determining `c_lambda(q)` explicitly for every collision partition requires
   decomposing the corresponding retained `Phi_lambda` divisor.
3. **Explicit relative Chebotarev errors.**  The equations and projective
   degrees bound the family complexity, so an effective Lang--Weil constant
   can be written.  The current note records only the `O_(N,lambda)` form and
   does not carry out that numerical estimate for every twist.
4. **Zero-dimensional strata.**  Their counts are exact Frobenius fixed-point
   counts rather than a power-saving asymptotic.  One must compute their
   finite etale algebras to obtain a closed residue-class formula.

These are effectivity and stratum-resolution problems.  None changes the
leading constant one or the rarity exponent for the full nonsurjective locus,
because that locus has one geometrically integral top-dimensional component.
