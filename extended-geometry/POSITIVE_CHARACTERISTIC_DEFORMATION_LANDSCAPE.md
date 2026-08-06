# Positive-characteristic deformation landscape

> **Status.** This is a research programme and an initial exact ledger, not a
> claim that every characteristic-zero construction has been spread out over
> `Spec Z`.  The proved inputs are the characteristic-\(p\) SIC theorem
> [`SIC2C4P`](TWO_PAIR_SIC_CHARACTERISTIC_P.md), the \(p\)-typical recovery
> and Frobenius-enriched boundary theorem
> [`HT1`](HASSE_TYPICAL_SEED_RECOVERY.md), and the characteristic-zero
> symmetric-monodromy theorem
> [`UM1`](../verified/UNIVERSAL_SYMMETRIC_MONODROMY.md).  The two elementary
> consequences proved below are the first calibration results for this
> research direction.

The organizing principle is:

\[
 \boxed{\text{reduction modulo }p\text{ is a deformation of the
 landscape, not merely a cheaper coefficient field}.}
\]

Three operations which coincide too easily in characteristic zero separate
in characteristic \(p\):

1. reducing an integral ordinary-derivative construction;
2. replacing ordinary derivatives by Hasse derivatives; and
3. retaining the integral or Witt lift and only then reducing its
   obstruction values.

The repository already contains examples where these three operations give,
respectively, a vacuous eventual-power theorem, a non-Mathieu Hasse kernel,
and a nontrivial \(p\)-adic valuation tower.

The
[characteristic-two plane theorem](../verified/HUQ_KURUVILLA_CHARACTERISTIC_TWO_AUDIT.md#8-the-preserved-coordinate-and-the-dimension-two-theorem)
is now the smallest-dimensional calibration in this landscape.  Mondello,
[arXiv:2608.02634v1](https://arxiv.org/abs/2608.02634), proves it over
(k=overline{mathbb F}_2): a polynomial coordinate change makes the
Huq--Kuruvilla threefold a skew product, and its preserved plane fibre is
etale of exact separable degree three.  The separately named repository
corollary `HKM2-ALLFIELDS` proves the same conclusion over every
characteristic-two field.  This shortens the characteristic-two
counterexample chain from dimension three to dimension two, but it does not
spread to odd characteristic: the identical integer plane formula has
nonconstant Jacobian.  Its ambient threefold normalization still supplies
the wild-radicial boundary model used below, while the plane hidden cubic
supplies the dimension-minimal noninvertibility theorem.

<!-- status-consumer: HKM2 6891c4426fa9c6ff -->

<!-- status-consumer: HKM2-ALLFIELDS bc2d0e1d0c37827c -->

The corresponding
[plane normalization theorem](../verified/HUQ_KURUVILLA_PLANE_BOUNDARY_NORMALIZATION.md)
now completes the first wild-ledger row: one retained sheet and one missing
Frobenius sheet over `Q=0`, with generic data
`(e,f_sep,f_insep,different,sheet loss)=(1,1,2,1,2)` and `S_3` monodromy.
The primitive cubic order has conductor `(P,T)`, but that conductor is not the
reconstruction boundary.  This supplies an exact warning that order
nonnormality and missing-source boundary must be tracked separately.  The
current Singular certificate plus separate identity replay is not an
independent second normalization implementation; that audit and two
independent monodromy routes remain open.

<!-- status-consumer: HKM2B1 c197aa4165670dc5 -->

## 1. Seven-axis survival ledger

| Axis | What survives without change | Characteristic-\(p\) deformation | Correct replacement or open datum |
|---|---|---|---|
| canonical boundary ledgers | finite normalization, valuations, residue extensions, Fitting ideals, and conductor ideals after a verified compatible base change | leading coefficients can vanish; horizontal primes can merge; vertical primes and radicial edges can appear; ordinary ramification indices no longer determine the different | spread out the complete finite normalization and record, prime by prime, horizontal, vertical, embedded, tame, wild, and radicial pieces |
| asymptotic conjectures | a fixed finite prefix survives at all sufficiently large primes | Frobenius can force every sufficiently high ordinary contraction to vanish; in ordinary SIC the single \(p\)-th pure moment already forces the conclusion | separate fixed-\(p\) eventual statements, uniform-in-\(p\) prefix statements, and Witt-level inverse systems |
| derivative data | \(D^{[1]}=d/dW\), affine covariance, and integral identities | ordinary derivatives kill \(k[W^p]\); ordinary Hessians can vanish or become inseparable; independently projectivized Hasse channels lose relative scale | use the common weighted tuple \(D^{[1]},D^{[p]},\ldots\), or the smaller intrinsic marked boundary edge when a map-level invariant is required |
| missing boundary primes | valuation degree sums remain the right exhaustion test | a root of multiplicity \(m\) gives a missing prime of index \(m\), wild when \(p\mid m\); in characteristic two the critical edge can be radicial of degree two | retain the extension, different, residue center, and conductor, not only the reduced boundary graph |
| stable moduli | stable normalization functoriality and the marked-edge reconstruction mechanism | the ordinary Hessian portrait can have positive-dimensional fibers and the derivative-only map can identify distinct seeds literally | Frobenius-enriched maps retain an \((N-3)\)-dimensional geometric family of stable classes in every characteristic; see Theorem 2.1 |
| moment and cotangent obstructions | exact integral matrices and complexes base-change derivedly | a scalar obstruction may vanish only because it is divisible by \(p\); ranks can drop; new Tor/Bockstein classes and new gauge directions can appear | record \(p\)-adic valuations for scalar moments and the derived reduction triangle for cotangent/Kuranishi complexes |
| arithmetic monodromy | good tame reductions can retain the characteristic-zero permutation group | Frobenius-additive terms can collapse \(S_N\) to an affine group; wild inertia and inseparable discriminant maps alter the branch generators | compute geometric monodromy, arithmetic Frobenius, tame inertia, and wild inertia separately; Theorem 3.1 is the first calibration |

The phrase “a boundary ledger survives” will mean more than preservation of
the same list of integers.  It means that a finite integral model has been
chosen and that normalization, the distinguished affine open, boundary
primes, target images, residue maps, and declared Fitting or conductor
modules commute with the base change under discussion.  Equality of a
degree sum after reduction is necessary but does not exclude a vertical or
embedded component.

## 2. Stable moduli survive, but the Hessian model does not

Fix an algebraically closed field \(k\) of characteristic \(p>0\).  A
normalized exact-double seed of degree at most \(N\) has the unique form

\[
 H(W)=W^2(1-W)Q(W),\qquad
 \deg Q\leq N-3,\qquad Q(1)=1.                    \tag{2.1}
\]

The exact-degree and exact-double conditions are the open conditions

\[
 [W^{N-3}]Q\ne0,\qquad Q(0)\ne0.                  \tag{2.2}
\]

### Theorem 2.1 — positive-characteristic stable-moduli lower bound

For every \(N\geq3\), the normalized exact-degree-\(N\) seed space contains a
nonempty open subscheme of dimension \(N-3\).  Applying the complementary
Frobenius-enriched Keller charts of `HT1` gives a family of determinant-one
polynomial maps \(\mathbb A^3_k\to\mathbb A^3_k\) of geometric degree \(N\)
whose geometric points are pairwise stably polynomially left--right
inequivalent whenever their normalized seeds differ.  Consequently the
family of stable classes has geometric dimension at least \(N-3\).

This is a lower bound on the dimension of the image in stable-class space.
It does not assert that a coarse moduli scheme for all Keller maps exists.
Over a finite ground field the open may have few or no rational points; the
statement is geometric, and rational points appear after finite extension.

#### Proof

The affine hyperplane \(Q(1)=1\) in
\(k[W]_{\leq N-3}\) has dimension \(N-3\).  For \(N=3\), it is the single
seed \(W^2(1-W)\).  For \(N\geq4\), neither \(Q(0)\) nor the leading
coefficient vanishes identically on the hyperplane, so their simultaneous
nonvanishing defines a nonempty open over the algebraic closure.

The two complementary source charts in
[`HASSE_TYPICAL_SEED_RECOVERY.md`](HASSE_TYPICAL_SEED_RECOVERY.md) cover
every normalized seed and produce a polynomial determinant-one map with
inverse pencil

\[
 H(W)-BCW+AC^2=0.                                  \tag{2.3}
\]

The full collision theorem in that note proves that the intrinsic marked
boundary edge

\[
 \left(\mathbb A^1_W;(W),(W-1),
 \operatorname{div}\frac{H}{W^2(W-1)}\right)       \tag{2.4}
\]

is preserved by stable left--right equivalence and reconstructs \(H\)
exactly.  Distinct geometric seed points therefore give distinct stable
classes. \(\square\)

The theorem illustrates the required change of invariant.  Ordinary
Hessian recovery is valid over the safe characteristic-zero model and over
large primes, but Frobenius directions can create positive-dimensional
fibers in small characteristic.  The stable modulus survives because the
full boundary edge sees the primitive-root divisor, including multiplicity
and radicial structure.

## 3. A first arithmetic-monodromy degeneration

The characteristic-zero theorem `UM1` gives \(S_N\) for every two-parameter
linear pencil \(H(W)-sW+t\).  Its proof uses ordinary critical points and
the absence of nontrivial finite étale covers of the affine line.  Both
inputs can fail in characteristic \(p\).

### Theorem 3.1 — the Frobenius pencil has affine monodromy

Let \(k\) be any field of characteristic \(p>0\), and set

\[
 E(W)=W^p-sW+t\in k(s,t)[W].                        \tag{3.1}
\]

Its geometric and arithmetic Galois groups in the natural degree-\(p\)
action are

\[
 \operatorname{AGL}_1(\mathbb F_p)
 =\mathbb F_p\rtimes\mathbb F_p^\times.             \tag{3.2}
\]

The order is \(p(p-1)\).  It equals \(S_p\) only for \(p=2,3\), and is a
proper solvable subgroup of \(S_p\) for every \(p\geq5\).

#### Proof

Adjoin \(\alpha\) with \(\alpha^{p-1}=s\), and put \(W=\alpha Y\).  Division
by \(\alpha^p\) turns (3.1) into the Artin--Schreier equation

\[
 Y^p-Y=-t/\alpha^p.                                 \tag{3.3}
\]

The Kummer extension \(k(s,t)(\alpha)/k(s,t)\) has degree \(p-1\).  Over it,
(3.3) has degree \(p\): the right side has a pole of order one at
\(t=\infty\), whereas a nonconstant Artin--Schreier coboundary \(Z^p-Z\)
has reduced pole order divisible by \(p\).  Its roots are \(Y+c\) for
\(c\in\mathbb F_p\), giving the translation subgroup.

Changing \(\alpha\) to \(a\alpha\), \(a\in\mathbb F_p^\times\), conjugates
translations by \(c\mapsto ac\).  The splitting field therefore has degree
\(p(p-1)\) and group (3.2).  All constants used by the action already lie
in the prime field, so extension of the constant field does not shrink or
enlarge the displayed generic group. \(\square\)

At \(s=0\), the polynomial becomes purely inseparable.  Thus the same
family contains both a separable affine-linear generic group and a radicial
boundary.  This is qualitatively different from reducing an \(S_p\)
branch-cycle list and deleting bad primes.

The next monodromy problem is to classify \(G_{\rm geom}\) for
\(H(W)-sW+t\) from the \(p\)-typical support of \(H\).  The first strata to
separate are:

1. \(H'\ne0\) with tame degree and a generically Morse tilt;
2. additive or linearized \(H\);
3. inseparable critical parameterizations with separable generic cover;
4. wild finite ramification; and
5. constant-field enlargement of the geometric group by arithmetic
   Frobenius.

## 4. Ordinary, Hasse, and Witt questions are different

The existing SIC calculation gives a sharp test case.

* **Ordinary reduction.** In characteristic \(p\),
  \(\mathcal E_{r,p}(f^p)=f(0,z)^p\).  The ordinary SIC conclusion follows
  from one pure moment and every multiplier vanishes from order \(p\)
  onward.  This is a complete theorem, but Frobenius makes the asymptotic
  question too weak to model characteristic zero.
* **Hasse replacement on the polynomial algebra.** The exact kernel uses
  the \(p\)-typical operators of orders \(1,p,p^2,\ldots\), but the kernel is
  not Mathieu already in one pair.  This is a different conjecture, not a
  repaired proof of ordinary SIC.
* **Fully divided-power multiplication.** Every positive-degree element is
  \(p\)-nilpotent, so eventual powers again become vacuous.
* **Witt lift.** Ordinary contraction over \(W(k)\) retains the integral
  mixed moments.  Reduction modulo \(p^a\) is governed by their exact
  valuations and can disappear and re-enter in nonradial families.

This research direction therefore attaches two quantifiers to every moment
problem:

\[
 \begin{array}{ll}
 \text{vertical:}&p\text{ fixed and }m\to\infty,\\
 \text{horizontal:}&m\leq M\text{ fixed and }p\to\infty.
 \end{array}                                         \tag{4.1}
\]

A fixed-\(p\) theorem cannot be used as evidence for a uniform horizontal
statement without an explicit comparison.

## 5. Wild missing primes and boundary survival

For the Frobenius-enriched weighted maps, `HT1` supplies the first complete
wild boundary calculation.  If \(\rho\ne0,1\) is a primitive root of
multiplicity \(m\), the local inverse equation at the second boundary has
the form

\[
 h_m\delta^m-B\rho C-BC\delta+AC^2
 O(\delta^{m+1})=0.                                 \tag{5.1}
\]

Its Weierstrass polynomial is Eisenstein at \(C\), so the missing boundary
prime has ramification index \(m\).  When \(p\mid m\), this is wild or
radicial; the index alone does not record the different.  In characteristic
two the critical boundary additionally carries the radicial edge

\[
 \mathbb A^1_W\longrightarrow\mathbb A^1_T,\qquad T=W^2. \tag{5.2}
\]

The next local ledger must add, for every missing prime \(E/Z\),

\[
 (e,f,\text{separable degree},\text{inseparable degree},
 v_E(\mathfrak D),\text{conductor},\text{sheet loss}). \tag{5.3}
\]

This is the minimum data capable of distinguishing a tame prime of index
\(m\), a wild separable prime of the same index, and a radicial edge.

For characteristic-zero boundary packages defined integrally, the initial
safe-prime test is:

1. choose one finite integral model of the normalization and affine open;
2. invert every denominator and every coefficient used to prove primeness,
   reducedness, degree, or noncontraction;
3. base-change the normalization and compare it with the normalization of
   the special fiber;
4. recompute the local degree sums including vertical primes;
5. compare Fitting ideals, differents, conductors, and nilpotency indices;
6. only then declare a ledger horizontal and surviving.

## 6. Lifting moment and cotangent obstruction classes

Scalar moment classes have a canonical integral lift.  If \(M_m\in\mathbb Z\)
is an exact mixed moment, then

\[
 M_m=0\pmod {p^a}\quad\Longleftrightarrow\quad
 v_p(M_m)\geq a.                                    \tag{6.1}
\]

Thus vanishing modulo \(p\) is only the first layer of an obstruction tower.
The radial SIC moments have monotone valuations, while a nonradial lift
already disappears and re-enters modulo \(11^2\).  Any modular obstruction
search should store the valuation, not a Boolean zero.

Cotangent and Kuranishi classes require a different lift.  Let
\(C^\bullet_{\mathbb Z_p}\) be an integral perfect complex.  The object to
reduce is

\[
 C^\bullet_{\mathbb F_p}
 =C^\bullet_{\mathbb Z_p}\otimes_{\mathbb Z_p}^{\mathbf L}\mathbb F_p,
                                                               \tag{6.2}
\]

not a matrix whose ranks happen to match.  Its base-change spectral sequence
contains the exact edge terms

\[
 0\to H^i(C_{\mathbb Z_p})\otimes\mathbb F_p
 \to H^i(C_{\mathbb F_p})
 \to \operatorname{Tor}_1^{\mathbb Z_p}
       (H^{i+1}(C_{\mathbb Z_p}),\mathbb F_p)\to0.   \tag{6.3}
\]

Accordingly, the lifting ledger distinguishes:

* reduction of an integral obstruction class;
* a new special-fiber class coming from \(p\)-torsion one degree higher;
* lifting the class as a module element;
* lifting a primitive or gauge which kills it;
* lifting the nonlinear Kuranishi section; and
* algebraizing a compatible Witt tower.

The boundary-obstruction package
[`BOUNDARY_OBSTRUCTION_THEORY.md`](BOUNDARY_OBSTRUCTION_THEORY.md) already
requires derived parameter base change and warns about rank drops and
vertical components.  The new task is to implement (6.2)--(6.3) for one
actual moment/cotangent complex, beginning with the relative canonical
symplectic control and then the first nontrivial SIC Kuranishi slice.

## 7. First work packages

### PCD.1 — spread-out boundary ledger

Choose one canonical characteristic-zero boundary theorem with integral
formulas.  Produce a machine-readable list of every coefficient inverted in
the proof, then recompute normalization, boundary primes, residue maps,
degree sums, Fitting ideals, and conductor at each exceptional prime.  The
output labels components horizontal, vertical, embedded, tame, wild, or
radicial.

### PCD.2 — Frobenius-vacuity registry

For every conjecture formulated by eventual powers or derivatives, record
whether Frobenius makes the premise impossible, the conclusion automatic,
or every positive-degree element nilpotent.  Rewrite only the vacuous rows
as a horizontal finite-prefix or Witt-level question.

### PCD.3 — Hasse/ordinary comparison

Retain the common weighted \(p\)-typical tuple, never separately
projectivized channels.  Determine which coefficient-level reconstructions
descend to map-intrinsic marked edges and which require an enriched map.

### PCD.4 — wild missing-boundary atlas

Starting from (5.1), compute the different and conductor for
\(m=p\), \(2p\), and \(p+1\), then allow collisions of several primitive
roots over one target point.  Compare separable wild and radicial rows with
the same \((e,f)\).  `HKM2B1` supplies the first completed plane row, with
generic data `(1,1,2,1,2)` and three reduced intersections between the
retained and missing components.  At each intersection the completed reduced
boundary is the node `k[[p,a]]/(pa)`; its conductor is `(p,a)`, while the
different remains `(p)`.  The general \(m=p,2p,p+1\) family and
bounded-degree algebraization of the compatible stable Witt tower remain
open.  `HKM2W1` obstructs both the exact plane representative and every
polynomially left--right equivalent plane representative.

The characteristic-divisible continuation now provides the reusable search
order.  Writing `N=p^n*m` with `(m,p)=1` separates the fierce component's
residue degrees as `(m,p^n)`; exact `N=6` controls in characteristics two and
three verify that the construction is not limited to prime powers.  The
canonical monomial gluing has a residual tame different for every `N>2`.
Balancing the gluing exponent removes that height-one different, but the
normal complement has a `G_m^2` core with valuation matrix

\[
 \begin{pmatrix}1&0\\N-2&N-1\end{pmatrix}
\]

and hence class group `Z/(N-1)`.  Thus a wild row must pass both a saturated
different gate and the complete
[class-trivial-core Smith gate](../plane-jc/BOUNDARY_LATTICE_PREFILTER.md#dual-torus-core-localization)
before point counts, Witt lifting, or coefficient algebraization are
informative.  These gates prove a characteristic-two exception only for the
complete one-boundary monomial-gluing architecture; other boundary
presentations remain open.
The exponent search itself is now complete for all `a>=0`: `a=0` has core
matrix `(N+1,N+1)`, a free unit, and class group `Z/(N+1)`; `a>=N` has a
wild index-`N` branch over `P=0`.  If `P^a` is replaced by an arbitrary
`C(P)`, every factor away from `P=0` creates an extra different divisor; on a
normalized noncollision prime its coefficient is `h*ord_D(R)`, not
automatically `h`.  For arbitrary `C(P,Q)`, the same support gate reduces the
search to `cP^aQ^b`.  The base-change Euler/transfer theorem and the uniform
full-source identity `Cl(U_(2,c))=Cl(D(xu))=(Z/2)^2` now close every `b>0` row.
For new families, the unequal-multiplicity form gives a cheap first sieve:
a principal reduced packet with generic unit rank one carries torsion exactly
when two fibre multiplicities share a prime.
For balanced multi-retained orders, the independent support identity forces
the retained polynomial to have the form `a0+T*B(T^p)`; therefore only cover
degrees congruent to one modulo `p` pass the divisorial different gate.  The
normalization then gives the exact root-count formula
`#(C_A-E_A)(F_q)=q^2+(n_q(A)-1)q` and
`chi_c(C_A-E_A)=deg(A)` after splitting `A`.  Thus every nonlinear
squarefree retained polynomial fails geometric affine-plane reconstruction.
The first degree-seven characteristic-three coefficient scan passes the
different gate in all six rows; four fail over `F_3`, and the two temporary
survivors fail over `F_27` with `810` rather than `729` open points.
Together with the linear class-group obstruction, this makes the
characteristic-two cubic exceptional throughout the balanced
squarefree-retained architecture.
Thus a next one-boundary candidate must change the boundary support or the
gluing architecture, not merely use a more complicated coefficient.

### PCD.5 — reduction of stable moduli

Use Theorem 2.1 as the base stratum.  Determine the quotient stack,
stabilizers, and rational point counts of the marked-edge invariant over
finite fields.  Compare them with the ordinary Hessian quotient and record
exactly where the latter ceases to be finite.

### PCD.6 — obstruction lifting

Freeze an integral basis and filtration for a bounded correction complex.
Compute Smith normal forms, the Bockstein maps, and the Tor term in (6.3).
Then test whether modular primitives lift through \(W_2,W_3,\ldots\), rather
than inferring a characteristic-zero primitive from one good-prime rank.

For plane boundary modules, run the Smith computation twice and keep the
matrices distinct: first on the codimension-one character valuations of the
normal reconstruction open, then on the integral correction complex used
for Bockstein/Tor.  The first is a geometric unit/class-group obstruction;
the second is a mixed-characteristic lifting obstruction.  Neither can
cancel a residual different divisor.

For the Huq--Kuruvilla--Mondello plane map, `HKM2W1` makes the first test
unbounded and terminal.  The complete cokernel is
`H^2_dR(F_2[x,y]) = xy F_2[x^2,y^2]`, and the integral Jacobian error has
class `xy(1+x^3(1+xy))^2`, with nonzero Cartier image.  Hence neither the map
nor any polynomially left--right equivalent plane representative has a
Keller lift through `W_2(F_2)=Z/4`.  After one identity stabilization,
however, every first obstruction is exact; for this map the truncated
geometric inverses of `1+2K` give compatible polynomial Keller lifts over
every finite `W_n(F_2)`.  Their inverse limit is restricted two-adic rather
than polynomial.  At `W_2`, the exact minimum stable degree is 18: the
functional `[x^13 y^4]+[x^14 y^5]` excludes every degree-at-most-17 lift and
the geometric correction attains degree 18.  Thus
`18 <= d_n <= 17(n-1)+1`.  A separate degree-25 construction at `W_3`
improves that level's upper bound from 35.  The unrestricted degree-18
coefficient system has 1,083 Boolean variables and 1,639 exact constraints;
Z3 proves it unsatisfiable.  At degree 19, a complete z-linear witness system
has 2,685 variables and 4,513 constraints; its SAT model has an odd constant
Jacobian modulo 8 under a direct sparse-polynomial replay and normalizes to
determinant one.  Hence `d_3=19` exactly.
This has no independent implementation or external human review.  The viable
continuation is quantitative and algebraic: later exact degrees, asymptotic
growth, bounded support, and algebraization of this stable tower.  The older
degree 25 remains sharp only in the canonical first-digit gauge; the optimal
lift necessarily uses a noncanonical kernel correction.
The preferred frozen optimal-degree witness has determinant one and support
`280+160` across the two Witt digits.  With its first digit fixed, the second
completion system has 1,485 variables, rank 681, and affine dimension 804.
It splits into independent `z0`, `z1`, and `z2` blocks and then into small
incidence components.  Exact componentwise minima are 92, 54, and 14, so
`s_2=160`; the determinant-one constraint has the same minimum.  The earlier
gauge gives `299+165`, showing cross-Witt support tradeoff.  Global support
minimality over all first digits is open.

The next digit separates degree optimality from tower compatibility.  The
preferred degree-19 `W_3` witness has a nonzero 48-term top de Rham class and
cannot extend to `W_4` at any degree.  A fixed degree-25 representative has
zero top class and exact extension degree 52: a two-coefficient functional
excludes degrees at most 51, while the third-coordinate correction `8zL`
attains 52.  Thus only `19 <= d_4 <= 52` is currently unrestricted.  Future
searches should carry the next determinant-digit cohomology class as a
frontier coordinate alongside present degree and support.
The zero-class locus is nonempty already in the degree-19 ansatz: a pinned
818-term `W_3` representative lies on it.  Its exact fixed extension degree
is nevertheless 52, so cohomological extendability and bounded-degree
extendability are separate frontier coordinates.  Exact unrestricted `d_4`
requires the joint bounded correction system, not just vanishing of the next
class.  The exact master--subproblem experiment has so far returned only
singleton codomain-hole duals through 64 cuts; a 599-hole master reaches the
600-second bound with status `unknown`, so this is evidence about the shape
of the obstruction space rather than an UNSAT theorem.  Factoring through
shared projected Jacobian minors makes the complete 4,340-hole quotient
computable, but its master also returns `unknown`; separately, its `z1` and
`z3` layers are SAT while `z0` and `z2` remain undecided at the same bound.

<!-- status-consumer: HKM2W1 904c57385ac0b0dd -->

<!-- status-consumer: HKM2W2 474e0d677133ee23 -->

<!-- status-consumer: HKM2W3 55d99efeee0298af -->

<!-- status-consumer: HKM2W4 6075dd4fbb9cb89c -->

### PCD.7 — arithmetic monodromy degeneration

Use Theorem 3.1 as the additive stratum.  For each \(p\)-typical support,
compute geometric monodromy over \(\overline{\mathbb F}_p(s,t)\), arithmetic
monodromy over \(\mathbb F_q(s,t)\), inertia filtrations at the discriminant
and infinity, and specialization maps from an integral characteristic-zero
model.

### PCD.8 — preserved-coordinate and Mondello-extension queue

For every new counterexample, search degree by degree for a polynomial source
coordinate \(c(x)\) and a target linear combination \(L\) with \(L(F)=c\).
Use linear algebra and Groebner reduction, and apply the detector to the
foundational, weighted, cancellation, quadratic-gauge, cubic BCW, and
positive-characteristic examples.  This is now a standard frontend beside
grading detection, LND searches, tangent quotients, and normalization
analysis.

For Mondello's plane map, the
[canonical follow-up queue](../verified/HUQ_KURUVILLA_CHARACTERISTIC_TWO_AUDIT.md#9-follow-up-experiments)
records the remaining independent normalization audit, two-route monodromy
verification, higher-Witt degree growth, and boundary-led searches for
odd-characteristic analogues.  The arbitrary-field corollary, full cokernel,
Cartier class, stabilization theorem, plane left--right lifting question, and
exact first stable Witt degree are closed.

## 8. Initial exact audit

Run

```bash
python3 scripts/verify_positive_characteristic_deformation_landscape.py
```

The dependency-free audit:

* verifies that the two normalization equations have rank two in every
  tested characteristic and degree;
* checks monomial recovery by the forced \(p\)-typical Hasse channels;
* verifies that ordinary differentiation misses \(W^{p^i}\) while the
  unique \(p^i\)-th Hasse channel detects it;
* constructs the affine group on \(\mathbb F_p\), checks its order and sharp
  two-transitivity, and compares it with \(S_p\); and
* enumerates the first tame and wild multiplicity rows for the missing-prime
  ledger.

These bounded loops are regressions for the proofs above.  They do not prove
normalization base change, a new cotangent lift, or a general monodromy
classification.

## 9. Claims not made

This programme does not currently claim:

* that every characteristic-zero canonical boundary ledger has good
  reduction outside the primes visible in its displayed formulas;
* that ordinary Hessian moduli commute with reduction;
* that Hasse derivatives alone define a map-intrinsic stable invariant;
* that equality of \((e,f)\) determines a wild boundary extension;
* that modular vanishing of a moment or cotangent class lifts to
  characteristic zero;
* that \(S_N\) monodromy survives at primes dividing the relevant degrees;
  or
* that a finite-field point count proves a geometric component or a
  characteristic-zero theorem.
