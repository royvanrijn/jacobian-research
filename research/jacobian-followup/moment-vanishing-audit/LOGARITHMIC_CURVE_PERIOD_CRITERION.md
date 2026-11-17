# An all-degree criterion for logarithmic curve periods

15 November 2026. Internal written proof extending the twice-punctured-sphere
estimate in [Quadratic W logarithmic periods](QUADRATIC_W_LOG_PERIODS.md).
No independent review, external review, formal verification or novelty claim
is made. The [exact checker](verify_logarithmic_curve_criterion.py) verifies
coordinate, connection and boundary examples only; it does not certify the
contour or comparison arguments below. Prior canonical proof bodies remain
unchanged.

## 1. Statement: arbitrary phase degree and moving punctures

Work over K=C(t), allowing a finite algebraic extension when needed to split
punctures. Restrict to a sufficiently small generic sector at t=infinity,
after a finite ramification and removal of finitely many exceptional base
values. Let D_t be a nonempty finite set of algebraic points on P1, distinct
on this sector, and let

    U_t=P1 minus D_t,
    Phi(t,x) in K(x), with every fibre pole in D_t.

Let omega be a rational relative one-form with at most simple poles, all in
D_t. Require its residues to be bounded as t tends to infinity on the
sector. This is a geometric condition on a differential, including its pole
at infinity; it is not just a condition on its coefficient in one chart.
Choose a nonzero algebraic normalization rho(t) with

    rho(t)=O(|t|^(-epsilon)) for some epsilon>0.

**Theorem.** In every zero-exponential block, any flat rapid-decay period

    F_Gamma(t)=integral_(Gamma_t) exp(Phi(t,x))*omega(t,x)

is O(1+log|t|). Consequently no flat period rho(t)*F_Gamma(t) is a
nonzero constant. The degree of Phi and the number of punctures are arbitrary
and fixed for the family. Multiple saddles, coincident actions, puncture
collisions at infinity and regular removed points are included.

Here an action is a critical value of Phi or its finite value at a regular
removed point, continued algebraically in t. A zero-exponential block has
no positive-power Puiseux part in its action; bounded constant actions are
included. Blocks with the same action are grouped before estimating them.
The statement does not bound the unnormalized nonzero-action blocks by a
logarithm: those may have exponential growth.

**Finite-cover corollary.** Let pi_t:C_t->P1 be any finite algebraic map of
fixed generic degree, with C_t smooth, and set V_t=pi_t^(-1)(U_t). If both
the phase and the form are pullbacks,

    Psi=pi_t^*Phi,        eta=pi_t^*omega,

the same no-constant conclusion holds for every normalized flat period on
V_t. The genus of C_t is unrestricted. Additional upstairs cycles can have
zero period; they cannot create a new constant period of these pulled-back
data. Both descent hypotheses are essential to this corollary.

## 2. Why logarithmic forms on P1 are controllable by residues

On a chosen affine chart write the distinct finite punctures as a_j(t).
A rational differential with only simple poles has the unique expression

    omega=sum_j r_j(t) dx/(x-a_j(t)).                 (1)

Its infinity residue is -sum_j r_j. There is no extra polynomial-differential
term: any nonzero polynomial times dx has a pole of order at least two at
infinity. There is no holomorphic differential on P1. Thus bounded residues
bound every coefficient r_j in this expression, irrespective of how the
punctures move.

This is invariant under parameter-dependent Mobius changes of coordinate.
Under a finite map, a simple logarithmic pole pulls back to a simple pole
whose residue is the ramification index times the old residue. An arbitrary
higher-pole differential does not satisfy (1).

## 3. Inverse-phase cuts and their exponential actions

First suppose Phi is nonconstant in x. For generic t it is a finite branched
map P1_x->P1_w. Cut the w-plane from each critical value and from every
finite value of a regular puncture toward the left half-plane. Away from
these cuts the inverse sheets are single-valued. Lifting the cuts and
gluing the boundaries gives adjacent-sheet differences at critical values,
descending tails at irregular ends, and small circle classes at regular
removed points. These span the relative rapid-decay homology: the remaining
inverse sheets are cut disks, so their interior has no additional one-cycle.
Use all local branches at a multiple critical point; no Morse assumption
is needed. A regular puncture is retained even when it is not a critical
point of Phi.

There are only finitely many sheets, endpoints and cut pieces, bounded in
terms of the fixed rational phase degree and puncture count. Puiseux
expansions label all action values and locations after finite ramification.
Choose a generic sector excluding the limiting Stokes directions of their
differences. Cut changes are then constant triangular changes of the
cycle basis. At a critical endpoint, the inverse differential has a
convergent fractional-power expansion in w-A(t); at a regular endpoint it
has a meromorphic expansion, with a logarithmic primitive when a simple
pole is present. Integrating these expansions against exp(w) gives the
usual power/logarithm factors times exp(A(t)). If endpoint scales merge,
keep their full branch star together before expanding. This groups equal
actions and does not create a different action. Descending tails carry no
additional finite endpoint action.

This inverse-sheet construction gives the action filtration of the period
space. A block whose action remains bounded can be represented by cuts
starting at its bounded endpoint values and regular-end detours, all lying
in Re(w)<=C for a constant C independent of large t. Unbounded action
blocks are separated first; exponentially small changes across a sector
do not change the power/logarithm estimate of a bounded-action block.
Distinct nonzero Puiseux exponential factors cannot supply a nonzero
constant solution: after grouping equal factors, uniqueness of their
sectorial asymptotic expansions separates them from the zero-exponential
block. We estimate that zero block directly next, without approximating a
coalescing saddle by a nonuniform Gaussian.

## 4. Uniform logarithmic length with many moving punctures

All critical locations, punctures, their nonzero separations and coefficient
transition scales are algebraic in t. Choose M large enough that all finite
locations are in |x|<=|t|^M, all distinct puncture separations exceed
|t|^(-M), and disks of a smaller power radius isolate their end models.
Increase M finitely many times as necessary. This is possible for a fixed
algebraic family: every relevant nonzero function has a Puiseux leading term.

Take cuts with a bounded number of horizontal or vertical pieces, and small
bounded-Re(w) detours when values coincide or would otherwise obstruct a
cut. A lifted straight piece lies on a real algebraic curve of uniformly
bounded degree B, given by Re(Phi)=constant or Im(Phi)=constant after
clearing denominators. The number of sheets and pieces is fixed. Such a
curve has length at most pi*B*R inside any disk of radius R: integrate the
total variation of its orthogonal projections over directions, use the
degree bound B on intersections with a generic line, and then the
one-dimensional coarea formula. Real line components obey the same bound.

For each puncture a_j and each dyadic annulus

    r<=|x-a_j|<=2r

in the truncated region, this gives

    integral |dx|/|x-a_j| <=2*pi*B.

The range of radii is between two fixed powers of |t|, so there are
O(1+log|t|) such annuli. Sum over the finitely many j in (1), with their
bounded residues. On the bounded-action cuts |exp(Phi)|<=exp(C); their
truncated contribution is therefore O(1+log|t|). This argument uses annuli
centered at each moving puncture, so collisions cannot introduce an
uncontrolled factor equal to an inverse puncture separation.

At an irregular end use y=x-a_j, or y=1/x at infinity. Beyond the chosen
transition scale the leading polar monomial of Phi dominates, and on a
lifted descending tail w=Phi is a coordinate satisfying

    |d log y/dw|<=C'/|w|.

The other logarithmic summands are bounded by a fixed multiple of
|d log y| there, because the end disk is smaller than every separation.
Integration of exp(Re(w))*|dw|/|w| bounds the tail. At a regular end whose
action is bounded, close with a small end circle on which Phi is uniformly
close to that action. Each logarithmic summand has bounded circular
length. This also covers a critical point approaching a regular puncture;
do not substitute the value of a singular amplitude at the critical point.

Thus the entire zero-exponential block is O(1+log|t|). Multiplication by
rho makes it tend to zero. Section 3 excludes a constant contributed by
the other action blocks. This proves the theorem for nonconstant Phi.

If Phi=phi(t) is independent of x, only ordinary puncture-circle periods
remain. By (1) every period is exp(phi(t)) times a fixed linear combination
of bounded algebraic residues. For bounded phi the normalized period tends
to zero. For unbounded algebraic phi it has a nonzero Puiseux exponential
factor, which cannot equal a nonzero constant times an algebraic function.
This proves the remaining case without an algebraic-independence assumption.

## 5. Finite covers: simultaneous descent is enough

Push a rapid-decay chain Gamma_t on V_t through pi_t. Finite maps are
proper, and exp(Psi)=exp(Phi) after composition, so the image is a
rapid-decay chain on U_t, with multiplicities and with boundary respected.
At a branch point, use a local coordinate x-a=y^e: integration of the
pulled-back form agrees with integration over the pushed-forward chain,
including the integer multiplicity e. Hence

    integral_(Gamma_t) exp(Psi)*eta
       =integral_(pi_*Gamma_t) exp(Phi)*omega.         (2)

Away from finitely many base values this pushforward commutes with flat
transport. The downstairs theorem applies to the right side. Kernel cycles
of pi_* have zero period for these data. This proves the corollary without
assuming that all upstairs cohomology comes from P1.

## 6. The exact seed-to-period bridge, including its limits

Here is the algebraic part isolated from the contour estimate. Let (H,nabla)
be a finite-dimensional differential module over K and v a seed. Suppose
a K-linear horizontal functional J:H->L exists in a differential field
extension, and J(v)=1. Let N be the cyclic submodule generated by
v,nabla v,nabla^2 v,... . Its first r independent generators form a basis
for N, and J takes these values to 1,0,...,0. Thus J|N is a rational
horizontal quotient N->K. Every rational relation among the generators has
zero coefficient of v; applying J proves this, without assuming v cyclic
in the whole H.

At an ordinary base point extend the functional on N to the fibre of H,
then solve the dual analytic connection. Horizontal uniqueness gives the
required extension on a simply connected ordinary neighborhood, and thence
on a generic sector avoiding the finite singular set. A rational extension
on all of H is neither needed nor generally possible. For example

    nabla e0=0,    nabla e1=e0/t,    J(e0)=1

extends analytically by J(e1)=log(t)+constant, while no rational primitive
of 1/t exists. Its residue at zero proves the rational obstruction.

For an exponential family of smooth affine curves use

    H=Omega^1(U/K)/(d_x+d_x Phi)(O(U)),
    nabla_0[f dx]=[(f_t+Phi_t*f) dx],
    nabla=nabla_0+(rho'/rho) id,    v=[omega].         (3)

Every derivative in (3) is taken at fixed x. In particular moving poles
of omega contribute their full derivatives; differentiating only Phi is
incorrect. The de Rham space is finite dimensional, and the normalized
connection is a scalar rational gauge of the actual Gauss-Manin connection.

The imported analytic comparison is [Hien--Roucairol,
arXiv:0704.1739v1](https://arxiv.org/html/0704.1739v1), Theorems 2.4 and
3.7 (=1.1), with Proposition 3.4. Its hypotheses concern smooth affine
total space, regular algebraically independent projection and phase, and
removal of a finite exceptional base set. Clear denominators by removing
their divisors from the total space and their bad base values. If Phi is
nonconstant on the generic fibre, it is algebraically independent of t.
A finite base change and smooth affine model give the same comparison
locally. The theorem represents analytic horizontal functionals by flat
rapid-decay cycles; it does not claim convergence of a formal Wick series.
The fibrewise constant phase was already handled separately above.

Combining this comparison with sections 1--5 proves a reusable exclusion:
no formal horizontal functional for (3) can take v to 1 when the stated
logarithmic/normalization hypotheses hold. To apply this to a Gaussian
family, one must first construct its exact formal J and prove its
integration-by-parts and tangent-derivative identities. The criterion does
not assert that every Gaussian polynomial admits such a presentation.

Equivalently, for a nonzero seed its minimal scalar relation

    nabla^r v+a_(r-1)nabla^(r-1)v+...+a_0 v=0

has a_0!=0: if a_0 were zero, assigning v to 1 and every positive derivative
to 0 would define the forbidden rational quotient on N, and the preceding
local extension/comparison would give a constant period. If v=0, the
order-zero relation already prohibits J(v)=1. This is an all-degree
criterion for the nonzero zeroth coefficient; no generic seed rank or
generic cyclicity assumption is allowed.

## 7. Essential hypotheses and exact counterexamples

* Higher-order poles cannot be allowed. On A1, take Phi=-x^2/(2t^2) and
  omega=dx, with the real rapid-decay contour for positive t. The period
  is sqrt(2*pi)*t, so rho=1/t gives a nonzero constant. The differential
  dx has a double pole at infinity and fails the logarithmic hypothesis.
* Growing residues cannot be allowed. On C*, Phi=x and omega=t dx/x.
  A fixed small positive circle has period 2*pi*i*t, and rho=1/t again
  makes it a nonzero constant. The logarithmic pole is simple, but its
  residue grows like t.
* Decaying normalization is essential. For Phi=x, omega=dx/x and rho=1,
  the same circle has the constant period 2*pi*i.
* Genus zero cannot be replaced by a residue-only condition on every
  curve. On a fixed elliptic curve choose a nonzero holomorphic form eta0
  and a cycle of nonzero period. The form t*eta0 has all residues zero,
  yet its period divided by t is a nonzero constant. Section 5 requires
  the form, as well as the phase, to descend; it does not include this case.

These counterexamples identify structural boundaries, not exceptional
numerical parameter points to discard.

## 8. Scope gained and remaining work

The old conic calculation is the two-puncture case omega=dx/x,
rho=1/(kappa*t). The new criterion permits any finite number of moving
punctures, any rational phase degree, and finite covers of arbitrary genus
with simultaneous phase/form descent. For example every rational phase
on P1 minus {0,1,infinity}, every bounded algebraic combination of
dx/x and dx/(x-1), and rho=1/t satisfy the period exclusion. This is a
new admitted period architecture; an exact Gaussian presentation must
still be supplied before claiming a new GMC family.

The formal cyclic-submodule bridge used in the existing notes passes this
scope audit. The contour/action-filtration proof here and its earlier
two-puncture specialization remain internally written analytic mathematics,
awaiting independent specialist review. The remaining broad task is to
find exact eliminated Gaussian families whose phase and full volume form
meet this geometric criterion, or whose data descend through a finite
cover. Arbitrary higher-genus forms and higher-pole amplitudes require a
different estimate. No least cubic prefix, global radical cofactor, or
completed witness calculation is reopened.

The [source audit](logarithmic_curve_criterion_source_audit.json) binds the
accepted input and receipt. The initial exact-check failure and input are
retained: rational cancellation alone did not normalize a product of
opposite exponential factors. Expanding exponential sums before rational
cancellation makes the unchanged full connection identity exactly zero;
the diagnostic also checks it with symbolic simplification. No geometric
hypothesis or identity was weakened.
