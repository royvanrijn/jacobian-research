# Integrable continuation densities and an octahedral Haar seed

15 November 2026. Internal written arguments, with selected exact algebraic
controls; independent specialist replay, external review and formal
verification remain pending. No novelty is claimed. See the
[checker](verify_integrable_period_density.py),
[receipt](integrable_period_density_verification.json) and
[source audit](integrable_period_density_source_audit.json).

The preceding [fixed-curve pole criterion](REAL_SEED_DENSITY_BOUNDS_AND_PERIOD_TRACE.md)
does not cover arbitrary moving periods. There is a useful extension that
does not require algebraic edge densities or a fixed compact curve: every
continuation branch must have integrable growth at every finite complex
singularity. A finite-rank differential equation can certify this condition.

This excludes every polynomial phase degree for a genuine octahedral
Haar seed, H=27a^2 b^2 c^2 on the unit sphere. Its density is a logarithmic
hypergeometric period, rather than an algebraic function. The theorem also
covers products of two independent beta variables with arbitrary positive
parameters, polynomial density tilts and finite affine mixtures.

This is a closure of seed classes, not a classification of all real seeds.
Both exceptional MZ problems, full non-torus-invariant thin phases and
image-interior marked safety remain **OPEN**.

## 1. The continuation hypothesis

Let gamma be a finite compact piecewise analytic path chain in the complex
t-plane. On its open pieces choose branches g(t) of multivalued analytic
functions. There is a finite set Sigma such that all these germs continue
along every path in C minus Sigma. Their continuations span a finite
dimensional space V of germs at an ordinary base point. Several such
spaces and path chains may be summed.

At every finite c in Sigma impose the following condition on **all**
continuation branches, not just the physical density on a real interval:
on each fixed local sector and sheet, for some delta_c>0 and integer N_c,

    |g(t)| <= C |t-c|^(delta_c-1) (1+|log(t-c)|)^N_c.       (1)

One may replace delta_c by min(delta_c,1). Finite dimensionality means
it suffices to check a spanning set; the constants may depend on the
branch. At ordinary points the same estimate holds with delta=1,N=0.
Endpoints in Sigma are allowed because (1) is integrable along the
finitely many incident analytic arcs. No condition at infinity is needed
for the argument below: all finite transports use compact contours, and
a polynomial phase has no finite preimages at infinity.

A useful sufficient certificate is a supplied regular-singular scalar
equation or system for the actual continuation space with finite singular
set, whose local solution functions all have Frobenius exponents of real
part strictly greater than -1. Resonance may introduce finite log powers.
One must check the actual scalar functions, including prefactors. A
physical branch can be integrable while another branch has a pole; that
does not satisfy the hypothesis. Holonomicity alone does not suffice.

## 2. All polynomial phase degrees are excluded at nonzero mass

Let p be any nonconstant complex polynomial and put

    A(s)=sum_gamma integral_gamma g(t)dt/(s-p(t)).           (2)

Under (1), a rational infinity germ A is identically zero. Therefore
eventual vanishing of the moments sum integral p(t)^m g(t)dt forces
vanishing of their zeroth moment as well. In particular a probability
law satisfying Section1 admits no nonconstant polynomial phase whose
positive moments all vanish. Positivity is only used to assert nonzero
total mass, not in the analytic argument.

Here is a contour proof that includes infinite branching of g. Subdivide
gamma at endpoints, singularities and relevant polynomial critical
points. Outside

    S = p(Sigma union endpoints(gamma)) union CritVal(p),

the moving points p^-1(s) are distinct and avoid the fixed density
singularities. Transport the finitely many contour pieces to avoid
these points, continuing their chosen g branches. For a fixed finite
parameter path this uses only finitely many sheets and compact paths.
Crossing a moving simple pole adds an integral over a small loop, hence
a residue proportional to g(xi)/p'(xi). This constructs continuation
along every parameter path outside S. Infinite global monodromy causes
no difficulty: no finite transport requires infinitely many pieces.

Fix s0 in S. Away from the finitely many t0 with p(t0)=s0, transported
contours may be fixed, giving holomorphic terms. Near such a t0 use a
holomorphic coordinate u with p(t)-s0=u^k. The density differential
has the bound C |u|^(delta-1)(1+|log u|)^N |du| with delta>0,
or the ordinary bound delta=1. Local paths are taken on finitely many
fixed sector rays; small connecting arcs obey the same estimate.
For r=|s-s0| in a sector away from these rays, splitting the path at
|u|=r^(1/k) bounds its Cauchy integral by

    O(r^(delta/k-1)(1+|log r|)^(N+1)) + O(1),              (3)

after decreasing delta to at most1 if necessary. The borderline
delta=k=1 has at most an additional logarithm, as included in (3).
For example use |s-s0-u^k|>=c(r+|u|^k) on the chosen rays and split
the resulting one-variable integral at r^(1/k). The finitely many
residue additions have the bound

    g(xi)/p'(xi)=O(r^(delta/k-1)(1+|log r|)^N).              (4)

If a transported contour winds around a fixed density puncture, keep
that finite winding on its chosen local sheets; (1) bounds its shrinking
arcs too. It does not create an integral with a worse power. Fixed
backbone portions are outside the disks and holomorphic in s. Thus every
continued germ has a bound o(r^-1) on a fixed local sector. This is the
only conclusion needed; no assertion of algebraicity is made.

If A is rational, a finite pole would contradict this subpolar estimate
on any such sector. There are no other finite singular points. A is
therefore polynomial, while its original infinity germ tends to0 by
compactness and integrability of the original chains. Hence A=0.
Eventual moment vanishing makes A a finite polynomial in 1/s, so the
same conclusion applies. The proof works for a finite sum even when
its individual summands are not rational. Multiplication of every g by
an arbitrary polynomial preserves (1), as do affine changes of t.

This follows the classical continuation-and-local-growth mechanism.
[Pakovich--Roytvarf--Yomdin, Theorem4.4 and Remarks2--3](https://arxiv.org/pdf/math/0312353v1)
describe continuation by sums of density branches and contour transport
for algebraic densities; Theorem4.7 characterizes rationality there by
trivial combinatorial monodromy. Their paper assumes algebraic g and
specific behavior on the contour. The infinite-branch, integrable-growth
argument (1)--(4) is written here, rather than silently attributed as
their stated theorem. Neither a vanishing jump at one vertex nor local
gluing is used as a substitute for the global argument.

## 3. Products of two beta variables

Let U,V be independent, with U~Beta(a,b), V~Beta(c,d), where all four
parameters are positive real numbers. The density of T=UV on (0,1) is

    w(t)=B(b,d)/(B(a,b)B(c,d)) * t^(c-1)*(1-t)^(b+d-1)
         * 2F1(c+d-a,b;b+d;1-t).                           (5)

Indeed the ordinary product-density change of variables gives

    w(t)=t^(c-1)/(B(a,b)B(c,d))
         * integral_t^1 x^(a-c-d)*(1-x)^(b-1)*(x-t)^(d-1)dx.

Set x=1-(1-t)v to obtain Euler's hypergeometric integral (5).
The normalization and positivity also follow directly from this
product construction. Euler's formula and the hypergeometric equation
are the standard identities in
[DLMF15.6.1](https://dlmf.nist.gov/15.6#E1) and
[DLMF15.10.1](https://dlmf.nist.gov/15.10#E1).

The complete two-dimensional solution space after the prefactor in (5)
has finite singularities only0 and1. At0 its exponent pair is

    c-1, a-1;

at1 it is

    b+d-1, 0.

All four exponents are strictly greater than -1. Integer differences
can add logarithms but cannot worsen these lowest powers; degeneracies
can reduce the continuation space. Thus **every** continuation branch
has (1), for all positive a,b,c,d. Section2 excludes every polynomial
phase degree. Polynomial density tilts with nonzero total mass, and
finite affine mixtures with nonzero total mass, have the same conclusion.
No inference that an arbitrary Haar density is a beta product is made.

## 4. A full all-degree octahedral seed family

Use the sphere with uniform probability area and coordinates
a^2+b^2+c^2=1. For its proper octahedral rotation group, define

    H=27a^2 b^2 c^2,  0<=H<=1.                             (6)

It is invariant under signed coordinate permutations, and is nonconstant.
The squared coordinates have Dirichlet(1/2,1/2,1/2) law. For completeness,
write independent standard Gaussians as their radius times a uniform
sphere direction; their squared values divided by their sum have that
Dirichlet law, by the gamma density and change of variables. Its beta
integral gives, for every m>=0,

    E(H^m)=27^m*(1/2)_m^3/(3/2)_(3m)
           =(1/2)_m^2/((5/6)_m*(7/6)_m).                 (7)

The second equality groups the 3m factors of (3/2)_(3m) by residue
modulo3. The final expression is exactly the moment sequence of the
product of independent Beta(1/2,1/3) and Beta(1/2,2/3) variables.
Both laws are probabilities on [0,1]; equality of all moments identifies
them, since polynomials are uniformly dense in the continuous functions
on this interval.

Putting these parameters in (5) and simplifying the beta constant gives

    w_H(t)=1/(3*sqrt(t))*2F1(1/3,2/3;1;1-t).               (8)

The constant is B(1/3,2/3)/(B(1/2,1/3)B(1/2,2/3))=1/3,
using Gamma(7/6)=Gamma(1/6)/6 and reflection at1/6. Conjugating the
hypergeometric equation by t^(1/2) gives the explicit certificate

    t^2(1-t) w'' + t(2-3t) w'
                    +(1/4-35t/36)w=0.                   (9)

Its finite indicial polynomials are (r+1/2)^2 at0 and r^2 at1.
Consequently all branches are bounded by
t^-1/2*(1+|log t|) at0 and by 1+|log(1-t)| at1 on local sectors.
This density is not algebraic: Euler's integral has a nonzero
logarithmic term at0. The product-density integrand in Section3, for these parameters,
has a nonzero positive 1/x leading part in the region t<<x<<1,
giving t^-1/2 log(1/t). An algebraic function has Puiseux expansions
without logarithms. Thus this is a genuine transcendental-period
application, outside the assumption that the density itself is algebraic.

In the original matrix ring C[x,y,z,w]/(xw-yz-1), take the left-torus
invariant regular sphere coordinates

    a=xz-yw,  b=i(xz+yw),  c=xw+yz.

They satisfy a^2+b^2+c^2=1, are real on SU2, and carry the standard
right rotation action. Choose the unitary realization of 2O aligned
with the coordinate octahedron. The function (6) is then in the
**original** A_(2O)^T, of matrix degree12, and its Haar law is (8).
Conjugating this realization only transports the same regular family.
The homogeneous degree12 form is nonzero on the rank-one matrix
(x,y,z,w)=(1,2,3,6), so is not divisible by xw-yz. It therefore survives
in the associated graded determinant ring; its powers have degrees12m
because that graded ring is a domain.

It follows that for every nonconstant p in C[t], the original phase
p(H) does not even have eventually vanishing pure moments. The only
pure phase in the entire subalgebra C[H] is0; nonzero constants are
immediate as well. This covers unbounded matrix phase degrees12deg p.
It asserts nothing about arbitrary other seeds or arbitrary 2O phases.
There is no corresponding 2I family in (6), since the icosahedral group
does not preserve the coordinate octahedron.

## 5. What remains to be opened

For an actual seed with a supplied finite-rank continuation system, the
new gate is its finite complex local exponents and actual continuation
space. If all branches have (1), the entire polynomial family is closed.
If this fails, isolate the nonintegrable branches and their continuation
maps: these are necessary possible sources of the pole1/s of a pure
phase resolvent. A bad exponent in an unnecessarily large annihilating
system is not enough; the physical branch might generate a smaller space.

Physical positivity and integrability on the real interval alone do
not suffice. For example -log(1-t)/t is positive and integrable on(0,1),
with a removable singularity at0 on that branch. Continuation around1
adds a nonzero multiple of1/t, so another branch has a finite pole at0.
This example fails Section1; no purity assertion about it is made.
Nor do regular holonomicity or an abstract period
representation without the local estimates. Inward-branch graph laws
for other Haar seeds, full non-T-invariant thin phases, image-interior
marked safety and the separate symmetrized height gate all remain open.
No expanding phase-degree, support or moment-prefix campaign is proposed.
