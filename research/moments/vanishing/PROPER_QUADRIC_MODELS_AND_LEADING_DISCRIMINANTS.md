# Proper quadric models and leading discriminants

15 November 2026. Internal written proofs with selected exact controls;
independent specialist replay, external review and formal verification
remain pending. No novelty is claimed. See the
[checker](verify_proper_quadric_gate.py),
[receipt](proper_quadric_verification.json) and
[source audit](proper_quadric_source_audit.json).

The [integrable continuation criterion](INTEGRABLE_PERIOD_DENSITIES_AND_OCTAHEDRAL_SEED.md)
can be supplied by geometry, without computing a separate density equation.
A proper map with isolated critical points and a holomorphic volume form
has only subpolar finite Cauchy-period singularities. For an affine quadric,
a smooth leading divisor constructs such a model explicitly.

Consequently every nonzero pure phase on either full exceptional quotient
must have a **singular leading divisor** on the projective rank-one quadric.
This is a necessary condition in all matrix degrees and includes phases
whose compact image has interior. It is not a classification of the
singular-leading locus or a decision of either exceptional MZ problem.
There are also two concrete all-degree seed closures: every polynomial
in the minimal real sphere seed is excluded on both 2O and 2I.

## 1. Proper maps with holomorphic volume forms

Let Y be a smooth complex algebraic variety of dimension N>=2, let
f:Y->A^1 be proper and nonconstant, with only finitely many critical
points, and let Omega be a holomorphic algebraic N-form on Y. Zeros
of Omega are allowed. Let M be a compact oriented real N-cycle in Y,
piecewise smooth, and suppose its volume period V=integral_M Omega is
nonzero. No reality of f on M is assumed.

For every nonconstant polynomial p, the moment sequence

    integral_M p(f)^m Omega, m>=0,                         (1)

is not eventually zero. Thus no nonzero element of C[f] is pure for
the normalized functional V^-1 integral_M(-)Omega. A nonzero constant
is immediate; zero is the trivial pure phase.

Here are the period and topology details used in this assertion.
For a regular fibre Y_t, write eta_t=Omega/df, meaning that
df wedge eta=Omega before restriction to Y_t. Over the complement of
the finite critical-value set, the compact fibre homology gives a
finite-rank local system. At a critical point, the local Milnor fibre
has finite middle homology. For a holomorphic top form, its local
Gelfand--Leray periods have convergent rational-power/log expansions
with every nonzero exponent strictly greater than -1. This is the
classical bound in [Malgrange, equation(6.2), printed420--421](https://numdam.org/item/10.24033/asens.1274.pdf),
using the positive exponents in(4.7.4). These are **local isolated-
singularity** statements, not a published general Haar theorem.

We need the bound also when a compact global cycle crosses the small
Milnor neighbourhood, rather than being a local vanishing cycle. It
can be seen in the following standard local decomposition. Choose
disjoint Milnor balls around all critical points of a fibre. Outside
these balls properness gives a compact regular fibration, with analytic
period contributions. Cutting that fibration along arcs to the critical
values reduces the remaining middle-dimensional chains to local
thimbles and Leray tubes. Indeed in a contractible Milnor ball B,
the boundary map H_N(B,F_t)->reduced H_(N-1)(F_t) is an isomorphism;
choosing a basis of local fibre cycles supplies the thimbles. Tubes
account for the boundary contributions in the complement of the
moving pole fibre. Subdivision and the relative homology exact sequence
glue these local pieces to the fixed regular portion. Components
entirely in a fibre contribute zero to an N-form pulled back to that
complex (N-1)-fold.

For completeness, relative local pieces do not hide an integer pole.
After embedded resolution of a critical fibre, in each local chart

    f-t0 = product u_i^m_i,
    pulled-back Omega = J(u) product u_i^k_i du_1...du_N,

with k_i>=0 and J holomorphic; coordinates with m_i=0 are harmless.
The nonnegative k_i follow from holomorphicity of Omega and the
Jacobian of a modification of smooth spaces. On a fixed sector,
integrating the relative form over a subdivided fibre chain gives
powers with lowest possible exponents

    (k_i+1)/m_i - 1 > -1,                                 (2)

and finite logarithmic powers when exponents coincide. To see the
bound, expand J in a convergent power series on smaller charts,
solve the monomial fibre relation for one coordinate, and integrate
the remaining monomials between fixed chart boundaries and their
Puiseux-moving boundaries. Each endpoint exponent is one of the
positive ratios (k_i+1)/m_i minus1; equal endpoint exponents give a
logarithm. Finite subdivisions and a common root cover suffice for
a fixed local sector. Terms from chart overlaps and the regular
portion are holomorphic or have the same bound. This also explains
the local analytic input to Malgrange's exponent estimate and why
relative chains crossing a node can add logarithms, not poles of
integer order. Only a positive delta is needed, not its optimum.

Now consider the infinity germ

    A_f(s)=integral_M Omega/(s-f).

It continues by transport of M in Y minus f^-1(s). Properness keeps
finite transports inside a compact part of Y over a bounded base
disk. The only finite singular parameter values are critical values.
The preceding decomposition expresses its local singular pieces as

    integral_arc g(t)dt/(s-t),   and residue additions 2*pi*i*g(s),

where g is a local Gelfand--Leray period; fixed regular pieces are
holomorphic. This identity is just integration of Omega=dt wedge eta
on a thimble and the residue formula on a tube. The local exponent
bound, or (2) for relative pieces, gives on each fixed sector

    A_f(s)=O(|s-t0|^(delta-1)(1+|log(s-t0)|)^K)+O(1),
    delta>0,                                              (3)

for every continuation branch. This is exactly the local Cauchy
estimate already proved in the integrable-continuation note.
The global regular portion cannot introduce another finite singularity.

For arbitrary p of degree d, partial fractions give, near infinity,

    A_p(s)=integral_M Omega/(s-p(f))
          =sum_(p(xi)=s) A_f(xi)/p'(xi).                    (4)

This is an identity of the principal infinity germs and remains true
under continuation, with the appropriate branches. At a finite
ramification p(xi)-p(c)~(xi-c)^k, (3) and division by p'(xi) give
power delta/k-1>-1. At an ordinary point of A_f use delta=1.
Thus every possible finite singularity of A_p is subpolar too.
If A_p is rational it has no finite poles; it tends to0 in its original
infinity germ, and hence is zero. Eventual vanishing of(1) makes A_p
a finite polynomial in1/s, contradicting its coefficient V at1/s.
This proves the assertion, for complex f as well as real seeds.

The global compact-cycle decomposition and its application here are
internal written arguments. The selected algebraic checker does not
certify Milnor topology, embedded resolution, or these analytic estimates.
Properness, holomorphicity on the whole finite model and isolated
critical points are substantive hypotheses. They are proved below
for the quadric model; they are not inferred for an arbitrary affine map.

## 2. A direct proper model of an affine quadric polynomial

Let Q_N={q(x_0,...,x_N)=1}, where q is a nondegenerate quadratic form,
and let F be a nonconstant regular function of its minimal polynomial
degree D. Its leading class is a nonzero section of O(D) on

    B={q=0} in P^N,

the smooth projective boundary of the projective quadric
barQ_N={q=u^2} in P^(N+1). Assume the zero divisor of this section on B
is smooth as a scheme. In particular repeated components and singular
crossings are excluded. Assume also D>=N.

Homogenize a degree-D representative to H(x,u), and take the closed graph

    Y={H(x,u)=t*u^D} in barQ_N times A^1_t.                 (5)

This is proper over A^1. It contains Q_N as the graph of F. At a boundary
point H=0,u=0, transversality on B permits smooth local coordinates
(u,v,z_1,...,z_(N-2)) on barQ_N with v=H in a local trivialization.
Equation(5) is simply v=t*u^D. Thus Y is smooth there, with coordinates
(u,t,z), and its projection to t is a submersion. Inside Q_N it is
isomorphic to the smooth affine quadric.

The invariant residue volume form Omega of Q_N has pole order N on B.
One can obtain this from adjunction K_barQ=O(-N), or directly from
the residue formula in a projective chart. Locally it is

    Omega=u^-N J(u,v,z) du wedge dv wedge dz,

with J a holomorphic unit. On the graph it becomes, up to sign,

    Omega_Y=u^(D-N) J(u,t*u^D,z) du wedge dt wedge dz.       (6)

It extends holomorphically precisely in the range used here, D>=N.
There is no finite boundary pole. Critical points of t on Y are absent
near its boundary. Any positive-dimensional irreducible component of
its critical locus would have constant t, would be proper inside that
fibre, and would lie in the affine open Q_N. A positive-dimensional
proper complex variety cannot be a closed subvariety of an affine
variety. Hence all critical points are isolated; algebraicity makes
their set finite. All hypotheses of Section1 are now established.

For q=sum x_i^2, normalized Omega on the real S^N is its positive
uniform volume form and has nonzero total period. Therefore

    smooth leading divisor and D>=N
       => only0 pure in the entire subalgebra C[F].        (7)

The same conclusion holds after a complex linear change of quadric
coordinates. F need not be real, and its compact image may have interior.
No implication is asserted when D<N or when the leading divisor is
singular. Polynomial composition usually makes the leading divisor
nonreduced; Section1 closes C[F] using the seed model, rather than
incorrectly applying the smooth-divisor test anew to p(F).

## 3. Every full exceptional pure phase has singular leading divisor

Use the determinant quadric SL2={xw-yz=1}, of complex dimension3.
Its real quaternionic core is SU2=S^3 and its residue volume normalizes
to Haar. The boundary is the smooth rank-one quadric

    B={xw-yz=0} in P^3 = P^1_left times P^1_right.

Every minimal-degree leading class of degree D is a section of O(D,D).
For F in A_Gamma=C[SL2]^Gamma, moments on the finite quotient equal
the normalized moments upstairs. The existing minimum nonconstant
matrix degrees are8 for2O and12 for2I, hence D>=3 automatically.
It follows from(7) that a nonzero pure exceptional phase **must** have
a singular leading zero divisor on B, in every original matrix degree.
This is an all-degree necessary gate on the full algebra, including
non-T-invariant and image-interior phases. It does not show that a
singular leading phase is pure, safe, or a failure.

The smooth-leading stratum is nonempty in infinitely many invariant
degrees. Write the retained binary invariants

    f_O=z^8+14z^4w^4+w^8,
    g_O=z^2w^2(z^4-w^4)^2,
    f_I=zw(z^10+11z^5w^5-w^10),
    g_I=z^20-228z^15w^5+494z^10w^10+228z^5w^15+w^20.

The pairs (f_O,g_O) and(f_I,g_I) have no common projective zero.
For every k>=1, let D=24k on2O or D=60k on2I. The right binary forms

    (f_O^(3k),g_O^(2k)), or (f_I^(5k),g_I^(3k)),

have degree D and no common zero. Tensor their span with the complete
degree-D left binary system. This gives an invariant base-point-free
linear system on B.

Here is a direct smoothness argument, without assuming invariants
separate all right tangent directions. Off the finite ramification
set of the nonconstant right ratio f^a/g^b, the value and two first
derivatives impose3 independent linear conditions on a section; the
point varies in dimension2. At a right ramification point, its value
and left derivative still impose2 independent conditions, while the
point varies only in the one-dimensional left P^1. In each case the
incidence of singular members has positive codimension. Properness
of B makes its projection closed. A general member is therefore smooth.
The Segre restriction map from degree-D polynomials in x,y,z,w onto
H^0(B,O(D,D)) is surjective; finite averaging preserves surjectivity
on right invariants. Thus these are actual invariant leading classes,
not sections with no regular lift.

Consequently there is a nonempty Zariski-open smooth-leading set in
each of those degree-D invariant coefficient spaces. Every seed F
in it, with arbitrary invariant lower-degree terms, has only0 pure
in C[F]. This concerns complex Haar moments; no arbitrary coefficient
ring or computed universal numerical cutoff is inferred.

## 4. Both minimal real sphere seeds are closed in every degree

For the complete original left-T-invariant algebras, use the unit
sphere quadric a^2+b^2+c^2=1. These coordinates have original matrix
degree2. Its projective null conic is parametrized by

    a=(z^2-w^2)/2, b=(z^2+w^2)/(2i), c=zw.                 (8)

For2O take the real harmonic degree4 seed

    h_O=a^4+b^4+c^4-(3/5)(a^2+b^2+c^2)^2.

Its null-conic restriction is f_O/8. It is squarefree, including at
the projective endpoint w=0, and so the leading divisor is smooth.
Its sphere degree4 is at least the dimension2 threshold.

For2I the binary form f_I has degree12 and is squarefree, including
at0 and infinity. A real harmonic degree6 lift is available explicitly.
Put q=a^2+b^2+c^2, let Delta be the ordinary three-variable Laplacian,
and set

    R=2a^5c-20a^3b^2c+10ab^4c+11c^6,
    h_I=R-q*DeltaR/22+q^2*Delta^2R/(22*36)
                         -q^3*Delta^3R/(22*36*42).        (9)

Its restriction under(8) is f_I, since the q terms vanish, and Delta h_I=0.
Restriction from harmonic forms of degree6 to binary forms of degree12
is an isomorphism compatible with the rotation action and real structure.
Thus h_I is invariant for the aligned icosahedral group that fixes f_I.
No general nonunitary finite-group embedding is being substituted for
the compact real core. Again the sphere degree6 is at least2, and its
squarefree restriction satisfies Section2.

Therefore only0 is pure in C[h_O] and C[h_I], for **all** complex
polynomial reparametrizations. Nonconstant polynomials do not even
have eventually vanishing pure moments. Their matrix phase degrees
are respectively8deg p and12deg p. This includes arbitrary real
affine changes of either seed.

The retained degree minima and one-dimensional binary invariant spaces
in degrees8 and12 identify every real minimal-degree nonconstant
original T-invariant seed, up to a real affine change, with these
respective seeds. Indeed the sphere coordinate ring is the direct sum
of harmonic spaces of degrees ell, which correspond to binary modules
V_(2ell), with original matrix degree2ell. Finite averaging takes their
Gamma invariants; the first nonconstant summand has dimension1 by
the retained Molien series. Its real part is the line generated above.
This closes that entire minimal-real-seed class in
both models, not all full phases of matrix degree8 or12 and not all
higher-degree real sphere seeds. The earlier product seed27a^2b^2c^2
has repeated null-conic zeros and is instead covered by its separate
hypergeometric-period argument; the two criteria are complementary.

## 5. Remaining problem

The full pure locus lies inside the singular-leading locus. Its
singularities may retain boundary poles or nonintegrable continuation
branches, so the proper holomorphic model above cannot be presumed
there. Some singular-leading phases are already safe; singularity
is a necessary condition, never a failure certificate.

Preserve full non-T-invariant thin phases, other actual Haar graph
laws, the full image-interior marked-safety problem and the separate
symmetrized height gate. The original matrix marker bounds, degree
gates, triangle conditions and Gaussian radial bridge remain in force.
Both exceptional MZ decisions remain **OPEN**. No expanding degree,
support, density or moment-prefix campaign is proposed.
