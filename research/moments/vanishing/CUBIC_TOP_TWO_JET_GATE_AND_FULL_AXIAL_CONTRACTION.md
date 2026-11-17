# The top two boundary jets and the full axial contraction-zero stratum

17 November 2026. Claim `SIC2B33TOPTWOGATE20261117`.
The geometric argument is an internal written proof importing the retained
proper holomorphic-period theorem. Its algebraic consequences below have
direct rational checks. External review and formal verification remain pending.

Let q=2zw+t^2 and F=H3+H2+H1 on q=1, with each Hj harmonic of degree j.
The constant is zero by the first moment. All assertions are field-valued
in characteristic zero. The full (3,3) Image problem and the three global
radical memberships remain OPEN; the overarching goal is ACTIVE.

**Top-two-jet gate.** If H3 is nonzero and F is pure, there is a point
on the projective conic q=0 where H3 has a multiple zero and H2 vanishes.
It is insufficient merely to require a repeated zero of H3 somewhere.

**Full axial consequence.** If B is nonzero semisimple axial and
C=Delta^2(H2 H3)=0, there is no pure phase, with arbitrary H1.
No finite consecutive prefix for this entire stratum is asserted.
The earlier [eigenvector theorem](AXIAL_CONTRACTION_ZERO_EIGENVECTOR_STRATUM.md)
still gives a separate first-eight-moment exclusion on its own stratum.

The [checker](verify_cubic_top_two_jet_gate.py) regenerates the boundary
jets and directly multiplies the two small rational units in
[the pole certificate](axial_contraction_zero_boundary_pole_units.json).
It does not certify the analytic period theorem or the global geometry.
No new Singular/msolve dependency enters the accepted algebraic proof.

## 1. A smooth proper graph with nontransverse boundary

Compactify the affine quadric to X={q=u^2} in P^3. Let B={u=0}.
Homogenize the phase as

    H=H3+u H2+u^2 H1.

The closed graph Y={H=tau u^3} in X times A1_tau is proper over A1.
Its affine open u!=0 is the smooth graph of F on q=1. At a point of B,
choose a coordinate s on the conic. In local coordinates (u,s), the
normal derivative of H at u=0 is H2, up to a tangential derivative
of H3. At a multiple zero of H3 that ambiguity is zero. Thus the graph
is smooth at every boundary point, and its map to tau is a submersion,
provided there is no point with

    H3=0,  d_B H3=0,  H2=0.                         (1)

Here are explicit charts, including the volume calculation. On w=-1/2
write z=s^2-u^2 and t=s. This parametrizes X, and on q=1 its affine
coordinates are (Z,T,W)=((s^2-u^2)/u,s/u,-1/(2u)). The residue form is

    Omega=dT wedge dW/(2W)=du wedge ds/(2u^2).

Write K=H-tau u^3. At the boundary, K_s=d_B H3 and K_u=H2.
If K_s is a unit, the graph has coordinates (u,tau), and

    Omega_Y = u/(2K_s) du wedge d tau.

If K_u is a unit, it has coordinates (s,tau), and

    Omega_Y = u/(2K_u) d tau wedge ds.

Both forms extend holomorphically. The analogous z=-1/2 chart covers
the missing point of the conic. Therefore absence of (1) constructs a
smooth proper Y with holomorphic extended volume and no boundary
critical point. This argument permits nontransverse intersections of
H3=0 with B when H2 does not vanish there; leading-divisor smoothness
is a stronger sufficient condition.

Every positive-dimensional irreducible component of the critical locus
would have constant tau in characteristic zero. It would be proper in
that fibre and contained in the affine quadric, which is impossible for
a positive-dimensional proper variety. Hence the critical set is finite.
The volume period on the real sphere is nonzero. All hypotheses of
`GPROPERHOLPERIOD20261115` now apply. They exclude a pure nonconstant F.
This proves the gate by contraposition. The same proper model excludes
eventual pure-moment vanishing for every nonconstant polynomial in F.

Complex coefficients need not be real. For other characteristic-zero
fields, embed the finitely generated coefficient field into C; the
moment equations and the existence obstruction are algebraic. This
argument establishes a necessary incidence condition, not a finite
ordinary ideal membership or an arbitrary coefficient-ring theorem.

## 2. The axial contraction kernel meets only the two poles

Normalize the nonzero axial quadratic to H2=t^2-zw. The complete
contraction-zero family is

    F=az+bt+cw+(t^2-zw)+dw^3+ew^2t+fz^2t+gz^3.     (2)

All seven coefficients remain free. On the boundary chart above,

    8 H3 = -d+2es+8fs^5+8gs^6,   H2=3s^2/2.

The quadratic's only projective zeros are s=0 and infinity. The cubic
has a multiple zero at s=0 exactly when d=e=0, and at infinity exactly
when f=g=0. Consequently every pure phase (2) has

    d=e=0  OR  f=g=0.                              (3)

For H3!=0 this follows from the top-two-jet gate. If H3=0, it holds
automatically, so no separate degree-three assumption is imposed in
the remaining algebra.

In either pattern dg=ef=0. The first two nontrivial moments give

    (2ac+b^2)/3+1/5=0,
    2(b^2-ac)/5+2/35=0,

and therefore

    ac=-16/105,   b^2=-31/105.                      (4)

In particular a,b,c are units at a field point. Use amplitude and the
actual stabilizing torus to put b=c=1, retaining the nonzero quadratic
amplitude lambda. Equations (4) become

    a=16/31,   lambda^2=-105/31.

This spends only available symmetries. Odd moments are divided by the
declared unit lambda, as in the retained complete linear normalization.
The two zero-coefficient patterns survive this transformation.

## 3. Two rational units finish the entire stratum

Use the exact normalized moments R4,R5,R6 from
`axial_contraction_zero.py`. On d=e=0 these are polynomials in (f,g).
The first two are independent affine linear equations, whose root is

    (f,g)=(91224/10571,-30240/10571).

On f=g=0 the variables are (d,e), and the corresponding root is

    (d,e)=(-29295/1408,11403/352).

For both roots the exact value of R6 is

    198023471767658880/121 != 0.

The certificate gives three multipliers on each chart, with respective
term counts (2,2,1), and the checker verifies

    A4 R4+A5 R5+A6 R6=1

over Q. No coefficient is inverted except the field units proved in
(4) and the nonzero rational constants. This excludes both patterns
(3), and hence every pure point of the full axial contraction-zero
stratum. Six moments exclude the two patterns; the geometry supplies
their necessity only for all-order pure phases. Do not conflate these
statements with an unconditional first-six or first-eight exclusion.

## 4. A larger closed pure-safe cone

Drop the eigenvector condition from the preceding closed container.
Require only that B belong to the closed harmonic-square quadratic
cone, and C=0. Its seven cubic quadratic-cone equations and three
quadratic contraction equations are exported in
`axial_contraction_zero_full_safe_container.json` in the full
fifteen-dimensional harmonic coefficient space.

Its nonzero non-nilpotent points are semisimple axial, now excluded
with arbitrary H1. Its nilpotent boundary has rank at most one and is
covered by the retained rank-one contraction-zero classification.
At B=0 the all-order odd classification applies. Every pure field
point of this larger closed container is therefore one-sided, with
sharp full-original-multiplier cutoff m>coordinate-degree Q.
The container is not entirely pure and has not been proved to contain
every pure phase. The old eigenvector theorem is preserved for its
stronger finite-prefix and purely algebraic assurance.

The large mixed-chart rational standard-basis attempts, degree-seventeen
partial checkpoint, degree-eighteen timeouts, M2 engine failures and the
later discriminant-augmented timeout remain preserved. Their failed
computational certificates are not promoted to successes. The new proof
supersedes the *all-order* axial contraction-zero question by a geometric
argument and two tiny units. The unrestricted mixed chart's first-eight
ordinary ideal/unit question remains uncomputed.

Next global work can use the necessary incidence (1), retaining the
entire flag locus and its projection. A repeated cubic root alone is
weaker. The other quadratic types, C!=0 pure points and all three
unrestricted global radical memberships remain open.
