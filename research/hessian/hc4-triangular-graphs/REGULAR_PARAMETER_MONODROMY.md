# A regular polynomial parameter preserves arbitrary class monodromy

14 November 2026. The cyclic restriction on Danielewski output fibres
is unnecessary. A nowhere-critical polynomial parameter is relatively
algebraically closed in the source rational function field. This keeps
the full monodromy of finite covers of that parameter, with no abelian
restriction, and forces the relevant geometric image components principal.

## 1. Nowhere-critical polynomial parameters are regular (HC4MYHX0)

Let f in C[z1,...,zn] be nonconstant, with no critical point on A^n.
Then C(f) is relatively algebraically closed in C(z1,...,zn).
In characteristic zero this is a regular field extension. Consequently
the source rational function field is linearly disjoint from every
finite Galois extension of C(f). Degrees and n>=1 are unrestricted.

First f cannot be composite. If f=H(g) for polynomial g and
deg H>1, choose a root rho of H'. A nonconstant polynomial g attains
rho over C, for instance on a line where its restriction is nonconstant.
Every point of g=rho would be critical for f, a contradiction.

The explicit [primary input](sources/regular-parameter-import.json) is
the classical closed-polynomial implication in
[Arzhantsev--Petravchuk, Theorem 1 (i) to (iii)](https://arxiv.org/html/math/0608157v2):
noncompositeness makes C[f] integrally closed in the ambient polynomial
algebra. We use that theorem, not a claim to reprove it here.

Here is the passage to function fields, also explained in the source's
Lemma 3. If h in C(z) is algebraic over C(f), choose a nonzero q in
C[f] clearing all denominators of a monic equation for h. Then g=q*h
is integral over C[f]. It is integral over C[z] as well, so normality
of C[z] puts g in C[z]. Closedness puts g in C[f], hence h in C(f).

For clarity the finite-cover consequence requires no further geometric
classification. Let k=C(f), L=C(z), and M/k finite Galois inside a
common algebraically closed field. Then M intersection L=k. The
extension LM/L is Galois. Restrict its group to M; the fixed field
of this image in M is precisely M intersection L, since the fixed
field of Gal(LM/L) in LM is L. Thus the image is the whole Gal(M/k),
and [LM:L]=[M:k]. This proves the asserted linear disjointness.

Every component of a Keller map is nowhere-critical. More generally,
if u is a polynomial on its target with no critical points, u composed
with K has none, by the invertible DK. The same conclusion therefore
holds for C(u(K)) in the source field. No claim about all special
fibres, or about the whole finite extension C(K) inside C(z), follows.

## 2. Nonabelian monodromy from one regular target parameter (HC4MYHX1)

Keep the closed stable hypersurface setup of
[HC4MYHW1](ABELIAN_CLASSGROUP_OUTPUT_OBSTRUCTION.md#2-a-class-group-and-geometric-fibre-output-criterion-hc4myhw1):
plane Keller K(z), j:A^(r+2)->A^(r+3), r>=0,
Pi composed with j=K, A=C[t1,t2], F=Frac(A), Omega an algebraic
closure, and ambient polynomial ring R over A. Require the same
common free nonnegative filtration, generic graded domain with only
scalar degree zero, geometrically integral affine height-one target
fibres, and normal geometric generic ambient ring with torsion-free
Weil divisor class group.

Replace the abelian-action condition by the following supplied certificate:
there is a polynomial target parameter u for which f=u(K) has no
critical point, and a finite Galois extension

\[
             M/C(u),\qquad Gal(Omega/FM)
             \text{ acts trivially on }Cl(R_{Omega}).     \tag{1}
\]

Thus the class action descends from a finite cover of that single
parameter. It may be nonabelian. A nowhere-critical u on the whole
target is a sufficient test for f; a target polynomial coordinate is
one particular such u. Then **K is invertible**, in every source
and observable degree.

Indeed put L=C(z), k=C(u), and H=Gal(Omega/L). HC4MYHX0 gives
L intersection M=k and says that H restricts onto the full Gal(M/k).
Here LM/L is finite Galois and its automorphisms extend to Omega,
which is also an algebraic closure of L. Since F is contained in L,
F intersection M=k too, and G=Gal(Omega/F) has the same full image.

For the divisor D0 of the geometric source component corresponding
to the fixed embedding L into Omega, H fixes its class alpha.
By (1), G and H have the same image in the class-group action.
Therefore alpha is fixed by all G. All d conjugate source components
have this same class, while their reduced union has principal equation
E=0. It follows that d*alpha=0. Torsion-freeness makes alpha zero,
so all components are principal.

The remaining argument is exactly the written HC4MYHW1 proof:
generic degree additivity yields their common component pencil;
primitive basis coefficients and affine height-one geometric integrality
give F[h] intersection R=A[h] globally; smoothness of the image and
DK invertible make the derivative of the minimal label relation a
nonzero constant, so d=1. The retained birational Keller endpoint
gives its polynomial inverse. Neither ambient factoriality nor abelian
class action is required. No individual fibre embedding is rectified.

The restriction in (1) matters. Monodromy depending genuinely on
both target coordinates need not descend to a regular one-parameter
subfield. Nor does regularity of one source polynomial make every
finite extension of the full two-variable target field disjoint.

## 3. Every polynomial Danielewski observable is covered (HC4MYHX2)

Let P in C[Z] be any nonconstant polynomial, of arbitrary degree
m>=1. In A^(r+3), r>=1, with unused affine variables allowed, set

\[
                         Pi=(T,xy-P(z)).                  \tag{2}
\]

Every polynomial closed stable embedding carrying a plane Keller
pair through (2) forces that pair invertible. In addition, for every
polynomial B(T,z), including a constant, and every r>=2,

\[
                         Pi=(T,T Y+xy-B(T,z))             \tag{3}
\]

has the same property. All source, z and T degrees are unrestricted.
Target shears also permit adding an arbitrary polynomial c(T) to
the second component of either pair. These are output observables;
they need not include the HC4 potential.

### Arbitrary root monodromy in (2)

Let V denote the second target coordinate. Assign x,y weight m,
z weight two, and unused variables weight one, with A=C[T,V]
in degree zero. The monic relation xy=P(z)+V gives the common
free A-basis of monomials not divisible by xy. Its generic graded
ring has equation xy=c_m*z^m, where c_m is the nonzero leading
coefficient of P; this is a domain. Every affine height-one fibre
xy=P(z)+V is geometrically integral: as a polynomial in y it is
primitive over kappa[x,z] and irreducible over its fraction field,
since x and the nonzero polynomial P(z)+V are coprime. This remains
true at all special target primes because P has constant nonzero
leading coefficient and positive z degree.

Over Omega the polynomial P(z)+V has m distinct roots a_i. A common
root with P' would require V to be one of finitely many constant
critical values, impossible for the transcendental generic V.
The ambient generic fibre is smooth for the same reason. Localizing
x gives a Laurent polynomial ring. As in HC4MYHW3, its prime
divisors D_i=(x,z-a_i) have the sole relation div(x)=sum_i D_i:
any further principal relation would be a unit c*x^k after localizing,
and therefore have all coefficients equal to k. Thus

\[
                       Cl(R_{Omega})=Z^m/Z(1,...,1),       \tag{4}
\]

which is torsion-free. Its Galois action permutes the roots and
factors through the finite splitting field M of P(z)+V over C(V).
There is no restriction on this group. Take u=V in HC4MYHX1:
the actual source parameter is the Keller component K2, so it is
nowhere-critical. All hypotheses hold and the conclusion follows.

For example P(z)=z^3-3z has splitting group S3 over C(V).
The polynomial P(z)+V is irreducible by its primitive degree-one
dependence on V. Its discriminant is 27*(4-V^2), not a square in
C(V), because it has simple zeros at V=2,-2. An irreducible cubic's
transitive root group is A3 or S3, and the Vandermonde square root
of the discriminant distinguishes them. Thus this example has a
nonabelian class action and is outside HC4MYHW1's abelian hypothesis;
the new regular-parameter proof covers it. For m>=2 its nonzero
generic class group also precludes every old generic-factorial or,
in HC4 size, smooth-infinity certificate in any ambient frame.

### The full B(T,z) family in (3)

Put m=max(1,deg_z B), give x,y weight m, z weight two, and Y
weight 2m. The same monic xy relation gives a common free A-basis.
In the generic graded ring T is nonzero, so the leading relation
T*Y+xy minus the top weighted part of B eliminates Y, giving a domain.
The actual generic ring is polynomial as well. At each height-one
target prime other than T, Y still eliminates. At T=0 the residue
field is C(V), and the fibre is xy=B(0,z)+V with an unused Y.
It is geometrically integral by the same primitive-linear argument,
even if B(0,z) is constant: B(0,z)+V cannot be identically zero.
HC4MYHW1 applies with class group zero, proving (3).

There is also a literal reduction of the T-dependent coefficients:
B(T,z)-B(0,z)=T*H(T,z) in C[T,z], and replacing Y by Y-H(T,z)
is a polynomial ambient automorphism with the opposite translation
as inverse. It turns (3) into (T,T*Y+xy-B(0,z)). No localization
at T or source Hessian invariance is involved.

### Actual embeddings and HC4 compatibility

For every P in (2), j(s,t,w)=(s,1,t+P(w),w) is closed with
inverse (T,y-P(z),z). Adding an unused v gives an HC4-size embedding.
The independent nondegenerate ambient potential T*z+y*v+x^2 pulls
back exactly to

\[
                     1+s*w+t*v+P(w)*v,                    \tag{5}
\]

whose actual Hessian determinant is one. Its gradient inverse sends
(a,b,c,d) to (c-P'(a)*b,d-P(a),a,b). These are positive compatible
examples for arbitrary P, not nonautomorphic Keller constructions.

For (3), j(s,t,w,v)=(T,Y,x,y,z)=(s,w,1,t-s*w+B(s,v),v)
is a genuine closed embedding with output (s,t) and literal inverse
(T,T*Y+xy-B(T,z),Y,z). No unit Hessian for an arbitrary independent
ambient potential is inferred from this embedding alone.

## 4. Exact boundaries

A regular polynomial parameter does not make all its special fibres
irreducible. For f=s+s^2*t,

\[
      (1-2st)f_s+4t^2 f_t=1,\qquad
      f^{-1}(0)=V(s)\ \sqcup\ V(1+st).                    \tag{6}
\]

The two components are A1 and Gm and are disjoint. The generic fibre
is Gm, via t=(f-s)/s^2, and C(f) is relatively algebraically closed
in C(s,t). This is a polynomial submersion example, not a stipulated
Keller pair. Only generic function-field regularity is used above.

Nor does one regular parameter control covers depending on both
target coordinates. The non-Keller map (a,b)->(T,V)=(a,b^2-a)
has C(T) relatively closed in C(a,b), yet sqrt(T+V)=b lies in that
source field. Its Jacobian 2b is not a unit. This illustrates the
specific one-parameter limit, without being a Keller counterexample.

The nonsingular-parameter and filtration certificates are hypotheses,
not conclusions for arbitrary observable extensions. In particular,
the full first-type family (T,xy-B(T,z)) with genuinely T-dependent
z coefficients is not classified here; its class monodromy can depend
on both target coordinates. Arbitrary nonintegral height-one fibres,
torsion class groups and ambient rings without the required common
filtration remain outside the theorem. Constant P in (2) leaves
generic Laurent units and is not included by this argument; it is
allowed in (3), where generic Y elimination changes the situation.

Both plane endpoints, their actual unit-Hessian source-orbit problems,
and unrestricted HC4/JC2 remain open. The
[checker](verify_regular_parameter.py) and
[receipt](receipt-regular-parameter.json) give exact finite controls;
the general result is the written proof with the explicit classical
closed-polynomial import. No formal verification, independent replay,
external review or novelty priority is asserted.
