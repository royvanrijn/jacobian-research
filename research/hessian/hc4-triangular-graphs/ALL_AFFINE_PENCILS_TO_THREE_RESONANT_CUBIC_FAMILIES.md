# All affine pencils reduce to the old tests and three resonant cubic families

15 November 2026. Claims **HC4MYJP0/1/2**. The previous simple-zero
quartic frontier no longer needs an arbitrary outer polynomial.
Every affine quadratic composite in four variables reduces its
universal full-source inverse question to the four retained tests
and the following three explicit cubic families:

    q=2xy+z^2,
    F_N=q+alpha x+e z+t(lambda q+N-c),
    N=0, x^2, or 2xz,
    lambda*c*e!=0, alpha arbitrary.                    (1)

No bound on the outer or source degree is assumed. All three
families in (1) remain open. In particular a passing critical-point
count or whole-value profile does not close them or identify them
with h0,h1. The original HC4/JC2 problem and both old endpoints
remain open.

The reduction uses two new geometric inputs: an actual quartic
plane cover of positive genus, and closure of every ternary-cone
multiplier, including the degenerate sections outside the earlier
JG genus-one criterion.

## 1. Precise statements and retained inputs

Work over C. In the affine family let

    h=q_A+H(g), g=q_B+2v^T A x+c0,
    A=A^T invertible, B=B^T, T=A^(-1)B, deg H=d>=2.

**JP0.** For b*nu!=0 and arbitrary e,beta, the generic curve

    f0=z^2+w^2+e z+H(nu z^2+2bw+beta)=kappa             (2)

has positive total compact genus G, and in fact
G>a_d=(d+epsilon-2)/2, epsilon=d mod 2, whenever either
d>=3 or d=2 and 1+4b^2 h_2!=0. Here h_2 is the quadratic
coefficient of H. Its actual t=g projection has degree four
and 3d simple finite ramifications. For d>=3 its true puncture
count is r=gcd(4,d), and

    G=(3d-4-r)/2+c_C,                                  (3)

where c_C counts components. For a nonresonant quadratic outer,
there are two punctures and G=c_C>=1. No connectedness is assumed.

**JP1.** Every full-source orbit of

    a(z)+ell(z)+t b(z), z in C^3,                       (4)

with a,b nondegenerate homogeneous quadratics and ell linear
is inverse-closed under unit Hessian, without any cyclicity
or generic-section hypothesis. This is inverse closure, not
nonexistence of every unit-Hessian member.

**JP2.** For the entire four-variable affine family above,
the universal full-source inverse question is equivalent to
the four retained tests together with all members of (1).
Equivalently any additional affine mechanism beyond the four
tests must lie in one of those cubic source-orbit families.
No source-degree restriction or native Hessian invariance is used.
This does not assert that the three families are pairwise distinct,
minimal or new universal Keller endpoint equivalences.

Retain [HG0](GLOBAL_VALUE_ORBITS.md): native critical schemes
must be finite Morse with common nonzero Hessian value under
unit-Hessian source reparameterization; hypothetical gradient
nonautomorphy also forces every value profile
L^3-NL+N_kappa L^2. The additive invariant and curve formula
are retained in [HC0](HODGE_TRACE_CLOSURE.md), and componentwise
Riemann--Hurwitz in
[IY0](CONTRACTIBLE_RADICAL_CURVES_AND_PENCIL_DEFECT.md).

Use the proved constant self-adjoint primary frames
[JD0](NILPOTENT_QUADRATIC_PENCILS_AND_ENDPOINT_ORBITS.md),
stationarity identities [JF](AFFINE_QUADRATIC_PENCIL_RESIDUES_AND_ENDPOINTS.md),
complete simple compatible groups and active pole budgets
[JH](SIMPLE_STATIONARY_POLES_AND_NONCYCLIC_KRYLOV_CLOSURE.md),
the earlier cone criterion [JG0](ELLIPTIC_CONE_TRACES_AND_COMPLETE_CYCLIC_CLOSURE.md),
and whole eigenvector fibre strata
[JL](JORDAN_POLAR_CURVES_AND_ALL_INVERTIBLE_PENCIL_CLOSURE.md).
[JO](ISOTROPIC_KERNELS_AND_SIMPLE_ZERO_QUARTIC_FRONTIER.md)
already reduces everything beyond the four tests to a simple
nondegenerate zero eigenline with nonzero offset.

## 2. The mixed zero/nonzero plane curve (JP0)

The literal global coordinate t=nu z^2+2bw+beta has inverse
w=(t-nu z^2-beta)/(2b). The plane potential is

    F(t,z)=z^2+e z+(t-nu z^2-beta)^2/(4b^2)+H(t).        (5)

Its equation is degree four in z with constant nonzero leading
coefficient, so t is finite free of degree four. Its polar is

    F_z=2z+e-nu z(t-nu z^2-beta)/b^2.

If e!=0, it is parametrized by z in Gm as

    t=nu z^2+beta+2b^2/nu+e b^2/(nu z).

On that polar w=b/nu+e b/(2nu z). For d>=3, F has exact pole
orders 2d at z=infinity and d at z=0. For d=2 the infinity
order is four and the zero order is two precisely when
1+4b^2 h_2!=0. The latter leading coefficient is
e^2 b^2(1+4b^2 h_2)/(4nu^2). Thus a generic value cuts the
polar transversely in exactly 3d points.

If e=0, the polar is the disjoint generic-value intersection
of the line z=0 and the parabola
t=nu z^2+beta+2b^2/nu. Their restrictions have exact degrees
d and 2d under the stated nonresonance condition. Their one
intersection and their finitely many higher-ramification values
are avoided by generic kappa. These are again exactly 3d
simple finite ramifications. In both cases F_zz is not
identically zero on any polar component, so ordinary polar
intersections give ramification index two. The finitely many
critical values of the polar restrictions and the singular
level values can be excluded directly; no Bertini theorem
or connectedness premise is needed.

For d>=3 put t=tau^(-4), z=tau^(-d)Z at infinity. The leading
equation is

    nu^2 Z^4+4b^2 h_d=0.

Its four roots are simple and nonzero. Every other term has
smaller pole order, including z^2 t since 2d+4<4d.
The fourth-root deck action has gcd(4,d) orbits. Hence r is
the actual puncture count and infinity ramification is 4-r.
Degree four accounts for every branch. Componentwise
Riemann--Hurwitz gives (3).

For d=2 instead put t=tau^(-2),z=tau^(-1)Z. The leading
equation is

    nu^2 Z^4-2nu Z^2+1+4b^2 h_2=0.

As a quadratic in Z^2 its discriminant is -16nu^2 b^2 h_2,
which is nonzero. Its constant is nonzero by nonresonance.
Thus it has four simple nonzero roots, paired by Z->-Z:
two true punctures, infinity ramification two. The six
finite ramifications give G=c_C>=1.

For d>=3, (3) implies

    G-a_d >= d-(r+epsilon)/2 > 0.

The d=2 conclusion also gives G>a_d=0. This proves JP0,
including every affine coefficient and true infinity branch.

## 3. An unmarked isotropic eigenvector: exact fibre strata

We need one calculation that applies to all remaining cone
and resonant multiplier degeneracies. Suppose a,b are
nondegenerate ternary forms and some isotropic eigenvector
for a^(-1)b annihilates ell. A constant Witt frame writes

    a=2xy+z^2,
    b=2lambda xy+E x^2+2D xz+nu z^2,
    ell=alpha x+e z, lambda*nu!=0.                     (6)

For a scalar c, put f=a+ell+t(b-c), and let

    t0=-1/lambda,
    K=a-b/lambda+ell+c/lambda
      =-(E/lambda)x^2-(2D/lambda)xz
          +(1-nu/lambda)z^2+alpha x+e z+c/lambda.

The y coefficient is 2x(1+lambda t). Over its nonzero
locus solve y, with polynomial inverse on that open stratum.
Its zero divisor in base (x,z,t) is the union of the two
planes x=0 and t=t0, with profile 2L^2-L.

On x=0 the closed level curve is

    z^2+e z+t(nu z^2-c)=kappa.

For generic kappa its z projection is A1 minus m points,
where m=1 if c=0 and m=2 if c!=0. At t=t0 the full closed
curve is K=kappa. Their intersection consists of n points,
where n=2 if nu!=lambda, n=1 if nu=lambda,e!=0, and n=0
if nu=lambda,e=0. These are distinct-point counts for
generic kappa, not intersection multiplicities.

Every closed point carries the free y line. Additivity gives
the exact whole generic profile

    E(f=kappa)=L^3-L^2+L E(K=kappa)-(m+n-1)L.           (7)

Write r_K for the rank of the binary quadratic part of K.
The profiles of K=kappa are L-1 for r_K=2, L for a nonzero
linear force on the quadratic kernel (including rank zero),
2L for rank one without kernel force, and 0 for rank zero
without force. Thus the latter two cases have a nonzero L^2
coefficient in (7). The first two require respectively

    N=m+n,                 N=m+n-1                    (8)

under hypothetical nonautomorphy.

Here are complete critical counts for the remaining two cases.
Every critical point satisfies b=c and x(1+lambda t)=0.
If c=0,e=0, then x=z=0 and, for every t!=t0, the x derivative
fixes y. This gives a full critical curve and forbids unit Hessian.
If c=0,e!=0, x=0 is impossible at a critical point.

If c!=0, an x=0 point has z=+-sqrt(c/nu). The z derivative
fixes t, and for t!=t0 the x derivative fixes y. There are
at most two such regular points. For t=t0,x=0 either the
remaining derivative is inconsistent or y is wholly free,
which is forbidden under the Morse assumption.

At t=t0,x!=0, the x,z equations are exactly dK=0. A rank-two
K has at most one such point, and b=c fixes y uniquely.
A forced K has no such point. No singular group is omitted.

Consequently if nu!=lambda, n=2, (8) is impossible: at c=0
there are at most 1 or 0 points, versus required 3 or 2;
at c!=0 there are at most 3 or 2, versus required 4 or 3.

If nu=lambda and e=0, a rank-two K has D!=0 and its unique
critical point has x=0. All possible regular x=0 points
also have t=t0. Thus there is no finite Morse point in this
case, versus required m>0. For a forced K, necessarily
rank zero with alpha!=0; when c!=0 it again has no critical
point versus required one, and c=0 was already a critical
curve. Every unforced lower rank has the L^2 obstruction.

If nu=lambda,e!=0 and c=0, n=1. Rank two has at most one
point versus required two; a forced K has no point versus
required one. Thus every unmarked cone case is inverse-closed.

For c!=0,nu=lambda,e!=0 the counts may pass. The only
possible ternary primary types are exactly scalar,
length-two Jordan plus one line, or length-three Jordan.
The JD constant primary frames give precisely (1).
We retain these families, rather than declaring a failed
obstruction to be a theorem or a new endpoint.

## 4. Every ternary-cone multiplier is inverse-closed (JP1)

For (4), ell=0 gives the full critical line z=0 and is
already impossible under unit Hessian. Suppose ell!=0.
If the JG0 geometric hypotheses hold, its actual conic
double cover has genus one and proves inverse closure.

Failure of the multiple-zero hypothesis gives an isotropic
common eigenvector annihilating ell, exactly the unmarked
case of Section 3 with c=0. Failure of the other hypothesis
also gives one. Indeed a_D=gamma ell_D^2 on the conic D={b=0}
implies a-gamma ell^2=mu b. Nondegeneracy of a makes mu!=0.
On the two-dimensional space ker ell choose nonzero p with
b(p)=0, possible over C. Then a(p)=0 and A p=mu B p,
so p is the required isotropic eigenvector and ell(p)=0.
This includes gamma=0. Section 3 closes it as well.

This exhausts all affine forces and all ternary Jordan types.
No cyclicity or generic-transversality assumption survives
in the conclusion. The exceptional cone and its vertex
are retained, and no arbitrary double-cover multiplicativity
or source-stability assertion is used. This proves JP1.

## 5. The all-degree simple-zero quartic outside quadratic resonance

By JO, beyond the four fixed tests normalize the potential as

    h=w^2+q_A3(z)+ell(z)+H(q_B3(z)+2bw+beta), b!=0,      (9)

where both ternary quadratics are nondegenerate. Write
S=A3^(-1)B3. The relative full operator is S direct-sum 0.
The stationary rational functions retain

    deg P=3, R(s)=beta-2b^2 s+O(s^(-2)),
    F=s-H'(R), W=1/(PF).                               (10)

If d>=3, W has decay at least order five. For a nonresonant
quadratic, W has decay order four. For a leading-resonant
quadratic with gamma=H'(beta)!=0 it has decay order three.
Only the double quadratic resonance can fail proper infinity
decay. The finite group analysis below never uses a compatible
multiple determinant residue as a finite critical mass.

Suppose first that S has no isotropic eigenvector annihilating
ell. Every finite stationary group satisfies the JH0 certificate.
Here is the fixed three-dimensional proof, independent of d.
A simple primary is active with a double R pole or is compatible
simple and counted in full. A repeated scalar two-primary has
nonisotropic force: otherwise its orthogonal force complement
contains an unmarked isotropic eigenvector. Its active budget
2(d-1)>=2 suffices. A scalar three-primary always has an
orthogonal isotropic eigenvector and therefore cannot occur here.
A length-three Jordan block is marked, has R pole order four
and budget 4(d-1)>=3. A length-two block plus a line is marked;
with a distinct eigenvalue the budgets are 3(d-1)>=2 and the
simple-line alternative; with the same eigenvalue its marked
cyclic plane has nondegenerate pairing and pole order three,
with 3(d-1)>=3. These exhaust all partitions in dimension three.
Together with (10), JH0 excludes unit Hessian unless the double
infinity resonance holds. This is a direct finite-boundary
check; JH1's invertibility assumption is not applied to the
singular full four-variable operator.

Now suppose an isotropic eigenvector annihilating ell exists,
with eigenvalue lambda!=0. Its constant Witt frame has inner

    g=2lambda xy+E x^2+2D xz+nu z^2+2bw+beta,
    q_A+ell=2xy+z^2+w^2+alpha x+e z.

The JL eigenvector surgery remains literal for this inner:
on x!=0 solve y from t=g. Its ternary quadric Y is

    E x^2+2D xz+(nu-lambda)z^2-lambda w^2
       -lambda alpha x-lambda e z+2bw+beta=P_kappa(t),
    P_kappa=t+lambda H(t)-lambda kappa.

At x=0 retain its whole binary cut Y0 and the actual level
f0=kappa from (2), with its free y line. Therefore

    E(h=kappa)=E(Y)-E(Y0)+L E(f0=kappa).                (11)

Since the w square has nonzero coefficient, Y0 has quadratic
rank one or two. Its only possible negative genus term is
-a_d L U, U=u+v. Every genus term of Y has a nonpositive
coefficient. The new plane has genus G>a_d by JP0, so the
L U coefficient in (11) remains strictly negative. This
excludes HG0's generic L^3-NL for all d>=3 and every
nonresonant quadratic. No invalid zero-eigenvalue substitution
into a nonzero-plane theorem or bare stabilization is used.

Thus the only remaining outer polynomials are quadratics with

    H(r)=H(beta)-(r-beta)^2/(4b^2)+gamma(r-beta).        (12)

## 6. Both quadratic infinity resonances at once

For (12), the literal global shear and its inverse are

    w=-b t-q_B3(z)/(4b),
    t=-(w+q_B3(z)/(4b))/b.

After the constant translation tau=t+gamma/2 the potential is

    q_A3(z)+ell(z)+tau(q_B3(z)-c)+constant,
    c=2b^2 gamma.                                     (13)

For gamma=0 this is JP1 and is inverse-closed. For gamma!=0,
if there is no unmarked isotropic eigenvector, the complete
finite JH0 certificate already closes it by (10). Otherwise
Section 3 applies to (13). All its cases are inverse-closed
except c!=0,nu=lambda,e!=0, and those have exactly the three
constant Jordan frames (1). This proves sufficiency in JP2.

All reverse inclusions are actual source-orbit inclusions,
not additional model assumptions. Given a member of (1),
choose b=1,gamma=c/2 and H(r)=-r^2/4+gamma r. Reverse (13)
to the potential w^2+q+ell+H(q_B+2w). Complete the nondegenerate
q+ell by a constant translation of z. This makes it a genuine
q_A+H(q_B+linear+constant) in four variables. Both q and q_B
are nondegenerate in the ternary part for every lambda!=0.
The old four tests already occur in the retained affine family.
Thus the universal reduction is exact. All maps and inverses
are polynomial with nonzero constant Jacobian. The literal
[HP0](QUADRATIC_ISOMETRY_SATURATION.md) and
[HM0](AMBIENT_QUADRIC_GRAPHS.md) transfers remain available.

## 7. Passing controls and the precise new limit

It is essential not to overstate the result. In (1), for all
three N choices there are two regular native critical points

    x=0, z=+-sqrt(c/lambda),
    t=-1/lambda-e/(2lambda z),

with y fixed uniquely. Their critical values are
c/lambda+-e sqrt(c/lambda), and both native Hessian
determinants are 4e^2 lambda^2. For N=0 or x^2 there
are exactly these two points, all Morse.

For N=2xz there is one further point

    t=-1/lambda, x=e lambda/2, z=alpha lambda/2,
    y=(4c-alpha^2 lambda^3-2alpha e lambda^2)
                                      /(4e lambda^2).

It is Morse, has the same Hessian determinant 4e^2 lambda^2,
and critical value c/lambda+alpha e lambda/2.
It stays a distinct critical point when values coincide.

The whole-value profiles pass HG0 as well, not just generically.
For N=0 or x^2 they are L^3-2L+n_kappa L^2; for N=2xz
they are L^3-3L+n_kappa L^2. To see all jumps directly in
Section 3's two-plane strata, C1 gains one free t line
at each of its two exceptional z values. This adds L^2
per regular critical point. In the length-three case the
rank-two K level degenerates from a hyperbola (L-1) to
two crossing lines (2L-1) at its unique critical value,
adding a further L^2. Its intersection C12 is always one
point because e!=0, including coincident-value cases.
For scalar/length-two K is always a forced affine line,
so there is no further jump. No source-orbit equivalence
with the old two-point models follows from these matches.

For example lambda=c=e=1,alpha=2rho,N=2xz gives, after
subtracting one from the potential, three critical values
-1,1,rho and common native determinant four. Their arbitrary
relative value parameter shows why a finite degree check
cannot settle these full source orbits. Full-source alignment,
an invariant stronger than these passing profiles, or actual
global embedding recovery is still required.

We rechecked the retained primary
[Riemann--Hurwitz formula](https://stacks.math.columbia.edu/tag/0C1B)
and [Soprunov residue context](https://arxiv.org/abs/math/0511684v2).
The related [Asok--Dubouloz--Ostvaer record](https://arxiv.org/abs/2112.08241v2)
is context only. No new cancellation, quadric classification,
generic connectedness or multivariable residue theorem is imported.
The finite checker controls the written symbolic proof; it is
not independent replay, formal verification or external review.
