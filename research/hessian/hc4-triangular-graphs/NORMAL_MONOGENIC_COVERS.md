# Normal monogenic covers and all monic two-parameter output equations

14 November 2026. A resultant and inertia argument removes the
one-parameter monodromy restriction for monic Danielewski equations.
The underlying statement concerns arbitrary normal finite monogenic
covers of complex affine space. Root transitivity is enough; preservation
of the entire splitting-field group is not asserted.

## 1. Normal finite monogenic covers stay integral (HC4MYHY0)

Let A=C[t1,...,tn], B=C[z1,...,zn], n>=1, and give B its
A-algebra structure through a Keller map K. Let

\[
       R(t,Z)\in A[Z]\quad\text{be monic and irreducible},
       \qquad X=Spec(A[Z]/(R))\quad\text{normal}.           \tag{1}
\]

Then R(K(z),Z) is irreducible over L=Frac(B). Equivalently this
finite cover stays integral under Keller base change. All polynomial
degrees and all n are unrestricted. The cover in (1) may be singular;
normality, global monicity and a polynomial source with constant units
are the hypotheses used below.

The [primary import record](sources/monic-cover-import.json) separates
the standard normality, ramification, purity and finite-cover inputs
from the written resultant argument. In particular, a Keller morphism
is etale, its base changes are etale, and etale pullback preserves
normality by [Stacks 025P](https://stacks.math.columbia.edu/tag/025P).
Thus

\[
                         X_B=Spec(B[Z]/(R(K,Z)))           \tag{2}
\]

is normal, though it has not yet been proved connected.

### A proper factor would have a constant cross-resultant

Suppose R(K,Z)=a(Z)c(Z) is a proper monic factorization over L.
Both factors lie in B[Z]: every coefficient is integral over B,
being a symmetric polynomial in integral roots, and B is integrally
closed. They have positive degrees. Characteristic zero and the
dominant target substitution make the polynomial separable over L,
so its generic factors are distinct.

A Noetherian normal scheme has disjoint irreducible components:
each local ring is a domain, and there are finitely many components.
Hence the two closed subsets a=0 and c=0 in (2) are disjoint.
Because the factors remain monic of the same degrees at every
source point, their resultant cannot vanish at any source point:
a zero would give a common complex root and an intersection. Thus

\[
                         Res_Z(a,c)=r\in C^*.             \tag{3}
\]

This uses all source points and polynomial source units. A resultant
nonzero only at the generic point, or a unit after inverting a target
polynomial, would not give (3).

### Integral roots prevent every cross-partition inertia collision

Put F=Frac(A) and embed L and a finite splitting field M of R over
F into a common algebraic closure. Label its distinct roots alpha_i.
The factor a selects a nonempty proper subset S of the roots; its
coefficients, and those of c, belong to M intersection L. Let
G=Gal(M/F) and G_S the subgroup preserving this subset.
Equation (3) is the field identity, up to the fixed resultant sign,

\[
                \prod_{i\in S,\ j\notin S}
                        (alpha_i-alpha_j)=r.              \tag{4}
\]

Every Galois conjugate subset sigma(S) has the same nonzero constant
cross-resultant, since sigma fixes C.

Take an arbitrary target prime divisor p and a prime q above it in
the finite normalization of A in M. The local ring at q is a DVR.
All alpha_i are integral there by the global monic equation (1), so
all differences in (4) have nonnegative valuation. If an inertia
element sent a root in S to a root outside S, those two roots would
have the same residue at q. One cross-difference would then have
positive valuation. Their product could not be the nonzero constant
r, a contradiction. Apply this to every conjugate subset as well.
Thus every such inertia group is contained in

\[
                        N=\bigcap_{sigma\in G}
                              sigma G_S sigma^{-1}.       \tag{5}
\]

This is the normal core. The argument includes all primes above all
finite target divisors, arbitrary discriminant multiplicities and
nontrivial residue field extensions. It only uses the defining fact
that inertia acts trivially on the residue field. It assumes neither
simple branch points nor transposition inertia.

### The normal quotient is trivial by purity

Let Y be the normalization of A^n in M^N. It is finite, integral
and normal. By (5), its Galois inertia is trivial at every target
codimension-one point. The ramification statements in
[Stacks 58.13](https://stacks.math.columbia.edu/tag/0BSD), followed by
[purity of branch locus, 0BMB](https://stacks.math.columbia.edu/tag/0BMB),
make Y->A^n finite etale. Here the base is regular, the cover is
finite hence quasi-finite, its local dimensions agree with those of
the base, and normality is built into the normalization. These are
the required purity hypotheses. Equivalently use the codimension-one
inertia formulation in Stacks Example 58.13.3.

A connected finite etale cover of complex A^n has degree one:
its analytification is a finite topological covering of contractible
C^n. Only the retained forward finite-etale-to-covering comparison
from Milne is used, not essential surjectivity of Riemann existence.
Thus M^N=F and N=G. Equation (5) then makes S G-invariant. Since
R was irreducible over F, G acts transitively on its roots, so a
nonempty proper invariant subset is impossible. This proves HC4MYHY0.

This result is irreducibility preservation for the given normal
monogenic cover. It does not make L/F trivial, or assert that the
source has the full Galois group of the splitting field M. The root
field stays linearly disjoint in the degree sense; the splitting field
need not be claimed linearly disjoint by this proof.

## 2. Every monic two-parameter Danielewski pair closes (HC4MYHY1)

Let P(T,z) have arbitrary coefficients in C[T], positive z degree
m and nonzero constant leading coefficient. Normalize that coefficient
to one by constant ambient and target rescalings. For r>=1, in
A^(r+3) with unused affine variables allowed, put

\[
                            Pi=(T,xy-P(T,z)).              \tag{6}
\]

If a polynomial closed stable embedding j:A^(r+2)->A^(r+3)
carries a plane Keller pair K through (6), then **K is invertible**.
Every T degree, z degree and embedding degree is unrestricted.
In particular the coefficients may genuinely depend on T, with
arbitrary nonabelian two-parameter root monodromy.

The root cover has equation P(T,Z)+V=0 over A=C[T,V]. Its total
space is A^2_(T,Z), by eliminating V; it is smooth, normal, connected
and finite over the target because P is monic. HC4MYHY0 makes

\[
                           P(K1,Z)+K2                     \tag{7}
\]

irreducible over L=C(source base). Thus the source stabilizer H
acts transitively on its m generic roots. No full-group equality is
needed.

The rest of the ambient algebra is explicit. Assign x,y weight m,
z weight two, and unused variables weight one; A has degree zero.
The monic relation xy=P(T,z)+V supplies a common free A-basis of
monomials not divisible by xy. Its geometric generic graded ring has
equation xy=z^m, a domain. Every affine height-one target fibre is
geometrically integral: the right side P(T,z)+V retains its positive
z degree and constant leading coefficient, so is a nonzero polynomial;
the equation is primitive linear in y with coprime coefficients x
and P(T,z)+V. The geometric generic ambient fibre is smooth, with
torsion-free class group

\[
                          C=Z^m/Z(1,...,1),                \tag{8}
\]

generated by the m root divisors, with Galois action by root permutation.
This is the same localization calculation as
[HC4MYHX2](REGULAR_PARAMETER_MONODROMY.md#3-every-polynomial-danielewski-observable-is-covered-hc4myhx2).

There is no nonzero class in (8) fixed by a transitive permutation
subgroup H. Indeed represent a fixed class by an integer vector v.
For h in H, h*v-v is a constant vector k_h*(1,...,1). Summing
coordinates gives m*k_h=0, hence k_h=0. Transitivity then makes
v constant, so its class is zero.

The divisor class of the geometric source component corresponding
to the fixed embedding L is H-fixed, hence zero. All its conjugates
are therefore principal too. The retained component-pencil argument,
global primitive saturation using the integral affine height-one
fibres, and smooth minimal-label derivative argument in
[HC4MYHW1](ABELIAN_CLASSGROUP_OUTPUT_OBSTRUCTION.md#2-a-class-group-and-geometric-fibre-output-criterion-hc4myhw1)
now force [L:F]=1 and the birational Keller inverse. This establishes
(6), without abelian or one-parameter monodromy and without ambient
factoriality or a source Hessian assumption.

For actual positive controls, j(s,t,w)=(s,1,t+P(s,w),w) is
closed with inverse (T,y-P(T,z),z). An unused v gives HC4 size.
The independent nondegenerate ambient quadratic T*v+y*z+x^2 pulls
back exactly to 1+s*v+t*w+w*P(s,w), whose actual Hessian determinant
is one for every P. Its gradient inverse sends (a,b,c,d) to

\[
 (s,t,w,v)=(d,\ c-P(d,b)-b P_z(d,b),\ b,\ a-b P_T(d,b)).   \tag{8a}
\]

This is a literal compatibility certificate, not a nonautomorphic
Keller construction. The previous all-degree
T*Y+xy-B(T,z) suspension theorem remains valid independently,
including nonmonic or constant B.

## 3. A nonautomorphic Keller normalization cannot be monogenic (HC4MYHY2)

Let K be complex Keller in any dimension, A=C[K], L=C(z),
F=Frac(A), and F subset E subset L a finite intermediate field.
Let A_E be the finite normalization of A in E. If it is globally
monogenic as an A-algebra,

\[
                              A_E=A[h],                   \tag{9}
\]

then E=F. To prove this, its integral generator h belongs to the
source polynomial ring: it is integral over A and hence over C[z],
which is normal. Its minimal polynomial R over F is monic with
coefficients in A, since A is integrally closed. The kernel of
A[Z]->A[h] is exactly (R), by monic division. Thus (9) is a normal
finite monogenic cover of the target. HC4MYHY0 says R(K,Z) is
irreducible over L. But it has the root h in L, so its degree is
one. This proves E=F.

In particular, for a hypothetical nonautomorphic plane Keller pair,
the entire finite normalization N of its target in its source field
requires at least two additional global algebra generators over A.
It cannot equal A[h] for any integral h. If it did, the preceding
argument would give birationality, and the retained plane birational
Keller theorem would give the polynomial inverse.

This concerns the finite normalization, not just the open source
polynomial ring. The generic primitive-element theorem does not
supply (9). Even an integral primitive element can generate a proper
nonnormal order whose normalization needs other generators. Every
such integral primitive order in the nontrivial plane Keller field
would have to be nonnormal; no hypothesis that it can be made
globally normal monogenic is asserted.

## 4. Controls and exact limits

Normality in Section 1 need not mean smoothness. The cover
Z^2=T1*T2 is normal with a singular point: its coordinate ring is
C[a^2,ab,b^2], the invariant ring of (a,b)->(-a,-b). An element
of its fraction field integral over that ring is integral over
C[a,b], hence polynomial, and is fixed by the involution; this proves
normality. It is thus a valid singular example of (1).

In contrast the monic primitive model

\[
                    q^2=T^2(T+1)                          \tag{10}
\]

is not normal. Its normalization is C[c,V], with T=c^2-1 and
q=T*c, and c=q/T generically. The element c is integral over the
model but is not in it: the two points c=1 and c=-1 have the same
T=q=0 and different c values. The normalized map to (T,V) is etale
near those two points, while the primitive generator q collides.
Thus a generic primitive element or smooth normalization does not
make an arbitrary monic equation normal. This example is not Keller
globally: its map (c,V)->(c^2-1,V) has Jacobian 2c.

The Keller condition is also essential. Pulling the smooth normal
cover Z^2=T back along (s,t)->(s^2,t) factors it as
(Z-s)(Z+s). The resultant is proportional to s, and its zero is
exactly where the source map ceases to be etale.

If that source is instead Gm, the same square map is etale and its
resultant 2s is a unit, but a nonconstant one. Its two pulled-back
components are then disjoint. This explains why an arbitrary etale
source cannot replace the stipulated affine-space source with scalar units.

Monicity is global. For T*Z^2=V the total space is smooth, but its
projection is not finite over T=V=0. Dividing by T creates roots
of Z^2=V/T with negative valuation at T=0. The nonnegative-valuation
argument in (4) then fails. No inversion of a nonconstant leading
coefficient is allowed here.

For permutation groups, a proper cyclic subgroup of S4 can already
act transitively on four roots and have zero invariant classes in
Z^4/Z(1,1,1,1). This is why the conclusion needed for Section 2
must not be reported as preservation of the full splitting group.

The remaining unsuspended xy-P(T,z) family with a nonconstant
leading coefficient, arbitrary nonnormal monic models, normal finite
algebras without global monogenic presentations, and the other
uncovered embedding classes remain distinct. Both plane endpoints
and unrestricted HC4/JC2 remain open. The
[checker](verify_monic_covers.py) and [receipt](receipt-monic-covers.json)
give finite algebraic controls alongside the written all-degree proof
and explicit primary imports. No formal verification, independent
replay, external review or novelty priority is asserted.
