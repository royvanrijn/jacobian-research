# Properness and relative Morse points: the exact HC4 boundary

14 November 2026. Properness and vertical nondegeneracy together close
a general HC4 relative-critical branch. Each condition separately still
contains the full plane Jacobian conjecture. The latter assertion is an
exact all-degree reduction, not a constructed counterexample.

## 1. One proper Morse component is enough (HC4MYHS0)

Let f in C[s,t,u,v] have nonzero constant full Hessian determinant. Put
E=span(partial_u,partial_v), Z=V(f_u,f_v) as a scheme, and let
pi:Z->A2_(s,t) be the retained projection. Write H_E for the vertical
2 by 2 Hessian block.

If **one nonempty irreducible component** Z0 of Z is proper over the
base and det H_E never vanishes on Z0, then grad f is a polynomial
automorphism. No condition is imposed on the other components. Indeed,
it already suffices that one component projects isomorphically to the
base. This strengthens the whole-scheme section criterion HC4MYHR0.

The scheme Z is smooth of pure dimension two whenever nonempty, since
its two Jacobian rows are independent rows of the invertible full Hessian.
Its irreducible components are therefore disjoint and open as well as
closed. At a point of Z, the tangent kernel of pi is the kernel of H_E.
Thus invertibility of H_E on Z0 makes pi|Z0 etale. Both source and base
are affine, so properness makes this map finite. A nonempty finite
etale map to connected A2 is surjective. Its complex points give a
finite covering of C2; C2 is contractible, so every connected covering
has degree one. The finite degree-one map to the normal base is an
algebraic isomorphism. Only this elementary direction of finite-cover
comparison is used; the construction of algebraic structures on arbitrary
topological covers is not needed.

It remains to justify the one-component section strengthening. Assume
nonautomorphy and use Tang's plane W. In projection rank zero, W=E,
and all components of Z project to points, never to base sections.
In rank two, Z is the irreducible polynomial graph w=P(z) from
HC4MYHR0. A component section is the whole section, so that retained
theorem closes the case by its binary envelope and known HC2 input.

In rank one use the retained normal form

\[
 f=D(t,U)+sA(t,U)+V C(t,U),\qquad J(A,C)\in C^*.
\]

The section component has U=h(s,t) with C(t,h)=0. Algebraicity over
C(t) forces h independent of s. Consequently the closed graph
U=h(t), an actual A1, is an irreducible component of C=0. The retained
affine-line Keller-component theorem HC4MYGR0 makes (A,C) invertible,
contradicting nonautomorphy. This exhausts all projection ranks.

In particular, if Z is nonempty, pi is proper and all its relative
critical points are Morse, then grad f is invertible. Moreover pi must
have degree one: gradient invertibility identifies Z with the affine
two-plane of target gradients with these two components zero, so Z is
connected; its finite etale covering degree is one. A putative proper
relative Morse covering of degree greater than one therefore forbids
unit full Hessian outright.

For a hypothetical nonautomorphic unit-Hessian f, every nonempty
relative-critical component over every constant two-plane must thus
either fail properness over its retained base or contain a point where
the vertical Hessian degenerates. This statement needs no ambient graph
embedding or cotangent alignment.

## 2. The everywhere-Morse case still contains every plane Keller map (HC4MYHS1)

Let K=(K1,K2):A2_z->A2 be any plane Keller map, and write J=DK.
In actual source coordinates (w1,w2,z1,z2) set

\[
                f_M(w,z)=w^T K(z)-\tfrac12 K(z)^T K(z).       \tag{1}
\]

Its full Hessian determinant is (det J)^2. Its gradient is invertible
if and only if K is, by the retained cotangent factorization. Take the
relative free directions to be z, with w the retained base. The two
relative equations and the vertical Hessian on their scheme are

\[
 \nabla_z f_M=J^T(w-K),\qquad Z_M=\{w=K(z)\}\simeq A2_z,
 \qquad H_E|_{Z_M}=-J^T J.                                 \tag{2}
\]

Every relative critical point is Morse, and projection to w is literally
K. It has finite reduced fibres and is etale, but is a finite morphism
only if K is proper. One must not replace finite fibres by finiteness
of the morphism. Projection is proper exactly when K is an automorphism:
a proper plane Keller map is finite etale over simply connected C2.
For a hypothetical nonautomorphic K, its omitted normalization lengths
and relative fibre count defects are exactly those already retained in
HC4MYOR2, now at the relative critical projection itself.

It follows that the assertion

> Every four-variable unit-Hessian polynomial having a nonempty relative
> critical scheme with everywhere nondegenerate vertical Hessian has an
> invertible gradient

is equivalent to **JC2**, not a weaker local substitute. JC2 implies
HC4 by retained Tang; conversely (1) applies the assertion to every K.
No example of a nonautomorphic plane Keller pair is supplied.

The full critical scheme of (1) is {w=0,K(z)=0}, and all its critical
values are zero. Thus even this common-value specialization does not
remove the plane problem. These potentials are polynomial source
substitutions of a quadratic using K; they are not asserted to have a
closed A4-in-A5 representation compatible with the fixed ambient Q.

## 3. An explicit finite free map from every unimodular matrix (HC4MYHS2)

The following general algebraic lemma shows why properness alone also
retains the full problem. Let k be any field, n>=1, and
M(z) in GL_n(k[z1,...,zn]). Write B=M inverse and take a positive integer
N strictly larger than every total degree of an entry of B. Define

\[
                   P(z)=M(z)(z_1^N,\ldots,z_n^N)^T.         \tag{3}
\]

Then P is **finite free of degree N^n**, and

\[
 \{z_1^{a_1}\cdots z_n^{a_n}:0\le a_i<N\}
\]

is an actual k[P]-module basis of k[z]. This is a global module identity,
not a generic degree count or a numerical properness estimate.

To prove it, introduce independent target variables p_i. The graph ideal
(p-P(z)) equals the ideal generated by

\[
        g_i=z_i^N-\sum_j B_{ij}(z)p_j,\qquad i=1,\ldots,n,  \tag{4}
\]

because B is polynomial and unimodular. Work over k[p], ordering the z
monomials by total degree with any monomial tie-breaker. Each g_i is
monic with leading monomial z_i^N; all other z degrees are less than N.
The leading monomials are pairwise coprime. Their overlap reductions
commute: replacing z_i^N and z_j^N in either order gives the same product
of their two lower right sides. Equivalently the pairwise monic
S-polynomials reduce to zero. Monic division terminates because z degree
strictly decreases, and it gives the unique normal monomials displayed
above. This proves the free basis over k[p]. Eliminating p identifies
the quotient with k[z], proves injectivity k[p]->k[z], and proves (3).
No flatness theorem or affine-space cancellation theorem is required.

The strict degree bound is necessary for a uniform statement. For n=2
take B=[[1,z1^N],[0,1]] and M=B inverse. At equality N=deg B, (3)
becomes (z1^N*(1-z2^N),z2^N), whose fibre over (0,1) contains the
entire line z2=1. The matrix is still unimodular, but the map is not finite.

Apply this with n=2, M=J^{-T}, B=J^T for any complex plane Keller K.
Choose N>=2 with N>deg J, and put

\[
 D_N(z)=-\frac{z_1^{N+1}+z_2^{N+1}}{N+1},\qquad
 f_N(w,z)=D_N(z)+w^T K(z).                                \tag{5}
\]

The source Hessian determinant is again (det J)^2, and gradient
invertibility is again exactly invertibility of K. Its relative critical
scheme is the polynomial graph

\[
                 w=P(z)=J^{-T}(z)(z_1^N,z_2^N)^T.          \tag{6}
\]

Projection to w is finite free of degree N^2 with the explicit monomial
basis, so every fibre has scheme length N^2. Its zero fibre is exactly
Spec C[z]/(z1^N,z2^N), because the graph equations differ by a unimodular
matrix. In particular DP(0)=0 for N>=2, and the vertical Hessian on
this relative critical point is zero. The projection is proper but
ramified, while the full four-variable Hessian remains a unit.

Therefore the assertion

> Every four-variable unit-Hessian polynomial with a nonempty finite
> relative critical projection has an invertible gradient

is also equivalent to JC2, even restricting to the explicitly finite
free family (5). A hypothetical nonautomorphic K gives a conditional
member of that family; no such K is constructed. The finite degree N^2
here is the degree of P, not the generic degree of K. Properness of P
does not establish properness of K.

## Consequence for the remaining graph programme

The component criterion is a usable all-degree obstruction for any
candidate source potential. But proving only everywhere-Morse relative
criticality, or proving only a finite free relative projection, does not
reduce the full plane conjecture: both already admit universal encodings
of every Keller pair. One needs the two properties on the same component,
or additional information from the genuine ambient embedding. Neither
universal encoding supplies that embedding automatically.

The [checker](verify_proper_morse_boundary.py) and
[receipt](receipt-proper-morse-boundary.json) keep exact identities and
small monic-division controls separate from the written arbitrary-n,N
free-basis lemma and the geometric component proof. Basic primary
sources are retained in the [import record](sources/proper-morse-boundary-import.json):
finite equals affine proper, and finite etale complex covers give finite
topological covers. The existing Tang, HC2/section, affine-line component
and cotangent inputs remain explicit.

The two plane endpoints, uncovered embeddings and unrestricted HC4/JC2
remain open. No independent replay, formal verification, external review,
nonautomorphic plane example or novelty priority is claimed.
