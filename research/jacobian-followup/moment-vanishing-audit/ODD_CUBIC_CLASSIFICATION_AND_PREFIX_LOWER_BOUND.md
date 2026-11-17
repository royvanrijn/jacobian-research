# Odd ternary cubics: classification, sharp mixed bound and prefix obstructions

Current quantitative refinement, 15 November 2026: the [exact Euler representative proof](EULER_REPRESENTATIVES_GLOBAL_297.md) gives global pure297, odd296, conic/full-W-linear258 and all-homogeneous258. Ordinary-node complements are297/271/244, with homogeneous nonisotropic-node101; cusp255/tangent284 are unchanged. Current intervals are global[20,297], odd[20,296] and all-homogeneous[14,258]. Earlier353/336 bounds below are historical valid coarser evidence. Leastness,20/14 sufficiency, expanded I297 cofactors and independent review remain open.

Current quantitative refinement, 15 November 2026: the [joint-degree actual-seed proof](JOINT_DEGREE_GLOBAL_353.md) gives global pure353, odd352, full W-linear336, rank-one140, mixed-parent255 and cusp255/tangent284. Current least-prefix intervals are global[20,353], odd[20,352] and all-homogeneous[14,336]. Earlier590/888 bounds and rank-one425/mixed-parent320 estimates below are historical valid coarser evidence. Leastness,20/14 sufficiency, expanded I353 cofactors and independent review remain open.

Current quantitative refinement, 15 November 2026: the [balanced logarithmic seed proof](BALANCED_LOGARITHMIC_GLOBAL_590.md) gives global pure590, all-homogeneous336, conic336 and full W-linear425. Current intervals are global/odd[20,590] and all-homogeneous[14,336]; smooth[20,94] and homogeneous-smooth[14,56] are unchanged. Earlier888/1855 bounds and frame-plus888 statements below retain their historical scopes. Leastness,20/14 sufficiency, expanded I590 cofactors and independent review remain open.

Current global refinement,15 November2026: the [rescaled singular seed proof](RESCALED_SINGULAR_GLOBAL_888.md)
now proves888 sufficient as a global pure-only prefix, without an extra
frame condition. Current intervals are global[20,888], odd[20,888] and
all-homogeneous[14,888]. Earlier1855 bounds and reservations about888
alone below are historical. The new proof remains internal and
independently unreviewed; leastness and expanded global cofactors are open.

15 November 2026. Internal written consequences of the retained complete
ternary cubic classification. The exact [checker](verify_odd_cubic_invariants.py)
and [receipt](odd_cubic_invariants_verification.json) establish the auxiliary
tensor, Jacobian and Wick identities below. They do not independently verify
the retained analytic theorems or replace the invariant-ring argument.
No external review, independent specialist replay or formal verification
is claimed. The [source audit](odd_cubic_invariants_source_audit.json)
records the precise import and the exploratory computation's boundary.

The new conclusions are:

- Every odd pure-null ternary polynomial of degree at most three has the
  positive circular-weight normal form (1) below. Its mixed cutoff is
  **m>deg Q**, sharply, including within genuine cubics. Its pure-null
  locus is irreducible of dimension six in the thirteen-dimensional odd
  coefficient space.
- There exists an odd ternary cubic over the algebraic numbers with its
  first **nineteen** pure Gaussian moments zero but not all moments zero.
  Thus the least global pure-only prefix is now in **[20,1855]**.
- There exists a homogeneous ternary cubic over the algebraic numbers
  with its first **thirteen** pure moments zero but not all moments zero.
  Thus its homogeneous least prefix is in **[14,1854]**.

The invariant arguments below are existence proofs. The subsequent
[certified smooth witness](SMOOTH_CUBIC_FOURTEENTH_MOMENT_WITNESS.md) now
supplies an isolated algebraic coefficient representation for the
homogeneous13-zero example, and proves its14th moment nonzero. The subsequent
[odd19-zero witness](ODD_CUBIC_TWENTIETH_MOMENT_WITNESS.md) now supplies
its isolated algebraic coefficients, nonzero20th and smooth leading form. Neither
argument proves sufficiency of twenty or fourteen moments. The upper bounds use the retained
[global1855 theorem](GLOBAL_TERNARY_CUBIC_PURE_PREFIX.md), with parity
removing its last odd moment on these two subclasses. The odd-only least
prefix is likewise in[20,1854]. The global least-prefix question stays open.

The subsequent [graded seed proof](GRADED_SEED_SMOOTH_CUBIC_PREFIXES.md)
now gives smooth-leading interval[20,94] and homogeneous-smooth[14,56].
The all-homogeneous interval[14,1854] is unchanged.

## Complete odd pure-null classification

Use the bilinear Gaussian covariance, with no complex conjugation. In
circular coordinates E(ZW)=E(T^2)=1 and every other coordinate pairing
zero. Let P(-x)=-P(x), deg P<=3. All positive pure moments vanish if and
only if some complex orthogonal circular frame gives

    P=Z*(a Z^2+b ZT+c T^2+d ZW+e).              (1)

Here all five coefficients are arbitrary, including zero. Equivalently,
apart from the zero polynomial which is included in (1),

    P=(ell.x)*(x^t A x+s),
    ell.ell=0, ell!=0, A=A^t, A ell is proportional to ell. (2)

The same isotropic line divides both the cubic and the linear part;
classifying only the homogeneous cubic part would not prove this.

For necessity, the retained [full cubic theorem](TERNARY_CUBIC_GMC_CLOSURE.md)
forces a W-linear frame. Oddness then leaves exactly

    P=W*(a+d Z^2+e ZT+f T^2)+G3(Z,T)+r Z+s T.  (3)

The rank-two coefficient case e^2-4df!=0 is excluded by the retained
W-linear conic theorem in the [rank-one/full-frame proof](RANK_ONE_WT2_CLASSIFICATION.md).
If e^2-4df=0 and f!=0, a linear Gaussian null rotation makes the
quadratic W coefficient f T^2. It preserves oddness. The same retained
rank-one classification then forces G3+rZ+sT=0: its necessary coefficient
b of WZ is already zero by parity. Thus P=W*(a+f T^2), which is (1)
after interchanging Z and W.

If f=0, the rank-one-or-zero condition forces e=0. Then (3) is in the
[cubic-height architecture](CUBIC_HEIGHT_AIRY_OBSTRUCTION.md), with
h(Z)=a+dZ^2. If a=0, its zero-branch theorem says h is divisible by Z^2
and every other coefficient of T^3,T^2,T,1 is divisible by Z. Parity
therefore leaves exactly (1). This includes h identically zero.

If a!=0, that theorem eliminates the T^3 term and the
[affine-W classification](CUBIC_AFFINE_W_CLASSIFICATION.md) applies.
Nonconstant h must have degree at most one, so d=0. On constant h=a,
the cubic degree restriction gives v=0, deg u<=1 and r=-Z*u^2/(2a).
Parity makes u constant. What remains is an isotropic linear polynomial,
also of form (1) in a suitable frame. These cases exhaust (3).

For sufficiency, the weights wt(Z)=1,wt(W)=-1,wt(T)=0 of the five
monomials in (1) are3,2,1,1,1. Wick contractions require total weight zero.
Every positive power of P has positive weight, proving pure nullity.
Every monomial of a multiplier Q has weight at least -deg Q, so

    E(Q P^m)=0 for m>deg Q.                      (4)

This coefficient one is sharp even with a genuine cubic: take P=ZT^2
and Q=W^D, D>=1. Then

    E(W^D P^D)=D! (2D-1)!! !=0.                 (5)

For dimension and irreducibility, let the isotropic factor line range
over the smooth projective conic. For a fixed line, (1) is a vector space
of dimension five in polynomial coefficient space. These spaces form a
vector bundle over that conic: factor divisibility and the eigenvector
condition are linear, of constant rank. The total incidence variety is
irreducible of dimension six, and its image is exactly the pure-null
locus. In fact this image is closed, since projection from the incidence
subset of coefficient space times the projective conic is proper.

The image has dimension at least six. Use

    Z'=Z+lambda*T-lambda^2*W/2, T'=T-lambda*W, W'=W,
    2 Z'W'+T'^2=2ZW+T^2,

and substitute the primes in (1). At lambda=0,e!=0 the five coefficient
directions are independent multiples of Z. The lambda derivative has
the nonzero linear term e*T, independent of them. The checker also gives
an exact nonzero six-by-six coefficient-Jacobian minor, equal to2 at
(a,b,c,d,e,lambda)=(1,2,3,4,1,0). Thus the dimension is exactly six.

## Ten independent invariants of odd cubic coefficients

Write the thirteen-dimensional odd representation as

    P=H(x)+(x.x)*(c.x)+ell.x,     Delta H=0.      (6)

The harmonic cubic H has seven independent coefficients; c and ell
have three each. This is a rational linear decomposition since
Delta((x.x)*(c.x))=10(c.x). Write T_H for the symmetric trilinear tensor
polarizing H, so T_H(v,v,v)=H(v), and put w=ell cross c. The following
ten polynomials in coefficient space are O(3,C)-invariants:

    I0=ell.ell, I1=ell.c, I2=c.c,
    I3=T_H(ell,ell,ell), I4=T_H(ell,ell,c),
    I5=T_H(ell,c,c), I6=T_H(c,c,c),
    I7=T_H(ell,ell,w)^2,
    I8=T_H(ell,c,w)^2, I9=T_H(c,c,w)^2.           (7)

Tensor contractions are SO(3)-invariant. A reflection changes the sign
of the pseudovector w; the squares remove that sign. Their coefficient
degrees are2,2,2,4,4,4,4,10,10,10. The checker verifies all thirty
infinitesimal invariance identities and all ten identities for a generating
reflection, universally over the rational polynomial ring.

Here is an explicit independence certificate, with zero-based indices.
Use harmonic coordinates

    H=h0*x^3+h1*x^2*y+h2*x*y^2+h3*y^3
      +h4*x^2*z+h5*x*y*z+h6*y^2*z
      -(3h0+h2)*x*z^2-(h1+3h3)*y*z^2-(h4+h6)*z^3/3.

At ell=(1,0,0), c=(0,1,0), take the ten Jacobian columns

    ell0,ell1,c1,h0,h1,h2,h3,h4,h5,h6.

The minor is 8*h4*h5*h6/6561; at h=(1,2,3,4,5,6,7) it is560/2187.
Consequently the invariant ring B=C[H,c,ell]^O(3) contains ten
algebraically independent elements. No claim that (7) generates B is needed.

## Why nine even moments cannot force all-order vanishing

We give the ring argument explicitly to avoid confusing the dimension of
the six-dimensional pure-null locus with the dimension of the quotient.

For a characteristic-zero orthogonal representation, a degree-preserving,
B-linear Reynolds projection C[V]->B exists. The general-representation
statement is [Maithani, Theorem4.2 and Example4.5, v2](https://arxiv.org/html/2412.18841v2).
Its algorithms for standard/dual copies are not being imported for (6).
One can also obtain this projection by averaging over compact O(3,R).

The projection gives the following elementary consequences. The ideal
C[V]*B_+ has finitely many homogeneous invariant generators b1,...,bt,
by Noetherianity. Averaging coefficients in any expression for a positive
homogeneous invariant and inducting on degree proves B=C[b1,...,bt].
Thus finite generation here requires no computed invariant generating set.

Every homogeneous positive-degree invariant vanishes on (1): the
orthogonal torus Z->tau*Z,W->tau^-1*W sends P to zero as tau->0,
because all five weights are positive. Invariance and continuity give
b(P)=b(0)=0. The classification proves this for every odd pure-null P.

Suppose that M2,M4,...,M18, nine homogeneous invariant polynomials,
forced all-order vanishing in V. Every bi would then vanish on their
common zero locus. Hilbert's Nullstellensatz gives

    bi^ni in (M2,M4,...,M18)*C[V].

Applying the B-linear Reynolds projection gives the same containment
inside B. Therefore B/J is finite dimensional over C, where
J=(M2,M4,...,M18)*B: all its algebra generators are nilpotent.

Choose a finite homogeneous basis of B/J and lift it to B. These lifts
generate B as a module over A=C[M2,M4,...,M18]. To see this, reduce a
homogeneous b modulo J. The remainder is a linear combination of the
chosen lifts; the difference is a sum M_(2j)*bj with bj of strictly
smaller degree. Induction on degree finishes the claim. Thus B is finite,
hence integral, over A. Their fraction fields have the same transcendence
degree, at most nine. This contradicts the ten independent invariants (7).

So an odd P has M2=...=M18=0 but is not pure-null. Odd Gaussian moments
are identically zero, so its first nineteen moments vanish. This proves
the global lower bound twenty. The nonzero later moment is even, but
its precise order has not been found. This invariant argument alone does not prove M20!=0. The subsequent
odd witness linked above now certifies that order by a separate construction.

The example can be chosen over Qbar. The global1855 theorem gives some
even k in[20,1854] with Mk(P)!=0. For one such k the polynomial system
M2=...=M18=0, u*Mk-1=0, defined over Q, is nonempty over C. The weak
Nullstellensatz over the algebraically closed field Qbar supplies a
Qbar point. This invariant argument is not an explicit coefficient construction; the
subsequent odd witness supplies the exact isolated-root representation.

## Seven invariants on the homogeneous subspace

In the homogeneous representation H+(x.x)*(c.x), set the auxiliary
vector ell=T_H(c,c,-). The Euclidean bilinear form identifies this
covector with a vector. The ten polynomials (7) remain invariants, with

    I5=I0, I6=I1, I9=0.

The remaining seven, indexed0,1,2,3,4,7,8, are independent. At c=(0,1,0),
h=(1,2,3,4,5,6,7), ell=(1,4,7/3), the checker uses the chain rule and
the seven columns c0,c1,c2,h0,h1,h2,h3. Its exact determinant is

    1431014859795988480/387420489 !=0.

The same finite-module argument with six moment polynomials
M2,...,M12 now gives a homogeneous non-pure-null cubic whose first
thirteen moments vanish. Algebraic coefficients follow as above. It
cannot be reducible or cuspidal, since the retained sharp class gates
are ten; this proof does not decide whether its leading curve is smooth
or irreducible nodal. This invariant argument alone does not improve the separate
smooth-leading lower bound. The subsequent certified witness linked above
proves a smooth example with first nonzero order14, raising that lower
bound to14 by a separate construction.

## Exploration and remaining obligations

The retained [Hilbert-series exploration](odd_cubic_hilbert_exploration.json)
counts SO(3) symmetric-power invariants through coefficient degree140
and keeps even degrees for O(3). Multiplying this truncated series by
the ten proposed moment degrees gives no observed negative coefficient.
This does not establish an exact Hilbert series, a system of parameters,
a regular sequence, or sufficiency of the first twenty moments. No such
conclusion is used in the proofs above.

Both lower-bound coefficient witnesses and their first nonzero orders
are now supplied by the separate certified constructions. The global
least prefix and expanded global radical cofactors,
and independent review remain open. The full qualitative cubic GMC
decision stays closed internally. Exceptional coefficient schemes,
global term-count minima and lower-power GVC are separate questions.
