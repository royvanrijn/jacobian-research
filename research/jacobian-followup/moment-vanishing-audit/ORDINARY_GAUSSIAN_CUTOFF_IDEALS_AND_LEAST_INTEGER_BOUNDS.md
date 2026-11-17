# Ordinary Gaussian cutoff ideals and least integer bounds

16 November 2026. General implementation of the retained ordinary-ideal
algorithms, with an explicit rational-residue stopping oracle. Internal
proof and exact controls; independent, formal and external assurance remain
false. No novelty or practical complexity claim is made.

The [producer](gaussian_ordinary_cutoff_producer.py),
[checker](verify_gaussian_ordinary_cutoff_producer.py),
[receipt](gaussian_ordinary_cutoff_producer_verification.json) and
[source audit](gaussian_ordinary_cutoff_producer_source_audit.json) retain
the implementation and its limits. The
[previous numerical producer](GENERIC_FIELD_CERTIFICATES_AND_NUMERICAL_GAUSSIAN_FAMILY_BOUNDS.md)
supplies upper bounds; it is not modified.

## 1. Exact objects computed

Let B=QQ[c_1,...,c_r], A=B/I0, and P(c,X) be any rational polynomial
parameter family in n Gaussian variables. Precision is a fixed rational
symmetric nonsingular matrix H; covariance is C=H^-1. Keep I0 ORDINARY.
For arbitrary nonnegative integers N,E, define

    J_(N,E)=I0+(E[X^alpha P^m]: m>N|alpha|+E, m>=1) in B.

The producer computes this FULL ordinary ideal, including nilpotents and
embedded structure, not just its radical. It also computes

    J_pure=I0+(E[P^m]:m>=1).

The former full affine filtration is K_N=J_(N,N); the pure linear cutoff
ideals are J_(N,0). There is no spatial dimension, phase degree, parameter
degree, parameter count or cutoff-index bound in the interface.

With the numerical family upper certificates, it computes the LEAST
nonnegative INTEGER N for either m>N(deg Q+1) or m>N deg Q:

* as an identity in the presented coefficient ring, when safe (and pure
  for the zero-intercept mode);
* uniformly on all safe field fibers, with the pure field fibers selected
  additionally in zero-intercept mode.

An empty target has least integer0. A nonpure safe ring phase has no
zero-intercept bound; an unsafe ring phase has no finite affine bound.
These are decisions about the specified shapes, not optimization over
all rational slopes or affine pairs.

After an ideal is known, the code computes its least generating total
index prefix |alpha|+m. For pure moments this is the usual least moment
prefix. It never uses two equal consecutive prefixes as a stopping test.

## 2. Faithful finite algebras without full Noether normalization

Given a proposed proper ideal I containing I0, compute its FULL rational
primary decomposition, including embedded components:

    I=intersection Q_i, p_i=radical Q_i.

For each p-primary Q choose a maximal algebraically independent subset U
of the ORIGINAL parameter variables modulo p, with |U|=dim(B/p), and let
V be the remaining variables. Set F=QQ(U) and

    T=F[V]/QF[V].

Every nonzero polynomial in QQ[U] avoids p. Because Q is primary, it acts
injectively on B/Q. Hence B/Q -> T is injective. The fraction field of
B/p is algebraic and finitely generated over F, so (B/p) localized over
F is zero-dimensional. Therefore T is a finite-dimensional F-algebra,
possibly nonreduced. A Groebner basis over F gives all standard monomials
and the complete multiplication table.

Finiteness of B/Q over QQ[U] is NOT asserted or needed. This localization
argument avoids constructing a full integral Noether normalization.
It is the classical maximal-independent-variable step used in primary
algorithms, with the primary injectivity proof kept explicit.

The official Macaulay2 documentation specifies
[independentSets](https://www.macaulay2.com/doc/Macaulay2-1.26.05/share/doc/Macaulay2/Macaulay2Doc/html/_independent__Sets.html)
and the full, embedded-component-preserving
[primaryDecomposition](https://macaulay2.com/doc/Macaulay2/share/doc/Macaulay2/PrimaryDecomposition/html/_primary__Decomposition.html).
The installed version is1.22; its actual operations and source files are
locked separately. The code checks the intersection, prime radicals,
independent-set dimensions, finite standard basis, ordinary relations and
multiplication table. Python independently checks the ordinary intersection
and finite algebra data, but primary status/independence retain their
trusted CAS boundary. Localizing an undecomposed nonprimary quotient is
not justified: it can erase embedded torsion.

## 3. A cone filter as one formal rational-residue period

Work in one such finite algebra T. Introduce independent scalar ray
parameters a_0,...,a_n, auxiliary variables v_0,...,v_n, and time u.
Extend F faithfully to F(a). The rational power series

    R_(N,E)(t,s)=t^(E+1)/((1-t) product_i(1-s_i t^N))
               =sum_(m>N|alpha|+E) t^m s^alpha

implements the STRICT forbidden cone, including N=0. Put

    Y(u)=CT_v E[exp(v_0 P+sum_i v_i X_i)
                    R_(N,E)(a_0u/v_0,a_i u/v_i)].        (1)

The coefficient of u^k is

    sum_(m+|alpha|=k, m>N|alpha|+E)
       E[X^alpha P^m] a_0^m a^alpha/(m! alpha!).         (2)

Distinct independent a-monomials prevent cancellation. Thus Y=0 iff EVERY
forbidden moment is zero in T. Setting every a_i=1 would be unsound.

For ordinary residue measure define the scalar rational amplitude

    g=R_(N,E)(a_0u/v_0,a_i u/v_i)/product_j v_j,
    f=(product_j v_j)(v_0-a_0u)
        product_i(v_i v_0^N-a_i a_0^N u^(N+1)).           (3)

Then g is in F(a)[u,v,f^-1]. In particular f^(E+2)g is polynomial.
At u=0, f is a nonzero scalar times a v-monomial. Consequently every
coefficient in the u-expansion of any polynomial amplitude divided by
f^h is a finite Laurent polynomial in v. Multiplication by the formal
exponential in (1) and taking its v-residue selects finitely many
exponential coefficients at each u-degree. Every Gaussian contraction
is finite. This justifies the period coefficientwise, including all
primitive amplitudes used below; no convergent integral is assumed.

## 4. An executable search for a scalar IBP certificate

On the finite free T-valued connection use

    nabla_Xi=partial_Xi-(HX)_i+v_i+v_0 P_i,
    nabla_v0=partial_v0+P,
    nabla_vi=partial_vi+X_i,
    nabla_u=partial_u.

All multiplications use the retained finite algebra table. These operators
commute. The localized connection over F(a) is holonomic; algebraic direct
image in X,v is holonomic. The cyclic class of g in its top integration
quotient has a nonzero SCALAR annihilator L in F(a)[u]<D_u>. Therefore
there exist rational amplitudes U_j in T[u,X,v,f^-1] with

    Lg=sum_j nabla_j U_j.                                (4)

This imports the same classical algebraic holonomic localization and
integration theorems as the sealed
[primary-thickening proof](PRIMARY_THICKENING_MOMENT_IDEAL_ALGORITHM.md)
and [Oaku source construction](WHOLE_GAUSSIAN_SOURCE_HOM_AND_EFFECTIVE_CUTOFFS.md).
No regularity, analytic solution interpretation or actual image kernel
is required. Nilpotents in T do not affect the finite-rank integrable
connection argument. Scalar coefficients of L are essential: a nonzero
T-coefficient may be a zero divisor.

The implementation FINDS (4), rather than leaving a general creative
telescoping routine unspecified. At box d=1,2,... it allows order and
u-degree of L at most d, and writes every primitive as a polynomial vector
of total degree <=d(deg f+1) divided by f^d. Clear by f^(E+2+d). Sparse
exact column reduction over F(a) finds a dependence containing a nonzero
L-column, retaining every coefficient in the original columns.

These boxes exhaust the needed identities: if a primitive has fixed pole
h and numerator degree b, multiply its numerator by f^(d-h) for d>=h.
Its degree is b+(d-h)deg f <=d(deg f+1) for sufficiently large d. A fixed
polynomial scalar L is also included eventually. Existence of (4) proves
unbounded search termination. No heuristic rank or finite degree cutoff
is substituted for that proof. A separate finite replayer reconstructs
ONLY the retained columns and verifies the original identity directly.

Gaussian IBP and total-residue differentiation show LY=0. For
L=sum c_ab u^a D_u^b define w=min(a-b) and

    b_origin(z)=sum_(a-b=w) c_ab (z)_b.

It is nonzero. Any nonzero T-valued formal solution has initial order
among its nonnegative integer roots, because a nonzero scalar in F(a)
acts invertibly on T. Exact denominator clearing and rational content
give the finite largest root R0. Check the ORIGINAL moments in (2) for
total index <=R0. All zero proves Y=0; otherwise a nonzero forbidden
moment is an explicit ordinary nonmembership witness.

For J_pure the smaller polynomial Gaussian connection suffices. Search
L G=sum nabla_X U for G=E exp(uP), using scalar L and polynomial primitives.
If p0=L(1), then G-1 is annihilated by L when p0=0, and otherwise by
(p0 D_u-p0') composed with L. The same origin-jet test applies.

## 5. Exact shortcuts and the genuine stopping loop

Several general cases avoid the rational search without weakening it.

* If P is spatially constant c, J_(N,E)=I0+(c^(E+1)), for every N;
  J_pure=I0+(c). This preserves the full affine thickness K_N=(c^(N+1)).
* At slope N=0, the complete ideal is the coefficient ideal of P^(E+1).
  Indeed the Gaussian pairing on polynomials of bounded degree is
  nondegenerate for every nonsingular H. Its finite Wick Gram matrix
  and inverse express every coefficient as a rational combination of
  original moments with that power. Higher powers are coefficientwise
  multiples. The code retains this matrix certificate.
* If P^ell=0 in T[X], only m<ell need checking. For a finite algebra of
  dimension s, a polynomial P is nilpotent iff P^s=0: all its coefficients
  then lie in the nilradical, whose power s is zero. The code only uses
  the directly verified sufficient equality P^s=0. For N>0 the remaining
  forbidden marker set is finite; pure mode has only the unit marker.
* If all active phase variables have pairwise zero covariance and every
  phase monomial has degree at least d0>0, a nonzero Wick contraction
  requires d0*m<=|alpha|. Hence all pure moments vanish and N*d0>=1
  certifies the specified zero tail. Coefficients may still have nilpotents.

For a proposed I, test EVERY primary component faithfully as above.
A failure returns an original forbidden moment outside I; adjoin it.
Every failure strictly enlarges the ordinary ideal in Noetherian B.
ACC forces termination. At the terminal positive zero test, all forbidden
moments lie in I and every added generator is forbidden, so the ideal is
exact. This is the stopping argument already proved abstractly in
MVLINEARIDEAL20261115, now with a complete implemented oracle.

From the finite retained generator indices obtain an upper total prefix.
Finite Groebner membership and binary search give its least generating
prefix. This does not assert stabilization of the descending cutoff
filtration as N grows.

## 6. Least integer bounds on rings and safe field families

The sealed family solver gives the full mixed-safe radical J and exact
normalization. The numerical producer supplies C,E (or E_ring), with
original certificates extending to every safe component point. For a safe
presented ring, ceil(max(C,E_ring)) is an affine-index upper bound.
If pure, ceil(C+E_ring) is a linear-index upper bound.

To find the least integer, test every N from0 up to that finite upper
bound using the exact ordinary zero oracle against I0. A failed test
retains a forbidden original moment, proving that smaller N fails.
The first successful test proves the bound for ALL markers and powers.

For all safe FIELD fibers use target J, and upper field constants C,E.
For their pure subset compute the full ordinary J_pure and use
radical(J+J_pure) as target. The full primary intersection is retained
when forming that radical. Use ceil(C+E) for this pure field target.
Thus least common integer bounds on the full bounded universal coefficient
space are computable by the implemented composition, without asserting
that large universal inputs have been executed.

All target ideals remain ordinary except where the task explicitly asks
for field points. In particular residue-field purity is not used to infer
ring purity. The retained P=eX over QQ[e]/(e^3) has affine least index2,
while its nonzero pure second moment proves NO zero-intercept linear bound.

## 7. Executed controls and the capped search

The focused controls retain complete affine-linear ordinary ideals,
least generating prefixes, a noncoordinate nilpotent primary algebra,
embedded-component detection, isotropic support and zero-ring behavior.
Least-bound controls include both positive answers and the purity/safety
obstructions. The pure scalar-IBP search executes on X^2 and detects a
nonzero actual moment through its certified jet.

The unrestricted cone-IBP ENUMERATION on P=0 in one variable reached its
20second cap and remains UNKNOWN, terminal, not restarted. No historical
source snapshot was taken at that cap; its retained record explicitly
says so. It is not used as a successful computation or a complexity claim.

A separate directly derived rational certificate validates the general
cone representation and identity replayer. At N=E=0, P=0 and n=1 put

    g=a0*u/[v0*(v0-a0*u)*(v1-a1*u)],
    h=1/(v0-a0*u)-1/v0,
    U0=-a0/[(v0-a0*u)(v1-a1*u)].

For L=D_u^2-a1^2*u*D_u-a1^2 the primitives are

    U_X=-a1*partial_u g,
    U_v1=-a1*partial_u g,
    U_v0=partial_u U0-a0*a1/(v0-a0*u).

Equation (4) replays by exact polynomial arithmetic after clearing f^2.
The origin polynomial is z(z-1), so the actual zero jets through1 prove
the filtered series is zero. The ambient source has other solutions;
the actual jets select the correct one. No actual image kernel is computed.
This 16-term supplied identity is distinct from restarting the capped
enumerator or claiming that it discovered the identity.

The default budget is30seconds; zero disables it. CAS subprocesses have
the remaining wall-clock timeout. Python search budgets are cooperative,
checked between exact columns/eliminations; a single exact algebra operation
can exceed the nominal budget. A cap returns UNKNOWN and errors fail closed.
The complete route is general but can be impractical.

## 8. What remains

This closes the executable ordinary cutoff-ideal and least INTEGER/prefix
gates at their stated scopes, subject to trusted primary computations
and the classical holonomic imports. It does not compute least real slopes,
classify every safe phase, or solve unrestricted positive GMC (already
false). The next general work is structural safe normal forms and actual
image/quotient geometry, and optimization of slope geometry where needed.
Avoid replacing those questions by deeper small-degree campaigns.

Sealed Oaku identities, torsion witnesses and specialization thresholds
remain intact. Sabbah/Maisonobe retain convergent rational-slope credit;
the arbitrary nonconvergent formal extension and all-ring descent remain
separately scoped with priority unassessed. No broad-goal completion or
assurance promotion is asserted.
