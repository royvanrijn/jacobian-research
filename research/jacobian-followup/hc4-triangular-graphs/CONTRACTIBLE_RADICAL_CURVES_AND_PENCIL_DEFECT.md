# Contractible radical curves and the exact pencil boundary defect

15 November 2026. Claims **HC4MYIY0/1/2**.
This closes the contractible singular-cusp boundary left open in
[IX2](UNBOUNDED_QUADRIC_CLASSES_AND_RADICAL_TRACE.md), under its actual
complete-preimage hypotheses. The unrestricted HC4/JC2 problem remains open.

## IY0: an exact generic-pencil degree identity

Let K:A2_z->A2_T be a complex plane Keller map with field degree d.
Suppose nonconstant polynomials u(T),v(z) satisfy

    u(K(z))=c v(z),       c in C*,

and both geometric generic fibres are integral.
Over an algebraic closure k of C(t), K induces a nonconstant map

    C_s={v=t/c} -> C_t={u=t}

between smooth affine curves, of degree d. Smoothness follows from generic
smoothness in characteristic zero. The degree is the field degree d:
C(t) acts on both function fields through the displayed identity, and
geometric integrality ensures the constant-field base change remains a
field on each side. This preserves the generic finite rank.

Let g_s,g_t be the genera of their smooth projective completions, and
r_s,r_t the respective numbers of removed points. Put

    nu_s=2g_s-2+r_s,       nu_t=2g_t-2+r_t.

The affine map is etale: DK is invertible and differentiation of the
pencil identity identifies its nonzero tangent map.
Its projective extension is a finite morphism of smooth projective curves.
No point in the source affine curve maps to a target puncture.

Let q be the number of source points above target punctures. They lie in
the source boundary, so q<=r_s. Their ramification contribution is
d r_t-q. All other ramification is at source boundary points mapping
to finite target points, because the affine map is etale. Write

    E=sum_(P in source boundary, image(P) finite) e_P >=0,   (1)

where e_P is the local ramification index.
Riemann-Hurwitz therefore gives the exact identity

    nu_s=d nu_t+E.                                        (2)

Indeed its total ramification is
(d r_t-q)+sum_(boundary over finite)(e_P-1); adding r_s gives (2),
since there are exactly r_s-q terms in the second sum.

This retains rather than discards finite-image boundary points. In
particular it does not repeat the excluded boundary identification
audited in [the leading-fibre import](sources/leading-fibre-import.json).
The retained characteristic-zero
[Riemann-Hurwitz formula](https://stacks.math.columbia.edu/tag/0C1B)
is the external foundation; the boundary accounting is supplied here.

Consequences in every degree:

- If nu_t>0, then d<=floor(nu_s/nu_t). HW0 excludes d=2,
  because every quadratic field extension in characteristic zero is
  abelian Galois. Thus nu_s<3nu_t forces d=1 and the Keller inverse.
- If the two generic fibres have the same genus and exactly one puncture,
  the projective map has a preimage of the target puncture. It must be the
  only source puncture. Hence E=0, nu_s=nu_t=2g-1!=0, and d=1.
  This includes genus zero; etaleness on the affine part is essential.
- If nu_s=nu_t>0, then again d=1 and E=0.
- If nu_s=nu_t=0, the identity alone gives no degree bound:
  t->t^m on Gm is finite etale for every m>=1.
- Etaleness alone does not imply finiteness of an affine curve map:
  t->t^m from Gm to A1 is etale on its source but omits zero.
  Its boundary point 0 maps to a finite target point with index m,
  so nu_s=0=m*(-1)+m exactly as (2) requires.

Once d=1, the retained HW0 coprime-pullback denominator argument gives
the actual polynomial Keller inverse. No properness of K is assumed.

## IY1: isomorphism over a complete contractible curve forces an inverse

Let K be complex plane Keller. Suppose A is a reduced irreducible closed
target curve, topologically contractible, and the *entire reduced*
inverse image C=K^(-1)(A) maps isomorphically onto A. Then K is a
polynomial automorphism.

This is a statement about the full inverse image, not one selected
branch or a generically injective normalization.

**New explicit external input.** We import the Lin-Zaidenberg theorem,
including its smooth Abhyankar-Moh-Suzuki case, in the exact form of
[Palka, arXiv:1405.5391v2, Theorem A](https://arxiv.org/html/1405.5391v2):
an irreducible topologically contractible plane curve is carried by
polynomial plane coordinates to x^p=y^q for coprime positive p,q.
The singular case has p,q>1; a smooth such curve is a coordinate line.
This theorem is not reproduced by the symbolic checker.

Apply it independently to A and C, which are isomorphic.
In the smooth case both curves can be made coordinate lines.
In the singular case they become monomial cusps with the same pair
of exponents up to interchange. Here is why the exponents agree:
an isomorphism of cusp rings extends to their normalizations C[t].
The unique singular point is t=0, so the normalized isomorphism is
t->a t, with a!=0. It preserves the valuation semigroup <p,q>.
Its smallest positive element and its smallest element not divisible
by that one recover the unordered coprime pair p,q.

After these independent source and target plane automorphisms, K
is still Keller. The target defining polynomial h pulls back to a
nonzero constant times the source defining polynomial g:

    h composed with K=c g.                                (3)

The full inverse-image assumption gives the unique prime support,
and the retained HW0 square-free pullback fact gives multiplicity one.
Thus no additional source component is omitted in (3).

In the smooth case h and g are coordinate polynomials, whose generic
fibres are A1, with genus zero and one puncture.
In the singular case both generic pencils, up to nonzero scalar levels,
are

    x^p-y^q=t,       gcd(p,q)=1,       p,q>1.                (4)

These geometric generic curves have one point at infinity and genus

    g=(p-1)(q-1)/2.

For clarity, this is not inferred from the topology of the singular
zero fibre alone. Project (4) to the y-line: the extension has degree p,
is totally ramified at its q simple finite roots, and has one point
above infinity of index p since gcd(p,q)=1. Eisenstein at one simple
root proves irreducibility. The Kummer valuation at infinity has
denominator p, so that one point exhausts the degree; no other
puncture is omitted. Riemann-Hurwitz gives
2g-2=-2p+(q+1)(p-1), which is the displayed genus.
The same calculation applies to both source and target pencils.

Now apply the equal-genus one-puncture consequence of IY0 to (3).
It forces d=1 and hence the polynomial inverse. In particular this
covers every monomial cusp x^p=y^q in every coprime degree, after
arbitrary polynomial plane rectifications.

The proof does not assume that an arbitrary isomorphism between plane
curves extends to an ambient automorphism. Lin-Zaidenberg rectifies
the two contractible curves separately; the preserved semigroup
identifies the exponents. Nor does it assert that every polynomially
parametric curve is contractible: identifying two normalization points
produces a genuine exception to that premise.

## IY2: the strengthened HC4 radical-trace obstruction

Retain the actual closed rank-one embedding and literal identity

    j(z,w)=(w1,w2,X,Y1,Y2),
    X^2+w1Y1+w2Y2=D(z)+w dot K(z),       J(K) in C*.

Let f=X(z,0). IX1 proves D=f^2 and that K restricted to V(f) is a
closed immersion.

**A componentwise saturation criterion without global trace descent.**
Let C be an irreducible component of the reduced f=0 curve and let
A=K(C). It is closed and isomorphic to C by the actual closed
immersion. If C is topologically contractible and the entire reduced
K^(-1)(A) is C, then IY1 forces the inverse. Equivalently it suffices
that an irreducible equation b of A have just one source prime factor
in b(K), namely the equation of C: HW0 supplies square-freeness.
This is an explicit full-preimage condition; it is not automatic from
injectivity of K on the radical curve.

**Every contractible descended component now closes.**
If f is nonconstant and f^2 belongs to C[K], IX2 gives f=a(K),
and K identifies the entire reduced a(K)=0 curve with the entire
reduced target a=0 curve. For any one irreducible target component A,
its full inverse image is its corresponding isomorphic source component.
Thus if even one such A is topologically contractible, IY1 forces
the inverse, regardless of all other components and their intersections.

Under hypothetical nonautomorphy IX2 already makes every target
component have affine normalization A1. Such a component with all
its finite singular points unibranch is topologically A1: its finite
normalization A1->A is proper and bijective, hence a homeomorphism.
The IY1 contradiction now proves a strictly stronger necessary condition:

    Every irreducible component of a=0 must have a singular point
    with at least two branches on that component itself.            (5)

Intersections between different target components do not satisfy (5).
Thus singular cusps, including every polynomially rectified monomial
cusp, no longer remain open within the descended-trace family.
An actual affine line is not silently substituted for its singular
contractible image; IY1 is the additional argument that closes it.

Writing b_P for the number of local branches at a singular point of
one component, its Euler characteristic is

    chi(A)=1-sum_P(b_P-1)<=0                               (6)

under hypothetical failure. This follows by gluing the finitely many
normalization points of A1. A cusp has b_P=1 and chi=1, so it is
excluded; a nodal rational curve can pass (5)-(6) and is not
declared closed.

**A further all-degree pencil test on the remaining components.**
Let b be the irreducible equation of a target component and set g=b(K).
The complete curve isomorphism makes g irreducible, up to a scalar:
its zero set is the one source component and HW0 gives multiplicity one.
Both b and g are noncomposite, since a nontrivial polynomial composition
would factor over C. The retained noncomposite-to-geometrically-integral
generic-fibre argument of IL0 applies over C; hence IY0 applies to
b composed with K=g without an extra generic-connectedness hypothesis.

For their respective generic fibre invariants,

    nu_g=d nu_b+E.

Consequently, if nu_b>0 and nu_g<3nu_b, the inverse follows.
In particular equal positive invariants close; equal genera and
one puncture on both sides close in every genus.
Hypothetical failure with nu_b>0 requires nu_g>=3nu_b, since HW0
excludes field degree two.
Matching these generic pencil invariants for the remaining multibranch
embeddings has not been established in this note. This inequality is
a certificate and a necessary boundary; its use does not assume that
every remaining component meets it.

## Assurance, sources and replay

The new Lin-Zaidenberg classification input is explicit and byte-pinned
in [the context](sources/contractible-pencil-context.json).
Palka's Theorem A is read from the primary research manuscript;
its proof and the original 1983 paper are not claimed independently
reproved here. The retained Jelonek statement in Theorem 2.3 provides
a consistent cross-check, not a new use of his theorem about the
entire nonproperness set. No subset of that set is identified with
the whole set.

IY0 reuses the explicit Riemann-Hurwitz foundation retained in HZ1.
The geometric-generic-fibre passage in IY2 reuses the field-agnostic
noncomposite argument of IL0; HW0 supplies divisor domination,
square-free pullbacks and denominator removal.
All prior proof/checker/receipt/source bytes remain retained.

Run

    .venv/bin/python research/jacobian-followup/hc4-triangular-graphs/verify_contractible_pencil.py

for [exact controls](receipt-contractible-pencil.json), and see
[scoped validation](validation-contractible-pencil.json).
Controls check the boundary-defect identity, monomial genus arithmetic,
a genuine nonlinear cusp-trace embedding, normalization-versus-embedding
distinctions and scope boundaries. The all-degree geometric arguments,
including the imported classification, are not a finite search or a
formal verification.

Multibranch descended components failing the pencil test, genuinely
non-descended unsaturated traces, higher helper ranks, unaligned
intrinsic ambient projections, actual endpoint representatives and
unrestricted HC4/JC2 remain open.
