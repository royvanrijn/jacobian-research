# A regular cubic approximate root and a quadratic remainder

PF2RULED171. 17 November 2026. Actual complex plane polynomial Keller maps,
within the complete marked quadratic source model.

The full polynomial part of the ruling-six square root is regular in
the source-image algebra. Every survivor has

    P=H3^2+R2, deg_w H3=3, deg_w R2<=2,

with H3,R2 regular source functions. Both remain outside the target
field. This reduces the remainder degree; it neither produces a new
target coordinate nor proves its integrality over the target.

## 1. The exact remaining coefficient

Use PF169's actual polynomial normal form

    w^2=N+M*z, N=u*Sigma^2, M=u*A,
    h*M divides r=b2-N*h^2-K^2,
    B0=r/(2hM), W1=(h*w+K)*z+B0*w,
    P=W1^2+C*w*z+D*z+E*w+F,                               (1)
    C=a1-2B0*K,
    D=b1-2B0*h*N-B0^2*M, E=a0, F=b0-B0^2*N.

All coefficients in(1) are in C[u]. The original fields and full
normalizations are those of PF165:

    B=C[P,Q], Ftarget=C(P,Q), Efield=C(u,w), L=C(x,t),
    u=t^2, w=t*(Sigma+A*x), z=2Sigma*x+A*x^2,
    dP wedge dQ=du wedge dw/M.                            (2)

We distinguish F in(1), a coefficient polynomial, from Ftarget.
No prime is removed or replaced.

[PF170](RULING_SIX_QUARTIC_TRACE_AND_CONSTANT_SQUARE_ROOT_CORRECTION.md)
establishes C/h regular at every M-root, every off-M h-root
where K is a unit, and every off-M h-root whose outer face consumes
four escaping roots. Its quartic trace argument includes all rational
outer slopes, not only actual A pole divisors.

We prove the missing global divisibility

    h divides C.                                          (3)

By PF168, all off-M h-roots have specialization degree two or four
and square degree-six outer initials. After the PF170 cases, only
off-M common zeros of h,K with outer escaping mass TWO remain.
Their specialization degree is two.

## 2. Exact orders at a remaining outer double cluster

Fix such a root c, put v=u-c and

    j=ord h>=1, k=ord K>=1.

The square outer initial has lowest Y power four. Comparing degrees
six, five and four gives an INTEGER slope

    ell=j-k, j>2k, ell>k,
    in(P*)=v^-(4ell-2k) C0*Y^4*(Y-c0)^2, c0!=0,
    Y=v^ell*w.                                            (4)

Here j>2k follows because the unit degree-two specialization term
lies STRICTLY above this face. The degree-six and degree-five
coefficients are h^2/M^2 and 2hK/M^2.

Suppose C!=0 and a=ord C<j. Since W1 is regular, its degree-three
square contribution has order at least k, apart from terms of still
higher order coming from N. The actual outer initial has no
degree-three term, so

    a>2k-ell.                                             (5)

Indeed otherwise C*w^3/M would lie on or below the face, and its
coefficient cannot cancel against W1^2: k>2k-ell.
All coefficients D,E,F are regular. Write d=ord D when D!=0,
and d=infinity otherwise.

The regular cubic W1 has an exact p-independent Laurent zero wo
at the outer simple root of its initial:

    wo=-K/h+o(v^-ell), ord wo=-ell,
    ord W1_w(wo)=k-ell.                                  (6)

This is ordinary constant-unit Hensel lifting of the cubic, after
the substitution w=v^-ell Y. Its initial is a nonzero constant
times Y^2*(Y-c0); the nonzero root is simple.

## 3. Avoiding an outer graph would require precise cancellation

At wo the two possible leading remainder orders in(1) are

    C*wo*z(wo): a-3ell,
    D*z(wo):   d-2ell.                                   (7)

Since a<j=ell+k<2ell, the first is strictly below -ell.
The terms E*wo and F are therefore higher than the smaller of(7).

Suppose the smaller terms of(7) do not cancel. Let

    nu=min(a-3ell,d-2ell)<-ell.

The corresponding leading coefficient at wo is a nonzero constant.
Set w=wo+v^beta Y, where

    beta=(nu-2k+2ell)/2.

Then beta+ell>0. For the cubic candidate this is exactly the positive
quantity (a+ell-2k)/2 from(5); for the quadratic candidate it is
(d+2ell-2k)/2>0. Thus the new displacement is smaller than wo, and
the full polynomial P*-p has, after scaling, initial

    c1*Y^2+c2, c1*c2!=0.                                 (8)

The cubic W1 linearizes at wo. Higher terms of W1 and of the
remainder have positive margins because beta+ell>0; differentiating
a degree-at-most-three remainder term raises its order at wo by
at least ell, before multiplication by the displacement.

The parameter term -p is higher because nu<0. At a simple root
of(8) the exact derivative order is

    ord P*_w=nu/2+k-ell<0.

Hence omega=-du/(M*P*_w) is holomorphic. Pass to v=s^2 if needed.
All-order constant-unit recursion is over C[p], the original source
escapes, and Q has polynomial residue in p. The original full
normalization therefore contains an omitted prime with whole closed
graph image, forbidden by the retained Chau component theorem.

Consequently avoiding this contradiction requires equal orders in(7)
AND cancellation of their leading constants. In particular

    a>=ell, d=a-ell<k.                                   (9)

## 4. That cancellation creates an inner graph instead

First suppose d=0. Then a=ell>k. On the LOWER slope k, the entire
initial of P*-p, scaled by v^(2k), is

    Y^2*((alpha*Y+b)^2+d0),                              (10)

where alpha!=0 is the leading coefficient of K/M, b=B0(c), and
d0=D(c)/M(c)!=0. The cubic C term is higher since a>k; the
degree-six terms are higher since j>2k. All other terms are higher.
The constant coefficient of the quadratic in(10) is nonzero because
P*(c,w) has degree two. Its discriminant is -4alpha^2*d0!=0.
It has two nonzero simple roots, giving holomorphic coframe and
forbidden graph images. Thus d cannot be zero.

Now 0<d<k. Since C(c)=h(c)=K(c)=D(c)=0 and the specialization
has degree two, its leading coefficient B0(c)^2 is nonzero.
The cubic W1 therefore has a second exact p-independent Laurent zero

    wi=-(B0(c)/(K/M)_k)*v^-k+...,
    ord W1_w(wi)=0.                                      (11)

For example substitute w=v^-k Y in W1: after scaling, its initial
is Y*(alpha*Y+B0(c)), with a simple nonzero root. The term h*w^3
is higher by j-2k>0.

At wi, D*z has order d-2k. It STRICTLY dominates every other
remainder term, because

    ord(C*wi*z)-(d-2k)=a-d-k=ell-k>0,
    ord(E*wi)-(d-2k)>=k-d>0,

and F is regular. Its leading coefficient is nonzero.
Use

    w=wi+v^(d/2-k) Y.

The displacement is smaller than wi because d>0. The ENTIRE
equation again has a quadratic-plus-nonzero-constant initial.
Its roots are simple and constant, and

    ord P*_w=d/2-k<0.

Its coframe is holomorphic. The same all-order/full-normalization
argument gives a forbidden closed graph. This contradicts(9).

Therefore a<j is impossible at every remaining root. Along with
PF170 this proves (3) globally.

## 5. The complete regular approximate root

Put e=C/(2h) in C[u] and

    H3=W1+e.

Equation(1) becomes exactly

    P=H3^2+R2,
    R2=(D-2eK)*z+(E-2eB0)*w+F-e^2.                      (12)

The w*z coefficient vanishes because C=2eh. Thus R2 has
ruling degree at most two and all its canonical coefficients are
polynomials in u. H3 is the full polynomial part of the formal
square root at w=infinity. Its correction terms B0 and e are
polynomials in u, so H3 is regular in the source-image algebra.
Its expansion in powers of w can still have M-denominators; no
coefficientwise regularity in the different algebra C[u,w] is claimed.

On the central whole source line t=0, H3=2Sigma(0)K(0)x+e(0)
is nonconstant affine. The source-intersection and one-line
arguments of PF165 imply H3 is outside Ftarget. If R2 belonged
to Ftarget, H3^2=P-R2 would belong to B; the squarefree pullback
argument of PF165 would then put H3 in B, again a contradiction.
Both functions remain outside the actual target field.

There is no claimed integral equation over B for H3: (12) has R2
in the SOURCE algebra. Nor is (12) a legitimate target reduction.
[PF188](RULING_SIX_SQUARE_DEFORMATION_AND_MARKED_QUADRATIC_EXCLUSION.md) excludes the remaining A faces and central quadratic residue in the
retained marked quadratic model. Automatic derivative and trace
integrality and JC2 remain open.

## 6. Proof and verification boundary

PF165 supplies the actual odd companion and full graph obstruction,
PF168 the square outer initials and even specializations, PF169 the
regular cubic with cubic remainder, and PF170 the four-sheet and
unit-K cases. The nested two-sheet argument above is the missing
coefficient descent. It retains the inner branch created by the
only possible outer cancellation.

Exact controls verify the entire outer/inner substitutions, both
Puiseux parities, derivative and coframe orders, cancellation and
non-cancellation cases, and the surface-ring identity(12). All-order
Hensel lifting and the strict order comparisons establish the
unrestricted result, with no coefficient or companion degree cutoff.
The checker passes 112 exact controls.

Replay:

    .venv/bin/python research/jacobian/affine-complements/verify_ruling_six_regular_cubic_root.py

Independent whole-proof replay, formal verification and external
human review remain unclaimed.
