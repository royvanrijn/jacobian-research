# Mixed spherical cubics with a conformal Euler grading

16 November 2026. Exact internal written classification of a supplied
architecture. This does not classify every pure spherical cubic.

## Result and the precise remaining question

Use circular coordinates z,w,t with q=2zw+t^2 and normalized sphere
functional h on q=1. Its Gaussian counterpart has E(zw)=E(t^2)=1.
Write the UNIQUE harmonic representative

    f=H1+H2+H3,    Delta Hj=0,    deg Hj=j,
    Delta=2 partial_z partial_w + partial_t^2.

Assume H2!=0 and H1+H3!=0. A **conformal Euler certificate** means a
linear vector field E on the ambient three-dimensional space satisfying

    E(q)=2q,                 E(f)=delta f.                 (1)

The second equality is an actual polynomial identity for this harmonic
representative. E need not be tangent to the unit sphere.

**SIC2B33EULER20261116.** If (1) holds and h(f^m)=0 for every m>=1,
then, after a complex orthogonal change of coordinates, f is either:

1. one-sided for the circular rotation z->lambda z,w->lambda^-1 w;
   every occurring charge is strictly positive, or every one is negative;
2. a member of the calibrated shear family

       f=c w+b zt+a z^3,   abc!=0,   b^2+2ac=0.            (2)

Every nonzero member of (2) is a scalar and orthogonal transform of

       f0=w+zt-z^3/2.                                    (3)

This is the Rodrigues orbit. Both alternatives satisfy full marked
spherical safety, hence full two-pair (3,3) Image safety. For a multiplier
of coordinate degree b, the uniform strict cutoff is m>3b. The coefficient
three is sharp on (3); one-sided phases already have the stronger m>b.

The classification and marked bound below have elementary algebraic
proofs. They do not require the analytic ternary Gaussian cubic theorem.
The [checker](verify_conformal_euler_cubic.py) gives finite exact checks,
not external, independent specialist or formal verification.

A sufficient exhaustion statement for the full (3,3) problem is now:

> Every pure phase outside the one-sided nullcone admits (1).

This is a proposed structural target, **not a theorem**. A pure semistable
phase without (1) would refute this proposed classification while still
possibly being Image-safe. Proving marked safety directly on that
complement would also suffice. No support census establishes exhaustion.

## 1. An exact linear recognizer on all sixteen coefficients

Every linear E with E(q)=2q is Eu+A, with Eu the ordinary Euler field
and A in so(q). A convenient basis for so(q) is

    H=z partial_z-w partial_w,
    Lplus=z partial_t-t partial_w,
    Lminus=w partial_t-t partial_z.

For a specified f, (1) is the LINEAR system in four unknowns

    (Eu+alpha H+beta Lplus+gamma Lminus)f-delta f=0.       (4)

Coefficient comparison decides it exactly. Inconsistency is not an
unsafe verdict. For example z+zt+z^3 is one-sided and pure, but admits
no certificate (4).

The [adapter](conformal_euler_cubic.py) also converts an arbitrary
bidegree-(3,3) biform to its unique harmonic representative. Put

    s=WZ+VY, z=WY, w=2VZ, t=WZ-VY.

Then q=s^2. Lifting a harmonic of degree j by s^(3-j) gives sixteen
biforms, a rational basis of the FULL sixteen-dimensional coefficient
space. The adapter inverts this fixed 16x16 matrix; it does not select
coefficient supports or require a sparse input. The scalar harmonic H0
is its first-moment obstruction and must be zero for purity.

For an ordinary polynomial g of degree<=3, harmonic reduction is also
explicit: with homogeneous pieces gj, set L=Delta g3/10 and
c=Delta g2/6. The sphere representative has

    H3=g3-qL, H2=g2-qc, H1=g1+L, H0=g0+c.

## 2. Why there are only four three-summand patterns

The field A preserves every homogeneous harmonic degree. A nonzero
element of so(3,C) is either semisimple or nilpotent. If A is nilpotent,
(j+A)Hj=delta Hj forces delta=j for every nonzero Hj. This is impossible
when two different harmonic degrees occur. Therefore A is semisimple
and an orthogonal change puts

    E=Eu+kappa H,    kappa!=0.

Each nonzero Hj is a SINGLE H-weight vector, of integral charge kj
between -j and j, and

    j+kappa kj=delta.                                   (5)

Each such weight space in the degree-j harmonics is one-dimensional.
If all three summands are nonzero, their charges form an arithmetic
progression with nonzero integral step d. Reverse all charges if
necessary so d>0. The bounds |kj|<=j leave exactly:

| (k1,k2,k3) | kappa | delta | Pure possibility |
|---|---|---|---|
| (-1,0,1) | -1 | 2 | Excluded by h(f^2)=b^2/5 with H2=b(t^2-zw). |
| (0,1,2) | -1 | 1 | Excluded by h(f^2)=c^2/3 with H1=ct. |
| (1,2,3) | -1 | 0 | f=cz+bz^2+az^3 is one-sided. |
| (-1,1,3) | -1/2 | 3/2 | Exactly (2). |

The exclusions assume, as in this table, all three components are
nonzero. They follow from harmonic orthogonality and rotational charge,
not a finite-prefix assumption about the general coefficient space.

For the last row, direct sphere averaging gives

    h(f^4)=(4/35)c^2(2ac+b^2).                          (6)

Thus purity forces (2). Section 4 proves all-order sufficiency.

## 3. Two nonzero harmonic summands cannot cancel opposite charges

Suppose exactly H2 and one of H1,H3 occur. By (5) they are scalar
multiples of two standard weight vectors U,V, of charges k,l.
If one charge is zero, its second moment is nonzero and the other
charge cannot cancel it. If both charges have the same strict sign,
the phase is one-sided and is automatically pure and safe.

If k>0>l, take r=|l|/gcd(k,|l|), s=k/gcd(k,|l|). The only charge-zero
term in (aU+bV)^(2(r+s)) is

    binomial(2(r+s),2r) a^(2r)b^(2s) (U^r V^s)^2.

In standard circular coordinates with w=conjugate(z) on the real
sphere, U,V have real associated-Legendre coefficients. Their product
U^r V^s has charge zero and is a NONZERO REAL polynomial in t on
[-1,1], after absorbing nonzero scalar normalizations into a,b.
The integral of its square is strictly positive. Therefore the displayed
moment is nonzero. Opposite-charge two-summand purity is impossible.

This proves exhaustion within (1), including all vanishing-component
boundaries required by the mixed-parity hypothesis.

## 4. The surviving cancellation is a sphere-preserving shear

Let L=t partial_w-z partial_t. It is a rotation derivation, annihilates
z and satisfies h(Lg)=0 for every polynomial g. Thus D=zL is locally
nilpotent and

    h(Dg)=h(L(zg))=0.

Its polynomial automorphism Phi_u=exp(uD) is

    z -> z,
    t -> t-u z^2,
    w -> w+u zt-u^2 z^3/2.                              (7)

It preserves q and h exactly, has inverse Phi_-u, and acts on the
original regular sphere ring. This is a concrete instance of the
retained [Haar replica theorem](HAAR_SHEARS_AND_GENERAL_NULLCONE_CLASSIFICATION.md).
No claim that it preserves the real integration cycle is needed.

For c!=0, (2) is exactly

    f=c Phi_(b/c)(w).

Consequently all positive pure moments vanish. For every marker Q,

    h(Q f^m)=c^m h(Phi_(-b/c)(Q) w^m).                  (8)

If deg Q<=b, the pullback in (8) has ordinary degree at most 3b and
therefore H-charge at most 3b. Charge invariance makes (8) zero for
m>3b. In the one-sided case the charge of f^m already dominates every
degree-b marker when m>b.

If abc!=0, a torus change z->lambda z,w->lambda^-1 w and an overall
scalar normalize b=c=1; the relation then gives a=-1/2. Hence this is
one scalar/orthogonal orbit, not an unclassified parameter family.

For f0 in (3) and Q=w^b,

    h(w^b f0^(3b))=(-1/2)^b (3b)!/(6b+1)!! !=0.         (9)

Indeed the unique highest-charge term in Phi_-1(w)^b is
(-1/2)^b z^(3b), and h(z^k w^k)=k!/(2k+1)!!.
This proves sharpness for every b>=1.

Balanced biform markers of bidegree (b,b) give all spherical markers
of degree at most b. For an unbalanced multiplier with dual degree a
and coordinate degree b, a>b vanishes by degree; otherwise the output
coefficient of z^gamma is the scalar contraction after multiplying by
a fixed dual monomial of degree b-a, divided by gamma!. The resulting
marker is balanced of degree (b,b). The
[full coefficient transfer](ODD_HOPF_CUBIC_IMAGE_TRANSFER.md)
therefore turns (8) into the full Image cutoff m>3b.

## 5. Identification of the historical Rodrigues form

The exact biform from the
[Rodrigues note](../../foundations/extended-geometry/TWO_PAIR_SIC_BIDEGREE33_RODRIGUES_SURVIVOR.md)
is, in the rational Hopf variables of Section 1,

    P=z^3+szt-s^2 w/2.

Its sphere restriction is

    f=z^3+zt-w/2,
    H3=z^3, H2=zt, H1=-w/2.

The normalized form f0 in (3) also makes the interaction explicit. Its
even part B=zt is pure, while its odd part O=w-z^3/2 is NOT pure:

    h(O^4)=-4/35,    6h(B^2 O^2)=4/35,    h(B^4)=0.

Thus the quadratic term repairs an actual nonzero odd fourth moment.
Purity of f does not imply purity of its odd part, so the closed odd
slice cannot simply be applied component by component.

Here b^2+2ac=1+2(1)(-1/2)=0. The parity mixture is thus not accidental:
the quadratic and cubic terms are the first and second corrections of
the polynomial shear (7). They carry charges 1 and 3, opposite the
linear charge -1. Their calibrated cancellation survives in every
moment and all original regular markers are controlled by (8).

## 6. A weighted polar bridge, and why ordinary parity was insufficient

For a phase satisfying E=Eu+kappa H, E(f)=delta f, take a marker Q
of ordinary degree e and charge k. The charge-zero part of f^m Q has
ONE ordinary degree

    N=delta m+e+kappa k.                                (10)

If N is not a nonnegative integer this part is absent. Otherwise Gaussian
polar decomposition gives the exact identity

    E_G[f^m Q]=E[R^N] h(f^m Q).                          (11)

Both functionals kill the nonzero-charge parts. Equation (11) transfers
Gaussian purity and marked safety even though f mixes ordinary parity.
For (3), kappa=-1/2 and delta=3/2; each charge-zero term in f^m has
ordinary degree 3m/2. This explains its simultaneous Gaussian and sphere
purity. The direct shear proof above is stronger evidence for this class
than an unnecessary import of the full Gaussian cubic theorem.

This grading is DIFFERENT from requiring both the Gaussian quadratic
and the phase to have weight two: here the phase eigenvalue is 3/2.
The earlier canonical homogeneous exceptional-lift coverage exclusion
does not exclude this three-variable architecture.

## Exact evidence and remaining work

Run from the repository root:

```sh
.venv/bin/python research/moments/vanishing/verify_conformal_euler_cubic.py
```

The [receipt](conformal_euler_cubic_verification.json) retains the full
Hopf-basis adapter, all four weight patterns, two-summand exclusions,
moment and shear identities, sharp marked values, and positive/negative
grading-recognizer controls. An unsupported control is itself pure and
safe, making the recognizer's boundary explicit.

The full problem remains OPEN. The next structural gate is semistable
pure phases WITHOUT a conformal Euler certificate. Neither the weighted
classification, the Rodrigues point's known local isolation, nor sparse
support exclusions prove this complementary locus empty. A proof of the
proposed exhaustion, an exhaustive new component description, or full
marked-moment control there would advance the original completion goal.
