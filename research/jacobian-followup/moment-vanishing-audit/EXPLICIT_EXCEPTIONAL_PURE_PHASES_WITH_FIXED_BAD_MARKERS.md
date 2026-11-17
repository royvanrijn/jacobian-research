# Explicit pure counterexamples on both exceptional quotients

17 November 2026. Written all-order counterexamples for the ORIGINAL
normalized Haar functional on SL2/2O and SL2/2I. Both full Mathieu--Zhao
statements are **FALSE**. The phase and marker below are fixed regular
invariants; neither depends on the moment order or the auxiliary prime.
External review remains deferred. Independent replay and formal verification
are not claimed. The finite checks accompanying this note verify constants
and normalizations; they do not supply the infinite-tail argument.

## 1. The two explicit original pairs

Write the original matrix and left root derivation as

    M = [[x,y],[z,w]],       xw-yz=1,
    D = z partial_x + w partial_y.

Use the standard invariant binary forms

    f8(X,Y)  = X^8 + 14 X^4 Y^4 + Y^8,                  Gamma=2O,
    f12(X,Y) = X Y (X^10 + 11 X^5 Y^5 - Y^10),          Gamma=2I.

In either row of this table put n=deg f, b=n/2, a=b-1, and set

    U=f(x,y),       F=f(z,w),
    R= f(z+F x,w+F y)/F,
    P=U^a R^b,      Q=F^a.                              (1)

The quotient defining R is a POLYNOMIAL. For
E=x partial_z+y partial_w its polynomial expansion is

    R=1 + sum_{k=1}^n F^(k-1) E^k F/k!.                (2)

All coefficients are integers. The two rows transform by the same right
Gamma action, F is invariant, and f is invariant; thus R,P,Q are ORIGINAL
regular functions on SL2/Gamma. No localization, extension of the marker
ring, finite-cover trace or change of the Haar functional is used in (1).

Explicitly, the fixed pairs are

    2O:  P=U8^3 R8^4,       Q=F8^3;
    2I:  P=U12^5 R12^6,     Q=F12^5.                   (3)

Let h be normalized Haar integration on SU2, restricted to these invariants,
equivalently normalized Haar integration on SU2/Gamma. We prove

    h(P^m)=0 for EVERY m>=1,                            (4)
    h(Q P^m)!=0 for infinitely many m.                 (5)

This is precisely a failure of the Mathieu property for ker h.

## 2. All-order purity by a REGULAR integration-by-parts identity

Give top coordinates weight +1 and bottom coordinates weight -1.
In the fraction field one has the finite identity

    R=F^(n-1) exp(F^-1 D) U.

For a fixed m>=1 put A=am and B=bm. Consequently

    P^m=sum_{j=0}^{nB}
         F^((n-1)B-j) U^A D^j(U^B)/j!.              (6)

The summand indexed by j has weight

    (n-2)(j-(n-1)B),

because nA=(n-2)B. These weights are distinct. Their sum is the regular
polynomial P^m, so they are its actual regular weight components; the
zero-weight term has no denominator at all. Haar invariance under the
compact left torus kills every other weight, giving

    h(P^m)=h(U^A D^j(U^B))/j!,      j=(n-1)B.

The original Haar functional satisfies h(DH)=0 for every regular H, by
left invariance and complexification of the Lie algebra. Repeated regular
integration by parts gives

    h(U^A D^j(U^B))=(-1)^j h(D^j(U^A) U^B)=0,

since deg_D(U^A)=nA=(n-2)B<j. This proves (4) in EVERY order, without
integrating a meromorphic function or transporting Haar through a rational
map.

## 3. The SAME fixed marker reduces to an integral polynomial

For Q=F^a, the unique zero-weight index in (6) is

    j=(n-1)B+b,        s=nB-j=b(m-1).

Thus

    h(Q P^m)=h(T_m),
    T_m=U^A F^-1 [t^j] f(x+tz,y+tw)^B
       =U^A F^-1 [t^s] f(z+tx,w+ty)^B.              (7)

Crucially, T_m belongs to Z[x,y,z,w]. In the last coefficient, at most
s of the B factors can contribute a positive power of t. Every other
factor contributes F. Its numerator is therefore divisible over Z by
F^(B-s)=F^b, which more than cancels the displayed F^-1. This also
justifies (7) without any meromorphic Haar integral.

The top and bottom degrees of T_m are equal, and both are

    N=nA+s=j-n=b(n-1)m-b.                             (8)

## 4. A rank-one congruence for the ORIGINAL normalized Haar functional

**Lemma.** Let p be an odd prime, and let T in Z[x,y,z,w] be homogeneous
of degree p-1 in each row. Then p h(T) is p-integral and

    p h(T) = [X^(p-1) Y^(p-1)] T(X,Y,X,Y) mod p.      (9)

**Proof.** A monomial has the form
x^i y^(p-1-i) z^k w^(p-1-k). Its Haar integral is zero unless
i+k=p-1. In that case the original SU2 model

    z=-conjugate(y),        w=conjugate(x)

and the uniform distribution of |x|^2 on [0,1] give

    h(monomial)=(-1)^k i! k!/p!.

Wilson's theorem gives i! (p-1-i)! = (-1)^(i+1) mod p.
Since k=p-1-i, multiplication by p makes the displayed Haar value
congruent to 1 mod p. These are exactly the monomials which contribute
to the coefficient on the right of (9). All denominators divide p!,
so p-integrality also follows. QED.

This lemma is a congruence of the original Haar functional. Evaluating
the polynomial on rank-one matrices is only an arithmetic device on its
coefficients; it does not replace the original SL2 phase space or its
volume by a singular degeneration.

## 5. Prime-indexed marked moments and the binary coefficient

Choose a prime of the form

    p=b(n-1)m-a.                                      (10)

Then N=p-1 in (8). At two equal rows, homogeneity gives

    T_m(X,Y,X,Y)=binom(nB,j) f(X,Y)^(A+B-1),
    A+B-1=(p-1)/b.

The polynomial divisibility in Section3 makes this evaluation valid even
on the zero set of f. Applying (9) yields

    p h(Q P^m)=binom(nB,j) H_p mod p,
    H_p=[X^(p-1)Y^(p-1)] f(X,Y)^((p-1)/b).           (11)

The binomial factor is NONZERO mod p. Indeed

    nB=p+q,       q=bm+a,
    j=p+n-1,      n-1<=q<p.

The coefficient identity (1+t)^p=1+t^p mod p gives

    binom(nB,j)=binom(q,n-1) mod p !=0.               (12)

The remaining issue is to prove H_p!=0 uniformly, not merely at a finite
list of primes. The explicit Klein relations do exactly that.

## 6. The Klein elliptic maps force H_p to be nonzero

We work over F_p with p>5 and p=1 mod b. Let C be the function field
of the superelliptic curve

    Y^b=f(t,1).

The polynomial f has simple roots in these characteristics. In particular
this is a field, t is a separating variable, and the differential dt/Y
is nonzero. Set

    O: G=X^2 Y^2 (X^4-Y^4)^2,           v=deg G=12,
    I: G=X^20+Y^20-228 X^15 Y^5
                    +494 X^10 Y^10+228 X^5 Y^15,     v=20,
    Cbin=(f_X G_Y-f_Y G_X)/(n v).

The identities, with these exact normalizations, are

    O: Cbin^2=f^3 G/36-3 G^3,
    I: Cbin^2=12 f^5+G^3/144.                       (13)

They are the retained [Klein-model identities](QUADRATIC_TRACE_MODELS_AND_HAAR_RESIDUES.md)
and are also directly checked by polynomial expansion in the accompanying
finite normalization checker.

In the function field C put

    g=G(t,1)/Y^(v/2),
    c=Cbin(t,1)/Y^((n+v-2)/2).

Equation (13) gives explicit maps to the nonsingular elliptic curves

    E_O: c^2=g/36-3g^3,
    E_I: c^2=12+g^3/144.                             (14)

Their pulled-back differential has an especially useful normalization:

    dg/c=-v dt/Y.                                    (15)

For completeness, Euler's homogeneous identity gives, on the affine
chart with second coordinate 1,

    Cbin=f'G/n-fG'/v.

Differentiating g and using Y'/Y=f'/(bf) proves (15). Its right side
is nonzero for p>5, so the function-field map from (14) is nonconstant
and separable. No assertion about a quotient graph or a resolution tree
is used here.

Here is the entire Cartier calculation needed. In a one-variable
function field of characteristic p, write a differential with respect
to a separating variable t as

    omega=sum_{i=0}^{p-1} a_i^p t^i dt;
    Cartier(omega)=a_(p-1) dt.

This definition is independent of the separating variable. It commutes
with separable pullback: use the same separating variable from the
smaller field in both fields, and the same displayed expansion proves
the assertion. The definition and coefficient calculation are standard;
see [Achter--Howe, Sections2.2 and3, pp5--6](https://www.math.colostate.edu/~achter/math/hwcm1017.pdf).
We use Cartier on differentials throughout, not an interchange of
Hasse--Witt and Cartier--Manin matrices.

Since p=1 mod b, we have

    dt/Y=f(t,1)^((p-1)/b) dt/Y^p.

The numerator degree is at most 2(p-1), so only its coefficient of
t^(p-1) survives Cartier. Therefore

    Cartier(dt/Y)=H_p dt/Y.                           (16)

All constants here are in F_p, so their p-th roots are unchanged.
For c^2=h(g) with the cubics in (14), the identical calculation gives

    Cartier(dg/c)=H_E dg/c,
    H_E=[g^(p-1)] h(g)^((p-1)/2).

Naturality and (15) imply H_p=H_E. These coefficients are explicit:

    O, p=1 mod4:
      H_p=binom((p-1)/2,(p-1)/4) (-1/12)^((p-1)/4);

    I, p=1 mod6:
      H_p=binom((p-1)/2,(p-1)/3)
                   12^((p-1)/6) 144^(-(p-1)/3).     (17)

Every factorial argument in these binomials is between 0 and p-1,
and all other factors are units for p>5. Thus H_p!=0 for EVERY prime
under consideration. This is a uniform algebraic calculation; neither
a conjecture about ordinary reductions nor a prime search is involved.

## 7. Infinitely many nonzero orders for ONE fixed original marker

Combining (11), (12) and (17) proves the stronger statement

    v_p(h(Q P^m))=-1
    whenever p=b(n-1)m-a is prime.                    (18)

In the two cases these progressions are

    O: p=28m-3,        p=25 mod28;
    I: p=66m-5,        p=61 mod66.                    (19)

They consist automatically of primes congruent to 1 mod4 and 1 mod6,
respectively. Since gcd(25,28)=gcd(61,66)=1, Dirichlet's theorem supplies
infinitely many primes in each progression, hence arbitrarily large
m with h(Q P^m)!=0. We use only the standard infinitude theorem; see
[Sutherland, Theorem18.1](https://ocw.mit.edu/courses/18-785-number-theory-i-fall-2021/mit18_785f21_lec18.pdf).
This proves (5) for the SAME P and SAME Q of (3), with no dependence
on p or m in either polynomial.

As normalization controls, direct original-polynomial integration gives

    O, m=2, p=53:
      h(QP^2)=-865594153113747456/298435517190385025,
      53 h(QP^2)=5 mod53;

    I, m=2, p=127:
      h(QP^2)=41958690403313837610721588134765625 /
        1532800285927159272199695730537217233992065276406,
      127 h(QP^2)=3 mod127.

The finite controls independently agree with the earlier exploratory
transvectant formula. The proof of the tail is (18)--(19), not those
two samples. See the [normalization checker](verify_exceptional_fixed_marker_counter.py)
and [its retained output](exceptional_fixed_marker_counter_checks.json).

## 8. What this says about the unit-invisible corrected obstruction

For each pair (3), the ACTUAL ordinary unit Cauchy transform is exactly

    C_P(z)=1/z

as its physical germ at infinity and analytic continuation. Yet the
ACTUAL ordinary marked transform for the displayed Q has infinitely
many nonzero coefficients, by (18). This is an ordinary moment statement,
not an exponential-generating-function substitution.

By the retained [intrinsic safety criterion](INTRINSIC_SOURCE_SPLITTING_AND_FULL_HAAR_SAFETY.md),
the intrinsic boundary extension is nonsplit and the actual corrected
Neron component at 0 is nonzero. By the proved
[all-genus unit-invisibility theorem](PURITY_UNIT_IDENTITY_AND_THE_CORRECTED_INVISIBLE_QUOTIENT.md),
this nonzero ORIGINAL rational component lies in the unit-invisible,
original-inertia-invariant obstruction space. Thus nothing inferred
solely from the simplicity of C_P can detect this component. These are
actual pure exceptional failures in that invisible space, rather than
abstract graph kernels or the previous nonpure controls with pairings
2 and -2.

This deduction consumes the existing extension calculation without
identifying an ambient resolution tree with the finite-fibre graph,
contracting any good component, or changing the volume. An explicit
graph presentation or a numerical value of its corrected pairing is
not needed for the counterexample and is not computed here.

## 9. Compatibility with capacity, automorphisms and family transfer

The rational Ga comparison which suggested (1) is informative but is
not used to transport Haar in the proof. Define

    R_-=F^-1 f(z-Fx,w-Fy),
    S=F^(n-1) U^b R_-^a.

Both R_- and S are original regular invariants. R_- has only nonpositive
weights, so S has only strict negative weights: its maximum is -na.
Thus S is pure SAFE. Nevertheless

    exp(F^-1 D) S=P.

This is a NONGLOBAL constant-volume rational Ga shift, along the one
remaining dominant-root direction. It fixes the regular marker Q, but
does not preserve the original Haar period. Indeed h(QS^m)=0 for all
m>=1, whereas (18) holds for P. Regularity of an inverse marker does
not repair the lost unit on the pole-free chart.

This seed has multiple negative weights. It is outside the proved
single-weight and positive-cone globalization theorems; none of those
theorems is contradicted. More explicitly:

    O: matrix degree(P)=280,  deg_D(P)=56,
       Ga-leading coefficient F^35, bottom degree280;
    I: matrix degree(P)=924,  deg_D(P)=132,
       Ga-leading coefficient F^77, bottom degree924.

The highest matrix part is

    F^((n-1)b) U^(a+b).

Its dominant left-root multiplicities are224 of280 and792 of924,
respectively, both strict. D annihilates this bottom kernel line.
The allowed pole order1 meets the two-integral-translate bounds
35>=ceil(56/2) and77>=ceil(132/2); the orbit-degree globalization
inequality fails as required. The phase has33 and73 distinct torus
weights, respectively, so it also lies above the retained sparse
capacity failure thresholds. These observations preserve, rather than
replace, the existing all-automorphism and capacity constraints.

The full [family-transfer theorem](FULL_MIXED_HAAR_FAMILIES_AND_TORSION_NORMALIZATION.md)
is retained. Its transfer equivalences and pure/safe ideals remain
valid; the proposed universal inclusion of the full pure locus in the
safe locus is now FALSE, already over C, at the explicit points (3).
No additional finite-degree or decidability assertion is being used to
substitute for that inclusion. These examples decide the two full
exceptional statements negatively by actual fixed-marker failure.

Earlier dated OPEN statements in retained notes describe the state
before this construction. All earlier proved results, replay inputs,
other research goals and the external-review deferral are preserved.

## 10. The retained two-model reduction now gives the full connected classification

Let G be any connected complex reductive group and H any closed reductive
subgroup, with no sphericity or component-group restriction. The retained
[two-model reduction](NONSPHERICAL_PAIRS_AND_EXCEPTIONAL_SL2.md), Theorem2,
now implies

    ker(normalized Reynolds on C[G/H]) is Mathieu--Zhao
        if and only if [G,G] is contained in H.        (20)

Indeed, if the derived group is contained in H, the quotient is a torus
and the torus theorem applies. Otherwise the reduction either already
provides a fixed pure/bad pair or gives a surjective equivariant morphism
q:G/H -> SL2/2O or SL2/2I. Normalized Reynolds satisfies h_X(q^*R)=h_Y(R).
Pulling back the corresponding FIXED P,Q of (3) preserves every pure
moment and the SAME infinitely many nonzero marked moments from (18).
Thus every non-toral pair fails the Mathieu property. The stronger
assertion of nonzero marked moments in EVERY order is not needed for
(20) and is not asserted for the exceptional pairs.

This is the previously recorded exact equivalence, now with both model
failures supplied. It is separate from, and preserves, the already proved
moment-nullcone classification. Other Gaussian/Image problems are not
closed by (20). No external review has been requested or performed.
