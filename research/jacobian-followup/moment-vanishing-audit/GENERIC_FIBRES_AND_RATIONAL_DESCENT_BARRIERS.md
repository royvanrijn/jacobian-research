# Generic fibres and rational descent barriers in every degree

15 November 2026. Internal written proofs with exact function-field,
elliptic-invariant and finite-group controls. No novelty, independent
specialist review or formal verification is claimed. The
[checker](verify_normalizer_fibre_geometry.py),
[receipt](normalizer_fibre_verification.json) and
[source audit](normalizer_fibre_source_audit.json) specify the evidence.
The general birational-group and differential-algebra arguments are
written proofs, not conclusions drawn from finite searches.

This note closes much larger construction classes than a degree or
support sweep. There are three results.

1. A primitive invariant phase on either exceptional SL_2 quotient cannot
   have generic fibre geometrically ruled over an elliptic curve. This
   follows from a general obstruction involving the second derived
   subgroup of a finite group acting birationally on a surface.
2. The **entire** birational orbit of Long's normalizer circuit, followed
   by **any** nonconstant rational reparameterization, fails to descend
   to either exceptional quotient. This includes all regular nonlinear
   changes of coordinates and all powers, in every degree.
3. The original Hopf circuit has rational generic fibre and does admit
   invariantization by birational coordinates for every finite subgroup.
   The resulting rational additive slice has unavoidable poles: no
   nonconstant rational function of it is regular. A general lemma for
   locally nilpotent derivations proves this last obstruction.

Thus birational invariantization is neither uniformly impossible nor
sufficient for moment descent. The remaining problem is a genuinely
regular invariant phase with an all-order pure identity and a fixed
infinite marked tail, or an all-degree positive safety theorem. Both
exceptional MZ decisions remain OPEN.

## 1. A general obstruction for elliptic ruled surfaces

Work in characteristic zero over an algebraically closed field k.
Write Gamma'=[Gamma,Gamma], Gamma''=[Gamma',Gamma']. Let Q_8 denote
the quaternion group, with its unique element of order2.

**Lemma 1.** PGL_2(F) has no subgroup isomorphic to Q_8 for any field
F containing k.

It suffices to extend to an algebraic closure of F. A nontrivial
projective involution has two fixed points. If it were the central
involution of Q_8, the whole group would preserve that unordered pair.
Its image would lie in G_m semidirect C_2, the normalizer of a torus.
Finite subgroups of this normalizer are cyclic or ordinary dihedral:
diagonal elements of finite order are roots of unity, and every
off-diagonal element is an involution. A noncyclic such group has more
than one involution. This excludes Q_8.

**Lemma 2.** If S is birational to E times P^1 with E a smooth
projective curve of genus1, then no finite group Gamma with
Q_8 contained in Gamma'' acts faithfully by birational automorphisms
on S.

Every birational self-map of S preserves its projection to E. Indeed,
composing it with that projection is constant on a general rational
fibre: a rational map P^1->E extends and is constant by
Riemann--Hurwitz. The inverse gives an inverse map on E. Thus there
is a homomorphism

    Bir(S) -> Aut(E), with kernel PGL_2(k(E)).           (1)

Choose an origin on E. Its full curve automorphism group is
E(k) semidirect Aut(E,0), and Aut(E,0) is cyclic of order2,4 or6.
In particular Aut(E) is metabelian, meaning its second derived
subgroup is trivial. The translation decomposition is recalled in
[Stanojkovski--Voll, Remark3.1](https://doi.org/10.1007/s00208-021-02193-8).
For the cyclic assertion, in a short Weierstrass model the origin-fixing
maps have the form (x,y)->(c^2x,c^3y), with c^4A=A and c^6B=B;
this gives exactly those possible orders in characteristic zero.

Consequently Gamma'' lies in the kernel of (1). Lemma1 excludes its
Q_8 subgroup. This proves the assertion. This is a birational result;
it does not depend on a chosen affine boundary or on a regular C^*
action on the surface.

The exact group controls give

    2O: derived-series orders 48,24,8,2,1,
        with (2O)''=Q_8;
    2I: perfect of order120, containing Q_8.             (2)

The checker constructs successive normal closures of generator
commutators. In both second derived groups it exhibits a subgroup of
order8 with element orders1,2,4 occurring1,1,6 times. This is Q_8,
not the ordinary dihedral group of order8.

## 2. Primitive pencils and reparameterization rigidity

Let X be an irreducible complex threefold with function field L, and
let p in L be nonconstant. Call p **primitive** when C(p) is
algebraically closed in L. Equivalently in characteristic zero, its
generic fibre is geometrically integral. Suppose that generic fibre
is geometrically ruled over a smooth genus1 curve E_t, t=p. Let
j(t) in C(t) be the j-invariant of that curve (or its Jacobian).

**Primitive-pencil criterion.** If the only Mobius map M with
j(M(t))=j(t) is the identity, then for every nonconstant R in C(T),

    Bir(X)_(R(p)) = Bir(X)_p.                           (3)

Here the subscript means fixing the indicated rational function
exactly, not merely its zero locus. Moreover this stabilizer has no
finite subgroup Gamma with Q_8 contained in Gamma''.

To prove (3), let sigma fix R(p). The element sigma(p) is algebraic
over C(p), since it satisfies R(sigma(p))=R(p). Primitivity puts it in
C(p). Applying the same argument to sigma^(-1) shows that sigma
preserves this subfield and acts on it by a Mobius map M. The induced
birational map of generic surfaces preserves their unique genus1
ruled bases, by the rational-fibre argument of Lemma2. Necessarily
j(M(t))=j(t), so M is the identity. The reverse inclusion in (3) is
immediate.

A finite subgroup of the stabilizer of p acts faithfully on its
generic fibre. Extend the base field C(t) to its algebraic closure;
geometric integrality preserves faithfulness. Lemma2 and (2) then
exclude both binary groups. This argument applies to arbitrary
birational changes of coordinates, including those with poles.
It asserts no preservation of Haar integration by such maps.

## 3. Exact elliptic ruled fibre of the normalizer circuit

Use the conventions of the
[nonspherical reduction](NONSPHERICAL_PAIRS_AND_EXCEPTIONAL_SL2.md):

    g=[[x,y],[z,w]], xw-yz=1,
    u=-2xy, v=2zw, tau=xw+yz, uv+tau^2=1,
    U=u^2, V=v^2, S=tau^2.

Long's regular Weyl-even normalizer witness is

    P_N=(1+U)*(V-S*(U+2)*(2-S*(U^2+2U+2))).             (4)

The pinned source supplies its all-order pure-null and bad marked
identities. We use its geometry here, without recomputing more moments.
Since UV=(1-S)^2, exact simplification gives

    P_N=(1+U)/U * (1-(1+U)^2 S)^2.                     (5)

Equation (5) is an equality of rational functions with the regular
function (4). It does not remove its regularity at U=0.

Set t=P_N, and on a dense open chart put

    b=(1+u^2)*tau,
    k=(1-b^2)/u.

Then

    E_t: k^2=t-(1-b^2)^2,                              (6)
    u=(1-b^2)/k, tau=b*k^2/t.                          (7)

The checker verifies both directions, including reduction modulo (6).
The right-T quotient field is therefore C(t)(E_t). The complete
SL_2 field is recovered by the independent right-torus coordinate x:

    y=-u/(2x), z=x*(1-tau)/u, w=(1+tau)/(2x).           (8)

These formulas verify determinant1 and the three definitions of
u,v,tau. Thus

    C(SL_2)=C(t)(E_t)(x).                              (9)

The quartic in (6) has four distinct roots for generic t. It remains
irreducible as a double-cover curve after algebraic closure of C(t),
so its smooth completion is a geometrically integral genus1 curve.
Equation (9) proves both primitivity of t and the elliptic ruled
birational type of the generic surface fibre. In particular (6) is
the base curve of that ruled surface, not the surface itself.

To compute j, multiply (6) by -1 and use the isomorphic quartic
y_E^2=b^4-2b^2+1-t. In the convention

    I=12ae-3bd+c^2,
    J=72ace+9bcd-27ad^2-27b^2e-2c^3,

the values are

    I=16-12t, J=144t-128,
    disc_polynomial=256 t^2(1-t),
    4I^3-J^2=27 disc_polynomial.                       (10)

The Jacobian is y^2=x^3-27Ix-27J, as recalled in
[Fisher, Section2, equation2](https://arxiv.org/pdf/2208.14977v1).
Consequently

    j(t)=64(4-3t)^3/(t^2(1-t)).                        (11)

The polynomial discriminant in (10) differs from the genus-one
discriminant normalized with an extra factor16 in that source; the
j computation uses the displayed Weierstrass equation and is
unambiguous. The checker independently computes the polynomial
discriminant and Weierstrass c4/discriminant ratio.

The divisor of poles of j is 2[0]+[1], and its zero divisor is
3[4/3]. Also j(infinity)=1728. Any Mobius symmetry of j must fix
0,1 and4/3 individually, since their multiplicities distinguish them.
It is therefore the identity. This checks the full Mobius symmetry
group, not just the possible sign change t->-t.

Applying Section2 gives the promised all-degree obstruction:

    For Gamma=2O or2I, every phi in Bir(SL_2) and every
    nonconstant R in C(T) satisfy
       R(phi^*P_N) not in C(SL_2)^Gamma.                (12)

If such an invariant existed, conjugating the right Gamma action by
phi would embed Gamma into the stabilizer of R(P_N), contrary to
(2)--(3). The result rules out invariant **rational** functions in this
construction class, and hence also regular ones. There is no remaining
degree, power, coordinate-change or unary rational parameter to test
within this class. Multivariable orbit norms and averages do not have
the form (12) and remain outside its scope.

## 4. Consequences for arbitrary invariant phases

The group obstruction also constrains phases other than P_N. On
SL_2, no primitive nonconstant rational function invariant under right
2O or2I has geometric generic fibre ruled over a genus1 curve. The
group already fixes that primitive parameter, so Lemma2 applies
without any hypothesis on its j-function.

For completeness, let P be any nonconstant **regular** invariant
function. Normalize C[P] in its relative algebraic closure inside
C(SL_2). That finite normal affine curve ring embeds in C[SL_2] by
normality. The curve is rational because SL_2 is rational; it has
only one point at infinity because C[SL_2] has only constant units.
Thus it is C[f], with f regular and primitive, and

    P=R(f), R in C[T].                                 (13)

This is the same normalization-curve argument used in the
[unary theorem](INTRINSIC_REYNOLDS_AND_FINITE_QUOTIENT_DESCENT.md),
now applied to the full relative algebraic closure. The finite group
acts on C[f] by affine transformations. After centering f, its image
is cyclic, of order dividing the exponent of the abelianization.

For perfect2I the image is trivial: f itself is invariant. Therefore
the primitive factor of **every** nonconstant regular invariant phase
on SL_2/2I has a generic fibre of a different birational type from
an elliptic ruled surface. This holds in every degree and does not
require pure nullity.

For2O the image is trivial or C_2. If it is trivial, the same
exclusion applies. If an elliptic ruled primitive factor occurs, the
only remaining possibility is a nontrivial involution f->2c-f, with
stabilizer2T on an individual generic fibre. Its genus1 j-function
must satisfy j(2c-t)=j(t). Such involutive pencils are not ruled out
by Lemma2. The j-function (11) has no such symmetry, explaining why
even this possibility is absent for (4).

These are necessary birational-type conditions for invariant functions.
They are not sufficient conditions for pure nullity, safety or failure.

## 5. Rational Hopf fibres really behave differently

The original Hopf witness, with radial variable specialized to1, is

    P_H=(1+u)*(v-(2+u)*tau^2)
       =(1+u)/u * (1-(1+u)^2*tau^2).                  (14)

Its bad marked identity is the original radial input in the pinned
Long draft. In contrast with (4), set b=(1+u)tau and t=P_H. Then

    u=(1-b^2)/(t-1+b^2),
    tau=b*(t-1+b^2)/t.                                (15)

Both directions are checked exactly. Together with (8), this gives

    C(SL_2)=C(P_H,b,x).                                (16)

The generic fibre here is rational, not elliptic ruled. Thus one
cannot turn the preceding obstruction into a classification of all
bad phases by asserting that every bad phase has elliptic ruled
fibre. The known Hopf witness itself would contradict that assertion.

There is even an invariant birational coordinate for every finite
right subgroup Gamma of SL_2. Let D=z partial_x+w partial_y and N=|Gamma|.
It commutes with right Gamma, and D(x/z)=1 in the function field. Put

    s_Gamma=(1/N) sum_(gamma in Gamma) (x/z)(g gamma).   (17)

This is invariant and D s_Gamma=1. For the first column (a,c) of gamma,
the determinant identity gives

    (xa+yc)/(za+wc)=x/z-c/(z*(za+wc)).                  (18)

Hence s_Gamma-x/z belongs to C(z,w), and

    C(SL_2)=C(z,w)(s_Gamma).                            (19)

The three generators in (16) and those in (19) are rational coordinate
systems. Sending (P_H,b,x) to (s_Gamma,z,w) defines a C-field
automorphism and therefore a birational self-map phi of SL_2 with

    phi^*P_H=s_Gamma in C(SL_2)^Gamma.                  (20)

This is an existence proof for birational invariantization, not a
preservation theorem for moments. The resulting function is not regular.

## 6. A general obstruction to regularizing a rational additive slice

**Slice lemma.** Let A be a C-domain, D a locally nilpotent derivation
of A, and h:A->C a linear functional with h(1)=1 and h(Da)=0 for all a.
If s in Frac(A) satisfies Ds=1, then

    A intersect C(s)=C.                                (21)

Indeed, if R(s) is regular, local nilpotence gives
D^n R(s)=R^(n)(s)=0 for some n. Since s is transcendental over C,
the rational function R has zero nth derivative and hence is a
polynomial. If its degree d were positive, D^d R(s) would be a
nonzero constant. This constant is in D(A), contradicting h(Da)=0
and h(1)=1. This proves (21), including rational reparameterizations
that are finite at infinity. No unitarity, normality, dimension or
degree bound is needed in this lemma.

For SL_2, D is the locally nilpotent left-unipotent derivation and
h(Da)=0 by Haar invariance. Applying (21) to (17) shows

    C[SL_2] intersect C(s_Gamma)=C,                     (22)

and a fortiori the same intersection with C[SL_2]^Gamma is C.
Thus no nonconstant rational unary transformation can regularize
the invariant rational phase in (20). This obstruction is stronger
than merely observing that the first slice x/z has a visible pole.

More generally, if a finite group commutes with D, averaging any
rational slice gives an invariant rational slice with the same
obstruction. For a nonzero locally nilpotent D a rational slice always
exists: choose a with D^n a!=0,D^(n+1)a=0 and use
s=D^(n-1)a/D^n a. In the homogeneous Haar setting this is an all-degree
obstruction to replacing a rational unipotent coordinate by a single
regular coordinate.

## 7. Closed routes and the remaining all-degree step

The entire birational-and-unary normalizer route is now closed by
(12), not just its previously tested averages, powers or finite degree
pieces. The rational Hopf route is possible as a field construction,
but its invariant slice cannot yield any nonconstant regular phase
by unary reparameterization, by (22).

Neither argument proves that multivariable regularization is impossible.
Combining the slice with functions of z,w can clear poles, but changes
the function and requires a new proof of all-order pure cancellation
and a genuinely infinite fixed marked tail. The intrinsic Reynolds
theorem applies to algebraic automorphisms and justified regular
constant-Jacobian maps. A birational coordinate isomorphism with poles
does not supply that period-transfer hypothesis. It must not be used
to declare (20) a regular MZ counterexample.

The next general avenues are therefore:

* a justified multivariable, regular and period-compatible construction,
  with a primitive fibration outside the exclusions of Section4;
* or an all-degree positive normal form placing every pure-null
  invariant function in the intrinsic marked-safe class.

Both exceptional MZ questions remain open. No new support enumeration
or degree-prefix campaign is proposed by these results. Independent
review of the general birational and slice arguments remains outstanding.
