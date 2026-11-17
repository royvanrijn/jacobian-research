# The Laurent gap-21 family: complete affine singularity strata

14 November 2026. PF2R7KCOLL1 completes the affine collision analysis
for the Laurent coefficient component K of PF2R7FAM1. PF2R7KCOLL2
combines it with the polynomial component and states the resulting
degree-six additional-image restrictions. Every parameter is covered;
no topology transport or global Keller realization is claimed.

Exact integer identities, quotient-ring identities, a characteristic-zero
point and finite-field certificates support the written proof below.
No printed discovery basis, capped run, algebraic-root approximation or
unverified parameter sample is a proof premise. This is an authored
written proof, without formalization or independent external review.

## 1. A polynomial parameter for the Laurent chart

Let x be the nonzero parameter of
[PF2R7FAM1](../jc2-r7-image-family/PROOF.md), and put z=x^3. Replace the
normalization parameter by t/x, and scale the target coordinates by
x^6,x^10. This gives polynomials P_z(t),Q_z(t) over Q[z], with leading
terms t^6,t^10. The reconstruction in [laurent.py](laurent.py) is exact:

    P_z(t)=x^6*p_x(t/x), Q_z(t)=x^10*q_x(t/x).

All powers of x on the right are multiples of three. Here z is nonzero.
Each z corresponds to three x parameters. These source and target
scalings preserve affine singularity types. For global image matching
they must be applied to every other target image in the same coordinates,
not separately to each curve.

Define C1=3z+2, C2=15z-4, W=25z^2+5z+4. The exact polynomials
F(z),T(z),D(z), of degrees23,27,50, are stored in
[exceptional_moduli.json](exceptional_moduli.json) and
[triple_point_certificate.json](triple_point_certificate.json).
Their definitions are bound to the reconstructed identities below.
The letter F here is a projection factor, not a Keller map.

All of z,C1,C2,W,F,T,D are squarefree and pairwise coprime. The complete
affine singularity classification is:

| z stratum | Number of z values | Affine singularities |
|---|---:|---|
| C1*C2*T*D nonzero, z nonzero | all other values | 12 ordinary nodes |
| D=0 | 50 | 1 ordinary triple point + 9 nodes |
| T=0 | 27 | 1 ordinary A3 tacnode + 10 nodes |
| C1=0 | 1 | A2 cusp + A5 two-branch tangency + 8 nodes |
| C2=0 | 1 | A2 cusp + 11 nodes |

Thus exactly79 nonzero z values, or237 Laurent x values, are nonnodal.
F=0 and W=0 are projection exceptions only; their members have12 nodes.

## 2. The unordered collision scheme and its apparent discriminant

For two normalization parameters s,t put v=s+t,w=s*t and form the
symmetric polynomial divided differences

    A(v,w,z)=(P_z(s)-P_z(t))/(s-t),
    B(v,w,z)=(Q_z(s)-Q_z(t))/(s-t).

The degrees in w are two and four. Their last two subresultants are
exactly

    -v^3*(k(v,z)*w+n(v,z))/12,
    v^10*H(v,z)/81,

where H has degree12 in v and constant leading coefficient81.
All three polynomials are reconstructed from P,Q. In particular the
factor v^10 is a leading-coefficient artifact, not ten finite pairs.
Up to nonzero rational factors the identities are

    Res_v(H,H_v)=z^8*C1^3*F^2*T,
    Res_v(H,k)=z^4*W^7*F^2.

For v nonzero, the w-leading coefficients of A and B are3v and5v.
Consequently each H root gives a finite common w. Away from F and W,
linear recovery gives w=-n/k. At v=0 the direct equations instead give

    w=-(15z^2+4z+4)/12,
    B(0,w,z)=constant*z^2*W.

The checker verifies the nonzero rational constant, rather than dropping
this chart. H(0,z)=-2z^2*W. Thus the zero-sum chart has a genuine pair
exactly at W=0.

At every F root a frozen modular gcd certificate proves that H and H_v
have exactly one monic linear gcd. This identifies their one double
v root. At this root k=n=0; the common incidence equation is the full
quadratic A in w. The checker reduces B modulo A and gets zero, and
verifies that A's discriminant and the norm of

    J=A_v*B_w-A_w*B_v

are units modulo F. Hence there are two distinct, reduced collision
pairs with the same v. The remaining H roots are simple. At W roots
another modular gcd proves gcd(H,k)=v, and direct substitution gives
J nonzero at the sole zero-sum collision. Every other pair is recovered
regularly. These calculations prove that F and W introduce no local
pair tangency; Section5 also excludes triple values there.

At a T root the resultant Res(H,H_v) has a simple zero as a function of
z, since the other displayed factors are nonzero and T is squarefree.
The Sylvester matrix therefore has corank exactly one. If its corank
were at least two, every determinant cofactor would vanish, and so would
the derivative of its determinant. With fixed polynomial degrees,
Sylvester corank is the gcd degree. Therefore H has exactly one double
root and all other roots simple. Linear w recovery is regular there.
Together with immersion and the absence of triples proved below, this
gives exactly one two-branch tangency of intersection multiplicity two,
an A3 tacnode.

## 3. Exhaust the possible triples with literal identities

Let a monic cubic encode an unordered triple of normalization parameters:

    c(t)=t^3-S*t^2+U*t-V.

For three distinct roots to have the same image it is necessary and
sufficient that both P and Q be constant modulo c. The coefficients
of t^2 and t in their remainders give four integer-primitive equations
f0,f1,f2,f3, reconstructed in [laurent.py](laurent.py). Their weighted
degrees in (V,U,S), of weights (3,2,1), are4,5,8,9.

No repeated cubic root occurs at an immersed parameter: a repeated root
of both P-p0 and Q-q0 would have P'=Q'=0. The complete nonimmersion
locus is C1*C2=0. Both retained nonimmersed classes have at most two
normalization branches at every affine image value, and their ordinary
cusps have multiplicity two. Thus they have no such cubic divisor,
even with repeated roots. It is enough to study the immersed parameters.

First keep S=0 explicitly. The two P equations there are

    48U^2+(180z^2+96z+48)U+225z^4+240z^3+144z^2+48z+16=0,
    -(24U+45z^2+24z+12)V-12Uz-15z^3-4z^2-4z=0.

If the coefficient of V vanishes, substitution gives the necessary
conditions

    C1*(75z^3+30z^2-12z+8)=0, z*C1*(5z+2)=0.

Their remaining factors are coprime. Otherwise solve for V and reduce
the two Q equations modulo the quadratic in U. The gcd of their two
resultants with that quadratic is a nonzero rational multiple of C1^2.
Thus S=0 supplies no immersed nonzero-z triple. The C1 class is already
known to have no triple, as explained above.

For S nonzero, f0 solves for V using only a constant times S. Substitute
it into the other equations. The first becomes a cubic in U with
constant leading coefficient384. Reduce the other two modulo it and
clear numeric denominators, obtaining three polynomials r0,r1,r2.
Every actual triple satisfies all three. The two long frozen integer
certificates are literal identities

    sum_i A_i*r_i = a*Lr(z)*S^2,
    sum_j B_j*f_j = b*Lw(z),

where a,b are nonzero integers and

    gcd(Lr,Lw)=constant*z^6*(3z+1)*D(z).

The default checker reconstructs both generator lists, performs the
sparse integer multiplications and compares every coefficient. It then
checks the exact univariate gcd. At z=-1/3 a separate small integer
unit certificate in the four original f_j excludes every solution.
Since S and z are nonzero, every remaining triple must satisfy D=0.
This necessity argument uses neither a printed elimination basis nor a
generic-rank assumption.

## 4. Construct and uniquely certify the triple at every D root

The polynomial D is primitive, irreducible over Q and has degree50.
The finite-field certificate proves irreducibility without factorization
in characteristic zero. Modulo109 its irreducible-factor degrees permit
only rational proper-factor degrees1 or49. A complete factorization
modulo151 permits neither. Every modular factor is independently checked
by its Frobenius identity and Rabin gcd tests, and their product is
checked against the reduction of D. Gauss's lemma then excludes every
proper rational factor.

Put K=Q[z]/(D). The point certificate gives S0,U0 in K and V0=N0/E0
with E0 nonzero in K. Exact polynomial reductions verify all four f_j
at this point after clearing powers of E0. The resulting cubic has
three distinct roots: D is coprime to C1*C2, so a repeated root would
contradict immersion. This constructs an unordered triple at every
complex conjugate root of D.

Uniqueness uses a different certificate. Consider all monomials in
(V,U,S) of weighted degree at most20. There are358. Form coefficient
rows of all allowed monomial multiples of the four f_j within that
bound. The certificate selects357 rows and357 columns. At the exact
finite-field specialization

    prime=101, z=21, D(21)=0 mod101,

its determinant is nonzero, as verified by exact modular elimination.
The leading coefficient of D is nonzero modulo101. Consequently this
minor is a nonzero element of K: if its integer determinant polynomial
were divisible by primitive D over Q, Gauss's lemma would make it
vanish at that modular root. This is a rank proof, not promotion of
a modular parameter sample to a curve theorem.

The characteristic-zero point supplies a nonzero evaluation vector on
these358 monomials, so the rank over K is at most357. The modular minor
makes it at least357. Therefore the kernel is exactly one-dimensional.
Any other geometric solution of the four f_j gives another kernel
vector with constant coordinate1. It must be the same vector; since
the monomial set includes V,U,S, every coordinate is the same. Thus
there is exactly one unordered triple at each D root.

In particular there is no value with four or more normalization branches:
four distinct parameters would give four different unordered cubics.
Nor can two distinct target triple values occur at the same parameter.
D is coprime to F,T,W,C1,C2, so all its collision pairs are reduced and
off-diagonal. Its unique triple point is ordinary, and the other nine
pairs give nine nodes.

## 5. Complete the counts and the nonimmersed rows

Away from F,W,C1,T, H has12 distinct roots, each with one recovered w.
The F and W controls in Section2 retain12 distinct reduced collision
points as well. Immersion, together with the exhaustive triple argument,
shows that every remaining image singularity is an ordinary node.
At T=0 there is one double collision and ten simple node pairs, with
no triple. At D=0 three of the12 simple pairs combine into the unique
ordinary triple, leaving nine nodes.

At z=-2/3 the scaled coefficients equal the retained nonimmersed class I.
At z=4/15 they equal class II after the source change t maps to -t.
Their complete affine types A2+A5+8A1 and A2+11A1 are imported from
PF2R7CUSP1, with exact coefficient comparison in the checker. This also
justifies the nonimmersed exclusion used in Section3. Every possible
finite singularity has now been covered, since a finite birational
polynomial normalization has a singular image value only at a multiple
normalization fibre or a nonimmersion parameter. This proves PF2R7KCOLL1.

## 6. Consequences for the two-family and Keller research frontiers

PF2R7COLL1 gives34 nonnodal coefficient points on polynomial L.
The Laurent normalization has237 such x parameters, with three points
in its intersection with L, all at C1=0. Thus the complete b=1
coefficient union has268 nonnodal coefficient points. This is not a
count of inequivalent affine-moduli classes; residual cube-root choices
are retained. Off L the inverse x in PF2R7FAM1 is unique; at the three
intersection points its d values are distinct. Thus no additional
identifications in the normalized coefficient image alter this count.
The two generic opens have twelve nodes, and every affine
exceptional stratum now has a known singularity type.

For a degree-six Keller candidate containing an immersed e2,f1 prime
with a Laurent D-root image, its one ordinary triple value supplies
three distinct omitted factors of ambient length two. They exhaust the
finite fibre. As in PF2R7COLL2, the original source omits that value and
no other omitted image can pass through it. This is a necessary matching
restriction, not a surjectivity contradiction or a global exclusion.
These combined counts and restrictions constitute PF2R7KCOLL2.

There is also a complete parameter restriction for the existing row-A
grouping H=C. Each normalization branch then supplies a length-two
ramified omitted point and a separate length-one unramified omitted
point, so a C triple value would require length9. Thus H=C excludes
all150 Laurent D-root x parameters. Combined with the13 polynomial-L
triple parameters already excluded in PF2OMIT11, this removes163
immersed coefficient points from that grouping. Their loci do not meet.
Within the immersed image family this branch-count test leaves the two
twelve-node opens and the99 tacnode coefficient points. These are
necessary eligible inputs, not realized covers. The remaining requirements
on B, including avoidance of Sing(C) in that grouping, are unchanged.

The remaining work is the actual global second-image and optional
unramified-image equations, their intersections and exceptional values,
their complement relations with the required S6 action, and the finite
normal cover with its source A2. Even with affine strata complete,
transporting a braid presentation requires a proper comparison including
infinity; no such topological transport is proved here. Whole normal r7,
other terminal families, (75,125) and JC2 remain open.

## 7. Discovery and assurance boundaries

Unsplit elimination and direct global lifts hit stated time caps. Reducing
the triple equations before elimination, and then using weights(3,2,1),
made the needed identities practical. An old `sat(...)[1]` idiom was
rejected after a membership check: this installed Singular returns an
ideal from `sat`, whereas `sat_with_exp` returns a list. Selecting the
first ideal generator had not preserved the ideal. Fractional export
parser errors were likewise rejected despite exit code zero. SymPy's
modular factor sorting also failed with the installed FLINT backend;
the final default does no modular factor discovery at all.

All large necessary identities replay with integer arithmetic. The
existence point, modular gcds and projection norms replay exactly in
polynomial quotient rings. The irreducibility and rank witnesses replay
by finite-field arithmetic. The geometric normalization, local-branch,
Sylvester-corank and rank-transfer arguments are written proof steps.
No formal verification or independent external review is claimed.

[Failed attempts](failed_attempts.json) and the hash-indexed optional
[discovery traces](discovery.json) preserve those computations separately
from the manifest inputs used by the default checker.
