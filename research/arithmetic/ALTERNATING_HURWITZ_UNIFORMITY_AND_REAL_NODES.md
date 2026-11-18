# Alternating Hurwitz spaces, symmetric pencils and real nodes

The first Hurwitz case beyond degree eight has one absolute full-alternating
braid orbit. A complete finite calculation proves this in degree ten. In every
degree the ordinary spin-lift sign fails to define a conjugacy-invariant
distinction for this passport. The symmetric coefficient locus has an exact
finite algebraic description, together with a composition obstruction and
uniform restrictions on real fibres. In degree eight it contains an infinite
rational family with the maximum four real collinear ordinary nodes.

Uniform rational existence in the admissible degrees 26, 50, 82, ... remains
open here, as do complex braid connectivity starting in degree fourteen and
unirationality of the Hurwitz spaces. The separate
[degree-twelve theorem](DEGREE_TWELVE_HURWITZ_ORBIT.md) now proves one
absolute full-A12 orbit. The [matching-walk theorem](REAL_PAIRED_LAURENT_MATCHING_WALKS.md)
now gives uniform real-algebraic existence, classifies the all-real chambers
and determines their complete real-root spectra. Rational points in those
chambers remain a separate arithmetic problem.

## Theorem 1. The degree-ten absolute Hurwitz orbit

Let C be the 630 double transpositions on {0,...,9}, let
P=(0 5 1 8 2 6 9 7 4)(3), and use composition pq=p composed with q. Among

\[
 \mathcal T_{10}=\{(a,b,c,d,e)\in C^5:abcde=P\},
\]

exactly 48,892,572 tuples occur, of which 17,006,112 are transitive. All the
transitive tuples form one B5 orbit and generate A10. The absolute Nielsen
class with distinguished infinity profile (9,1) therefore has 1,889,568
points and one braid orbit.

**Proof, including the finite reduction.** Every class element is obtained
uniquely by choosing four symbols and one of three pairings. Thus
|C|=3 binomial(10,4)=630. For every product m of two class elements, store
all ordered pairs (a,b) with ab=m. There are 64,144 distinct such products.
For each m, each c in C, and each stored pair (a,b), the remaining pairs
are exactly those satisfying

\[
 de=c\,m^{-1}P.
\]

This enumerates every five-tuple once. It gives the displayed total count.
A tuple is nontransitive exactly when every entry fixes 3. Indeed any orbit
of its generated group is invariant under P, whose only orbits are {3}
and the other nine symbols. If some generator moves 3, its group orbit
meets both P-orbits and is the entire ten-point set. The calculation finds
31,886,460 tuples fixing 3 and hence 17,006,112 transitive tuples.

The Hurwitz moves on adjacent involutions are

\[
 (a,b)\longmapsto(b,bab),\qquad
 (a,b)\longmapsto(aba,a).
\]

They are mutually inverse and preserve both product and generated subgroup.
Breadth-first traversal of every such move, starting at the literal tuple

\[
 ((0\,1)(2\,3),(1\,4)(2\,3),(0\,5)(1\,6),
                       (1\,2)(6\,7),(1\,8)(6\,9)),
\]

visits all 17,006,112 transitive tuples. Its generation of A10 is the b=5
case of the
[alternating-target construction](../meridional-rank/representations/ALTERNATING_TARGET_PRIME_FAMILIES.tex).
Since the generated subgroup is constant on an orbit, every transitive
tuple generates A10. This argument uses no classification of primitive
groups containing double transpositions.

An arbitrary product of type (9,1) can be relabeled to P. Its centralizer
in S10 is <P>, of order nine. This centralizer acts freely on the generating
tuples: a permutation centralizing all their entries centralizes A10,
and that centralizer in S10 is trivial. For example, its value at 0 must
be fixed by the A9 point stabilizer and hence must be 0; transitivity then
forces every symbol to be fixed. The full braid twist effects simultaneous
conjugation by P or P^{-1}, so this quotient cannot split the orbit.
Consequently the absolute class has 17,006,112/9=1,889,568 points. ∎

The [checker](verify_degree_ten_hurwitz_orbit.py) contains the complete
C++17 algorithm as an embedded string and compiles it with assertions
enabled. Permutations and tuples are encoded injectively in 64-bit
integers. Every enumerated tuple, every neighbor lookup, the domain counts,
and both directions of every braid move are checked exactly. Wrong-product
and unconjugated-swap controls are rejected. The
[receipt](degree_ten_hurwitz_orbit_verification.json) pins the source,
embedded algorithm, literal seed inputs, compiler and counts.

Thus the [rational degree-ten pencil](RATIONAL_DEGREE_TEN_ALTERNATING_PENCIL.md)
lies in the same absolute Hurwitz component as the original knot tuple.
There is a rational point on this component, supplied by that explicit
cover. In the usual characteristic-zero absolute Hurwitz moduli problem,
the unique geometric component, stable under Galois, is defined over Q.
Neither the finite calculation nor existence of this rational point proves
unirationality or classifies rational points. Inner classes and their
field-of-definition refinements are not classified here.

## Proposition 2. The ordinary spin sign does not distinguish this passport

For n>=8, the double-transposition class C in A_n admits no
conjugacy-equivariant section to its inverse image in the ordinary spin
double cover of the natural permutation representation. In particular a
product of arbitrarily chosen lifts does not define a sign invariant of
tuples of this passport. In the corresponding reduced quotient, the
spin-central kernel is killed.

**Proof.** Work in the Clifford algebra with e_i^2=1 and
e_i e_j=-e_j e_i for i!=j. Let

\[
 d=(1\,2)(3\,4),\qquad h=(1\,3)(2\,4).
\]

These commute in A_n. Lifts in Spin(n) are

\[
 D=\frac{(e_1-e_2)(e_3-e_4)}2,\qquad
 J=\frac{(e_1-e_3)(e_2-e_4)}2.
\]

Direct multiplication gives D^2=J^2=-1 and JD=-DJ. Their actions on
the vector basis by conjugation are exactly d and h. Hence
JDJ^{-1}=-D. Both lifts of d have order four and are interchanged by a lift
of an element centralizing d. An equivariant section would have to assign
to d a lift fixed by this conjugation, which is impossible. Changing the
chosen lift of one tuple entry changes the product sign. Even choosing
the unique odd-order lift of the long cycle does not remove that ambiguity.

The commutator of lifts of the commuting pair h,d is -1. The reduced
cover construction kills just such central commutators with an inertia
class element. It therefore kills the spin-central -1 in this case.
This proves the statement about this double cover, without asserting a
classification of every possible reduced Schur multiplier. ∎

The exact Clifford calculation is independently recomputed by the
[spin and real-geometry checker](verify_alternating_real_and_spin.py).
For the lifting-invariant framework, Section 7.1 of Ellenberg,
Venkatesh and Westerland,
[Homological stability for Hurwitz spaces ... II](https://math.stanford.edu/~akshay/research/EVW2arxiv.pdf),
requires precisely a conjugacy-equivariant lift.
Their fixed-group stabilization theorem, Theorem 7.6.1, assumes sufficiently
large multiplicity in every inertia class. Here the alternating group grows
with n, while the long-cycle multiplicity stays one. That theorem does not
supply the missing uniform transitivity assertion. The
[source ledger](alternating_hurwitz_uniformity_source_ledger.json) records
the exact retained bytes and inspected locations.

## Proposition 3. A composition obstruction

Every characteristic-zero rational function with exactly two poles of
coprime orders (n-1,1) is indecomposable over an algebraic closure.

**Proof.** Suppose f=g composed with q, with both degrees greater than one.
The poles of f are the preimages under q of the distinct poles of g.
If g has one pole, its order is deg(g), which divides each pole order of f.
Coprimality forces deg(g)=1. If g has two poles, each must have exactly
one q-preimage, since f has only two poles. Both are totally ramified
for q, so deg(q) divides each pole order of f. Again this gives degree one.
Three or more poles of g would give at least three poles of f. ∎

In particular composition of smaller Laurent examples cannot provide the
desired degree recursion. This is an obstruction to that construction
method, rather than an obstruction to rational points in higher degrees.

## Proposition 4. An exact symmetric coefficient problem

Let n=2b>=10, with b odd. Restrict to monic rational even H with H(0)=0
and f(z)=(H(z)+a)/z, a!=0, having b distinct finite branch values, each
with two simple critical points. After a simultaneous source/target
quadratic twist preserving rational coefficients, the zero branch may
be normalized as follows:

\[
 H(z)+a=N(z^2),\quad N(u)=(u-1)^2G(u),\quad
 G\in\mathbb Q[u]\text{ monic},\quad \deg G=b-2.
\]

Put

\[
 S(u)=(3u+1)G(u)+2u(u-1)G'(u).                                  \tag{2}
\]

Conversely this normalized model has the stated paired-simple-critical
passport if and only if all of the following hold:

\[
 a=G(0)\ne0,\quad S(1)\ne0,\quad S\text{ squarefree},
\]
\[
 \det(Y-M\mid\mathbb Q[u]/S)=Q(Y)^2,                            \tag{3}
\]

where M acts by multiplication by N(u)^2/u, and Q is a monic squarefree
polynomial of degree (b-1)/2 with Q(0)!=0. Here division by u is valid
because S(0)=a. These conditions imply a/(n-1) is a square in Q.

If n-1=m^2 with m odd, the discriminant is the literal square

\[
 \operatorname{Disc}_z(H(z)-sz+a)
      =\left(m^{\,n-1}s\,Q(s^2)\right)^2.                       \tag{4}
\]

For n=26 this is a problem in the eleven coefficients of monic G of degree
eleven: S has degree twelve and (3) requires the square of a squarefree
degree-six polynomial. No solution of that coefficient problem is asserted.

**Proof.** Since f is odd, its finite branch set is stable under s->-s.
Its odd cardinality forces zero to be a branch value. The two critical
points over zero are r,-r; the only repeated root of the rational
polynomial H(sqrt(u))+a is therefore u=v=r^2, of multiplicity two.
Uniqueness implies v is rational and nonzero. Replacing its numerator
N_0(u) by v^{-b}N_0(vu) normalizes v to one. Over the algebraic closure
this comes from the source change z->sqrt(v)z and the target change
s->v^{-b+1/2}s; the resulting coefficients are rational. It need not be
a Q-isomorphism of covers. Factorization then gives the displayed N.

Direct differentiation gives

\[
 2uN'(u)-N(u)=(u-1)S(u),\qquad
 z^2f'(z)=(z^2-1)S(z^2).                                       \tag{5}
\]

The leading coefficient of S is n-1, while S(0)=a and S(1)=4G(1).
Thus the three first conditions are exactly simplicity and nonvanishing
of the critical roots. At a root u of S, the two critical points
z=+-sqrt(u) have opposite values whose square is M(u)=N(u)^2/u.
The characteristic polynomial in (3) records these squared values with
multiplicity. Requiring it to be a square with squarefree nonzero Q
says exactly that two distinct u-roots give each nonzero squared value,
and no u-root has zero value. This gives precisely two simple critical
points over each of the two corresponding values. The pair z=+-1
accounts for zero.

The even-degree resultant identity gives, with monic H and odd b,

\[
 \operatorname{Disc}_z(H(z)-sz+a)
   =(n-1)^{n-1}s^2
                 \det(s^2-M\mid\mathbb Q[u]/S).                 \tag{6}
\]

Indeed its critical-point factors are s-f(r), and their scalar factor
is (-1)^{b+1}(n-1)^{n-1}. Equations (3) and n-1=m^2 prove (4).
Finally the norm of u in the even-dimensional algebra Q[u]/S is
a/(n-1). Consequently

\[
 \operatorname{Norm}(M)=\frac{n-1}{a}\operatorname{Norm}(N)^2.
\]

Equation (3) and Q(0)!=0 make this a nonzero rational square, forcing
a/(n-1) to be a square. In the admissible degrees a itself is a square. ∎

This criterion is restricted to the even-H locus. It neither classifies
general pencils nor proves that every cover has such a symmetry.
For the real version, positivity of all roots of S, together with (3),
implies that every critical point and branch value is real.

A model satisfying these conditions has full geometric A_n: it is
transitive because H(z)-sz+a is irreducible over Qbar(s); the (n-1)-cycle
at infinity forces primitivity. The primitive double-transposition theorem
for n>=9, imported as Jones, Proposition 6.1(3), then gives A_n, and
the even inertia gives the reverse inclusion. The exact imported statement
and inspected source are in the
[existing priority ledger](alternating_hurwitz_priority_ledger.json).
This conditional group conclusion does not construct a degree-26 model.

## Lemma 5. No class-power obstruction in the admissible degrees

For odd m, each of the two A_{m^2+1} classes of cycle type (m^2,1)
is stable under every power coprime to m. The double-transposition class
also has that property. Thus the standard class-power rationality test
does not obstruct rational covers in these degrees.

**Proof.** Label the long cycle by Z/m^2, with one additional fixed symbol.
For a unit a modulo m^2, multiplication j->aj conjugates it to its a-th
power. Partition Z/m^2 by residues modulo m. Let epsilon be the sign of
multiplication by a on Z/m. Its permutation of the m blocks has sign
epsilon^m=epsilon. Within each block the map is an affine map of Z/m
whose translation has sign +1, since m is odd, and whose linear part has
sign epsilon. The product of these m signs is epsilon^m=epsilon.
The total sign is +1. Thus the conjugator belongs to A_{m^2+1}.
The class does split: the centralizer of a long odd cycle with one fixed
symbol consists of its powers and is entirely even. The calculation works
inside either split class. For an involution every coprime power is itself. ∎

This removes one necessary obstruction only. It does not establish descent
of individual covers, rational points on all components, or a rational
parametrization.

## Theorem 6. Uniform restrictions on real-root spectra

Let n>=8 be even, H real of degree n with H(0)=0, and a real nonzero.
Suppose f=(H+a)/z has exactly n simple finite critical points, paired two
at each of n/2 distinct finite branch values.

For every real regular value s the number R(s) of real roots of
H(z)-sz+a satisfies

\[
                         R(s)\equiv2\pmod4.                    \tag{7}
\]

If H is even, R(s)=R(-s). In that symmetric locus,

\[
 R(s)\le
 \begin{cases}
 n-2,&n\equiv0\pmod4,\\
 n-4,&n\equiv2\pmod4.
 \end{cases}                                                  \tag{8}
\]

In particular symmetric admissible-degree models, including degree ten,
have no totally real regular fibre. This does not prevent all their
critical points and branch values from being real.

**Proof.** For sufficiently large |s| there is one real root near the
simple pole and one near the odd-order pole at infinity; the other
n-2 roots are nonreal. Hence R(s)=2 on both outer intervals.
At a real branch value the two critical points are either real or a
nonreal conjugate pair. In the latter case all four nearby roots remain
nonreal, so R does not change. In the former case the two simple real
folds each change R by +-2, with total change -4,0,or4. This proves (7).
Oddness of f gives R(s)=R(-s), and (7) proves the first bound in (8).

Suppose n=4q+2 and H is even. There are 2q+1 branch values, an odd number,
so zero is a branch value and its critical points are r,-r.
If r is real, oddness gives opposite second derivatives at this pair.
For small nonzero real s, exactly two of the four local roots near this
pair are real. The remaining roots number n-4, so R(s)<=n-2; by (7)
it is at most n-4. If r is nonreal, reality and the uniqueness of the
critical pair over zero imply bar(r)=-r. The four local roots remain
nonreal for small s, again giving R(s)<=n-4.

There are at most q positive real branch values. Starting from R=2
on the outer positive interval, reaching n=4q+2 would require exactly
q jumps of +4 while moving toward zero. Every positive branch value
would then already have been used, so R would stay n near zero,
contradicting the preceding paragraph. Thus R<n for every positive
regular s; (7) gives R<=n-4. Reflecting s proves the negative case. ∎

The ordered jump sequence at real branch values is locally constant
while critical points remain simple and branch values remain distinct.
It therefore gives a necessary invariant of connected components of the
real coefficient locus. No completeness of this invariant is asserted.

## Theorem 7. An infinite rational family with four real nodes

For 1<x<2 put

\[
 d_x=\frac{47+7x^2}{168},\quad
 e_x=\frac{379+357x^2+196x^3}{21168},\quad
 a_x=-\frac{(7x^2-3)^2}{87808},
\]
\[
 H_x(z)=z^8-z^6+d_xz^4-e_xz^2.
\]

Every real x in this interval has eight distinct real critical points,
four distinct real paired branch values, and four real transverse
ordinary nodes on t=a_x in its discriminant curve. For

\[
                     x_k=1+\frac{11}{1+121k},\quad k=1,2,\ldots,
\]

the coefficients are rational, the geometric monodromy is A8, the
arithmetic monodromy is S8, and the Galois-closure constant field is
Q(sqrt(-7)). These covers are pairwise geometrically nonisomorphic.

**Proof.** This is the
[existing rational coefficient family](RATIONAL_ALTERNATING_PENCIL_REALIZATION.md)
with c=-1,t=-x. Its discriminant is -7^7 Q_x(s)^2, where
Q_x=s^4+B_xs^2+C_x and

\[
 B_x=-\frac{6912-14112x^2+5488x^3+8575x^4-7203x^5+
                     31213x^6+16807x^7}{533655864},
\]
\[
 C_x=\frac{(7x^2-3)^2(4-7x+7x^2)^2(48+56x+35x^2)^3}
                   {379718108242115328}.
\]

The critical polynomial in u=z^2 is
T_x(u)=7u^4-5u^3+3d_xu^2-e_xu-a_x. Its discriminant is

\[
 \frac{x^2(x-2)^6(x+1)^2(48+56x+35x^2)}{1259712}>0.
\]

At x=1 the signs of T_x at 0,1/100,1/10,2/5,4/5 are +,-,+,-,+.
Its four roots are therefore positive and distinct. As x moves in [1,2),
no collision is possible and no root can cross zero since a_x!=0.
All four roots remain positive. They lift to eight distinct real
critical points. Also

\[
 B_x^2-4C_x
       =-\frac{x^6(x-2)^5(x+1)^3}{10806531246}>0,\qquad C_x>0.
\]

At x=1, B_x<0. It cannot change sign while B_x^2>4C_x>0.
The two roots in s^2 are therefore positive and distinct. This gives
four distinct real branch values. Each contains two simple critical
points by the full discriminant identity.

Let f_x^+ denote the original c=1,t=x family and f_x the present one.
The identity f_x(z)=i f_x^+(iz) shows geometric equivalence over C.
For x=x_k all denominators are units modulo 11 and x_k=1 modulo 11.
The original family has the same good reduction at s=2 as its x=1 seed,
whose exact irreducible (3,5) factorization proves full geometric A8.
Thus the present twist has that group too. The discriminant has
nonsquare constant class -7, so its arithmetic group is S8 and its
constant field is Q(sqrt(-7)).

Finally the invariant of the rational function with its two unequal poles

\[
                  I=\frac{b_1 b_7^2}{b_5^3}=e_x
\]

is strictly increasing for x>0, because
e'_x=(714x+588x^2)/21168>0. Here b_j denotes the coefficient of z^j
in f_x, so b_1=-e_x,b_5=-1,b_7=1. To verify invariance, any isomorphism
of covers preserving the pole of order seven and the pole of order one
has source change z->lambda z and target affine change s->mu s+nu.
For two odd normalized models the constant term forces nu=0, and
the displayed coefficient ratio cancels lambda and mu. More generally
infinity is the unique branch value of inertia index six, so a
cover isomorphism must preserve it and its two pole orders.
The distinct x_k therefore give distinct geometric covers.

At a simple critical point r the discriminant normalization
r->(H'_x(r),rH'_x(r)-H_x(r)) has tangent slope r. The paired real
points are distinct and nonzero, giving two real transverse branches,
each transverse to t=a_x. Thus each value is an ordinary node and
consumes two intersections. The four nodes exhaust the degree-eight
intersection bound. ∎

Write 0<beta_x<alpha_x for the positive roots of Q_x. Exact Sturm
counts at x=1 give R(0)=2, R(1/150)=6, R(1)=2, with the middle
sample between beta_1 and alpha_1. The branch ordering and fold jumps
stay constant on 1<x<2. Hence the whole family has the exact spectrum

\[
 R(s)=
 \begin{cases}
 6,&\beta_x<|s|<\alpha_x,\\
 2,&|s|<\beta_x\text{ or }|s|>\alpha_x,
 \end{cases}
\]

away from branch values. For comparison, the known degree-ten real model
has R=6 for |s|<alpha_10, excluding its branch values, and R=2 outside;
the intermediate pair +-beta_10 and the zero branch are neutral jumps.
Its exact Sturm controls are R(1)=R(1040000)=6 and R(1100000)=2.
Both spectra are {2,6}, with different interval patterns.

## Verification, dependencies and remaining cases

Replay the two exact checkers from the repository root:

~~~sh
.venv/bin/python research/arithmetic/verify_degree_ten_hurwitz_orbit.py
.venv/bin/python research/arithmetic/verify_alternating_real_and_spin.py
~~~

The first requires g++ and about 0.8 GB RAM. Its receipt pins the exact
compiler and assertion-enabled flags. The second uses the pinned SymPy
1.14.0 environment and its
[receipt](alternating_real_and_spin_verification.json) records the full
symbolic family identity, exact Clifford arithmetic, rational Sturm
counts, modulo-11 family congruence, odd-square parity controls and
the generic eleven-coefficient degree-26 derivative identity.
The all-degree assertions use the written proofs, rather than the
finite odd-square controls or the degree-26 shape test.

The finite orbit theorem imports generation of the knot seed, not a
Hurwitz connectivity theorem. The rational real family imports the
previous seed monodromy proof and invariant argument. The conditional
full-group conclusion in Proposition 4 imports the precisely retained
primitive double-transposition theorem. All verification and checking
reported here are internal; no formal verification, external human
review or independent external replay is claimed.

The smallest unresolved complex Hurwitz case is n=14,b=7. The next
admissible rational literal-square degree is n=26,b=13. Degree twelve
is the first remaining rational real-node construction after degrees
eight and ten if a constant nonsquare discriminant factor is allowed;
its leading square class is -11. A rational degree-twelve literal square
is impossible by the leading-coefficient test. The exact normalized
coefficient problem above is a uniform reduction in the symmetric
locus, not a uniform existence proof. All-real chamber classification is
supplied by the matching-walk theorem; higher rational points remain open.
The [unmarked reduction](UNMARKED_POLYNOMIAL_HURWITZ_REDUCTION.md) now
supplies uniform marking and attachment control, leaving two unmarked
polynomial connectivity assertions as a sufficient all-degree criterion.

