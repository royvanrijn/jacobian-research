# Etale centralizers and Galois completion

15 November 2026. **QET1.** For a finitely generated actual scalar
centralizer with individual conjugates, its map to the space of
commuting images is etale. Weyl units force its image to retain the
generic point of every base divisor. A birational or Galois centralizer
field must therefore be the polynomial field, and QFC1 gives finite
integral full completion.

Equivalently, a nontrivial Galois centralizer field is an all-partner
obstruction to simultaneous individual conjugates. A nonpolynomial
finitely generated centralizer admitting those conjugates must have
non-Galois field degree at least three. These are statements about the
actual centralizer field, not leading-symbol fields or the full quantum
division indices. No such field certificate is supplied for a new Long
pair. Actual parameter poles of individual mates can be removed
finitely, and an already flat generic tuple descends integrally
without a finiteness hypothesis.

## Hypotheses and conclusions

Retain the algebra and classical-centralizer assumptions of
[QFC1](FINITE_CENTRALIZERS_AND_COMPLETION.md): E is an h-torsion-free,
h-separated K[[h]]-algebra with commutative reduction; the commuting Fi
have independent reductions fi and full joint Poisson centralizer K[fi].
Actual integral Ti satisfy [Fj,Ti]=h*delta_ij. One may also start
with actual generic mates: the parameter-pole lemma below normalizes
them integrally before applying this theorem. Write
\[
 L=K((h)),\quad P=L[F_1,\ldots,F_m],\quad
 B=\operatorname{Cent}_{E[h^{-1}]}(F_1,\ldots,F_m),\quad
 Q=\operatorname{Frac}P,\quad H=\operatorname{Frac}B.        \tag{1}
\]
QFC1 makes B commutative, and its commuting translation derivations
nabla_i=-ad(Ti)/h extend partial_(Fi). Their joint constants in B are L.

Assume additionally that B is finitely generated as an L-algebra,
H/Q is finite, and
\[
                            P\cap B^\times=L^\times.      \tag{2}
\]
For E the scalar finite-PBW algebra of rank m, the finite-field
hypothesis and (2) are automatic once B is finitely generated.
The dimension input is the classical commutative-subalgebra bound
trdeg C<=m in A_m. It follows from
[Bavula, Corollary1.6](https://arxiv.org/html/math/0401263):
extend scalars faithfully to an algebraic closure, apply that bound,
and descend dimension. Since B contains P of transcendence degree m,
its finitely generated fraction field is finite over Q.
For units, Bernstein degree is additive on products, so a unit of
A_m(L) has degree zero and is scalar.

**Centralizer geometry theorem.** Under these hypotheses:

1. B is smooth over L and P->B is etale. The open image of Spec B
   in Spec P contains every height-one point.
2. In the finite normalization C of P in H, above each base prime
   divisor at least one unramified prime meets Spec B. Thus every
   geometric codimension-one inertia generator has a fixed sheet.
3. If H=Q, or more generally H/Q is Galois, then B=P.
4. The integral centralizer is K[[h]][F]. Polynomial curvature
   homotopy gives a finite full completion; corrections can be h-small
   if the prescribed classical individual conjugates are flat.

In particular, a proposed fixed commuting scalar pair with finitely
generated centralizer and nontrivial Galois field H/Q cannot admit all
the required individual conjugates. No bound on their orders, poles
or degrees is involved. Nontrivial quadratic fields are included.

For a projective corner or another actual algebra, H/Q finiteness
and the unit condition (2) require proof; scalar Weyl units are not
being assumed for an arbitrary endomorphism algebra. Without finite
generation of B the etale and normalization argument below is not
asserted.

## Algebraic translation extensions are differentially simple

We first prove a commutative algebra fact, independently of the Weyl
realization. Suppose B is a domain containing L[t1,...,tm], Frac B is
algebraic over L(t), and B has derivations extending all coordinate
partials. Then every ideal stable under those derivations is zero
or the unit ideal.

Indeed, if 0!=b lies in such an ideal I, an algebraic equation with
nonzero constant term, after clearing scalar rational denominators,
gives
\[
                 a_db^d+\cdots+a_1b+a_0=0,\qquad
                 a_j\in L[t],\quad a_0\ne0.               \tag{3}
\]
Thus a0 belongs to I intersect L[t]. This intersection is stable
under every partial derivative. Finite Taylor coefficient extraction
makes every nonzero such polynomial ideal the unit ideal, as proved
in QFC1. Hence I=B. No finite-generation assumption was needed for
this simplicity statement.

## Finite generation makes the centralizer etale

Now impose finite generation. Omega_(B/P) is a finite B-module and
has connections
\[
         \nabla_i(b\,dc)=(\nabla_i b)\,dc+b\,d(\nabla_i c). \tag{4}
\]
Separability of H/Q gives Omega_(B/P) tensor_B H=0.
Since Omega is finite, its annihilator is a nonzero ideal of B.
The connections preserve that annihilator: differentiate a*omega=0.
Differential simplicity forces the annihilator to be B, so
\[
                     \Omega_{B/P}=0.                     \tag{5}
\]
The exact sequence of differentials and the derivations give mutually
inverse maps
\[
 B^m \xrightarrow{e_i\mapsto dF_i}\Omega_{B/L}
       \xrightarrow{db\mapsto(\nabla_i b)_i} B^m.           \tag{6}
\]
Surjectivity of the first map is (5), and the composite is the identity.
At every closed point the residue field is finite separable over L,
so (6) makes the cotangent-space dimension m. The integral finite-type
L-algebra B has dimension m, and every closed local ring has dimension
m. Thus every such point is regular; the Jacobian criterion in
characteristic zero makes B smooth over L. The differential of the
map to affine m-space is an isomorphism by (6), so that map is etale.
This proves flatness as well as unramifiedness; vanishing relative
differentials alone would not have sufficed.

## Weyl units retain all base divisors

An etale morphism is open. Let U be the image of Spec B in Spec P.
If the generic point of an irreducible divisor q(F)=0 were absent
from U, its whole divisor would be absent: an open set containing
any point of an irreducible divisor contains its generic point.
Then q(F) would vanish at no prime of B, hence be a unit. By (2),
q would have to be constant, a contradiction. Therefore
\[
       \operatorname{Spec}P\setminus U
                         \text{ has codimension at least two}. \tag{7}
\]
This asserts no rational point over a base field; it retains geometric
points over every generic divisor, possibly with nontrivial separable
residue-field extensions.

Let C be the integral closure of P in H. It is finite over P, since P
is a polynomial ring over a field and H/Q is finite. Smoothness makes
B normal, so C is contained in B. Zariski's Main Theorem identifies
\[
          \operatorname{Spec}B\hookrightarrow\operatorname{Spec}C
                  \longrightarrow\operatorname{Spec}P     \tag{8}
\]
as an open immersion followed by a finite map. We use the
[normalization form of Zariski's Main Theorem, Stacks02LR](https://stacks.math.columbia.edu/tag/02LR).
In this normal dominant setting its normalization is precisely C.
The result is an external geometric input, not a finite experiment.

At each height-one base prime, (7) gives a point in Spec B over that
prime. There (8) identifies the local normalization with the etale
source, so this is an unramified prime of C. In characteristic zero
all residue extensions are separable and inertia is tame. Over a
geometric generic divisor the cycle lengths of an inertia generator
are the ramification indices, repeated by residue degree. Hence
\[
       \text{each codimension-one inertia generator has a fixed point}.
                                                               \tag{9}
\]
This does not supply one sheet fixed by the whole monodromy group.

## Birational and Galois fields collapse to polynomiality

If H=Q, then C=P because P is normal. Equation (8) makes B an affine
open subring in the same function field. Every b in B is regular
at all height-one points of P by (7). Since P is a factorial normal
domain,
\[
                  B\subset\bigcap_{\operatorname{ht}\mathfrak p=1}
                                      P_{\mathfrak p}=P.  \tag{10}
\]
Thus B=P. This is the birational case; no module-finiteness hypothesis
on B/P was needed.

Suppose next that H/Q is Galois. At a base prime divisor the Galois
group acts transitively on the primes of C above it, preserving
ramification indices. Since at least one is unramified, all are.
This standard DVR fact is recorded in
[Stacks09E3](https://stacks.math.columbia.edu/tag/09E3).
Consequently the finite map Spec C->Spec P is unramified in
codimension one. The source is normal and the target regular;
[purity of the branch locus, Stacks0BMB](https://stacks.math.columbia.edu/tag/0BMB)
makes this finite map etale everywhere.

Finite etale covers of affine space over characteristic-zero fields
are constant, using the same explicit external covering theorem as
QFC1, [Milne, Remark3.3(c), Theorem3.4 and the fundamental group sequence](https://www.jmilne.org/math/CourseNotes/LECc.pdf#page=27).
Thus C is C0 tensor_L P for a finite etale L-algebra C0. Because C
is a domain, C0 is a field. It lies inside B, and its elements are
annihilated by all the nabla_i, being separably algebraic over L.
The joint constants in B are L, so C0=L, H=Q, and the birational
argument gives B=P.

For scalar Weyl algebras one can also exclude C0 directly: an element
algebraic over L in A_m(L) has Bernstein degree zero and is scalar.
The formal-horizontal-constant argument permits the more general
algebra hypotheses in (1)--(2).

Return now to the integral algebra. QFC1's valuation argument gives
\[
 \operatorname{Cent}_{E}(F)=K[[h]][F],\qquad
                   h^{-1}[T_i,T_j]\in K[[h]][F].          \tag{11}
\]
Its finite polynomial de Rham homotopy kills all mixed curvature
inside E, and is h-small when the classical tuple is flat.
Thus Galois or birational field data furnish a second positive
completion certificate after the individual conjugates are found.

## Degree and inertia restrictions without a degree cutoff

If all hypotheses hold and B is not polynomial, then
\[
                         [H:Q]\ge3,\qquad H/Q\text{ is not Galois}.
                                                               \tag{12}
\]
Every quadratic extension in characteristic zero is Galois; degree
one was treated separately. A cubic nonpolynomial survivor must have
normal closure S3 and may have generic inertia type (2,1), but not
(3). Two transpositions can each fix a different sheet and generate
transitive S3, so no common branch or algebraic rank-one extraction
follows from (9).

The lower degree three is sharp in rank at least three.
The retained rescaled cotangent examples of
[QNI1](IDENTITY_LIFTS_AND_UNBOUNDED_INVERSES.md) have B=L[x1,x2,x3]
and P=L[f_(h^N)] with identity classical reduction and all full
relations. [QDI1's cubic chart](QUANTUM_DIVISION_AND_COMPLETION.md)
proves their position-field extension has degree three. With its
unnormalized output parameters A,B,C and source parameter w, the
equation is
\[
                        Cw^3-2w^2+Bw-2A=0.               \tag{13}
\]
Its discriminant is
4B^2-4CB^3+72ABC-64A-108A^2C^2.
As a quadratic in A its discriminant is
64*(4-3BC)^3, nonzero, so it is squarefree in A over L(B,C)
and not a square. Therefore the normal closure is S3. These
centralizers are etale but nonfinite over their position-image rings,
as QFC1's cotangent boundary shows. This is a reused higher-rank
control, not a rank-two or Long construction.

The degree in (12) is [Frac B:L(F_i)]. It is distinct from QDI1's
left/right indices of whole noncommutative Weyl division rings and
their completed degrees. No equality between them is inferred.

## Two controls that protect the hypotheses

First, a Galois cover with nonconstant units really escapes. Set
\[
 B_0=K[[h]][t_1,t_2,z,z^{-1},(z+1)^{-1}]
                                  /(z^2-1-h t_1).         \tag{14}
\]
It is a domain, h-torsion-free and h-separated by its faithful
formal branch z=(1+h*t1)^(1/2). Reduction at h=0, after inverting
z+1, is K[t1,t2]. Its commuting derivations are
\[
 \nabla_1t_1=1,\quad \nabla_1z=h/(2z),\quad
 \nabla_2t_2=1,\quad \nabla_1t_2=\nabla_2t_1=\nabla_2z=0.  \tag{15}
\]
Form the Ore algebra with [pi,b]=h*nabla_i(b) and [p1,p2]=0.
Positions ti and individual mates -pi already form a flat full tuple.
Their centralizer is B0, generically quadratic Galois over L[t1,t2].
But 1+h*t1=z^2 is a nonconstant unit. This is an abstract Ore
boundary algebra, not a finite scalar Weyl algebra. The h-adic
branch has coefficient binomial(1/2,n)*t1^n at h^n, nonzero for
every n>=1, so its spatial degrees are unbounded. It demonstrates that (2)
cannot be dropped, even with a polynomial completed centralizer,
correct classical graph and all full relations.

Second, leading-symbol birationality cannot replace actual field
birationality. In the commutative polynomial ring K[x,y], take
\[
                u=x^2y+y,\qquad v=xy,\qquad
                         v x^2-u x+v=0.                  \tag{16}
\]
The leading forms x^2*y and x*y generate the full rational field:
x=u_lead/v_lead and y=v_lead^2/u_lead. The actual extension has
degree two, since the quadratic discriminant u^2-4v^2 has two
simple prime factors. Its Jacobian is y*(x^2-1), so it is a
non-Keller control. Replacing QFC1's full module-finite leading
centralizer certificate by a leading field-degree calculation
would lose the crucial condition.

A scalar quantum control is F1=x+h^2*x^2, F2=y. Their classical
joint centralizer is polynomial, while their actual joint centralizer
is L[x,y]. Its field over L[F1,F2] is quadratic Galois, with involution
x->-h^-2-x. Hence simultaneous finite individual mates are impossible.
The rational mate -(1+2h^2*x)^-1*p_x exists in the division algebra;
it inverts precisely a missing divisor. A direct leading-order
check also forces a finite mate to have momentum order at most one,
with linear coefficient -(1+2h^2*x)^-1, which is not polynomial.
This recovers an older
[QCR1 boundary](CENTRALIZER_NORMAL_FORM.md) by the new all-partner
geometric criterion.

## Remove parameter poles inside the actual generic algebra

This lemma needs only QFC1's algebra and polynomial classical
centralizer hypotheses, not finite generation, a field-degree
condition or the unit condition. Suppose actual generic Ti in E[h^-1]
satisfy [Fj,Ti]=h*delta_ij. There are unique polynomial polar parts
\[
 T_i=U_i+a_i(F),\qquad U_i\in E,\qquad
       a_i\in h^{-1}K[h^{-1}][t_1,\ldots,t_m].             \tag{17}
\]
To construct them, if Ti has pole order N>=1, multiply by h^N.
Its reduction Poisson-commutes with all fj, since
[Fj,h^N Ti]/h is divisible by h. Write that reduction uniquely
as a polynomial q(fi), subtract h^-N*q(F), and the pole order
drops by at least one. After at most N steps Ti is integral.
There are finitely many polynomial corrections, so this remains
in the actual finite algebra. Uniqueness follows because a
nonzero polynomial in F with only negative h coefficients cannot
lie in E: its first coefficient is nonzero in the independent fi.

Subtracting ai(F) does not change the translation derivations on
B. Thus the geometric theorem above can start with actual generic
individual mates; one first applies (17). Their classical images
can be matched to any prescribed classical individual conjugates:
each difference lies in K[f], and its polynomial lift can be
subtracted finitely.

There is an extra full-tuple statement with no centralizer-finiteness
hypothesis. Direct commutation gives
\[
 h^{-1}[U_i,U_j]=h^{-1}[T_i,T_j]
                         +\partial_i a_j-\partial_j a_i. \tag{18}
\]
If the original curvature divided by h is integral, then the last
term is integral too. But its F-coordinate coefficients have only
negative h powers. By the same independence argument it is zero.
In particular, a generic flat tuple stays flat after (17), with
no additional curvature correction.

If its integral classical tuple is to equal a prescribed flat
tuple t_i^0, the differences are polynomials
\[
             \bar U_i=t_i^0+q_i(f),\qquad
                         \partial_i q_j-\partial_j q_i=0. \tag{19}
\]
Closedness follows from the two classical flatness relations.
Replacing Ui by Ui-qi(F) preserves every quantum relation exactly
and attains the prescribed full classical tuple. Hence an actual
generic full tuple with the regular Long positions already suffices
for a finite integral full Long lift. This is a parameter-pole
descent theorem, not a bounded comparison between arbitrary lifts.

Only h is inverted in E[h^-1]. A rational coefficient such as
(1+2h^2*x)^-1, an operator denominator or an arbitrary formal series
does not become an actual finite-PBW element by this lemma. The
quadratic scalar rational mate above remains outside its scope.
The finite controls normalize explicit generic tuples with
momentum-dependent positions, two polar orders and a prescribed
classical shift. They check both flatness preservation and exact
preservation of a nonzero integral curvature.

## Evidence and remaining work

The checker verifies the retained QFC1/QDI1 hashes, the quadratic
and cubic field equations, the nonsquare discriminant witnesses,
permutation inertia controls, the formal quadratic branch and its
unbounded coefficient degrees, the compatible Ore derivations, the
birational-leading countercontrol and the scalar quadratic rational
mate. Universal differential simplicity, etaleness, boundary
retention, normalization, purity and all-rank completion are written
proofs using the stated external geometric results. Finite controls
do not prove finite generation, determine a Long centralizer field
or establish all-representative nonlifting. No independent theorem
replay, formal verification, external review or novelty priority
is asserted.

For a new actual scalar Long pair, determine the full centralizer,
its finite generation and actual fraction field before trying
another frame. A birational field plus individual conjugates now
suffices for finite full completion; a nontrivial Galois field
forbids those conjugates. A non-Galois extension of degree at least
three remains possible and must still pass every global spectral,
curvature and finite-carrier condition.

Run the declared exact controls from the repository root:

```sh
.venv/bin/python research/dixmier/dixmier-formal-lifting/check_etale_centralizers.py
```
