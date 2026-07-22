# The degree-eighteen triple intersection

Work over a field of characteristic zero.  This note studies the component
`C_(3,4)` over the exact collision stratum `E_(6,6,6)`.  It proves the
pairwise and ordered-triple transverse intersection algebras on a nonempty
open of the ordered-root chart.  The reconstruction of the entire completed
three-branch component ring, and hence its total conductor equalizer, remains
open.

## 1. Three normalization branches

At a collision root of multiplicity `m`, a normalization branch chooses an
allocation

\[
2i+3j=m,
\]

where `i` is its multiplicity in `Q` and `j` its multiplicity in `R`.  For a
sixfold root the two choices are `(3,0)` and `(0,2)`.  Consequently

\[
[U^3V^4](U^3+V^2)^3=3.
\]

Thus the normalization of `C_(3,4)` has three points over the generic
`(6,6,6)` collision polynomial.  They are indexed by the root assigned to
`Q`; the other two roots are assigned to `R`.

An exact rational witness is

\[
M=(W+3/4)^6(W-1)^6(W-2)^6.                       \tag{1}
\]

It satisfies the normalized incidence equation

\[
\Phi=M(1)-M(0)-M'(0)=0
\]

and has

\[
D=M'(1)-M'(0)=729/64\ne0.
\]

For the three root coordinates `(-3/4,1,2)`, the derivatives of `Phi` are

\[
-243/2,\qquad -2187/32,\qquad -2187/128.          \tag{2}
\]

In particular, the incidence cuts the three-dimensional root-position space
smoothly.

## 2. Tangent geometry

Restrict the seven coefficient parameters of a monic cubic `Q` and monic
quartic `R` by `dPhi=0`, then differentiate the normalized coefficient map

\[
(Q,R)\longmapsto \rho(Q^2R^3)
\]

on each of the three allocation sheets at (1).  Exact rational linear
algebra gives:

\[
\begin{array}{c|c}
\text{space} & \text{dimension}\\ \hline
\text{image tangent of each sheet} & 6\\
\text{span of any two sheet tangents} & 7\\
\text{intersection of any two sheet tangents} & 5\\
\text{span of all three sheet tangents} & 8\\
\text{common tangent of all three sheets} & 5.
\end{array}                                      \tag{3}
\]

The triple collision stratum has dimension `3-1=2`.  Hence its triple
intersection has three excess common tangent directions, exactly one for
each sixfold root.

There is a useful distinction between the pairwise and triple reduced
supports.  Two allocation sheets agree at the third sixfold root.  Their
common quadratic factor may split, so their reduced intersection moves in
the collision stratum

\[
E_{(6,6,3,3)},
\]

which has dimension `4-1=3`.  The five-dimensional pairwise common tangent
therefore consists of three reduced directions and two excess directions.
All three allocation sheets can remain equal only when every one of the
three roots stays sixfold, so the reduced triple support is `E_(6,6,6)`.

## 3. Pairwise transverse algebra

For each of the three pairs, take a linear slice transverse to its
three-dimensional `E_(6,6,3,3)` support.  The simultaneous lifting map has
two excess common tangent directions.  After quotienting arbitrary
second-order tangent corrections, its scalar quadratic obstruction has rank
two and is a zero-dimensional complete intersection of two quadrics.

Its Groebner initial ideal has four standard monomials.  On the other hand,
the stronger polynomial equality separates at the two exchanged sixfold
roots into two copies of the universal block

\[
Q^2=T^3.
\]

Depressing the cubic as in the degree-twelve calculation gives the exact
block ideal

\[
(a,v,2u-3b,b^2),
\]

so the stronger-equality pairwise slice is already a tensor product of two
dual numbers and has length four.  The quadratic upper bound and this closed
length-four subsystem agree.  Therefore, for all three pairs,

\[
\boxed{
D_{ij}^{\mathrm{tr}}
 \simeq k[\epsilon_i,\epsilon_j]/
 (\epsilon_i^2,\epsilon_j^2),
\qquad \operatorname{length}D_{ij}^{\mathrm{tr}}=4.}       \tag{4}
\]

The word *transverse* in (4) refers to the actual pairwise support
`E_(6,6,3,3)`, not only to the smaller triple stratum.

## 4. The triple second-jet obstruction

Let `V_i` be the `17 x 6` differential matrix on the `i`-th sheet.  A triple
of source tangent vectors has the same target image precisely when it lies in
the kernel of

\[
L=
\begin{pmatrix}
V_1&-V_2&0\\
V_1&0&-V_3
\end{pmatrix}.                                  \tag{5}
\]

At (1), `L` has rank thirteen and a five-dimensional kernel.  Quotient its
two-dimensional collision-stratum tangent plane and retain transverse
coordinates `x,y,z`.  The second derivatives of the two differences in (5),
modulo arbitrary accelerations on the three sheets, give a three-dimensional
space of scalar quadrics.  One exact basis, up to nonzero scalar multiples,
is

\[
\begin{aligned}
g_1={}&1178793x^2+16765056xy-25147584xz-3407147008y^2\\
     &\quad+4698820608yz-1856240640z^2,\\
g_2={}&(4y-3z)^2,\\
g_3={}&z(8y+z).
\end{aligned}                                    \tag{6}
\]

In graded reverse lexicographic order their ideal has Groebner basis, after
rescaling,

\[
\begin{aligned}
z^3,&\\
x^2+\frac{128}{9}xy-\frac{64}{3}xz+\frac{2560}{27}z^2,&\\
y^2+\frac34z^2,&\\
yz+\frac18z^2.&
\end{aligned}                                    \tag{7}
\]

The initial monomials are `z^3,x^2,y^2,yz`.  Hence the standard monomials are

\[
1,z,z^2,y,x,xz,xz^2,xy,                         \tag{8}
\]

and the quadratic tangent-cone quotient has length eight.  It follows that
the full affine-difference triple intersection has transverse length at most
eight.

## 5. The length-eight sandwich

Impose the stronger equations

\[
M_1=M_2=M_3.
\]

Hensel separation at the three distinct roots of (1) makes the local
factorization equations independent.  At each root one sheet supplies the
monic cubic factor and the other two sheets supply the same monic quadratic
factor.  The local equation is the universal `Q^2=T^3` block, hence one dual
number.  The three roots therefore give the closed transverse subsystem

\[
D_0=k[\epsilon_1,\epsilon_2,\epsilon_3]/
       (\epsilon_1^2,\epsilon_2^2,\epsilon_3^2),  \tag{9}
\]

of length eight.  Equation (2) removes one smooth root-position direction
and introduces no relation among the three nilpotent directions.

Equality of normalized seeds is weaker than equality of the `M_i`: it only
requires every difference `M_i-M_j` to be affine in `W`.  Thus (9) is a
closed subscheme of the full triple intersection and supplies a length-eight
lower bound.  The second-jet calculation (6)--(8) supplies the matching upper
bound.  The induced surjection of transverse Artin rings is an isomorphism.

Therefore, on a nonempty open of the ordered `(6,6,6)` chart,

\[
\boxed{
D_{123}^{\mathrm{tr}}
 \simeq k[\epsilon_1,\epsilon_2,\epsilon_3]/
 (\epsilon_1^2,\epsilon_2^2,\epsilon_3^2),
\qquad \operatorname{length}D_{123}^{\mathrm{tr}}=8.}       \tag{10}
\]

The cubic socle class

\[
\epsilon_1\epsilon_2\epsilon_3\ne0              \tag{11}
\]

is the first genuinely three-way infinitesimal datum.  It disappears from
every pairwise transverse slice.

## 6. What remains: the total conductor equalizer

Let `B_1,B_2,B_3` be the three completed regular normalization sheets.  The
desired completed component ring should be recovered as the simultaneous
coefficient equalizer

\[
A\longrightarrow B_1\oplus B_2\oplus B_3.
\]

For two minimal primes, the ordinary fiber-product exact sequence used in
degree twelve makes this automatic.  For three minimal primes, pairwise
congruence need not be sufficient: exactness in the middle of the Cech-type
complex

\[
0\longrightarrow A\longrightarrow
\bigoplus_i B_i\longrightarrow
\bigoplus_{i<j}D_{ij}\longrightarrow D_{123}     \tag{12}
\]

requires a distributivity or effective-descent argument.  Establishing (12)
for the completed coefficient ring is the genuinely new remaining theorem.

This can be reduced to one explicit ideal identity.  Write
`p_1,p_2,p_3` for the three minimal primes of `A`.  After gluing residues on
the first two branches, a residue on the third branch can also be glued
precisely when

\[
(\mathfrak p_1+\mathfrak p_3)\cap
(\mathfrak p_2+\mathfrak p_3)
=\mathfrak p_3+(\mathfrak p_1\cap\mathfrak p_2). \tag{13}
\]

Indeed, pairwise compatibility places the remaining difference in the
left-hand side; changing the first two glued lifts changes it by the
right-hand side.  Thus the only obstruction to the middle exactness of (12)
is the distributivity-defect module

\[
\mathcal H=
\frac{(\mathfrak p_1+\mathfrak p_3)\cap
      (\mathfrak p_2+\mathfrak p_3)}
     {\mathfrak p_3+(\mathfrak p_1\cap\mathfrak p_2)}.       \tag{14}
\]

Away from `E_(6,6,6)` at most two of these three sheets remain, so
`mathcal H` is supported on the triple stratum.  Computing `mathcal H` in
the completed coefficient ring is the sharply formulated next step.

If (12) is exact, the conductor component on the `i`-th normalization sheet
is

\[
\mathfrak c_i=
\ker(B_i\to D_{ij})\cap\ker(B_i\to D_{ik}),       \tag{15}
\]

and the pairwise conductor maps together with the triple compatibility map
determine the entire completed local ring.

The tangent calculation makes the resulting conductor formula concrete.  On
a fixed sheet `B_i`, the two pairwise common tangent spaces coincide with the
projection of the five-dimensional triple common tangent space.  Hence the
two pairwise quotient kernels have the same linear generator.  Conditional
on the formal straightening implicit in (12), choose parameters

\[
B_i=k[[t_1,t_2,e_1,e_2,e_3,w_i]],
\]

where `t_1,t_2` parameterize `E_(6,6,6)` and the `e_r` are the three transfer
directions.  For `{i,j,k}={1,2,3}`, the expected pairwise kernels on the
`i`-th sheet are

\[
I_{ij}=(w_i,e_i^2,e_j^2),\qquad
I_{ik}=(w_i,e_i^2,e_k^2).                         \tag{16}
\]

After slicing by `t_1,t_2`, their sum recovers the proved transverse triple
quotient:

\[
B_i/(I_{ij}+I_{ik})
\simeq k[[t_1,t_2,e_1,e_2,e_3]]/(e_1^2,e_2^2,e_3^2).
\]

Their intersection gives the predicted conductor component

\[
\boxed{
\mathfrak c_i=I_{ij}\cap I_{ik}
=(w_i,e_i^2,e_j^2e_k^2).}                         \tag{17}
\]

Thus the equalizer theorem would immediately produce an explicit conductor,
including the quartic equation where the two pairwise conductor components
meet on one normalization sheet.

The completed component ring itself is reduced; the nilpotents in (4) and
(10) belong to its scheme-theoretic branch-intersection quotients.  The
nonreduced pairwise conductor already shows that the component is not
seminormal along the adjacent generic `E_(6,6,3,3)` locus.  Formula (12),
rather than another nilpotent search, is what is needed for a complete
pointwise conductor description at the triple stratum.

## 7. Finite-jet evidence for the equalizer

The defect module (14) has an important finiteness property.  It is
annihilated by the sum of the three minimal primes:

\[
(\mathfrak p_1+\mathfrak p_2+\mathfrak p_3)\mathcal H=0.   \tag{18}
\]

For example, if `x` represents a defect class and
`x in p_2+p_3`, then

\[
\mathfrak p_1x\subset
\mathfrak p_1\mathfrak p_2+\mathfrak p_3
\subset(\mathfrak p_1\cap\mathfrak p_2)+\mathfrak p_3;
\]

the other two primes are identical.  Thus `mathcal H` is a module over the
triple-intersection ring.  On the transverse slice (10), its possible
nilpotent degrees stop at degree three.

This makes low formal jets unusually informative.  At the rational witness
(1), construct the normalized coefficient map on each smooth branch modulo
the fourth power of its maximal ideal.  Inside the direct sum of the three
branch jets, take the subalgebra generated by the seventeen normalized seed
coefficients.  Compute its three branch-prime kernels and apply (14).  Exact
rational linear algebra gives:

\[
\begin{array}{c|c|c|c|c}
\text{jet order}&\dim B_i&\dim A&\dim\mathfrak p_i&
\dim(\text{numerator})/\dim(\text{denominator})\\ \hline
1&7&9&2&3/3\\
2&28&45&17&25/25\\
3&84&158&74&108/108.
\end{array}                                                     \tag{19}
\]

The distributivity defect is therefore zero in every order through the
cubic socle degree.

The successive Hilbert dimensions of the generated coefficient algebra are

\[
1,\quad8,\quad36,\quad113.                              \tag{20}
\]

They agree exactly with the Cech prediction

\[
\boxed{
\operatorname{Hilb}_A(z)=
\frac{3}{(1-z)^6}
-\frac{3(1+z)^2}{(1-z)^3}
+\frac{(1+z)^3}{(1-z)^2}.}                         \tag{21}
\]

Here the three terms are respectively the three smooth six-dimensional
branches, the three expected pairwise rings

\[
k[[t_1,t_2,s,\epsilon_i,\epsilon_j]]/
 (\epsilon_i^2,\epsilon_j^2),
\]

and the triple ring

\[
k[[t_1,t_2,\epsilon_1,\epsilon_2,\epsilon_3]]/
 (\epsilon_1^2,\epsilon_2^2,\epsilon_3^2).
\]

This calculation includes the class
`epsilon_1 epsilon_2 epsilon_3`; there is no hidden conductor relation at the
first order where genuinely triple multiplication occurs.

The remaining proof obligation is now narrow.  One must show that formation
of (14) commutes with a transverse base change to `E_(6,6,6)`, or directly
prove that the two collision-stratum parameters are regular on the relevant
kernel and cokernel modules.  Then (18), the transverse nilpotency bound, and
the order-three vanishing in (19) imply `mathcal H=0` by Nakayama.  Without
that base-change lemma, (19)--(21) are exact finite-jet evidence rather than
an all-order completed-ring theorem.

Run the exploratory calculation with

```bash
python scripts/explore_degree18_conductor_equalizer.py --order 1
python scripts/explore_degree18_conductor_equalizer.py --order 2
python scripts/explore_degree18_conductor_equalizer.py --order 3
```

The cubic calculation is intentionally not part of the default verification
target because exact row reduction takes several minutes.

## 8. Executable transverse-intersection certificate

Run

```bash
python scripts/verify_degree18_triple_intersection.py
```

The script checks the branch points, the admissible witness, every tangent
rank in (3), all three pairwise transverse length-four calculations, the
triple obstruction (6)--(8), the universal sixfold block, and the matching
length-eight upper and lower bounds.
