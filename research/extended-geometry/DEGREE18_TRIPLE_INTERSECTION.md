# The degree-eighteen triple intersection

Work over a field of characteristic zero.  This note studies the component
`C_(3,4)` over the exact collision stratum `E_(6,6,6)`.  It proves the
pairwise and ordered-triple transverse intersection algebras on a nonempty
open of the ordered-root chart.  An excess-conormal calculation then proves
the completed three-branch conductor equalizer on a possibly smaller
nonempty open.

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
length-four subsystem agree.  Equivalently, the
[Kuranishi nilpotence cutoff theorem](DEFECT_SYMBOL_APOLARITY.md#2-kuranishi-nilpotence-cutoff-theorem)
supplies the upper bound and the stronger-equality subsystem supplies the
matching lower bound.  Therefore, for all three pairs,

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
the Kuranishi nilpotence cutoff theorem gives transverse length at most
eight and fourth-power-zero transverse maximal ideal for the full
affine-difference triple intersection.

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
bound through the cutoff theorem.  The induced surjection of transverse
Artin rings is an isomorphism.  This is the three-variable quadratic
sandwich.

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

## 6. The total conductor equalizer

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

requires a distributivity or effective-descent argument.  Section 7.2 proves
this exactness generically by an excess-intersection calculation.

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
`mathcal H` is supported on the triple stratum.  Section 7.2 proves that it
vanishes on a nonempty open.

By (12), the conductor component on the `i`-th normalization sheet
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
two pairwise quotient kernels have the same linear generator.  Apply the
simultaneous formal straightening proved in Section 7.3 and choose parameters

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

Thus the equalizer theorem gives the intrinsic conductor formula (15).  The
simultaneous straightening displayed in (16) gives the monomial form (17),
including the quartic equation where the two pairwise conductor components
meet on one normalization sheet.

The completed component ring itself is reduced; the nilpotents in (4) and
(10) belong to its scheme-theoretic branch-intersection quotients.  The
nonreduced pairwise conductor already shows that the component is not
seminormal along the adjacent generic `E_(6,6,3,3)` locus.  Formula (12),
rather than another nilpotent search, is what is needed for a complete
pointwise conductor description at the triple stratum.

## 7. Finite jets and the excess-conormal equalizer proof

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

### 7.1 A relative transverse-jet refinement

The constant-rank route is better served by jets relative to the ordered-root
stratum than by the absolute jets in (19).  Write `c_d` for the normalized
coefficient of `W^d`.  At the witness (1), the restrictions of `c_2,c_3` to
`E_(6,6,6)` have independent differentials.  In the tangent basis obtained by
solving `dPhi=0` for the first root coordinate, their exact Jacobian is

\[
\begin{pmatrix}
-14&-11/4\\
56/3&121/24
\end{pmatrix},
\qquad \det=-77/4.                              \tag{22}
\]

After shrinking the ordered-root chart, `c_2,c_3` are therefore etale base
coordinates.  Cut all three normalization sheets by the same equations
`c_2=c_2(M)` and `c_3=c_3(M)`, and truncate the resulting four transverse
variables in degree four.  The exact calculation then gives

\[
\begin{array}{c|c|c|c|c}
\text{jet order}&\dim B_i&\dim A&\dim\mathfrak p_i&
\dim(\text{numerator})/\dim(\text{denominator})\\ \hline
1&5&7&2&3/3\\
2&15&28&13&19/19\\
3&35&77&42&61/61.
\end{array}                                                     \tag{23}
\]

Thus the relative transverse defect also vanishes through the cubic socle
degree.  Its successive algebra dimensions are `1,6,21,49`, exactly the
coefficients through degree three of the base-sliced Cech prediction

\[
\frac{3}{(1-z)^4}
-\frac{3(1+z)^2}{1-z}
+(1+z)^3.                                         \tag{24}
\]

This relative calculation supports a clean generic-open argument.  Over the
ordered-root coordinate ring, the degree-at-most-three branch jets are finite
free after localizing the branch-coordinate Jacobians.  The generated
algebra, its three branch kernels, every sum and intersection, and the final
defect are kernels, images, and cokernels of finite matrices.  Localizing at
the nonzero witness minors makes all these ranks constant.  Formation of the
finite modules then commutes with residue-field base change, and (23) makes
the universal *truncated* defect zero on a nonempty open.

There is still one logically separate comparison to prove.  Let
`mathcal H^{<=3}` denote the defect formed after relative transverse
truncation.  One needs, at the geometric generic point, an injection

\[
\mathcal H\longrightarrow \mathcal H^{\le3},       \tag{25}
\]

or an equivalent strictness statement for the transverse filtration on the
Cech complex.  Neither (18) nor the relation `mathfrak n^4=0` in the triple
quotient alone implies (25): quotienting can enlarge intersections and can
also enlarge the denominator in (14).  Sufficient replacements for (25)
include strictness of the filtered kernels and images through degree three,
vanishing of the relevant Rees-homology torsion, or the corresponding
Tor-independence statements for the intersection presentations.  Once this
comparison is established, (18), `mathfrak n^4=0`, and (23) give
`mathcal H=0` immediately.

The Tor formulation is explicit.  Let `S` be the completed regular local ring
of the ambient normalized-coefficient space and let `I_i` be the ideal of the
`i`-th branch, so `B_i=S/I_i`.  Put

\[
C_{12}=S/(I_1\cap I_2),\qquad D_{12}=S/(I_1+I_2).
\]

Tensor the two-branch Mayer--Vietoris sequence

\[
0\longrightarrow C_{12}\longrightarrow B_1\oplus B_2
\longrightarrow D_{12}\longrightarrow0
\]

over `S` with `B_3`.  The kernel of

\[
C_{12}\otimes_S B_3\longrightarrow
(B_1\otimes_S B_3)\oplus(B_2\otimes_S B_3)
\]

is exactly (14).  The Tor long exact sequence therefore gives

\[
\boxed{
\mathcal H\simeq\operatorname{coker}\!\left(
\operatorname{Tor}_1^S(B_1,B_3)\oplus
\operatorname{Tor}_1^S(B_2,B_3)
\longrightarrow
\operatorname{Tor}_1^S(D_{12},B_3)
\right).}                                             \tag{26}
\]

Thus one may replace the filtered comparison (25) by generic surjectivity of
the Tor map in (26).  This pinpoints which Tor-independence statement is
needed; asking for all pairwise Tor groups to vanish would be unnecessarily
strong and is incompatible with the nontransverse pairwise intersections.

### 7.2 Excess conormals prove the generic equalizer

We use the following elementary Koszul lemma.  Let `S` be regular local, let
`I=(f_1,...,f_a)` be generated by a regular sequence, and put `R=S/(I+J)`.
Suppose that the image of `I` in `S/J` is a complete-intersection ideal of
height `c`.  After an invertible change of the `f_i` over `S/J`, their images
have the form

\[
g_1,\ldots,g_c,0,\ldots,0,
\]

where the `g_i` are a regular sequence.  Tensoring the Koszul resolution of
`S/I` with `S/J` consequently gives

\[
\operatorname{Tor}_q^S(S/I,S/J)
\simeq \bigwedge^q_R E,
\qquad \operatorname{rank}_R E=a-c.              \tag{27}
\]

Here `E` is the kernel of the excess conormal map.  The construction is
functorial for inclusions of complete-intersection ideals.  Therefore the map
in (26) is the map of excess conormal modules

\[
E_{13}|_{D_{123}}\oplus E_{23}|_{D_{123}}
\longrightarrow E_{(12),3}.                       \tag{28}
\]

Use the common etale coordinates `t_1=c_2,t_2=c_3` from (22) at the rational
witness and write the full ambient ring as
`S=QQ[[t_1,t_2,y_1,...,y_15]]`.  The transverse slice in the computation only
sets `t_1,t_2` to their witness values; it computes the residue map of the
full conormal modules, rather than base-changing the defect module.  The
branch ideals have height eleven.  By (4), their pairwise sums have height
fourteen and complete-intersection quotient

\[
QQ[[t_1,t_2,s,\epsilon_i,\epsilon_j]]/
 (\epsilon_i^2,\epsilon_j^2).
\]

By (10), the triple sum has height fifteen and complete-intersection quotient

\[
QQ[[t_1,t_2,\epsilon_1,\epsilon_2,\epsilon_3]]/
 (\epsilon_1^2,\epsilon_2^2,\epsilon_3^2).
\]

Thus (27) applies to every Tor module in (26).  The pair excess modules have
rank eight over their pairwise rings, while `E_(12),3` has rank ten over the
triple ring.

The fibers of these conormal modules and all maps between them depend only on
the branch two-jets.  Exact rational linear algebra at (1) gives

\[
\begin{array}{c|c}
\text{conormal fiber}&\text{dimension}\\ \hline
I_i/\mathfrak m I_i&11\\
(I_i+I_j)/\mathfrak m(I_i+I_j)&14\\
(I_1+I_2+I_3)/\mathfrak m(I_1+I_2+I_3)&15\\
E_{ij}\otimes QQ&8\\
E_{(12),3}\otimes QQ&10.
\end{array}                                             \tag{29}
\]

The map (28) has rank ten; one exact maximal minor is

\[
\frac{2401}{256}\ne0.                                  \tag{30}
\]

It is therefore surjective on the residue field.  Since its target is free of
rank ten over the local triple-intersection ring, Nakayama makes (28)
surjective over that entire Artin ring.  Equations (26)--(28) give
`mathcal H=0` at the witness.

The branch, pair, and triple complete-intersection conditions and the
nonvanishing minor (30) are open conditions on the ordered-root chart.
Shrinking once more gives `mathcal H=0` throughout a nonempty open.  Hence
(12) is exact there, proving the generic completed three-sheet conductor
equalizer.  Notice that this argument bypasses the unproved filtered
comparison (25); the relative cubic calculation (23) remains an independent
check of the resulting theorem.

### 7.3 Simultaneous straightening and the quartic conductor

It remains to justify the simultaneous coordinates used in (16).  Work over
the geometric generic point of `E_(6,6,6)` on one branch `B_i`, and write

\[
I=\ker(B_i\to D_{ij}),\qquad
J=\ker(B_i\to D_{ik}),\qquad K=I+J.
\]

The ideals `I,J` are height-three complete intersections and `K` is a
height-four complete intersection.  Under the compatible identifications
(4) and (10), the image `K/J` is generated by the square of the free pairwise
parameter and is a nonzerodivisor.  Applying (27) inside the regular ring
`B_i` therefore shows that

\[
F=\operatorname{Tor}_1^{B_i}(B_i/I,B_i/J)
  =(I\cap J)/IJ
\]

is free of rank two over `B_i/K`.  The two-jet conormal calculation at (1)
gives

\[
\dim(I/\mathfrak mI,\ J/\mathfrak mJ,\ K/\mathfrak mK)=(3,3,4), 
                                                               \tag{31}
\]

and the two projections of `F\otimes k` to the minimal-generator spaces of
`I` and `J` both have rank two.  Their common core has one linear and one
quadratic generator.

Choose `h_1,h_2` in `I cap J` lifting a basis of `F` and complete their images
to minimal generating sets with `a` in `I` and `b` in `J`.  The projection
ranks and Nakayama give

\[
I=(h_1,h_2,a),\qquad J=(h_1,h_2,b).                \tag{32}
\]

Since `K=(h_1,h_2,a,b)` has height four, these four elements are a regular
sequence.  Modulo `(h_1,h_2)`, the usual intersection identity for two
members of a regular sequence gives

\[
I\cap J=(h_1,h_2,ab).                              \tag{33}
\]

The labeled pair and triple algebras (4) and (10) identify the common linear
generator with `w_i`, the common quadratic generator with `e_i^2`, and the
two remaining quadratic generators with `e_j^2,e_k^2`.  Thus (32) is exactly
(16), and (33) proves the quartic conductor formula (17) on a nonempty open.
Together with Section 7.2, this gives both the generic three-sheet equalizer
and its explicit conductor on every normalization sheet.

Run the exploratory calculation with

```bash
python scripts/explore_degree18_conductor_equalizer.py --order 1
python scripts/explore_degree18_conductor_equalizer.py --order 2
python scripts/explore_degree18_conductor_equalizer.py --order 3
python scripts/explore_degree18_conductor_equalizer.py --transverse --order 1
python scripts/explore_degree18_conductor_equalizer.py --transverse --order 2
python scripts/explore_degree18_conductor_equalizer.py --transverse --order 3
python scripts/verify_degree18_conductor_equalizer.py
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
