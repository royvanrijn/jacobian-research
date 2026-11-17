# Ferrand norms and the first higher transfer blocks

Work over a field `k` of characteristic zero.  This note tests the proposed
identification of the degree-`h` square/cube transfer block with the Ferrand
norm of a dual-number thickening of the universal degree-`h` collision
algebra.

The answer is sharp:

\[
\boxed{\mathcal N_2\simeq Z_2,}
\qquad
\boxed{\mathcal N_3\not\simeq Z_3.}                 \tag{1}
\]

The positive quadratic result is an isomorphism over the full monic
quadratic base, not just at the double-root fiber.  The cubic failure already
occurs at the triple root, even though both algebras are finite flat of rank
eight and have the same Hilbert vector and socle dimension there.

There is nevertheless a more precise positive replacement.  If
\(\mathcal I=(X,Y,T)\) is the nilpotent ideal of the universal cubic
transfer block, then

\[
\boxed{
\mathcal N_3\simeq
\operatorname{gr}_{\mathcal I}(Z_3).}               \tag{2}
\]

The cofactor-sign coordinates for (2) are

\[
 \alpha=X,\qquad\beta=-Y,\qquad\tau=T.
\]

Thus the cubic Ferrand norm is the full relative normal cone of the transfer
block along its reduced monic-cubic base.  The transfer block retains a
nontrivial filtered extension in which one square-zero tangent class acquires
a cubic square.

The next bounded test has the same outcome over the full monic-quartic
base.  If \(\mathcal I_4=(A,B,C,D)\), then

\[
\boxed{
\mathcal N_4
\simeq\operatorname{gr}_{\mathcal I_4}(Z_4).}
\]

Both sides are finite locally free of rank sixteen.  The proof uses saturated
Rees elimination, the canonical cofactor map, and an exact integral
calculation on the only collision stratum with residual affine moduli.  An
all-`h` normal-cone theorem is not proved here.

There is now also an exact bounded theorem through `h=17`.  A Hensel
factorization lemma reduces every geometric collision fiber to tensor
products of maximally collided fibers.  At the maximally collided polynomial
\(S=Z^h\), exact saturated Rees calculations give cone length \(2^h\) for
every \(1\leq h\leq11\); at `h=12`, the raw finite filtered transfer algebra
has length \(4096\), hence so does its associated graded.  Since the
canonical cofactor map is surjective onto the rank-\(2^h\) Ferrand norm, it
is an isomorphism.  Degrees `13` through `17` use good-prime upper-bound
certificate described below.  Thus

\[
\boxed{
\mathcal N_h\simeq\operatorname{gr}_{\mathcal I_h}(Z_h)
\quad(1\leq h\leq17).}
\]

The remaining all-degree gap is no longer a list of collision strata: it is
the single maximally-collided colength identity for arbitrary `h`.

## 1. The construction being tested

Let

\[
 A_h=k[s_1,\ldots,s_h],\qquad
 S(Z)=Z^h+s_1Z^{h-1}+\cdots+s_h,
\]

and let

\[
 B_h=A_h[z]/(S(z)),\qquad
 C_h=B_h[\epsilon]/(\epsilon^2).                    \tag{3}
\]

The algebra `B_h` is finite free of rank `h` over `A_h`, and `C_h` is a
rank-two algebra over `B_h`.  Ferrand's norm of the underlying `B_h`-module
is

\[
 \mathcal N_h
 :=N_{B_h/A_h}(C_h)
 =\Gamma_{A_h}^h(C_h)
   \otimes_{\Gamma_{A_h}^h(B_h)}A_h,                \tag{4}
\]

where the map

\[
 \Gamma_{A_h}^h(B_h)\longrightarrow A_h
\]

is induced by the determinant norm `Nm_(B_h/A_h)`.  Because `C_h` is an
algebra, the norm module (4) has a natural commutative algebra structure and
the universal normic polynomial law is multiplicative.

This is precisely Ferrand's functor, not merely the ordinary norm of selected
elements of `C_h`.  The construction and its arbitrary-base-change property
are reviewed in
[Gille--Neher--Ruether](https://arxiv.org/abs/2401.15051); the original
construction is
[Ferrand, *Un foncteur norme*](https://smf.emath.fr/publications/un-foncteur-norme).

In characteristic zero,

\[
 \Gamma_A^h(M)\simeq (M^{\otimes_A h})^{S_h},       \tag{5}
\]

so (4) is amenable to an exact symmetric-tensor computation.  If `[i_1...i_h]`
denotes the orbit sum of a tensor-basis word, the relations imposed by (4)
are the full polarizations of

\[
 [bm,\ldots,bm]
 =\operatorname{Nm}_{B_h/A_h}(b)[m,\ldots,m].       \tag{6}
\]

Equation (6), followed by componentwise tensor multiplication, is the
algorithm used below.

Over the discriminant complement `B_h/A_h` splits etale-locally.  There (4)
becomes a tensor product of `h` dual-number algebras and has rank `2^h`.
Ferrand base change makes the collision fibers canonical.  What is not
formal is that this canonical degeneration must equal the square/cube
factorization degeneration.

## 2. The quadratic norm

Put

\[
 A=k[p,q],\qquad B=A[z]/(z^2+pz+q),\qquad
 C=B[\epsilon]/(\epsilon^2).                        \tag{7}
\]

Use the `A`-basis

\[
 c_0=1,\quad c_1=z,\quad c_2=\epsilon,\quad
 c_3=z\epsilon
\]

of `C`.  In the symmetric-tensor basis of `Gamma_A^2(C)`, write `[ij]` for
`c_i tensor c_j+c_j tensor c_i` when `i<j`, and for `c_i tensor c_i` when
`i=j`.

Tensoring over `Gamma_A^2(B)` with the determinant norm gives the six monic
relations

\[
\begin{aligned}
[01]&=-p[00],& [11]&=q[00],\\
[12]&=-[03]-p[02],& [13]&=q[02],\\
[23]&=-p[22],& [33]&=q[22].                       \tag{8}
\end{aligned}
\]

Consequently \(\mathcal N_2\) is free over `A` on

\[
 [00],\quad [02],\quad [03],\quad [22].            \tag{9}
\]

Set

\[
 \alpha=[02],\qquad \beta=[03],\qquad \delta=[22].
\]

Componentwise tensor multiplication and (8) give

\[
\begin{aligned}
\alpha^2&=2\delta,&
\alpha\beta&=-p\delta,&
\beta^2&=2q\delta,\\
\alpha\delta&=0,&
\beta\delta&=0,&
\delta^2&=0.                                      \tag{10}
\end{aligned}
\]

Since `2` is invertible, put

\[
 X=\alpha,\qquad Y=-\beta.
\]

Then (10) is exactly

\[
 X^3=0,\qquad 2XY=pX^2,\qquad Y^2=qX^2.            \tag{11}
\]

Therefore

\[
\boxed{
\mathcal N_2
\simeq
k[p,q,X,Y]/(X^3,\,2XY-pX^2,\,Y^2-qX^2)
=Z_2.}                                             \tag{12}
\]

At `p=q=0`, this specializes to

\[
 k[X,Y]/(X^3,XY,Y^2)
 \simeq k[X]/(X^3)\times_k k[Y]/(Y^2).             \tag{13}
\]

The fiber-product description explains its two-dimensional socle:
`X^2` is the socle class on the length-three curvilinear branch and `Y` is
the socle class on the dual-number branch.

This result does not revive the archived norm-polynomial map from the moving
divided-power conductor ribbon.  That map selected the coefficients of two
particular norm polynomials and killed `X^2` at the double root.  Formula
(12) instead computes the entire Ferrand norm algebra (4), whose extra
polarized class retains `X^2`.  The functors and the induced maps are
different.

## 3. The universal cubic Ferrand norm

Now put

\[
 A=k[p,q,r],\qquad
 B=A[z]/(z^3+pz^2+qz+r),\qquad
 C=B[\epsilon]/(\epsilon^2).                        \tag{14}
\]

Order the `A`-basis of `C` as

\[
 1,z,z^2,\epsilon,z\epsilon,z^2\epsilon
\]

and use orbit-sum notation `[ijk]` in `Gamma_A^3(C)`.  Direct reduction of
the normic relations (6) leaves the free basis

\[
\begin{gathered}
[000],\ [003],\ [013],\ [033],\\
[113],\ [133],\ [233],\ [333].                    \tag{15}
\end{gathered}
\]

Define the three algebra generators

\[
 \alpha=[003],\qquad
 \beta=[013],\qquad
 \tau=[113].                                       \tag{16}
\]

The remaining five basis classes are generated by their products.  Exact
tensor multiplication gives the following relative presentation:

\[
\boxed{
\mathcal N_3
\simeq A[\alpha,\beta,\tau]/I_{\mathrm{Fer}},}       \tag{17}
\]

where

\[
\begin{aligned}
I_{\mathrm{Fer}}=\big(&
\alpha^4,\\
&3\alpha^2\beta+2p\alpha^3,\\
&3\alpha^2\tau-q\alpha^3,\\
&\beta^2+2\alpha\tau+2p\alpha\beta
 +(p^2-q)\alpha^2,\\
&2\beta\tau-2q\alpha\beta-(pq-r)\alpha^2,\\
&\tau^2+2r\alpha\beta+pr\alpha^2
\big).                                             \tag{18}
\end{aligned}
\]

With variable order `tau,beta,alpha`, (18) is a relative Groebner basis.
Its standard monomials can be taken as

\[
 1,\ \alpha,\ \alpha^2,\ \alpha^3,\
 \beta,\ \alpha\beta,\ \tau,\ \alpha\tau.          \tag{19}
\]

Thus

\[
\boxed{\operatorname{rank}_A\mathcal N_3=8.}        \tag{20}
\]

This verifies all three formal attractions of the conjecture in the first
new case: the Ferrand algebra has the expected rank, commutes with base
change, and becomes the threefold tensor product of dual numbers over the
split etale locus.

## 4. Independent computation of the transfer block `Z_3`

Let

\[
 S=Z^3+pZ^2+qZ+r,\qquad
 V=S^2+XZ^2+YZ+T.                                  \tag{21}
\]

Take `U` monic of degree nine.  The coefficient equations in degrees
seventeen through nine of

\[
 U^2=V^3                                             \tag{22}
\]

eliminate the nine nonleading coefficients of `U` triangularly.  The
remaining ideal has the relative Groebner basis

\[
\begin{aligned}
T^3={}&6T^2r^2+6X^2pr^3-12XYr^3,\\
T^2X={}&2T^2pr+2X^2p^2r^2-4XYpr^2,\\
TX^2={}&2T^2q+2X^2pqr-4XYqr,\\
X^2Y={}&4T^2p+4X^2p^2r-8XYpr,\\
X^3={}&6T^2+6X^2pr-12XYr,\\
2TY={}&-X^2pq+X^2r+2XYq,\\
Y^2={}&-2TX-X^2p^2+X^2q+2XYp.                    \tag{23}
\end{aligned}
\]

The standard monomials are

\[
 1,\ T,\ T^2,\ X,\ TX,\ X^2,\ Y,\ XY,              \tag{24}
\]

so

\[
\boxed{\operatorname{rank}_A Z_3=8.}               \tag{25}
\]

The same calculation shows that discarding the linear and constant
coefficients of `U^2-V^3` does not enlarge the ideal.  Thus (23) is both the
strong square/cube block and the affine-difference block relevant to
normalized seeds.

## 5. The triple-root obstruction

Base change (18) to `p=q=r=0`.  The cubic Ferrand fiber is

\[
\boxed{
(\mathcal N_3)_0
=k[\alpha,\beta,\tau]/
(\alpha^4,\alpha^2\beta,\alpha^2\tau,
 \beta^2+2\alpha\tau,\beta\tau,\tau^2).}            \tag{26}
\]

It has basis (19), Hilbert vector `(1,3,3,1)`, and a
two-dimensional socle.  Crucially,

\[
 \tau\notin\mathfrak m^2,\qquad \tau^2=0.           \tag{27}
\]

So (26) has a square-zero element with nonzero tangent class.

On the other hand, (23) specializes to

\[
\boxed{
(Z_3)_0=
k[X,Y,T]/
\begin{pmatrix}
T^3,T^2X,TX^2,X^2Y,TY,\\
X^3-6T^2,\ Y^2+2TX
\end{pmatrix}.}                                    \tag{28}
\]

This algebra also has length eight, Hilbert vector `(1,3,3,1)`, and
two-dimensional socle.  Nevertheless it has no square-zero element outside
\(\mathfrak m^2\).

Indeed, let an element have tangent part

\[
 aX+bY+cT.
\]

Modulo `mathfrak m^3`, its square is

\[
 a^2X^2+2abXY+2(ac-b^2)TX.                         \tag{29}
\]

The three displayed quadratic classes are independent.  If the square
vanishes, (29) forces `a=b=0`.  The element is therefore `cT+h` with
\(h\in\mathfrak m^2\).  But

\[
 T\mathfrak m^2=0,\qquad \mathfrak m^4=0,\qquad
 T^2=X^3/6\ne0.                                    \tag{30}
\]

Hence `(cT+h)^2=c^2T^2`, which forces `c=0`.  The existence of a square-zero
tangent class is an isomorphism invariant, so

\[
\boxed{(\mathcal N_3)_0\not\simeq(Z_3)_0.}          \tag{31}
\]

By base change, (31) rules out an isomorphism of the universal algebras:

\[
\boxed{\mathcal N_3\not\simeq Z_3.}                \tag{32}
\]

## 6. What the Ferrand norm does capture

Let \(\mathcal I=(X,Y,T)\) in the universal algebra (23).  For every
relation `f`, form its Rees homogenization

\[
 \widetilde f(h)
 =h^{-\nu_{\mathcal I}(f)}f(hX,hY,hT).
\]

Saturate the ideal of these homogenizations by `h` and then set `h=0`.
The resulting relative normal-cone ideal is

\[
\begin{aligned}
I_{\mathrm{cone}}=\big(&
X^4,\ 2pX^3-3X^2Y,\ 3X^2T-qX^3,\\
&Y^2+2TX-2pXY+(p^2-q)X^2,\\
&2TY-2qXY+(pq-r)X^2,\\
&T^2-2rXY+prX^2
\big).                                             \tag{33}
\end{aligned}
\]

Under the cofactor-sign substitution

\[
\alpha\longmapsto X,\qquad
\beta\longmapsto -Y,\qquad
\tau\longmapsto T,                                \tag{34}
\]

the six relations (18) become exactly (33).  Hence the replacement theorem
holds over the full monic-cubic base:

\[
\boxed{
\mathcal N_3
\simeq\operatorname{gr}_{\mathcal I}(Z_3).}        \tag{35}
\]

The saturation is essential: merely taking the lowest-degree parts of the
seven displayed generators in (23) misses consequences such as `X^4`.
The checker computes both saturated ideals and verifies containment in both
directions.

At `p=q=r=0`, equation (35) becomes

\[
\boxed{
(\mathcal N_3)_0
\simeq\operatorname{gr}_{\mathfrak m}(Z_3)_0.}     \tag{36}
\]

The special cone presentation is

\[
k[X,Y,T]/
(X^4,X^2Y,X^2T,Y^2+2XT,YT,T^2),
\]

which is (26) under (34).  The mismatch between the cone and the transfer
fiber is the filtered extension

\[
 T^2=X^3/6.                                        \tag{37}
\]

The associated graded kills `T^2`; the transfer block remembers it as the
top cubic class.  This also explains why rank, tangent dimension, Hilbert
vector, and socle dimension all fail to distinguish the two fibers.

For `h=2`, the collided transfer presentation is already homogeneous in its
maximal-ideal filtration, so its tangent cone equals the algebra itself.
That is consistent with the exact norm identification (12).  The first
filtered correction appears at `h=3`.

## 7. The quadruple-root test

The next computation tests whether (35) is accidental.  Put

\[
B_0=k[z]/(z^4),\qquad
C_0=B_0[\epsilon]/(\epsilon^2).
\]

In `Gamma_k^4(C_0)`, order the basis of `C_0` as

\[
1,z,z^2,z^3,\epsilon,z\epsilon,z^2\epsilon,z^3\epsilon.
\]

The fully polarized Ferrand quotient has dimension sixteen.  Four tangent
generators are

\[
\alpha=[0004],\qquad
\beta=[0014],\qquad
\gamma=[0114],\qquad
\delta=[1114].                                    \tag{38}
\]

Independently eliminate the twelve nonleading coefficients of `U` from the
quadruple-root equation

\[
U^2=(Z^8+AZ^3+BZ^2+CZ+D)^3.
\]

Rees saturation of the resulting transfer ideal at
\(\mathfrak m=(A,B,C,D)\) gives

\[
\begin{aligned}
I_4^{\mathrm{cone}}=\big(&
BC+AD,\ C^2+2BD,\ CD,\ D^2,\\
&AB^2+A^2C,\ ABD,\ B^3-3A^2D,\ B^2D,\\
&A^3B,\ A^3C,\ A^3D,\ A^5
\big).                                             \tag{39}
\end{aligned}
\]

The standard monomials may be taken as

\[
\begin{gathered}
1,A,B,C,D,A^2,AB,B^2,AC,AD,BD,\\
A^3,A^2B,A^2C,A^2D,A^4.
\end{gathered}                                     \tag{40}
\]

Direct multiplication in the Ferrand quotient shows that (39) vanishes
under

\[
\alpha=A,\qquad
\beta=-B,\qquad
\gamma=C,\qquad
\delta=-D.                                        \tag{41}
\]

The sixteen monomials (40) map to a basis—the determinant in the orbit-sum
basis is `4608`.  Therefore

\[
\boxed{
N_{B_0/k}(C_0)
\simeq\operatorname{gr}_{\mathfrak m}(Z_4)_0.}     \tag{42}
\]

The quartic transfer fiber again contains filtered data absent from the
norm.  In particular, its relations include

\[
A^2B=2D^2,\qquad A^3=12CD,
\]

whereas the normal cone has `D^2=CD=0`.  Thus the same mechanism seen in
(37) persists at the next collision.

## 8. The relative quartic theorem

Now let

\[
S=Z^4+pZ^3+qZ^2+rZ+s
\]

over \(A_4=k[p,q,r,s]\), and let
\(\mathcal I_4=(A,B,C,D)\) be the transverse ideal in the universal
quartic transfer block.  Triangular elimination of the twelve nonleading
coefficients of `U`, followed by saturated Rees homogenization, produces a
finite relative cone.  Over the fraction field of \(A_4\), its standard
monomials are (40), with graded ranks `(1,4,6,4,1)`.  Its base change to
`p=q=r=s=0` is exactly (39).

There is a canonical map from this cone to the Ferrand norm.  Over a
splitting ring, write

\[
S(Z)=\prod_{j=1}^4(Z-r_j).
\]

The four degree-one norm classes are the signed coefficients of

\[
\sum_{j=1}^4\epsilon_j\frac{S(Z)}{Z-r_j}.
\]

They give the map `(alpha,beta,gamma,delta)=(A,-B,C,-D)`.  On the
discriminant complement this is the usual identification of both algebras
with a tensor product of four dual-number algebras.  Since Ferrand's norm is
projective over \(A_4\) for the monogenic algebra \(B_4\), it is
torsion-free.  Consequently all relative cone relations, which vanish after
inverting the discriminant, already vanish in the norm.  Thus the generic
map extends over the whole base.

It remains to prove that the map does not lose generators on the
discriminant.  The geometric root types of a quartic are

\[
1+1+1+1,\quad 2+1+1,\quad 2+2,\quad 3+1,\quad 4.
\]

The split type is formal and the type `4` is (42).  Up to an affine change
of `Z`, the only nonsplit type with a modulus lies in the family

\[
S_\lambda(Z)=Z^2(Z-1)(Z-\lambda)
=Z^4-(1+\lambda)Z^3+\lambda Z^2.
\]

Exact Rees saturation over \(k[\lambda]\) gives a monic Gröbner basis whose
leading monomials are

\[
BC,\ C^2,\ CD,\ D^2,\ AB^2,\ ABD,\ B^3,\ B^2D,\
A^3B,\ A^3C,\ A^3D,\ A^5.
\]

Hence its cone is free over \(k[\lambda]\) on the sixteen monomials (40).
Independently, the fully polarized Ferrand quotient for \(S_\lambda\) is
free on its sixteen orbit classes, and the images of (40) have determinant

\[
-4608.
\]

Therefore the cofactor map is an isomorphism throughout this family.
The specializations `lambda=0` and `lambda=1` cover the types `3+1` and
`2+2`; the open values cover `2+1+1`.  Together with the split and
quadruple-root cases, the map is an isomorphism on every geometric fiber.
Both algebras are finite over \(A_4\), and the norm is projective, so the
fiberwise isomorphism criterion gives

\[
\boxed{
\mathcal N_4
\simeq\operatorname{gr}_{\mathcal I_4}(Z_4).}      \tag{43}
\]

This is more than a collection of sample specializations: the parameter
\(\lambda\) covers the complete moduli-bearing collision stratum, and its
boundary covers the two remaining intermediate root partitions.

## 9. The factorization-algebra reformulation

Ferrand's projectivity theorem supplies more structure than the rank count.
For a finite locally free rank-`h` algebra `B/A`, the norm of the free
rank-two `B`-module decomposes canonically by the two multidegrees:

\[
 N_{B/A}(B^{\oplus2})
 \simeq
 \bigoplus_{i=0}^h P^{(h-i,i)}(B).                 \tag{44}
\]

Here \(P^{(h-i,i)}(B)\) is Ferrand's algebra representing factorizations of
the determinant norm into multiplicative polynomial laws of degrees `h-i`
and `i`.  When \(B=A[z]/(S)\) is monogenic, it equivalently represents
factorizations

\[
 S(Z)=S_{h-i}(Z)S_i(Z)
\]

into monic factors of the indicated degrees.  Ferrand proves that this
summand is finite locally free of rank

\[
 \binom hi.
\]

See §§4.1--4.3 of
[Ferrand's original paper](https://www.numdam.org/item/10.24033/bsmf.2319.pdf),
especially (4.1.2), Lemma 4.3.3, and Theorem 4.3.4.

Since the dual-number algebra is \(B^{\oplus2}\) as a `B`-module, (44)
applies to \(\mathcal N_h\).  It explains the graded ranks seen in every
calculation:

\[
\begin{array}{c|c}
h&\text{ranks by epsilon/cone degree}\\ \hline
2&(1,2,1)\\
3&(1,3,3,1)\\
4&(1,4,6,4,1).
\end{array}
\]

It also isolates a sharper all-`h` target.  It would suffice to construct
isomorphisms

\[
 \operatorname{gr}_{\mathcal I_h}^{\,i}(Z_h)
 \simeq P^{(h-i,i)}(B_h)                           \tag{45}
\]

for every `i`, and prove that these identifications respect multiplication.
Over the split-etale locus, (45) is the elementary identification between
`i`-subsets of the `h` roots on both sides.  In coordinates, the degree-one
map is canonical: if

\[
S(Z)=\prod_{j=1}^h(Z-r_j),
\]

then the cone coordinates are the signed coefficients of the cofactor
polynomial

\[
 \sum_{j=1}^h\epsilon_j\frac{S(Z)}{Z-r_j}.          \tag{46}
\]

Formula (46) gives `(alpha,beta)=(X,-Y)` for `h=2`,
`(alpha,beta,tau)=(X,-Y,T)` for `h=3`, and
`(alpha,beta,gamma,delta)=(A,-B,C,-D)` for `h=4`.
Thus “cofactor-sign coordinates” are intrinsic, not a low-degree fit.

## 10. Reduction to the maximally collided cone

Two structural observations turn the all-`h` comparison into one local
colength problem.

First, the cofactor classes generate the whole Ferrand algebra, including
at ramified fibers.  More explicitly, the degree-one summand in (44) is

\[
P^{(h-1,1)}(B_h)\simeq B_h.
\]

On a splitting algebra, write \(L(a)=\sum_j a(r_j)\epsilon_j\).  A product
\(L(z^{a_1})\cdots L(z^{a_i})\), restricted to an `i`-element subset of the
roots, is the monomial symmetric function with exponent multiset
\(\{a_1,\ldots,a_i\}\), up to the usual nonzero multiplicity.  These
monomial symmetric functions span the coordinate algebra of the universal
degree-`i` factor of `S`.  Ferrand's factorization description and fiberwise
Nakayama therefore show that the multiplication maps

\[
\operatorname{Sym}^i_{A_h}
P^{(h-1,1)}(B_h)\longrightarrow P^{(h-i,i)}(B_h)
\]

are surjective.  Consequently \(\mathcal N_h\) is generated by its cofactor
degree-one classes.

Second, suppose a geometric fiber has root partition
\(h=m_1+\cdots+m_r\).  The coprime factors supported at the distinct roots
lift uniquely under infinitesimal deformation.  Applying this Hensel
factorization simultaneously to `V` and `U` in \(U^2=V^3\) identifies the
completed filtered transfer algebra with the completed tensor product of
the `r` maximally collided transfer algebras.  Passing to the normal cone
gives

\[
(\operatorname{gr}_{\mathcal I_h}Z_h)_S
\simeq
\bigotimes_{\nu=1}^r K_{m_\nu}^{\max},
\]

where \(K_m^{\max}\) is the cone at \(S=Z^m\).  The norm has the matching
tensor decomposition because
\(B_S\simeq\prod_\nu k[z]/(z^{m_\nu})\).

There is a canonical graded cofactor map

\[
\phi_h:\operatorname{gr}_{\mathcal I_h}(Z_h)\longrightarrow\mathcal N_h.
\]

It is an isomorphism on the split-etale locus.  Every cone relation therefore
vanishes in the norm: the norm is projective, hence torsion-free, over the
universal monic-polynomial base.  The preceding generation lemma makes
\(\phi_h\) surjective on every fiber.  It follows that the single equality

\[
\dim_k K_h^{\max}=2^h
\]

implies \(\phi_h\) is an isomorphism at the maximally collided fiber.  If
this equality is known for every \(m\leq h\), the Hensel tensor
decomposition proves that \(\phi_h\) is an isomorphism on every geometric
fiber, hence over the full universal base.

The maximally collided transfer ideal admits a compact exact construction.
Write

\[
V=Z^{2h}+\sum_{j=0}^{h-1}X_jZ^{h-1-j},\qquad
F(w)=1+\sum_{j=0}^{h-1}X_jw^{h+1+j}.
\]

After \(w=Z^{-1}\), the high coefficient equations uniquely force the
reversed polynomial for `U` to be the degree-`3h` truncation of
\(F(w)^{3/2}\).  The coefficients of

\[
\big(\operatorname{trunc}_{\leq3h}F^{3/2}\big)^2-F^3
\]

in degrees \(3h+1,\ldots,6h\) generate the exact residual transfer ideal.
Rees homogenization, saturation, and specialization of these relations give

\[
\operatorname{Hilb}(K_h^{\max})
=\left(\binom h0,\binom h1,\ldots,\binom hh\right),
\qquad
\dim K_h^{\max}=2^h
\quad(1\leq h\leq11).
\]

All operations are over the rationals and the dimensions are exact Gröbner
colengths, not numerical samples.  Combining them with the structural
reduction proves the relative normal-cone theorem for every \(h\leq11\).
For `h=12`, direct Gröbner reduction of the unhomogenized filtered transfer
algebra gives length \(4096\).  A finite-dimensional filtered algebra and
its associated graded have the same total dimension, so the same
equal-length argument proves the theorem at `h=12` without constructing the
slower Rees basis.

Reducing the primitive integral residual equations modulo `3` and using
Singular's `slimgb` engine gives finite fibers of exact lengths

\[
2^h\qquad(11\leq h\leq16).
\]

Fiber length is upper semicontinuous in this integral family, so these
special fibers give the same upper bounds for the characteristic-zero
generic fiber.  The cofactor surjection onto the rank-\(2^h\) Ferrand norm
gives the opposite lower bounds.  Hence equality follows in characteristic
zero without assuming that a modular Gröbner basis lifts.  Independent
certificates at prime `32003` for `h=13,14` and prime `5` for `h=15` agree
with this uniform characteristic-`3` run.  The explicit divided-Frobenius
presentation and its intrinsic tail-weight order give the next exact length
\(131072=2^{17}\).  More importantly, no proof of the colength identity for
arbitrary `h` is known here.

The computations also identify a plausible uniform presentation.  Put

\[
X(w)=x_0+x_1w+\cdots+x_{h-1}w^{h-1}
\]

and let \(J_h^{\mathrm{cur}}\) be generated by

\[
[w^n]X(w)^k,\qquad
2\leq k\leq h+1,\quad n\geq h+1-k.                \tag{47}
\]

For `h=3`, (47) gives the six relations in (36); for `h=4`, it reduces to
(39).  Exact bidirectional ideal reduction now verifies

\[
I(K_h^{\max})=J_h^{\mathrm{cur}}
\qquad(2\leq h\leq11).
\]

There is a useful intrinsic rewrite.  If

\[
R(Z)=Z^{h-1}X(Z^{-1}),\qquad S(Z)=Z^h,
\]

then (47) says exactly

\[
S^{k-1}\mid R^k\qquad(2\leq k\leq h+1).           \tag{48}
\]

Indeed, \([w^n]X(w)^k\) contributes the power
\(Z^{k(h-1)-n}\), which lies below degree \(h(k-1)\) precisely when
\(n\geq h+1-k\).  Only the `h` border coefficients

\[
h+1-k\leq n\leq2h-k
\]

are needed for each `k`: higher coefficients of \(X^k\) follow inductively
from those of \(X^{k-1}\) using \(X^k=X\cdot X^{k-1}\).  This border-window
observation makes the ideal comparison through `h=11` nearly as cheap as
the colength calculation.

The identification of (47) with the fusion algebra is exact, not merely a
pattern match.  Specialize Definition 2.2 of
[Feigin--Feigin, *Q-characters of the tensor products in the
\(\mathfrak{sl}_2\) case*](https://arxiv.org/abs/math/0201111)
to \(a_1=\cdots=a_h=2\).  Their relation
\[
e^{(h)}(w)^k\ \text{is divisible by}\ w^{h(k-1)}
\]
and then their index-reversing involution give exactly
\([w^n]X(w)^k=0\) for \(n\geq h+1-k\).  Their Theorem 2.1 and Proposition
2.1 identify this quotient with the collision limit of `h` dual-number
algebras and prove its dimension \(2^h\).  Its degree-`k` piece has the
monomial basis

\[
x_{i_1}\cdots x_{i_k},
\qquad
0\leq i_1\leq\cdots\leq i_k\leq h-k,
\]

and hence dimension \(\binom hk\).  Chari--Venkatesh give a later
current-algebra treatment and compatible monomial bases
([arXiv:1305.2523](https://arxiv.org/abs/1305.2523)).
These results close the target-presentation and colength questions.  They do
not by themselves prove the transfer theorem: the precise remaining
symbolic step is to derive every current-power relation (47) from the
saturated square/cube Rees equations for arbitrary `h`.  Proving that
inclusion would immediately give the upper bound \(2^h\), while the
cofactor surjection already gives the opposite bound.

There is also a sharper form of the remaining saturation problem.  Put
\[
\Phi(t,w)=\bigl(1+t w^{h+1}X(w)\bigr)^{3/2},
\qquad
G(t,w)=\operatorname{trunc}_{\leq3h}\Phi(t,w).
\]
Because \(G+\Phi\) has constant term `2`, it is a unit in the formal
`w`-series ring.  Thus the saturated residual equations
\(G^2-(1+t w^{h+1}X)^3\) are equivalently the tail equations
\[
[w^m]\Phi(t,w)=0,\qquad m>3h.                    \tag{48a}
\]
For the border generator \(a_{k,n}=[w^n]X^k\), take
\(m=k(h+1)+n\).  Its tail equation contains
\[
\binom{3/2}{k}t^k a_{k,n}.
\]
Every lower-power term has index
\[
n_j=(k-j)(h+1)+n>2h-j,
\]
so it lies beyond the `j`th border window and is already propagated from
lower-power border relations; every higher-power term carries an additional
power of `t`.  Exact saturation checks at `h=4,5` confirm that these tail
equations specialize to (47).  A uniform triangular tail-saturation lemma
formalizing this cancellation is therefore sufficient, and is now the only
unproved algebraic step.

### A divided-Frobenius good-prime reduction

There is a second, potentially shorter route that avoids proving the
characteristic-zero tail lemma directly.  It comes from the unexpectedly
simple primitive fiber in characteristic `3`.

Write
\[
X(w)^2=q(w)+p(w),\qquad
\deg q\leq h-2,
\]
and set \(P_r=[w^r]X(w)^2\), with \(P_r=0\) outside
\(0\leq r\leq2h-2\).  The integral divided cube
\[
\Delta_3(X)=
\frac{X(w)^3-\sum_{i=0}^{h-1}x_i^3w^{3i}}{3}
\]
has integer coefficients.  Direct expansion of the square/cube residual
gives
\[
\begin{aligned}
64\bigl(G^2-F^3\bigr)
={}&-48w^{2h+2}p\\
&+8w^{3h+3}X(q-8p)
+9w^{4h+4}q^2,                                  \tag{48b}
\end{aligned}
\]
where \(F=1+w^{h+1}X\) and
\(G=\operatorname{trunc}_{\leq3h}F^{3/2}\).
Since
\[
X(q-8p)=X^3-9Xp,
\]
primitive reduction of each coefficient modulo `3` has the following
uniform presentation:
\[
\boxed{
I_h^{(3)}
=\left(
P_{h-1},P_h,\;
x_i^3\ (0\leq i<h),\;
[w^r]\Delta_3(X)+P_{h+1+r}\ (3\nmid r)
\right),}                                       \tag{48c}
\]
where \(0\leq r\leq3h-3\).  There are exactly `3h` displayed generators.
They are homogeneous for the positive weight
\[
\operatorname{wt}(x_i)=h+1+i.
\]

Equation (48c) is an all-`h` identity, not a guessed fit: it follows
coefficientwise from (48b).  The checker independently constructs both
ideals and performs bidirectional reduction through `h=10`.  Exact
Gröbner calculations give
\[
\dim_{\mathbf F_3}\mathbf F_3[x_0,\ldots,x_{h-1}]/I_h^{(3)}
=2^h
\qquad(2\leq h\leq17).
\]

This produces a sharper closure criterion:

> **Divided-Frobenius PBW lemma.**  The algebra in (48c) has length \(2^h\)
> for every \(h\).

One proof of this lemma would settle the characteristic-zero theorem for
all `h`.  Indeed, the primitive integral transfer family has generic fiber
equal to the rational transfer algebra, so the characteristic-`3` length
is an upper bound for its characteristic-zero length.  The Ferrand
cofactor surjection gives the opposite bound \(2^h\).

There are now four concrete attacks, in decreasing order of apparent
leverage:

1. **Divided-Frobenius PBW.**  Split
   \(X(w)=A(w^3)+wB(w^3)+w^2C(w^3)\), use the cubic polarization formula
   for \(\Delta_3\), and construct a recursive normal-form set of size
   \(2^h\).  The residue-class split makes this a three-block recurrence,
   not an unstructured Gröbner calculation.  Concretely, for
   \(a=A(w^3), b=wB(w^3), c=w^2C(w^3)\), the Witt addition polynomial gives
   \[
   \begin{aligned}
   \Delta_3(a+b+c)={}&
   \Delta_3(a)+\Delta_3(b)+\Delta_3(c)\\
   &+a^2b+a^2c+ab^2+b^2c+ac^2+bc^2+2abc
   \quad\text{in }\mathbf F_3.
   \end{aligned}
   \]
   Meanwhile \(X^2\) splits into the same three residue classes.  Thus a
   finite-state induction only needs to remember the two adjacent boundary
   coefficients \(P_{h-1},P_h\).  Producing that boundary-state transition
   table, and proving that its total count doubles when `h` increases by
   one, is the most concrete next proof task.
2. **Divided-generator Ferrand comparison.**  Adjoin the genuine
   divided-degree Ferrand generators, compare the resulting algebra with
   the integral norm, and eliminate those extra generators back to (48c).
   A naive comparison using only cofactor products cannot work in
   characteristic `3`: multiplication into degree `i` carries the
   multiplicity `i!`, and already for `h=3` the ordinary characteristic-`3`
   fusion monomials span rank `7` rather than `8`.  The divided generators
   are therefore essential, not cosmetic.  If the enlarged presentation
   eliminates to (48c), Ferrand projectivity supplies length \(2^h\).
3. **Tail-syzygy PBW.**  Prove that the border generators (47) form the
   lowest-degree standard basis of the tail equations (48a).  The main
   ambiguity is the wrap term: cancelling a lower-degree tail coefficient
   can create a higher-degree normal term.
4. **Apolar duality.**  Embed the dual of each cone-degree piece into the
   Feigin--Feigin symmetric-polynomial model.  Diagonal divisibility would
   give the upper bounds \(\binom hk\) degree by degree.

The first two attacks exploit the new good-prime formula and would prove
the whole theorem at once.  The latter two remain useful characteristic-zero
fallbacks.

The first step of the missing induction is visible without elimination.
Write

\[
V=S^2+tR,\qquad
U=S^3+\frac32tSR+t^2Q.
\]

The coefficient of \(t^2\) in \(U^2-V^3\) is

\[
2S^3Q-\frac34S^2R^2.
\]

Polynomiality of `Q` therefore produces \(S\mid R^2\), the `k=2` case of
(48).  Rees saturation experimentally propagates this to
\(S^{k-1}\mid R^k\) at each subsequent order.  A uniform proof must
formalize that propagation; neither the Ferrand references nor the
current-algebra references state this square/cube saturation lemma.

## 11. Consequences and corrected direction

The proposed universal identity

\[
 Z_h\stackrel?=\operatorname{Nm}_{B_h/A_h}
 \big(B_h[\epsilon]/(\epsilon^2)\big)               \tag{49}
\]

is true for `h=1,2` and false for `h=3`.  The failure is not caused by a rank
jump or bad base change.  Both sides are canonical finite-flat rank-eight
degenerations in the cubic test, but they are different degenerations.

The corrected conjecture is:

> For the nilpotent ideal \(\mathcal I_h\) of the universal degree-`h`
> transfer block, the Ferrand norm is the relative normal cone
> \[
> N_{B_h/A_h}(B_h[\epsilon]/(\epsilon^2))
> \stackrel?{\simeq}\operatorname{gr}_{\mathcal I_h}(Z_h),
> \]
> in cofactor-sign coordinates.

This is proved relatively over the full base for every `h<=17`.  The cases
`h<=4` have explicit relative presentations, and the cases `5<=h<=11`
follow from the exact maximally-collided colength certificates and the
Hensel reduction above; `h=12` follows from the equal raw filtered
colength, and `h=13,14,15,16,17` from good-prime upper bounds.  The all-`h`
theorem remains open.  Its exact remaining lemma is

\[
\dim_k K_h^{\max}=2^h\qquad\text{for every }h,
\]

equivalently (45), including multiplication.  Any proof must control Rees
saturation at the maximally collided root; split-etale agreement and rank
`2^h` do not by themselves control that closure.

## Executable certificate

Run

```bash
python scripts/verify_ferrand_norm_transfer_blocks.py
```

The checker:

1. constructs the Ferrand quotient (4) from the fully polarized normic
   relations;
2. verifies the universal quadratic presentation and its isomorphism with
   `Z_2`;
3. computes the universal cubic norm presentation (18) and its free
   rank-eight basis;
4. independently eliminates `U` from `U^2=V^3` to recover (23);
5. proves the square-zero tangent obstruction (31); and
6. verifies the saturated relative cubic normal-cone isomorphism (35);
7. independently eliminates the quadruple-root `Z_4` equations;
8. constructs the length-sixteen quartic Ferrand fiber and verifies (42);
   and
9. computes the relative quartic Rees cone, its binomial generic basis, its
   quadruple-root base change, and the monic double-root-family basis.

The longer integral Ferrand calculation on the moduli-bearing quartic
collision family is optional:

```bash
python scripts/verify_ferrand_norm_transfer_blocks.py --quartic-strata
```

It certifies the \(k[\lambda]\)-free Ferrand basis and determinant `-4608`
used in the proof of (43).

The maximally collided certificates are independent and run with:

```bash
python scripts/verify_maximally_collided_transfer_cones.py
```

For each `5<=h<=11`, this constructs the residual transfer ideal from the
truncated \(F^{3/2}\) recurrence, performs saturated Rees elimination in
Singular, and verifies the exact binomial Hilbert vector and colength
\(2^h\).  Passing
`--min-degree 1` also rechecks the four lower maximally collided fibers.

The next bounded case uses the faster equal-total-length route:

```bash
python scripts/verify_maximally_collided_transfer_cones.py \
  --method raw --min-degree 12 --max-degree 12
```

It verifies length \(4096\) for the finite filtered `h=12` transfer algebra,
which is also the total length of its associated graded cone.

The remaining bounded cases use one uniform characteristic-`3`
upper-bound certificate:

```bash
python scripts/verify_maximally_collided_transfer_cones.py \
  --method modular --basis-engine slimgb --prime 3 \
  --min-degree 11 --max-degree 16
```

It returns \(2^h\) in every displayed degree, ending with
\(65536=2^{16}\).  These values, combined with the characteristic-zero
Ferrand lower bounds, prove equality in characteristic zero by upper
semicontinuity.

The weighted divided-Frobenius calculation reaches the next degree:

```bash
python scripts/verify_maximally_collided_transfer_cones.py \
  --method divided-three --basis-engine slimgb \
  --min-degree 17 --max-degree 17 --timeout 240
```

It uses \(\operatorname{wt}(x_i)=h+1+i\) and returns
\(131072=2^{17}\).

The divided-Frobenius presentation (48c) is checked independently with:

```bash
python scripts/verify_maximally_collided_transfer_cones.py \
  --method modular --basis-engine slimgb --prime 3 \
  --min-degree 2 --max-degree 10 \
  --compare-characteristic-three-model
```

For each degree, this reduces the primitive residual ideal and the explicit
divided-Frobenius ideal against one another before checking the length.

The optional command

```bash
python scripts/verify_maximally_collided_transfer_cones.py \
  --min-degree 2 --max-degree 11 --compare-current-power-ideal
```

also verifies equality with (47) by reducing each Gröbner basis against the
other.
