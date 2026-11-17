# Two-pair Image in bidegree (3,3): partial classification and exact reductions

The full sixteen-parameter problem remains open. This note gives an exact
coordinate model for that full family, covers both parity strata, and classifies
mixed slices with explicit finite certificates. These slices and their orthogonal
translates have not been proved to exhaust the mixed family. The least balanced
degree of an Image failure therefore remains either three or four; the
[quartic counterexample](TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md) supplies the
upper bound.

Work over C. For a balanced-degree-three polynomial F in two contraction pairs,
the target is the following implication, for every polynomial multiplier Q:

\[
 \mathcal E_2(F^m)=0\ (m\geq1)
 \quad\Longrightarrow\quad
 \mathcal E_2(QF^m)=0\ \text{for all sufficiently large }m,
 \qquad
 \mathcal E_2(\xi^\alpha z^\beta)=\partial_z^\alpha z^\beta.
 \tag{1}
\]

This is an Image/Mathieu assertion. Membership in the orthogonal nullcone is a
sufficient condition, but the safe nonlinear shears below show that it is not a
necessary condition. An argument proving only a nullcone criterion would miss
admissible pure cubics.

## The full coefficient space and its spherical functional

Set

\[
 R=\xi_1z_1+\xi_2z_2,\quad X=\xi_1z_2,\quad
 Y=2\xi_2z_1,\quad T=\xi_1z_1-\xi_2z_2,
 \qquad T^2+2XY=R^2.
\]

Restriction to R=1 identifies the balanced space V_3 with

\[
 A_{\leq3},\qquad A=\mathbb C[x,y,t]/(q-1),\qquad q=t^2+2xy.
\]

Indeed, eliminate every xy using xy=(1-t^2)/2. The resulting basis is

\[
 1,t,t^2,t^3;\quad
 x,xt,xt^2,x^2,x^2t,x^3;\quad
 y,yt,yt^2,y^2,y^2t,y^3.
 \tag{2}
\]

To see independence, order monomials so that xy is the leading term of q-1;
the monomials in (2) are distinct unreduced monomials. Every member of (2)
has a balanced-degree-three lift, by homogenizing it with R. Conversely the
sixteen original monomials
xi_1^i xi_2^(3-i) z_1^j z_2^(3-j), 0<=i,j<=3, map to this basis by an
invertible rational matrix. In the displayed basis order and lexicographic
(i,j) order its determinant is 2^(-24). Thus this identification discards no
coefficient stratum.

Let L be the normalized sphere functional. Its exact algebraic definition is

\[
 L(x^a y^b t^c)=0\quad\text{unless }a=b\text{ and }c\text{ is even},
 \qquad
 L(x^a y^a t^{2k})=\frac{a!(2k-1)!!}{(2a+2k+1)!!}.
 \tag{3}
\]

Here (-1)!!=1. Equivalently, take the phase constant term after substituting
y=(1-t^2)/(2x), and then integrate t uniformly on [-1,1]. This formula
annihilates (q-1)A and is SO_3(C)-invariant. For every balanced polynomial B of
bidegree (n,n), with sphere restriction b,

\[
 \mathcal E_2(B)=(n+1)!L(b).
 \tag{4}
\]

One can derive (4) without analytic integration: the formal beta identity
\(\int_0^1 s^i(1-s)^j ds=i!j!/(i+j+1)!\), followed by phase extraction,
gives precisely the monomial contraction. Thus (1) for balanced multipliers
is exactly the eventual mixed-moment assertion for L.

Antipodal parity gives a second useful full coordinate model:

\[
 p=H+E\pmod{q-1},\qquad H\text{ homogeneous cubic},\quad
 E\text{ homogeneous quadratic}.
 \tag{5}
\]

This representation is unique: homogenize odd linear terms by q, and even
constants by q; comparing degrees in (q-1)(a+ell) proves injectivity. The
parameter dimensions are ten and six. These are two linear subspaces of the
full space, not two components whose separate classification settles their
mixed sums. In particular the summands in L((H+E)^m) have different total
degrees. The Gaussian-to-sphere factors in the odd-stratum proof cannot be
pulled out of this sum.

### Why balanced multipliers suffice

Suppose Q has bidegree (a,b). If a>b its product with F^m contracts to zero.
Otherwise r=b-a and P_m=E_2(QF^m) is homogeneous of coordinate degree r. For
each multi-index gamma with |gamma|=r,

\[
 \mathcal E_2(\xi^\gamma QF^m)=\partial_z^\gamma P_m.
 \tag{6}
\]

The multiplier xi^gamma Q has bidegree (b,b). Vanishing of all these scalar
moments forces every coefficient of P_m to vanish, in characteristic zero.
A balanced multiplier of degree b restricts to a sphere polynomial of degree
at most b. Hence a sphere bound m>C deg g implies the original Image bound
m>C deg_z Q. Split an arbitrary Q into its finitely many bidegrees to obtain
the same statement for arbitrary, unbalanced multipliers.

## An exact frame criterion, with the completeness obligation left open

For the canonical lift p=H+E in (5), consider a nonzero null vector v. There
exists a cubic representative of p which is affine linear in the coordinate
in direction v if and only if

\[
 q(v)=0,\qquad E(v)=0,\qquad \nabla_q H(v)\in\mathbb C v.
 \tag{7}
\]

Here the gradient uses the bilinear form of q. Representatives of degree at
most three differ by (q-1)(a+ell), with a constant a and a linear form ell.
At a null vector, adding q ell changes the gradient of H by 2 ell(v)v and
adding a q does not change E(v). Consequently (7) permits a representative
whose cubic gradient vanishes at v and whose quadratic part vanishes there.
For a cubic polynomial P with homogeneous pieces P_3,P_2,

\[
 D_v^2P(z)=2B(\nabla_q P_3(v),z)+2P_2(v).
\]

This is zero exactly when P is affine linear in the v coordinate. The
converse proves necessity of (7) for such a representative.

This criterion is effective without choosing an orthogonal frame. Parametrize
the null conic by v(s,r)=(s^2,-r^2/2,sr) and form the binary sextic
f(s,r)=H(v(s,r)) and binary quartic g(s,r)=E(v(s,r)). Then (7) is equivalent
to a common projective zero of f_s,f_r,g. Euler's identity gives f=0 at such
a zero, and the derivative condition says that the zero of f on the conic
is repeated. If f is identically zero, retain both zero derivatives and
apply the same criterion; no division by f or by a discriminant is allowed.

**Unproved obligation:** all-order purity has not been shown to imply (7).
Even if this frame-extraction step is proved, the resulting affine-linear
sphere family still requires classification with its actual spherical
weights. Neither step follows from the Gaussian affine-frame theorem.

## Safe mixed cubics outside the orthogonal nullcone

Let delta=x partial_t-t partial_y. It is an infinitesimal orthogonal action,
delta(q)=0, and delta(x)=0. Sphere invariance gives L(delta h)=0; therefore
L(x delta h)=L(delta(xh))=0. The derivation x delta is locally nilpotent,
so its exponential preserves L by a finite polynomial expansion on every
input. In particular the automorphism phi=exp(-x delta) is

\[
 \phi(x)=x,\qquad\phi(t)=t-x^2,\qquad
 \phi(y)=y+xt-\tfrac12x^3=:p_*.
 \tag{8}
\]

Its inverse replaces x delta by -x delta in the exponential. Both maps
preserve q. Thus L(p_*^m)=L(y^m)=0 for every m>=1, and

\[
 L(gp_*^m)=L(\phi^{-1}(g)y^m)=0\qquad(m>3\deg g).
 \tag{9}
\]

The reason is that deg(phi^(-1)g)<=3 deg g, and its positive phase weight
cannot exceed its degree. Equation (6) gives the Image cutoff
m>3 deg_z Q for the balanced cubic

\[
 F_*=R^2Y+RXT-\tfrac12 X^3.
\]

This coefficient three is necessary even on this safe family. For D>=1,
only the highest positive-weight monomial in phi^(-1)(y^D) can balance
y^(3D), so

\[
 L(y^D p_*^{3D})=
 (-\tfrac12)^D L((xy)^{3D})
 =\frac{(-1)^D 2^{2D}((3D)!)^2}{(6D+1)!}\ne0.
 \tag{10}
\]

Consequently E_2(Y^D F_*^(3D)) is (10) times (10D+1)!, and is nonzero.
There is also an explicit invariant proving that F_* is outside the
orthogonal nullcone (the set whose orbit closure contains zero). Its odd
homogeneous component is H_*=qy-x^3/2, and

\[
 L(H_*^4)=L((y-x^3/2)^4)=-2L(x^3y^3)=-4/35\ne0.
 \tag{11}
\]

The map (H,E) to L(H^4) is an invariant polynomial of positive degree, so it
vanishes on that nullcone. Equation (11) also shows that purity of a mixed
sum does not force purity of its two parity summands. This does not make
F_* an Image counterexample; (9) proves its safety.

## The odd parity stratum

Work over \(\mathbb C\), put \(q=t^2+2xy\), and let
\(A=\mathbb C[x,y,t]/(q-1)\). Let \(L\) be the normalized sphere
functional, with
\[
 L(x^a y^b t^c)=0\quad\hbox{unless }a=b\hbox{ and }c\hbox{ is even},
 \qquad
 L(x^a y^a t^{2c})
 =\frac{a!(2c-1)!!}{(2a+2c+1)!!}.
\]
We use \((-1)!!=1\). These formulas follow by first extracting phase
weight zero and then integrating
\(2^{-a-1}(1-t^2)^a t^{2c}\) on \([-1,1]\).

The odd part of \(A_{\leq3}\) has dimension ten. Every member has a
unique representation as the restriction of a homogeneous ternary
cubic:
\[
 p=C_3+L_1\pmod{q-1}
 \quad\longleftrightarrow\quad H=C_3+qL_1 .
\]
This map is well-defined despite the ambiguity in \(C_3,L_1\). If a
homogeneous cubic vanishes in \(A\), it is divisible by \(q-1\); the
degree bound makes its quotient linear, and comparing homogeneous
degrees one and three makes that quotient zero. Thus the map is
injective, and replacing each linear monomial by \(q\) times itself
gives surjectivity.

Let \(G\) be the Gaussian Wick functional with
\(G(x^a y^b t^c)=\delta_{ab}a!(c-1)!!\) for even \(c\), and zero for
odd \(c\). Monomial inspection proves, for every homogeneous polynomial
\(K\) of even degree \(n\),
\[
 G(K)=(n+1)!!L(K).
 \tag{O1}
\]
For odd \(n\), both functionals vanish by antipodal symmetry. Therefore
\[
 G(H^{2r})=(6r+1)!!\,L(p^{2r}),\qquad r\geq1,
 \tag{O2}
\]
and both odd-power sequences vanish automatically. This transfer uses
the homogeneous lift \(H\); applying it directly to an arbitrary mixed
parity phase would not be valid.

The internal Gaussian classification is the theorem "Normal forms for
pure ternary cubics", label \(\texttt{thm:normalforms}\), in
[the Gaussian-boundary manuscript](../../../unpublished/gaussian-boundary/manuscript.tex).
It states that a complex polynomial of degree at most three with all
positive Gaussian moments zero is, after a complex orthogonal linear
change of coordinates and a nonzero scalar multiple, one-sided in
circular weight or the exceptional phase
\[
 W+ZT-\tfrac12Z^3.
\]
Its proof uses the frame-extraction proposition
\(\texttt{prop:frameextract}\) and the affine-linear classification
\(\texttt{prop:wlinear}\), proved through the manuscript's included
sections. This is a dependency on the existing internal written
classification; equations (O1)--(O2) do not independently verify that
classification.

For homogeneous \(H\), the exceptional alternative is impossible:
orthogonal linear substitutions and nonzero scalar multiplication
preserve homogeneity, whereas the displayed exceptional polynomial has
nonzero terms in three distinct degrees. Consequently (O2) and that
classification give the following entire-stratum result.

**Odd-stratum proposition, dependent on the internal Gaussian
classification.** Every odd \(p\in A_{\leq3}\) with \(L(p^m)=0\) for all
\(m\geq1\) becomes, after a complex orthogonal change of coordinates,
\[
 H=x(a x^2+bxt+ct^2+dxy)
\]
for some \(a,b,c,d\in\mathbb C\). Conversely every such form has all
positive sphere moments zero. The nonzero monomials have strict
positive circular weights \(3,2,1,1\). The zero polynomial is included.
If the orthogonal substitution has determinant \(-1\), compose with
\(-I_3\); this changes a homogeneous cubic by a scalar and gives
determinant \(+1\).

This conclusion is stronger than safety on the odd stratum. It must
not be extended to mixed parity phases. For example, the sphere
automorphism
\(\exp(a x\delta)\), \(\delta=x\partial_t-t\partial_y\), sends \(y\) to
\[
 y-axt-\tfrac12a^2x^3.
\]
This is a mixed parity pure family and does not contradict the
homogeneous argument.

### Transfer back to arbitrary Image multipliers

Write the four original bilinears as
\[
 R=\xi_1z_1+\xi_2z_2,\quad
 X=\xi_1z_2,\quad Y=2\xi_2z_1,\quad
 T=\xi_1z_1-\xi_2z_2,
 \qquad T^2+2XY=R^2.
\]
For odd sphere restriction, the unique balanced-degree-three lift is
the homogeneous cubic \(F=H(X,Y,T)\). The normalization action of
\(\mathrm{SO}_3(\mathbb C)\) lifts through
\(\mathrm{SL}_2(\mathbb C)\to\mathrm{SO}_3(\mathbb C)\) to a linear
substitution of the original contraction pairs, preserving \(R\).
Concretely this is conjugation of the traceless matrix
\[
 \begin{pmatrix}T/2&X\\Y/2&-T/2\end{pmatrix},
 \qquad \det=-\tfrac14(T^2+2XY).
\]
The contraction \(\mathcal E_2\) is equivariant under this simultaneous
linear substitution of \(z\) and contragredient substitution of \(\xi\);
this follows directly from the chain rule
\(\mathcal E_2(\xi^\alpha z^\beta)=\partial_z^\alpha z^\beta\).
Bidegrees of arbitrary multipliers are preserved by these linear maps.

Choose torus weights
\[
 \mathrm{wt}(\xi_1)=\mathrm{wt}(z_2)=1,\qquad
 \mathrm{wt}(\xi_2)=\mathrm{wt}(z_1)=-1.
\]
Then \(X,Y,T,R\) have weights \(2,-2,0,0\), and each monomial of the
normalized \(F\) has weight at least two. Let \(Q\) have bidegree
\((a,b)\). If \(a>b\), its product with \(F^m\) contracts to zero
because every term has too many derivatives. If \(b\geq a\), its
contraction has coordinate degree \(b-a\) and therefore weights at most
\(b-a\). On the other hand, every input weight in \(QF^m\) is at least
\(2m-a-b\). Equivariance under the torus now forces
\[
 \mathcal E_2(QF^m)=0
 \quad\text{if}\quad
 2m-a-b>b-a,
 \quad\text{equivalently }m>b.
 \tag{O3}
\]
For a general polynomial \(Q\), split it into its finitely many
bihomogeneous pieces. Equation (O3) gives the uniform sufficient cutoff
\[
 \mathcal E_2(QF^m)=0\qquad(m>\deg_z Q).
\]
The same degree bound applies before the normalization. Thus the odd
stratum has the Image conclusion for completely arbitrary polynomial
multipliers, not merely balanced ones.

## The even parity stratum: a proof using three moments

Every even member of \(A_{\leq3}\), where
\(A=\mathbb C[x,y,t]/(t^2+2xy-1)\), is the restriction of a unique
homogeneous quadratic \(E\). A constant \(c\) is homogenized as \(cq\).
The injectivity follows by comparing degrees in a proposed identity
\(E=(q-1)c\). Hence this stratum has dimension six.

Let \(B\) be the bilinear form with \(q(v)=B(v,v)\), and write
\(E(v)=B(v,Av)\) for the unique \(B\)-selfadjoint operator \(A\). The
first three sphere moments are
\[
 \begin{split}
 L(E)&=\tfrac13\operatorname{tr}A,\\
 L(E^2)&=\tfrac1{15}
 \left((\operatorname{tr}A)^2+2\operatorname{tr}A^2\right),\\
 L(E^3)&=\tfrac1{105}
 \left((\operatorname{tr}A)^3+
 6\operatorname{tr}A\,\operatorname{tr}A^2+
 8\operatorname{tr}A^3\right).
 \end{split}
 \tag{E1}
\]
These identities hold for every selfadjoint \(A\), including defective
ones. One derivation first proves the Gaussian identities by Wick
pairings, then divides by \(3!!,5!!,7!!\), respectively. Alternatively
they are the first three coefficients of the formal identity
\[
 G(e^{sE})=\det(I-2sA)^{-1/2};
\]
this identity follows by diagonalizing a dense open set of symmetric
matrices, multiplying the one-variable formal Gaussian identity, and
then extending each polynomial coefficient to all symmetric matrices.
Thus no spectral diagonalization assumption is imposed on \(A\).

If the first three sphere moments vanish, equation (E1) successively
gives
\[
 \operatorname{tr}A=\operatorname{tr}A^2=
 \operatorname{tr}A^3=0.
\]
Newton's identities in dimension three now give the characteristic
polynomial \(\lambda^3\). Cayley--Hamilton proves \(A^3=0\).

Assume \(A\neq0\). If \(A^2\neq0\), choose \(v\) with \(u=A^2v\neq0\).
Then
\[
 Au=0,\qquad B(u,u)=B(v,A^4v)=0.
\]
If \(A^2=0\), choose \(v\) with \(u=Av\neq0\). Then again
\[
 Au=0,\qquad B(u,u)=B(v,A^2v)=0.
\]
In both cases \(\ker A\) contains a nonzero null vector \(u\). Complete
\(u\) to a Witt basis and identify it with the \(y\)-axis. Because
\(A\) is selfadjoint,
\[
 D_uE(v)=2B(u,Av)=2B(Au,v)=0.
\]
Thus \(E\) is independent of \(y\), and has the form
\[
 E=ax^2+bxt+ct^2.
\]
In these Witt coordinates \(\operatorname{tr}A=c\). Therefore \(c=0\).
The resulting form
\[
 E=ax^2+bxt
 \tag{E2}
\]
has only strict positive circular weights two and one. Its all-order
pure sphere moments vanish by phase extraction, and every fixed
sphere multiplier has eventual mixed vanishing. The case \(A=0\)
gives \(E=0\) and the same conclusion.

Consequently, for the full even parity stratum, the following
conditions are equivalent: the first three pure sphere moments
vanish; all positive pure sphere moments vanish; and the phase is
orthogonally equivalent to (E2). No Gaussian ternary-cubic
classification is used in this argument.

In original balanced bidegree \((3,3)\), the even phase lifts as
\(F=R(aX^2+bXT)\). The factor \(R\) has torus weight zero, so the
same weight argument as (O3) proves
\[
 \mathcal E_2(QF^m)=0\qquad(m>\deg_z Q)
\]
for every polynomial multiplier \(Q\).

## A complete seven-parameter phase slice

Let \(L\) denote the normalized spherical moment functional, and consider
\[
 p=y(c_0+c_1x+c_2x^2)+d_0+d_1x+d_2x^2+d_3x^3.
 \tag{S1}
\]
This slice has dimension seven. Its complete classification is a proper
subfamily result and supplies no stratification of all sixteen parameters.

**Proposition.** For (S1), the following conditions are equivalent:

1. \(L(p^m)=0\) for \(1\le m\le6\).
2. \(p\) has strictly one-sided phase support in the displayed coordinates,
   including the zero polynomial.
3. \(L(p^m)=0\) for every \(m\ge1\).

More precisely, in the \(c_0\ne0\) branch, \(p=c_0y\). In the
\(c_0=0\) branch, \(c_1=d_0=0\), leaving only positive weights. The fixed
bound six is sharp for this slice.

**Proof.** The monomial rule
\[
 L(x^iy^j)=\delta_{i,j}\frac{i!}{(2i+1)!!}
 \tag{S2}
\]
follows by phase averaging and \(xy=(1-t^2)/2\):
\[
 L((xy)^i)=2^{-i}\int_0^1(1-t^2)^i\,dt
          =\frac{i!}{(2i+1)!!}.
\]
Thus for \(p=yC(x)+D(x)\), the exact identity is
\[
 L(p^m)=\sum_{k=0}^m\binom mk\frac{k!}{(2k+1)!!}
             [x^k]C(x)^kD(x)^{m-k}.
 \tag{S3}
\]

First suppose \(c_0\ne0\). Divide \(p\) by \(c_0\), write
\(C=1+ax+bx^2\), and use \(d_i\) for the rescaled coefficients of
\(D\). The first moment gives \(d_0=-a/3\). Successively setting
moments two, three, and four to zero gives
\[
\begin{aligned}
 d_1&=-\frac{a^2}{30}-\frac{2b}{5},\\
 d_2&=\frac{a(31a^2-108b)}{1890},\\
 d_3&=-\frac{289a^4}{28350}
       +\frac{47a^2b}{1050}-\frac{2b^2}{75}.
\end{aligned}
 \tag{S4}
\]
The coefficient of each new variable \(d_1,d_2,d_3\) in its equation is,
respectively, \(2/3,2/5,8/35\), so no division by a parameter is hidden
in (S4). After these substitutions, the next two moments are
\[
\begin{aligned}
 L(p^5)&=-\frac{4a}{9823275}
 (2219a^4-11874a^2b+13176b^2),\\
 L(p^6)&=-\frac8{13408770375}
 (482933a^6-1979487a^4b+218781a^2b^2+1571724b^3).
\end{aligned}
 \tag{S5}
\]
If \(a\ne0\), put \(z=b/a^2\). Vanishing in (S5) forces
\[
 P(z)=2219-11874z+13176z^2=0,
 \quad Q(z)=482933-1979487z+218781z^2+1571724z^3=0.
\]
These two equations are incompatible: the exact certificate is
\[
 A(z)P(z)+B(z)Q(z)=76274195239375,
 \tag{S6}
\]
where
\[
\begin{aligned}
 A(z)&=6676468225812z^2+2727590293062z-7662133538547,\\
 B(z)&=35364219296-55969842888z.
\end{aligned}
\]
It follows that \(a=0\). Then (S5) becomes
\(L(p^6)=-32b^3/34125\), so \(b=0\). Equations (S4) and the first
moment give \(d_0=d_1=d_2=d_3=0\). Undoing the scaling gives \(p=c_0y\).

If \(c_0=0\), the only phase-zero terms are \(c_1xy+d_0\); all other
terms have positive phase weight. They cannot contribute to any pure
moment. The first moment gives \(d_0=-c_1/3\); after this substitution,
the second moment is \(c_1^2/45\). Therefore \(c_1=d_0=0\), leaving
strictly positive phase weights. In either branch strict one-sidedness
proves all pure moments vanish, and it proves eventual mixed vanishing
for each fixed polynomial multiplier by the same weight argument.

Finally, take
\[
 p_*=y+x^2y-\frac25x-\frac2{75}x^3.
 \tag{S7}
\]
This is (S4) with \(a=0,b=1\), so
\[
 L(p_*^m)=0\quad(1\le m\le5),
 \qquad L(p_*^6)=-\frac{32}{34125}\ne0.
\]
Thus five moments do not suffice, even inside this slice. The proof above
shows that six do. This is an exact cutoff for (S1), not a cutoff for the
full cubic family. \(\square\)

The symbolic certificate is replayed by
`research/foundations/scripts/verify_two_pair_image_cubic_phase_slice.py`
with the existing `sympy==1.14.0` environment. It checks the literal
triangular identities, both exceptional branches, the nonzero Bezout
constant, rejection of an altered Bezout coefficient, and (S7). It uses
no floating-point arithmetic. The retained JSON records these equations
and their exact scope; the proof's support argument supplies the
all-order implication after the finite equations have classified (S1).


## A complete ten-parameter affine slice

Consider

\[
 p=ay+(v_0+v_1x)t^2+(c+dx+ex^2)t+\sum_{j=0}^3r_jx^j.
 \tag{A1}
\]

The ten coefficients vary independently, including \(a=0\). This slice
omits \(yt,yt^2,t^3\) from the thirteen-parameter affine family.

**Proposition.** For (A1), vanishing of the first eight pure sphere moments
is equivalent to either of the following disjoint conditions:

1. \(a=0\) and \(v_0=c=r_0=0\).
2. \(a\ne0\), \(v_0=v_1=e=0\), and
   \[
   \sum_{j=0}^3r_jx^j=-\frac{x(c+dx)^2}{2a}.
   \tag{A2}
   \]

In both cases all pure moments vanish, and every fixed multiplier has eventual
mixed vanishing. The sphere cutoff is \(m>3\deg g\); the Image cutoff is
\(m>3\deg_z Q\). No optimality of the finite prefix eight is asserted.

**Proof.** First take \(a=0\). Phase averaging removes every term of positive
\(x\)-weight, leaving the univariate polynomial \(v_0t^2+ct+r_0\).
The first moment gives \(r_0=-v_0/3\). The second and then the third give

\[
 L(p^2)=\frac{4v_0^2}{45}+\frac{c^2}{3},\qquad
 L(p^3)\big|_{c^2=-4v_0^2/15}=-\frac{256v_0^3}{4725}.
 \tag{A3}
\]

Thus \(v_0=c=r_0=0\), leaving strict positive weights. This also proves the
converse and mixed vanishing in this branch.

If \(a\ne0\), divide by \(a\) and rename coefficients, so \(a=1\). The
orthogonal null rotation

\[
 x=x',\qquad t=t'+kx',\qquad y=y'-kt'-k^2x'/2
 \tag{A4}
\]

preserves the slice, changes the constant \(t'\)-coefficient from \(c\) to
\(c-k\), and preserves \(v_0,v_1\). Take \(k=c\). Hence it suffices to
classify the branch \(a=1,c=0\).

Formula (3) computes each pure moment exactly. The first four moments solve
successively for \(r_0,r_1,r_2,r_3\), with respective nonzero constant
coefficients \(1,2/3,2/5,8/35\). The solutions are

\[
\begin{aligned}
 r_0&=-v_0/3,\\
 r_1&=-(2v_0^2+3v_1)/15,\\
 r_2&=-4v_0(31v_0^2+54v_1)/945,\\
 r_3&=-(14175d^2+4624v_0^4+10152v_0^2v_1+3024v_1^2)/28350.
\end{aligned}
 \tag{A5}
\]

After substitution, write \(L(p^j)=\kappa_jN_j\), \(5\leq j\leq8\), where

\[
 (\kappa_5,\kappa_6,\kappa_7,\kappa_8)=
 \left(\frac8{9823275},-\frac8{13408770375},
 \frac{32}{273648375},-\frac{256}{1465387048125}\right).
 \tag{A6}
\]

The checker and its retained certificate give the four literal polynomials
\(N_j\). If \(v_0=0\), their specializations are

\[
 N_8=9144576v_1^4,\qquad N_6\big|_{v_1=0}=-58046625e^2.
 \tag{A7}
\]

Thus \(v_1=e=0\), and (A5) leaves only \(r_3=-d^2/2\), as required.

If \(v_0\ne0\), put

\[
 b=v_1/v_0^2,\qquad D=d^2/v_0^4,\qquad E=de/v_0^5,\qquad F=e^2/v_0^6.
\]

The four equations \(N_j/v_0^j=0\) are linear in \(D,E,F\). Their augmented
matrix, with columns \(D,E,F,1\), is

\[
 \begin{pmatrix}
 155925&155925&0&B_5\\
 77395500-116093250b&77395500&-58046625&B_6\\
 226800-297675b&226800+127575b&-212625&B_7\\
 4819500-216877500b&4819500-65063250b&-75907125&B_8
 \end{pmatrix},
 \tag{A8}
\]

where

\[
\begin{aligned}
 B_5&=35504+94992b+52704b^2,\\
 B_6&=30907712+63343584b+3500496b^2-12573792b^3,\\
 B_7&=53504+17328b-211968b^2-108864b^3,\\
 B_8&=7986496-27544704b-90351504b^2-16632864b^3+9144576b^4.
\end{aligned}
\]

The coefficient columns satisfy \(C_D-C_E-2bC_F=0\), so the coefficient
matrix has rank at most two. Consistency forces all three-by-three augmented
minors to vanish. Columns \(D,F,1\) and row triples \((5,6,7)\), \((5,6,8)\)
give, after removing nonzero rational factors,

\[
\begin{aligned}
 f(b)&=2485161b^3+6647562b^2+4908168b+1097984=0,\\
 g(b)&=10216206b^4+24356619b^3-25372278b^2-41158800b-11711872=0.
\end{aligned}
 \tag{A9}
\]

But

\[
 \operatorname{Res}_b(f,g)
 =-3179537325694882741702679455938048000000000\ne0.
 \tag{A10}
\]

The receipt also retains rational polynomials \(U,V\) with the checked identity
\(Uf+Vg=1\). This excludes \(v_0\ne0\), including all possible degeneracies
of the linear equations. Undoing (A4) and the nonzero scalar normalization
recovers (A2).

For sufficiency in the nonzero branch, let \(h(x)=(c+dx)/a\) and use

\[
 \phi(x)=x,\quad \phi(t)=t-xh(x),\quad
 \phi(y)=y+h(x)t-xh(x)^2/2.
\]

This is \(\exp(-h(x)\delta)\), with \(\delta=x\partial_t-t\partial_y\).
Since \(\delta h=0\), infinitesimal sphere invariance gives
\(L(h\delta g)=L(\delta(hg))=0\). Local nilpotence proves that the
polynomial automorphism preserves \(L\). Condition (A2) says \(p=a\phi(y)\).
Its inverse increases degrees by at most three, proving the stated all-order
and mixed conclusions exactly as in (9). This completes the proposition.

The [affine-slice checker](../scripts/verify_two_pair_image_cubic_affine_slice.py)
derives the moments directly from the sphere monomial rule and verifies the
normalizations, all exceptional branches, (A5)–(A10), and the Bezout identity.
This result does not prove that every pure cubic enters (A1).


## What still needs full-family coverage

The frame criterion above identifies a precise remaining reduction, rather than
a proof of purity. If it holds, every phase can be written in a suitable frame as

\[
 p=y(a+bt+ct^2)+G(x,t),\qquad \deg G\leq3.
 \tag{12}
\]

To obtain this thirteen-parameter form, eliminate xy in a representative which
is affine linear in y. Its coefficients are independent in basis (2). The
ten-parameter slice omits yt, yt^2, and t^3. The following two obligations
therefore remain:

1. Prove that an arbitrary all-order pure point in the original sixteen-parameter
   family satisfies the null-conic incidence (7), or treat the points which do
   not satisfy it by another argument.
2. Classify the complete thirteen-parameter family (12), including all vanishing
   and repeated-factor branches, with the spherical functional (3).

These are actual missing implications. The parity results, the ten-parameter
slice, its orthogonal translates, and finite checks of further examples do not
supply them. No universal cubic moment cutoff is claimed. The one-profile
quartic mechanism and radial padding do not address either obligation.

## Reproduction and assurance

Use the existing environment pinned in
[research/requirements.txt](../../requirements.txt), from the repository root:

```sh
.venv/bin/python research/foundations/scripts/verify_two_pair_image_cubic_structure.py
.venv/bin/python research/foundations/scripts/verify_two_pair_image_cubic_phase_slice.py
.venv/bin/python research/foundations/scripts/verify_two_pair_image_cubic_affine_slice.py
```

The [structure checker](../scripts/verify_two_pair_image_cubic_structure.py)
verifies the two invertible sixteen-dimensional coordinate maps, the complete
quadratic trace identities, the exact binary frame identities and affine
representative correction, and the sphere shear and its inverse. It also runs
bounded monomial and power regressions with failure controls. The all-order
conclusions use the displayed invariant and locally nilpotent derivation
arguments, not those finite power loops. Its
[receipt](two_pair_image_cubic_structure_verification.json) distinguishes these
scopes.

The [seven-parameter checker](../scripts/verify_two_pair_image_cubic_phase_slice.py)
verifies exact symbolic moment equations and the literal Bezout certificate,
including both parameter branches and the sharp six-moment control. Its
[receipt](two_pair_image_cubic_phase_slice_verification.json) retains the exact
inputs and identities. The
[ten-parameter checker](../scripts/verify_two_pair_image_cubic_affine_slice.py)
verifies the normalization, eight exact symbolic moment equations, the
rank-deficient coefficient matrix and coprime augmented minors, and all
exceptional branches. Its
[receipt](two_pair_image_cubic_affine_slice_verification.json) retains the
eliminations and rational Bezout coefficients.

All three runs are internal exact replays with SymPy 1.14.0. Additional internal
agents reviewed the coordinate/multiplier/frame argument, the parity proofs,
and the seven-parameter certificate. A separate internal calculation recovered
the ten-parameter obstruction polynomials in different coordinates. This is not
external mathematical review or formal verification. The odd-parity result
retains its explicit dependency on the existing internal Gaussian
classification; none of these checkers certifies that theorem or closes the
full mixed family.
