# Linear coefficient pencils and units on marked components

Work geometrically in characteristic zero. The preceding
[coefficient graphs](COEFFICIENT_GRAPHS.md) leave q-dependent graphs
outside the strict weight test. This note treats every graph that is
affine-linear in q,r, with arbitrary dependence on b, and gives separate
source obstructions for components selected by common factors.

The classical input is characteristic-zero symmetric sectional monodromy:
[Kadets, Proposition1.1](https://arxiv.org/html/1809.07293v4).
The cofactor-unit fact is
[Fonollosa, Lemma6.2](https://arxiv.org/html/2609.33491v1);
we give an explicit inverse and derive the applications below.
The [source pin](sources/linear-pencil-manifest.json) distinguishes these
prior inputs from the coefficient-graph deductions.
The argument does not use the positive-characteristic classifications.

[Exact field helper](linear_pencil.py) · [checker](verify_linear_pencil.py) ·
[model](linear_pencil_model.json) · [receipt](linear_pencil_receipt.json).

## Statements and scope

Let F=k(b), let q,r be independent, and write
\[
 P(T)=qH(T)+rJ(T)+L(T),\quad
 \deg H=n,\quad\deg J=n-1,\quad\deg L\leq n-2,                 \tag{1}
\]
with H,J monic and L nonzero. In the coefficient-graph application
\[
 h_j=a_j(b)q+d_j(b)r+c_j(b),\quad
 H=T^n+\sum a_jT^j,\quad J=T^{n-1}+\sum d_jT^j,\quad
 L=T-b+\sum c_jT^j.                                         \tag{2}
\]
The coefficients may be rational in b for the field assertions.
For statements about actual polynomial inverse-image varieties,
all graph coefficients and factorizations are global polynomials.

**KMFI77 (the moving pencil has symmetric or imprimitive root action).**
Remove the common factor D=gcd(H,J,L), and put h=H/D,j=J/D,l=L/D,
g=deg D and N=n-g. The moving polynomial qh+rj+l is geometrically
irreducible. Let
\[
 e=[\overline F(T):\overline F(h/l,j/l)].
\]
Then e divides N, and d=N/e>=2 is the degree of the image plane curve.
If e=1 the moving root group is S_N. If e>1 its root action is
imprimitive, with d blocks of size e and full S_d on the blocks.
No uniform full-wreath equality is asserted in the imprimitive case.
Consequently a primitive moving root action is always symmetric.
In particular an irreducible prime-degree pencil, including degree11,
cannot have proper primitive PSL2(11) root monodromy.
Static factors D are recorded separately, not silently discarded.

When D=1 and e=1, a general ordinary tangent gives exactly one double
root, with leading coefficient q a unit. Thus the derivative-marked
field K(t,a), a^(n-2)=P'(t), has the same full G_n, arithmetic action
and constants as KMFI69--70. This does not supply an affine source
for every graph in(2).

**KMFI78 (complete root classification without the offset).**
For arbitrary a_j(b),d_j(b) and c_j=0, the full root group is:
\[
 \begin{cases}
 S_n,&(H(b),J(b))\ne(0,0),\\
 1+S_{n-1},&H(b)=J(b)=0.
 \end{cases}                                                \tag{3}
\]
The second line fixes the rational root T=b and acts symmetrically
on the remaining roots. There are no other cases, for all n>=3.
In the first case the full marked and arithmetic groups are G_n
and the cyclotomic groups of KMFI70.
Thus q-dependent linear graphs without an offset cannot produce a
proper transitive root group, even outside the earlier weight gate.
For prime-degree irreducible graphs the same exclusion holds with
an arbitrary offset c_j(b), by KMFI77.
A common factor of a higher-degree pencil may have another root
group;(3) is not a theorem about every possible static factor.

**KMFI79 (a base factor forbids the actual moving affine component).**
Suppose (2) is polynomial over k[b], n>=4, and D is a global primitive
representative of the complete gcd of H,J,L, of positive degree. On the actual marked component
belonging to the moving factor, D_h(x,p) is a nonconstant unit,
where p=1+xy and D_h is the degree-g homogenization.
This component has a rational function field, but is not stably
ordinary affine space. Its marked degree over the target is N(n-2).
This excludes that actual component as an affine-space Keller source;
it is not an absolute exclusion of its abstract monodromy action
or of every different birational field model.

**KMFI80 (monic factors and constant leading coefficient).**
More generally pull the original marked-root map back along any
polynomial target morphism Gamma=A^s -> A^n.
If an actual component belongs to a generic irreducible factor
G(T) in k[Gamma][T] that is monic of degree >=2, then p=1+xy is
a nonconstant unit on that component. Hence it is not stably
ordinary affine space. More generally this remains true if a polynomial
unimodular change of the projective root coordinates makes the selected
factor monic: the corresponding denominator linear form in x,p is
the nonconstant unit.

In particular, if the leading coefficient q is a nonzero constant
on Gamma, no component is a nonlinear Keller self-map of ordinary
affine space: nonlinear root factors have this unit obstruction;
linear root factors have cyclic marked extensions, excluded by
KMFI2 unless their degree is1. This closes the literal constant-leading
base-change route, including any fixed monic Klein curve polynomial
with independent target spectators. It does not exclude a different
birational target/source presentation of a PSL2(11) action.

The main objective remains open, including exceptional degree11 PSL2(11).

## 1. The moving plane curve

Divide the three sections in(1) by their common factor. Their
homogenizations to degree N have no common zero on P1: the finite
ones were removed, and h is monic of degree N at infinity.
They define a morphism
\[
       \nu:\mathbb P^1_T\longrightarrow C\subset\mathbb P^2,
                     \qquad [h:j:l].
\]
The image is integral and is not a line, since the three polynomials
have distinct degrees N,N-1,at most N-2. Over closed constants let
the map to the normalization of C have degree e.
The pullback of a line has degree N, so N=e*deg C=e*d, and d>=2.

A general line qX+rY+Z=0 cuts C in d points. Their sectional group
is S_d by the cited characteristic-zero theorem. The e preimages
of each point form a block of roots of qh+rj+l. The incidence over
P1 is irreducible (a projective-line bundle before restriction to the
dual chart), so the moving polynomial is irreducible and the action
is transitive. If e>1 the blocks are proper; if e=1, the same
sectional cover is the root cover, giving S_N.
The degree calculation proves the prime-degree consequence without
a bounded list of subgroups or a sampled coefficient search.

For completeness, sectional S_d can also be seen from ordered
secants, whose incidence is irreducible, giving double transitivity,
and a general ordinary tangent, giving one transposition.
In characteristic zero the Gauss map is generically one-to-one:
the tangent point is recovered by intersecting the generic tangent
line with its derivative along the dual curve. The nonzero
three-section Wronskian gives ordinary contact generically.
Thus a general tangent meets the curve at only one double point.
Its dual locus is not a line, so it meets the chart with constant
coefficient1 and q!=0. In the birational case its roots supply the
single colliding pair required by the Kummer proof of KMFI69.
The generic derivative-product identity and arithmetic bounds
then give exactly its marked group and cyclotomic constants.

An exact way to certify e in finite examples is to take, over
F(T), the gcd in S of
\[
 h(S)l(T)-h(T)l(S),\qquad j(S)l(T)-j(T)l(S).                  \tag{4}
\]
Its degree is e. Removing D first is essential: otherwise base
points introduce spurious common factors. The helper implements(4).

## 2. Why no offset gives only the two groups in(3)

Here L=T-b. A common factor is therefore either1 or T-b.

If H(b),J(b) are both nonzero, H/L and J/L have respective
rational-map degrees n and n-1. The intermediate degree e divides
both, so e=1.

If H(b)=0 and J(b)!=0, then H/L is a polynomial of degree n-1,
whereas J/L has a unique finite pole, at b, of order1.
In an intermediate field containing a nonconstant polynomial,
the place under T=infinity is totally ramified in F(T).
The simple finite pole lies above a different place (the polynomial
is finite there). All points above that place must be poles of J/L;
only b is available. Its ramification index is therefore e, and
must divide the pole order1. Hence e=1.
Interchanging H,J gives the case H(b)!=0,J(b)=0, using the
polynomial J/L of degree n-2.

Finally, if both vanish, divide by T-b. The denominator becomes1,
and H/(T-b),J/(T-b) are polynomials of coprime degrees n-1,n-2.
The intermediate degree again divides both and equals1.
The moving polynomial has group S_(n-1), and T=b is a single
fixed root. This proves all cases, including exceptional equalities
in the earlier weight test.

## 3. A nonconstant cofactor unit on every moving component

The marked-root construction satisfies the global polynomial identities
\[
 P_h(x,p)=0,\quad
 (\partial_XP_h)(x,p)=p,\quad
 (\partial_ZP_h)(x,p)=-x,\qquad p=1+xy.                       \tag{5}
\]
The first two are the marking identities; the third follows from
Euler's relation on the dense p!=0 chart and then polynomiality.

For any global factorization P_h=G_h H_h, on a component with
G_h(x,p)=0, (5) gives the explicit unit inverse
\[
 H_h(x,p)\left((\partial_XG_h)(x,p)
                    +y(\partial_ZG_h)(x,p)\right)=1.          \tag{6}
\]
This is the attributed cofactor-unit fact with its polynomial inverse.

Apply it to G=qh+rj+l, H=D. Let m=n-2 and work in the generic
moving root field. Solving the root equation gives
\[
 q=-(rj+l)/h,\qquad
 P'(t)=D(t)(rW(t)+V(t)),                                    \tag{7}
\]
where W=j'-j*h'/h and V=l'-l*h'/h.
W is nonzero because h,j have distinct degrees.
Thus P'(t) has a simple linear r-prime over F(t).
The marking a^m=P'(t) has full degree m by its valuation1.
Equation(7) solves r rationally from t,a, then q, so the marked
field is exactly F(t,a), a rational field of dimension3 over k.

On the actual component, the unit from(6) is
\[
                         U=D(t)/a^g.
\]
Since D divides L, 1<=g<=m. If g<m, a deck rotation of order m
changes U, so U cannot be a constant. If g=m, U=1/(rW+V),
which is nonconstant because W!=0.
The unit persists after independent polynomial stabilization,
excluding the actual source component from being affine space.
No assumption about rationality being sufficient is used.

The special common factor T-b has unit x-bp on its moving component.
Equation(6) specializes to the explicit inverse
\[
 (x-bp)^{-1}
      =(\partial_XG_h)(x,p)(1-by)-y,
\]
where now G_h denotes the moving degree-(n-1) factor.
Its generic value (t-b)/a is nonconstant. On the fixed-root
component p has inverse1-by, recovering KMFI76 when its derivative
is nonconstant.

## 4. Monic factors force a different nonconstant unit

Let A=k[Gamma] be a polynomial ring, and let G in A[T] be monic
and generically irreducible of degree e>=2, with P=G H.
Let S be an actual component belonging to G.
It is smooth, hence normal, because it is a component of an etale
base change over smooth Gamma.

Generically t=x/p satisfies the monic equation G(t)=0.
It is integral over O(S) and thus regular on S by normality.
Consequently
\[
                              p(1-yt)=1
\]
on all of S. This proves p is a global unit, not just a nonzero
generic rational function.

It cannot be a constant. Otherwise a=1/p is a nonzero constant,
and the marking gives
\[
                  G'(t)H(t)=P'(t)=a^m\in k^*.
\]
This identity holds already in R=A[T]/(G), since R embeds into
the selected root field and then into the component field.
Thus G' is a unit in R. The finite free algebra R/A is etale.
It is connected by generic irreducibility, so simple connectedness
of A^s in characteristic zero forces its rank e to be1,
contradicting e>=2. Hence p is a nonconstant unit.

The same argument is invariant under polynomial root frames.
Let M in GL2(A), with determinant delta in k*, and put
(X',Z')=M(X,Z), (x',p')=M(x,p). The pair x',p' is unimodular
since x,p is: p-xy=1. If the selected factor becomes monic in
u=X'/Z', then u=x'/p' is integral and regular on S. A polynomial
Bezout relation for x',p' then makes p' a unit.
For the transformed whole binary polynomial the marking gradient is
\[
       \nabla P'_h(x',p')=M^{-T}(p,-x)
                            =\delta^{-1}(p',-x').
\]
If p' were constant, its dehomogenized derivative would be the
constant delta^(-1)*p'^(-m), again making the monic root algebra
finite etale and contradicting its nonlinear rank. Thus p' is a
nonconstant unit. Polynomial translations, inversions and their
unimodular parameter-dependent combinations cannot remove the
obstruction by merely changing the root frame.

If q is a nonzero constant, P is monic after a scalar division.
Every generic irreducible factor is monic over A: its coefficients
are integral over the integrally closed A.
Nonlinear factors fall under the preceding argument.
For a linear factor t is in A and the marked component is a cyclic
extension of K of degree dividing m. If that degree exceeds1,
KMFI2 excludes an affine-space Keller model. Degree1 Keller maps
are automorphisms. This proves the constant-leading assertion with
all components and independent stabilizations accounted for.

This obstruction concerns the actual inverse-image component.
Replacing the affine target/source by a different birational model
can change which functions remain regular; an abstract action
exclusion is not inferred from this unit alone.

## 5. Boundary controls

There really are composite-degree imprimitive pencils. Over k(b),
\[
 H=T^4,\quad J=T^3-2bT^2,\quad
 L=-T^2/(4b)+T-b
\]
factor through phi=T^2/(T-2b):
\[
        P/(T-2b)^2=q\,\phi^2+r\,\phi-1/(4b).
\]
The root group is D4=C2 wr S2, not S4.
Over the splitting field of the two phi-roots u1,u2, the inverse
quadratic discriminants u_i(u_i-8b) have independent square classes:
the divisor u_i=0 has value1 only on the corresponding one.
Thus the order is exactly8. This example has a rational-in-b
offset, and its b=-1 slice is a polynomial two-parameter example.
It shows why the general statement is symmetric-or-imprimitive and
why the prime-degree and no-offset hypotheses must be respected.

There are also explicit rational marked components with monic factors
of any degree. Let f be monic of degree d>=2, f(0)=0,f'(0)=1, and
let c be its T^(d-1) coefficient. Put n=d+3,m=d+1 and
\[
 P=(f(T)-b)\left(1+rT^2+q(T^3-cT^2)\right).                 \tag{8}
\]
Its top coefficients are exactly q,r, its linear coefficient is1
and its constant is-b, so it is a polynomial affine-linear graph
of the original family. On the f(t)=b component,
\[
 a^m=f'(t)\left(1+rt^2+q(t^3-ct^2)\right).
\]
Solving q shows its function field is k(t,a,r), rational3.
Its actual source nevertheless has the nonconstant unit p by Section4.

If f has polynomial monodromy G_f on d roots, this relative marked
field has geometric group C_m wr G_f: over the root splitting field,
the cofactor linear forms at distinct nonzero roots are distinct,
and each gives an individual value-one Kummer divisor.
All m rotations of every root block are therefore independent.
This supplies a precise relative group calculation, not a new
ordinary-affine realization. In particular inserting a monic
PSL2(11) curve factor this way does not solve the affine-source problem.

Exact checks retain base factors, fiber-gcd degrees, the D4 composition,
independent root classes, rational field reconstructions, homogeneous
derivative identities and polynomial unit certificates.
The unbounded conclusions have the written proofs above.
No formal verification, independent whole-proof review, full
action classification or absolute exceptional PSL2(11) decision is claimed.
