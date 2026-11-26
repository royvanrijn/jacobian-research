# Hyperelliptic octic conductors and four-period rigidity

Let \((P,Q)\) be an actual complex polynomial Keller map with
\(P\) monic of positive degree in the source coordinate \(x\).
Write \(K=\mathbb C(y)\). Suppose its generic coordinate order
has the hyperelliptic presentation
\[
K[P,Q]=K[u^2,u g(u^2)]\subset K[u],\qquad
u=x-h(y),\qquad \deg_T g=4.
\tag{1}
\]
The same hypothesis can be made over an algebraic closure of \(K\)
and after an affine change of \(x\) there. Then this Keller map is
impossible. This excludes the hyperelliptic branch of global
conductor degree eight, including changes of ambient embedding
and possible drops in the specialized conductor degree. Equivalently,
an open degree-eight case must have generic order with least positive
source degree at least three, even after algebraic extension of the
parameter field.

## 1. The centre and the four reflected pairs

The
[exact coordinate-order conductor theorem](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md)
identifies the specialized order on each \(y=c\), its conductor and
its collision support. Its conductor is a nonzero nonconstant
polynomial on every source fibre of a noninvertible Keller map.
To see the degree-two characterization, centre and rescale any
quadratic element of the generic order to give \(T=u^2\). Every
element has even part in \(K[T]\), so its odd part lies in the order
as well. The odd parts are \(u\) times an ideal of \(K[T]\), because
the order contains \(T\). That ideal is nonzero since the order
has fraction field \(K(u)\). Its generator is \(g\), giving (1)
over the chosen field. Degree eight of the conductor gives
\(\deg_T g=4\) by the calculation below. An element of degree one
would make the whole order the normalization. Thus the residual
degree claim follows from the exclusion of (1).

The conductor of (1) is \(g(u^2)\,K[u]\). For membership in the
conductor, both a polynomial and its product by \(u\) must lie in
the order, which is sufficient since \(1,u\) generate \(K[u]\)
over \(K[T]\). The two tests require both coefficients in the
\(K[T]\)-basis \(1,u\) to be divisible by \(g\). The centroid of the eight roots, counted with
multiplicity, is \(h=-\kappa_7/(8\kappa_8)\), where \(\kappa_i\) are the
coefficients of the exact source conductor \(\kappa\). Thus \(h\)
lies in \(K\), and scaling the affine coordinate has no effect on
the presentation in the centred variable over the algebraic
closure. The even and odd parts descend over \(K\), so (1) follows
over \(K\) with \(g(u^2)=\kappa(u+h,y)\).

If \(h\) had a finite pole, choose a conductor root branch approaching
a finite root of the nonconstant specialized conductor. Such a branch
exists by Weierstrass preparation, or by continuity of the finite
roots near a small circle about the chosen root. Its mate \(2h-x\)
has the same value of the monic \(P\) but tends to infinity, which
monicity forbids. Hence \(h\in\mathbb C[y]\). In the polynomial
frame \((u,y)\), the Jacobian is unchanged and
\[
g(T)=\sum_{j=0}^4 a_j(y)T^j,\qquad a_j\in\mathbb C[y],\quad a_4\ne0.
\]
Its constant coefficient \(d=a_0\) is a nonzero constant. An identically
zero constant coefficient would give \(P_u=Q_u=0\) at \(u=0\) from (1).
At a finite zero, a reflected root pair approaches \(u=0\). Dividing
the collision identities by \(2u\) and taking the limit again gives
\(P_u=Q_u=0\). This excludes every finite zero.

The four roots \(A_i\) of \(g\) are generically simple. At a
persistent repeated root the odd parts of \(P,Q\) vanish to order
at least two at both \(r\) and \(-r\), with \(r^2=A_i\ne0\).
The \(u\)-derivatives at the mates are therefore opposite, and
differentiating the equality of the two values makes their
\(y\)-derivatives equal. The Jacobian would change sign between
the two mates. We may work with four simple pairs over a parameter
open and then use the exact conductor again at exceptional
parameters.

## 2. The two generator types for the exact order

The pole semigroup of (1) is \(\langle2,9\rangle\). A pair of actual
generators can be reduced over \(K\) by interchanging them or
subtracting a monomial in the lower-degree generator from the
higher whenever its degree is a multiple of the lower degree.
Each step is an ambient automorphism and lowers the degree sum.
The reduced degrees satisfy \(n>m>\gcd(n,m)\); neither generator
can be constant because the normalization has fraction field
\(K(u)\) and the order has no element of degree one. The same
reduction is available independently at each complex fibre if its
exact order is (1).

The one-place semigroup theorem of
[Assi and García-Sánchez, Propositions 2 and 4 and §6](sources/assi-garcia-sanchez-1407.0490v1.pdf)
gives a characteristic sequence \(r_0=n,\ r_1=m,\ldots,r_h\).
To check its hypotheses, extend the ground field to an algebraic
closure. The normalization is the \(u\)-line, which supplies a
single place at infinity. The irreducible relation, monic in the
coordinate of degree \(m\), has degree \(n\) in that coordinate.
The coefficient at index \(i\) has degree at most
\(\lfloor im/n\rfloor<i\) in the other coordinate: take elementary
symmetric functions over the field of that coordinate and bound
their poles at infinity. The norm of the coordinate of degree
\(m\) has degree \(m\), and a constant shear removes the first
coefficient. Define the decreasing gcds
\(d_i=\gcd(r_0,\ldots,r_{i-1})\), \(e_i=d_i/d_{i+1}\), with
\(d_{h+1}=1\). The theorem supplies the conductor formula
\(\mu=\sum_{i=1}^h(e_i-1)r_i-r_0+1=8\), and the \(r_i\) generate
the semigroup.

Here is the short specialization of that formula. For a prefix
ending before \(r_i\) put
\[
\mu_i=\frac{\sum_{j<i}(e_j-1)r_j-r_0+d_i}{d_i}.
\]
The first prefix has \(\mu_2=(n/d_2-1)(m/d_2-1)\ge2\), since
the two scaled degrees are distinct coprime integers greater than
one. Each extension has
\(d_{i+1}\mu_{i+1}\ge d_i\mu_i\), using
\((e_i-1)r_i\ge d_i-d_{i+1}\). Some \(r_i\) must be 2. If \(r_1=2\),
the gcd is already 1 and the formula gives \(n=9\).
If the first 2 occurs later and all preceding degrees are even,
the indispensable generator 9 must follow it, reducing the gcd
from 2 to 1. The formula then gives \(8=2\mu_h+8\), contrary to
\(\mu_h\ge2\). In the remaining case 2 is last and \(d_h\) is odd
and at least 3. Then \(8=d_h\mu_h+d_h-1\), so \(d_h=3\) and
\(\mu_h=2\). The prefix inequality forces \(d_h=d_2,\ h=2\), and
the first scaled degrees are 3 and 2. This leaves exactly
\[
(r_0,r_1)=(9,2),\qquad (r_0,r_1,r_2)=(9,6,2).
\tag{2}
\]

In the \((9,2)\) case a degree-two generator is affine in
\(T=u^2\), and the degree-nine generator has odd part a nonzero
constant multiple of \(V=u g(T)\). A triangular shear yields
\((T,V)\). In the \((9,6,2)\) case let \(U,W\) be the reduced
generators of degrees 9 and 6. Proposition 2 supplies the monic
approximate root of degree 3 in \(W\), with source degree 2 and
hence affine in \(T\). The weights of the summands in \(W^2,W,1\)
are distinct modulo 9. Cancellation with \(W^3\) at weight 18
therefore permits only a quadratic in \(U\) in the constant
coefficient, an affine coefficient in \(U\) for \(W\), and a
constant coefficient for \(W^2\). After scaling the quadratic term,
the resulting identity has the form
\[
U^2+\alpha UW+\beta U+\psi(W)=\eta T+\theta,\qquad \eta\ne0,
\]
where \(\psi\) is cubic. An affine shear completes the square
in \(U\). The other terms are even in \(u\), since \(W\) has
degree 6. Since the odd part of \(U\) is a nonzero multiple of
\(V\), evenness of the square then forces its even part to vanish.
Scale \(W=q(T)\) to a monic cubic and shift it by a constant to
depress the cubic in \(W\) on the right. With \(G=g/a_4\) monic, this gives
\[
T G(T)^2=q(T)^3+kq(T)+c+eT,\qquad e\ne0 .
\tag{3}
\]
These changes are affine or triangular with nonzero determinant.
Over a parameter open where the reductions are valid, their
determinants and coefficients are rational in \(y\).

## 3. The degree-six generator case

Write \(q=T^3+pT^2+sT+r\) and \(G=T^4+A T^3+B T^2+C T+D\).
Comparison in (3) at degrees 8 through 5 determines \(A,B,C,D\).
The coefficients at degrees 4 and 2 after removing the cubic
term by \(kq\) give, up to nonzero rational factors,
\[
(p^2-4s)(p^3-4ps+16r),\qquad
(p^2-4s)(25p^5-128p^3s+408p^2r+112ps^2-480sr).
\]
On \(p^2=4s\) the coefficient \(e\) vanishes. Otherwise the first
factorization gives \(r=-p(p^2-4s)/16\), and the second remaining
factor specializes to \(-p(p^2-4s)^2/2\). Thus \(p=r=0\) and
\(s\ne0\). The exact coefficient identities are reproduced by
[the octic checker](verify_hyperelliptic_octic_conductors.py).
In this case
\[
q=T^3+sT,\qquad g=a_4(T^4+\tfrac32sT^2+\tfrac38s^2),\qquad a_2^2=6d a_4 .
\tag{4}
\]
At a zero of \(a_4\), (4) forces the conductor to specialize to
the constant \(d\), excluded by the exact order theorem. Thus
\(a_4,s,g\) are constant.

Define the physical action of the \(i\)-th matched pair by
\[
I_i(y)=\int_{-r_i}^{r_i}P(u,y)\,\partial_u Q(u,y)\,du,\qquad r_i^2=A_i.
\]
On the simple-root open the endpoint equalities give \(I_i'=-2Jr_i\).
In detail the boundary contributions in the derivative combine
into the common endpoint value of \(P\) times the derivative of the
difference in the endpoint values of \(Q\). The integral left
after integration by parts has integrand
\(P_y Q_u-P_uQ_y=-J\). Under an ambient affine or triangular change
this action scales by the determinant: the difference of the
corresponding one-forms is closed on the affine target plane,
hence exact, and both endpoints are the same target point.
Thus all physical actions are given by the canonical actions in
the coordinates \((q,V)\) times one common nonzero rational factor
\(\Lambda(y)\), with orientation incorporated into \(\Lambda\).
Integration by parts gives the canonical actions divided by \(r_i\) as
\[
-4\sum_j a_j\left(\frac{3 A_i^{j+3}}{2j+7}
+\frac{s A_i^{j+1}}{2j+3}\right).
\]
All indices with \(a_j\ne0\) are even, and the four roots occur
in pairs \(A,-A\). The displayed function has opposite values
at them. It cannot, after multiplication by \(\Lambda'\), equal
the same nonzero constant \(-2J\) for each root. This excludes
the generic degree-six generator case.

## 4. The degree-two generator case

For \((T,V)\), integration by parts gives
\[
I_i=\Lambda r_i C(A_i),\qquad
C(T)=-4\sum_{j=0}^4 a_j\,T^{j+1}/(2j+3),
\tag{5}
\]
for a common nonzero rational \(\Lambda\). The canonical vector
is nonzero whenever \(g\) has four simple nonzero roots:
vanishing would make \(\sum_j a_j T^j/(2j+3)\) proportional to
\(g\), impossible by comparing the constant and leading terms.

Differentiate the canonical action using its integral
\(-\int_{-r_i}^{r_i}2u^2g(u^2)\,du\) multiplied by \(\Lambda\).
The moving endpoint terms vanish since \(g(A_i)=0\). Combining
with \(I_i'=-2Jr_i\) gives
\[
\sum_{j=0}^4(\Lambda a_j)'T^{j+1}/(2j+3)-J/2
=(\alpha T+\beta)g(T),\qquad
\alpha=(\Lambda a_4)'/(11a_4),\quad\beta=-J/(2d).
\]
Equivalently
\[
(\Lambda a_j)'/(2j+3)-\alpha a_j=\beta a_{j+1},
\qquad j=0,1,2,3. \tag{6}
\]

Suppose \(a_4\) has a finite zero of order \(m>0\).
Let \(k\) be the largest index with \(a_k\) a unit there.
The specialized conductor is nonconstant, so \(1\le k\le3\).
Write \(n=\operatorname{ord}\Lambda\) and \(w_j=\operatorname{ord}a_j\)
for nonzero \(a_j\). The term of order \(n+w_j-1\) on the left
of (6) has coefficient proportional to
\[
\frac{n+w_j}{2j+3}-\frac{n+m}{11}.
\]
For \(n<1\), the unit coefficients at \(j=0,k\) would require
both these numbers to vanish, impossible for \(m>0\).
For \(n>1\), the equation at \(j=k-1\) has positive-order left
side and unit right side; the same contradiction holds if
\(a_{k-1}=0\). For \(n=1\), the equation at \(j=k\) forces
\(m=(8-2k)/(2k+3)\), never a positive integer for \(k=1,2,3\).
This proves \(a_4\) constant.

Now (6) at \(j=0\) makes \(\Lambda'\) polynomial, so \(\Lambda\)
is polynomial. It cannot be constant, since (6) would force
all \(a_j\) with \(j>0\) to vanish. If \(\deg \Lambda=L\ge1\),
induction in (6) gives \(\deg a_j=j(L-1)\). Indeed at step \(j<4\)
the leading coefficient on the left is nonzero by
\[
j(L-1)/(2j+3)+L\bigl(1/(2j+3)-1/11\bigr)>0 .
\]
Since \(a_4\) is constant, \(L=1\) and all \(a_j\) are constant.
With \(\xi=\Lambda'/(11\beta)\ne0\), the recurrence gives
\[
g(T)/d=1+\tfrac83(\xi T)+\tfrac{16}{5}(\xi T)^2
+\tfrac{64}{35}(\xi T)^3+\tfrac{128}{315}(\xi T)^4. \tag{7}
\]

Finally let \(y_0\) be the zero of \(\Lambda\).
The specialized order at \(y_0\) is still \(\mathbb C[T,V]\).
The four fixed distinct pairs remain matched by specialization.
Their four independent equality conditions on \(\mathbb C[u]\)
cut out \(\mathbb C[T,V]\), by interpolation and the description of
the odd part. The exact specialized conductor has degree eight,
so the specialized order also has codimension four and equals
this ring. Its generator degrees can again be classified by (2).
The second type would require (4) for this same \(g\), whereas (7)
has nonzero cubic coefficient, preserved under any nonzero scaling
of the centred variable. Hence the first type applies, and its
action vector is nonzero by (5) and the nonzero determinant of
the specialized target change. But (5) continues polynomially
in \(y\) for these fixed endpoints, so all physical actions
vanish at \(y_0\). This contradiction excludes (1).
