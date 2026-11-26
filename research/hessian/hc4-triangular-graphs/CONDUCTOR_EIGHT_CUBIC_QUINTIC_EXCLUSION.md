# The cubic–quintic conductor-eight exclusion

Let \(F=(P,Q)\) be an actual complex polynomial plane Keller map
in a source frame \((x,y)\) with \(P\) monic of positive degree in \(x\).
Suppose its exact coordinate conductor has global \(x\)-degree eight.
Then the generic order \(K[P,Q]\subset K[x]\), \(K=\mathbb C(y)\), has
pole semigroup \(\langle4,5,6\rangle\), and its reduced one-place
characteristic sequence is \((6,4,5)\). Thus the remaining embedding
has generators of source degrees six and four and an approximate
root of degree five. This last type is treated in the separate
[global octic conductor proof](GLOBAL_OCTIC_CONDUCTOR_EXCLUSION.md).

The proof first reduces the singular rational orders directly,
using the one-place characteristic sequence. It then excludes the
cubic–quintic pair by a parameter-dependent target action equation
and a valuation argument through the finite source fibres.

## 1. The characteristic sequences of the exact order

The
[coordinate-order conductor theorem](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md)
gives the exact specialization on each line \(y=c\). Its defect is
half the specialized conductor degree. In particular the generic
order is birational of defect four, and each specialized defect is
between one and four in the hypothetical noninvertible case.

In any two-generator polynomial order with fraction field \(K(x)\),
reduce the two degrees by subtracting an appropriate monomial in the
lower-degree generator whenever it divides the higher degree. These
are triangular target automorphisms, over \(K\) if the starting pair is
defined over \(K\). The process ends at degrees \(n>m>\gcd(n,m)\)
when the order has positive defect. After extending \(K\) algebraically,
the plane curve has one place at infinity, since its normalization is
the \(x\)-line. Its equation is monic of degree \(n\) in the lower-degree
generator. The coefficient of index \(i\) in that equation has degree
at most \(\lfloor im/n\rfloor<i\) in the other generator, by the pole
bound on the elementary symmetric functions over that generator.
A constant translation kills the first coefficient. The norm of
the reduced generator of degree \(m\) has degree \(m\) in the
higher-degree generator.

Thus the standing hypotheses of
[Assi–García-Sánchez, Propositions 2 and 4 and §6](sources/assi-garcia-sanchez-1407.0490v1.pdf)
apply. The characteristic degrees \(r_0=n,r_1=m,\ldots,r_h\)
generate a symmetric numerical semigroup. Its number of gaps equals
the defect of the order, by the filtration by \(x\)-degree. Therefore
the semigroup conductor is eight. Let \(d_i=\gcd(r_0,\ldots,r_{i-1})\),
\(e_i=d_i/d_{i+1}\), \(d_{h+1}=1\); these gcds strictly decrease.
If \(h=1\), then \((n-1)(m-1)=8\), giving \((9,2)\) or \((5,3)\).
For \(h\ge2\) set
\[
\mu_i=\frac{\sum_{j<i}(e_j-1)r_j-r_0+d_i}{d_i}\quad (i\ge2).
\]
The first prefix has \(\mu_2=(n/d_2-1)(m/d_2-1)\ge2\).
Since \((e_i-1)r_i\ge d_i-d_{i+1}\), the quantities \(d_i\mu_i\)
are nondecreasing. For \(h\ge3\) one would have \(d_2\ge4\) and
\[
8=d_h\mu_h+(d_h-1)(r_h-1)\ \ge\ 2d_2+1,
\]
as \(d_h\ge2\) and no degree-one element belongs to the order.
For \(h=2\) the same formula leaves just
\((d_2,\mu_2,r_2)=(2,2,5),(3,2,2)\), giving \((6,4,5)\) and
\((9,6,2)\), respectively. These are necessary sequences.
The two sequences with a degree-two element in the order fall
under the
[hyperelliptic octic exclusion](HYPERELLIPTIC_OCTIC_CONDUCTORS_AND_FOUR_PERIOD_RIGIDITY.md).
The calculation uses the singular generic order itself.

## 2. Polynomial parameter compatibility for the cubic–quintic pair

Suppose the sequence is \((5,3)\). The reductions above, followed
by a source translation \(u=x-h(y)\) with \(h\in K\) and affine target
normalizations, give an actual target automorphism over \(K\) with
determinant \(\lambda\in K^*\) and outputs
\[
U=u^3+a_0u,\qquad V=u^5+b_0u^4+d_0u^2+e_0u .
\]
Normalize \(U\) to be monic, depressed and without constant term.
Subtract a multiple of \(U\) and a constant from the monic \(V\).
The reciprocal determinant is the common scale from the canonical
actions to the physical actions.

Here is the differential constraint for any such parameter-dependent
target automorphism \(\Phi\). Its target vector field
\(\partial_y\Phi\circ\Phi^{-1}\) is polynomial in \(U,V\), with
divergence \(c=\lambda'/\lambda\), by differentiating the target
Jacobian determinant. After subtracting \((0,cV)\), this field has
zero divergence and is of the form \((-H_V,H_U)\) for a polynomial
\(H\), by polynomial integration in characteristic zero. The source
Jacobian \(J=\operatorname{Jac}_{x,y}(P,Q)\) therefore gives
\[
U_u V_y-U_y V_u=\lambda J+c\, U_uV+\partial_u H(U,V).
\tag{1}
\]
Only the value \(H(U,V)\) is needed. Its degree in \(u\) is at most
eight. In pole semigroup \(\langle3,5\rangle\), the degrees at most
eight give the basis \(1,U,V,U^2,UV\). If \(\delta\) is the coefficient
of \(UV\) in \(H(U,V)\), comparison at degrees seven and six in
(1) gives \(\delta=-3c/8\) and \(b_0'=b_0 c/8\).

If \(b_0\ne0\), write \(\rho=b_0\). Otherwise choose an eighth root
\(\rho\) of \(\lambda\), working over \(L=K(\rho)\).
Make the longitudinal source change \(u=\rho z\) and rescale the
two target outputs by \(\rho^{-3},\rho^{-5}\), respectively. In
the resulting frame
\[
U=z^3+az,\qquad V=z^5+bz^4+dz^2+ez,\qquad b\in\{0,1\}.
\]
The target determinant is now the nonzero constant
\(\lambda_*=\lambda/\rho^8\), and the original source Jacobian in
\((z,y)\) is \(\rho J\). The analogue of (1) has \(c=0\) and
\(j=\lambda_*J\rho\) in place of \(\lambda J\).

Let \(B\) denote the coefficient of \(V\) in \(H(U,V)\) and put
\(m=d-4ba/3\). The coefficient equations and their integrals give
\[
\begin{split}
m'&=5B/3,\\
C&=e+10a^2/9-4bm/5\ \in\mathbb C,\\
S_0&=\tfrac53 a^3-9Ca-\tfrac{36}{5}bam-\tfrac{27}{5}m^2
        \ \in\mathbb C,\\
G&=(a^2/3-3C/5)m-6bm^2/25+4ba^3/27,\qquad G'=j.
\end{split}
\tag{2}
\]
For example the degree-one coefficient is
\(24(2ab-3d)B=4a'(5a^2+9e)\), which after inserting the preceding
equations is \(S_0'=0\). The constant coefficient gives \(G'=j\).
The
[exact checker](verify_conductor_eight_cubic_quintic.py)
verifies these identities before division by any coefficient.
The constants of the finite extension \(L/\mathbb C(y)\) with its
extended derivation are still \(\mathbb C\). In particular (2)
also holds on the loci where one of the coefficients vanishes.

## 3. Finite pairs and coefficient poles

For distinct matching \(z_\pm\) write \(s=z_++z_-\).
The cubic divided difference is \(s^2-z_+z_-+a\).
With this zero, the quintic divided difference gives
\[
s^4+bs^3+a s^2+(2ab-d)s-(a^2+e)=0,\qquad
(z_+-z_-)^2=-3s^2-4a.
\tag{3}
\]
Every actual finite source parameter \(y=c\) has at least one
generic matching pair whose \(x\)-coordinates are both finite
and have distinct limits there, after a local extension. Indeed
the specialized conductor is a nonzero nonconstant polynomial.
A generic conductor root branch can be chosen tending to a finite
root of the specialized conductor, by continuity of polynomial
roots near a small circle about that root. The exact generic
conductor identifies a mate. It also stays finite, since its
value of the original monic \(P\) equals that of the finite root.
If both limits coincided, the divided differences of the original
\(P,Q\) would give \(P_x=Q_x=0\) there. Thus at least one physical
gap in the \(x\)-coordinate is a unit at every finite parameter
place. The argument applies over any specified place of \(L\) above
\(c\), by further extending the local field to the paired roots.

Suppose \(a\) or \(m\) has a pole at a place of \(L\) over a finite
parameter with ramification index \(l\). Completing the square
in (2) shows
\[
\left(m+\tfrac23ba\right)^2
=\tfrac{25}{81}a^3+\tfrac49 b^2 a^2-\tfrac53 Ca-\tfrac5{27}S_0.
\]
Thus both have poles, of orders \(2k,3k\) for some \(k>0\); \(G\)
has pole order \(7k\). By \(G'=\lambda_*J\rho\), the pole order
of \(\rho\) is \(7k+l\). In a local uniformizer the leading
coefficients \(A,M\) of \(a,m\) satisfy \(M^2=25A^3/81\).
The Newton polygon of the quartic (3) forces every pair sum
to have order \(-k\), with leading coefficient \(s_0\) satisfying
\[
s_0^4+A s_0^2-Ms_0+A^2/9=0.
\]
This polynomial has no common root with \(3s_0^2+4A\).
At a putative common root the equation would give
\(Ms_0=5A^2/9\), and squaring contradicts \(M^2=25A^3/81\)
and \(s_0^2=-4A/3\). Hence the gap in \(z\) has order \(-k\)
for every pair. Its physical gap \(\rho(z_+-z_-)\) diverges for
every pair, violating the unit-gap condition.

Thus \(a,m\), and hence \(e,G\), are regular at every finite
parameter place. If \(\rho\) had a zero there, all \(z\)-gaps
are finite by (3), so all physical gaps tend to zero, which is
also impossible. If \(\rho\) has no finite poles either, then
\(\rho^8\in K\) is a nonzero constant, hence \(\rho\) is constant.
Now \(a,m,G\) are polynomials in \(y\), with \(G'\ne0\) constant.
If \(a\) is constant, (2) makes \(m,G\) constant. If \(a\) is
nonconstant, the same leading-order calculation at infinity
makes the degree of \(G\) equal to \(7\deg(a)/2>1\).
Both contradict \(G'\ne0\) constant.

When \(b_0\ne0\), the field \(L=K\) has \(l=1\) at every finite
parameter, so \(G'\) is regular there and \(\rho\) has no finite
pole. This closes that case.

## 4. The algebraic scale and its only possible degeneracy

It remains to use \(b=0\) and \(\rho^8=\lambda\), with \(m=d\).
At a finite pole of \(\rho\), now \(l\mid8\) and \(l>1\), since
\(G\) is regular there and \(G'\) has a pole. Inertia acts on
\(\rho\) by a primitive \(l\)-th root of unity \(\zeta\). It acts
on \(a,d,e\) by \(\zeta^{-2},\zeta^{-3},\zeta^{-4}\), since the
unscaled coefficients are rational in \(y\). In particular the
limit of \(d\) is zero. Write \(A\) for the limit of \(a\).
The unit physical gap at this place forces a normalized pair gap
to limit to zero, so (3), \(d\to0\), and (2) give
\(C=5A^2/9\). If \(l=4\) or 8, then \(A=0\) by inertia.
Whenever \(A=0\), (2) now gives \(C=S_0=0\).

In the remaining case \(l=2,\ A\ne0\), write
\[
d^2=\tfrac{25}{81}a^3-\tfrac53 Ca-\tfrac5{27}S_0 .
\tag{4}
\]
The polynomial on the right has exactly a double root at \(A\):
its derivative vanishes by \(C=5A^2/9\), and its second derivative
is nonzero. Thus \(v(d)=v(a-A)\). The right side is even because
\(a\) is inertia invariant. The left side is odd because \(d\)
changes sign under inertia. These valuations are finite since
\(G'\ne0\) excludes the corresponding constant solution.
This contradiction shows that any finite pole of \(\rho\)
requires \(C=S_0=0\).

For \(C=S_0=0\) equation (4) gives, with \(\tau=d/a\),
\[
a=\tfrac{81}{25}\tau^2,\qquad d=\tfrac{81}{25}\tau^3,\qquad
e=-10a^2/9.
\]
The collision equations (3) become fixed equations after scaling
\(s\) by \(\tau\), and every gap is a nonzero constant multiple
of \(\tau\). The absence of a zero gap in the fixed scaled
equations follows from the same coprimality calculation in
Section 3, now an identity of the complete equations at \(b=0\).
Every physical gap is consequently a nonzero constant multiple
of \(w=\rho\tau=d_0/a_0\in K\). The unit-gap condition forces
\(w\) to have neither zeros nor poles at finite \(y\), hence
\(w\in\mathbb C^*\). Now \(G\) is a nonzero constant multiple of
\(\tau^7\). The identity \(G'=J\rho=Jw/\tau\) makes \(\tau^8\)
a nonconstant affine polynomial in \(y\). Therefore
\(\lambda^{-1}=\rho^{-8}=w^{-8}\tau^8\) has a finite zero \(y_0\).

The unscaled reduced pair in \(u=x-h(y)\) is fixed:
\(a_0,d_0,e_0\) are constants as \(w\) is constant. The common
centre \(h\) also has no finite pole, since a pole would make both
endpoints of every matched pair escape at that parameter. Hence
\(h\) is polynomial. Write \(A_*\) for the fixed reduced order in
\(\mathbb C[u]\). Since \(P,Q\) belong to \(K\otimes A_*\) and
are polynomials in \(u,y\), they belong to \(\mathbb C[y]\otimes A_*\);
one can test this coefficientwise by projection to
\(\mathbb C[u]/A_*\). At any specialization their generated order
is contained in \(A_*\) and has colength at most four by the exact
specialized conductor. Since \(\operatorname{codim} A_*=4\), it
equals \(A_*\), in particular at \(y_0\).

For a matched pair the physical action is the integral of
\(P\,dQ\) between the two source points on a line \(y=\text{constant}\).
Its derivative is \(-J\) times the difference of the endpoints,
by differentiating and integrating by parts, using the matched
endpoint values. Under a target automorphism the action scales
by its determinant, since the remaining polynomial one-form
is exact on the target plane. Thus each fixed canonical pair
here has a nonzero action, because its generic physical action
has nonzero derivative. At \(y_0\) all physical actions vanish
by the generic identity with scale \(\lambda^{-1}\), which
continues polynomially for the fixed endpoints in \(u\).
On the other hand the specialized generators of \(A_*\) can
be reduced at \(y_0\) to degrees three and five by Section 1.
In these degrees they are affine changes of the fixed
canonical pair: in degree at most five the filtered order has
basis \(1,U,V\). The specialized physical actions are therefore
nonzero multiples of the fixed canonical actions, a contradiction.

This excludes the entire cubic–quintic branch, including
contacts, coincident collision values and source-fibre degree
drops. In global conductor degree eight it suffices to exclude
sequence \((6,4,5)\), as in the linked global proof.
