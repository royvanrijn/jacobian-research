# Quadratic outer maps on sextic tangent-conic pencils

**Lemma KellerSexticConicLeadingDiscriminantExclusion.**
Let \(A,B,C\) be ordinary homogeneous sextic forms in
\(k[x,y,z]\), where the field has characteristic zero.
Assume \(AC=B^2\), \(A,B\) are algebraically independent,
\(k(A,B)\) is relatively algebraically closed in \(k(x,y,z)\),
and the degree-six part of \(k(A,B)\cap k[x,y,z]\) is
\(\langle A,B,C\rangle\), with no positive homogeneous
parts of degrees below twelve other than six. Then there
is no Keller completion of \((A,B,C)\).

**Theorem SexticKellerTangentConicQuadraticOuterNonpencilSquarefreeExclusion.**
Over the algebraic closure of a characteristic-zero field, consider
the quadratic outer case of the sextic tangent-conic reduction.
Write the highest tuple as \(g(p^2,pq,q^2)\), for coprime
quadratics \(p,q\) and the quadratic common factor \(g\).
If \(g\) is squarefree and lies outside \(\langle p,q\rangle\),
there is no Keller completion. Consequently the only quadratic
outer cases on a tangent-conic pencil are, up to the indicated
source and target normalizations, the following leading tuples:
\[
 p(p^2,pq,q^2)\ (p=xz+y^2,\ q=x^2),\qquad
 q(p^2,pq,q^2)\ (p=xz+y^2,\ q=xy),\qquad
 x^2(p^2,pq,q^2)\ (p=xz+y^2,\ q\text{ binary}).
\]

## Descent of the conic relation

Suppose a completion \(F\) of the lemma existed. Its components
are algebraically independent by the Jacobian criterion.
The polynomial \(D=F_1F_3-F_2^2\) has degree below twelve
by the leading relation. For every constant choice of
\(c_0,c_1,c_2,c_3\), the polynomial
\(D_0=D-\sum_i c_iF_i-c_0\), where the sum is over \(1,2,3\),
is nonzero by independence, and
\[
                   \det J(F_1,F_2,D_0)=(F_1-c_3)\det JF .
\]
If \(D_0\) has positive degree \(d\), the highest
determinant coefficient of expected degree \(9+d>6\)
must vanish. Thus its leading part is algebraic over
\(k(A,B)\) and hence belongs to that field. The degree
hypothesis forces \(d=6\) and that leading part to belong
to \(\langle A,B,C\rangle\). Subtracting the corresponding
linear combination of \(F_i\) lowers the degree.
Once \(D_0\) is constant its bracket is zero, contrary to
the displayed identity. This proves the lemma.

## The leading field of a nonpencil squarefree factor

Work over the algebraic closure. The independent quadratic
outer coordinates span the binary quadratics, so the chosen
target form of the tuple is valid.
Set \(T=q/p\) and \(A=gp^2\). The two leading components
\(A,AT\) are independent since \(T\) is nonconstant of
degree zero and \(A\) has positive degree. Their fraction
field is \(k(A,T)\).

If \(g\) is irreducible, the
[tangent-conic regularity reduction](SEXTIC_KELLER_INDEPENDENT_LEADING_PENCIL_RESTRICTIONS.md#regularity-restrictions-in-the-tangent-conic-cases)
gives
\[
 p=xz+y^2,\qquad q=x^2,\qquad g=\mu p+\nu x^2+\rho xy ,
 \qquad \mu\ne0 .
\]
The nonpencil condition is \(\rho\ne0\). On \(x\ne0\)
put \(r=y/x\). Then \(p=x^2/T\) and
\[
                   A=x^6(\mu/T+\nu+\rho r)/T^2.
\]
Thus \(r\), and thereafter \(y,z\), can be recovered
rationally from \(A,T,x\). The full field is
\(k(A,T)(x)\), proving relative closedness.

If \(g\) is reducible and squarefree, it is a product
of two distinct lines. By the
[intersection partition restriction](SEXTIC_KELLER_INDEPENDENT_LEADING_PENCIL_RESTRICTIONS.md),
on a generic smooth pencil conic there is a point
of intersection multiplicity at least three with \(g\).
One of the two lines must be tangent there, and the
other must pass through that point. In fact the tangent
point is fixed in the pencil. To see this restrict the
pencil to the fixed tangent line. Each generic member
restricts to a binary square, so its linear space of
restrictions has dimension one: two independent
binary squares have linear combinations of rank two.
Consequently the common double root is fixed.

Place that basepoint at \([0:0:1]\), let \(x=0\) be the
common tangent, and choose the smooth member \(p=xz+y^2\)
in a linear frame. Subtracting a scalar multiple of \(p\)
from the other member removes its first-order jet at the
basepoint, making \(q\) binary in \(x,y\). The binary
pencil change preserves the Veronese factorization up to
a constant target change. We have \(g=x\ell\) with
\(\ell\) binary and independent of \(x\). Since \(g\)
is outside the pencil, \(q\) is not proportional to \(g\).
Put \(r=y/x\), \(g_0(r)=g(1,r)\), \(q_0(r)=q(1,r)\).
Over \(k(A,T)\) the generic fibre is given by
\[
                         x^6=\frac{AT^2}{g_0(r)q_0(r)^2}.
 \tag{1}
\]
In fact \(p=q/T\) and \(z=(p-y^2)/x\).

The gcd of six and the pole orders on the right of (1)
is one. At the finite root of \(g_0\), its order is
\(1+2m\) for \(m\in\{0,1,2\}\), where \(m\) is the order
of \(q_0\) there. The only order among these with
nontrivial gcd with six is three. In that case \(q_0\)
has a second distinct finite root giving order two,
since otherwise \(q\) would be proportional to \(x\ell\).
The degree-six equation is therefore irreducible even
over the algebraic closure of the generic constants.
Indeed, the field obtained by adjoining a root of
\(X^6-f\) to a field containing the sixth roots of unity
is cyclic of degree dividing six. A proper degree makes
\(f\) a square or a cube there in the base rational
function field, which the orders exclude. This proves
geometric integrality of the generic function field
and hence relative closedness of \(k(A,T)\).

## Homogeneous polynomial constants

In both nonpencil cases the rescaling of the source
scales \(A\) by the sixth power and fixes \(T\).
Any homogeneous element of \(k(A,T)\) of positive
ordinary degree less than twelve therefore has degree
six and form \(A f(T)\) with \(f\in k(T)\).
Indeed the action of the sixth roots of unity first
forces the degree to be a multiple of six. Dividing
by the corresponding power of \(A\) then gives a
rational function invariant under every rescaling
of \(A\), which belongs to \(k(T)\).

Suppose \(A f(T)\) is polynomial. At a finite value
\(T=\lambda\ne0\), there is a prime divisor of
\(q-\lambda p\) on which \(A\) is a unit. In the
irreducible case \(g\) is outside the pencil, so it
does not divide this member. In the reducible case
\(g\) is binary, while the member has a nonzero \(xz\)
term, so it has a nonbinary factor. In either case
the member is coprime to \(p\). A pole of \(f\) there
would give a pole in the source, which is impossible.

At \(T=0\), in the irreducible case \(q=x^2\) and \(A\)
does not vanish generically along \(x=0\). In the
squarefree binary case the same argument applies
if a factor of \(q\) does not divide \(g\).
If both factors of \(q\) belong to \(g\), then \(q\)
must be the square of just one of them, since \(q\)
is not proportional to \(g\). Along that factor
\(T\) vanishes to order two whereas \(A\) has
order one. This again excludes a pole of \(f\).
Finally at \(p=0\), \(T\) has a simple pole and \(A\)
vanishes to order two, because \(p\) is irreducible
and coprime to \(g,q\). Hence \(f\) is a polynomial
of degree at most two, and \(A f(T)\) is in
\(\langle A,gpq,gq^2\rangle\).
The conic-relation lemma now excludes a completion.

If \(g\) is irreducible and belongs to the pencil,
the same tangent-conic restriction gives \(q=x^2\),
and a pencil shear takes the smooth generator to
\(g\), with normalization \(p=xz+y^2\). If \(g\) is
squarefree reducible and belongs to the pencil, in
the preceding frame it is proportional to \(q=x\ell\).
Use a linear change preserving \(p=xz+y^2\) to
replace \(\ell\) by \(y\). Indeed a shift of \(y\)
by a multiple of \(x\) is compensated by a linear
shift of \(z\), and a rescaling normalizes the
coefficient. Finally a nonreduced quadratic is the
square of a line. The same intersection restriction
forces that line to be tangent at a fixed basepoint,
so the preceding tangent frame gives \(g=x^2\).
These are the displayed remaining tuples.

## The first correction when the common factor is the reduced member

**Theorem SexticKellerTangentConicReducedCommonFactorFirstCorrection.**
Let \(P=xz+y^2\) and \(Q=xy\). For a characteristic-zero
three-variable Keller map with highest tuple
\((P^2Q,PQ^2,Q^3)\), its degree-five terms can be absorbed
into a reference
\[
        (\widetilde g\,\widetilde p^{\,2},
          \widetilde g\,\widetilde p\,\widetilde q,
          \widetilde g\,\widetilde q^{\,2})
\]
with
\(\widetilde p=P+M\), \(\widetilde q=Q+\sigma-L\),
\(\widetilde g=Q+\sigma+2L\), for linear forms
\(M,\sigma,L\) over the algebraic closure.

### The first discriminant orders

Work over an algebraically closed field \(k\).
The pair \(P,Q\) has relatively closed field, since \(x\) is
a rational free variable and \(y,z\) can be solved when
\(x\ne0\). It is onto even at \(Q=0\), using \(y=0,x\ne0\)
and choosing \(z\). Polynomial contraction therefore gives
constant ring \(k[P,Q]\).

Use ordinary rescaling in \(\epsilon\), so the determinant
is \(\kappa\epsilon^{15}\). For the rescaled tuple write
\(D=F_1F_3-F_2^2\), with coefficients \(D_i\) of
\(\epsilon^i\) of ordinary degree \(12-i\).
The bracket identity \(\det J(F_1,F_2,D)=
\kappa\epsilon^{15} F_1\) shows that the first nonzero
\(D_i\) below order fifteen belongs to \(k[P,Q]\):
the leading pair has the same rational constant field
as \(P,Q\). In particular \(D_1=0\) by its odd ordinary
degree. Thus the degree-five error \(E\) is tangent to
the conic. Since \(P,Q\) are coprime, contraction with
\((Q^2,-2PQ,P^2)\) gives
\[
                     E=(2P U_0,\ Q U_0+P V,\ 2Q V)
\]
for cubic forms \(U_0,V\).

At order two \(D_2\in k[P,Q]\) and, on \(Q=0\), its
value is \(-P^2V^2\). As it has degree ten, this
value must be a scalar times \(P^5\). On the component
\(y=0\), \(P=xz\), whose cube is not a square in
\(k(x,z)\). Thus the scalar is zero, and \(V\) vanishes
on both components. Since \(Q\) is reduced, \(V/Q\)
is a linear form. Put \(\sigma=2V/(3Q)\), and begin
with reference
\[
                 (P^2 q,\ Pq^2,\ q^3),\qquad
                 q=Q+\epsilon\sigma .
\]
In the scaled tuple the residual order-one
coefficient is \((2P U,qU,0)\), with
\(U=U_0-P\sigma/2\), where \(Q\) may be replaced by
\(q\) in an order-one coefficient at the expense of
order-two terms. Use this convention, so the remainders
after that order-one correction have polynomial coefficients
\(A_j,B_j,C_j\) at orders \(j\ge2\).

For local calculations use \(P,q\) and one free coordinate
along each component of \(Q=0\). At \(x=0\) take that
coordinate to be \(z\), and at \(y=0\) take it to be \(x\).
These charts are regular at generic points because
\(dP,dQ\) are independent there. We can use \(q\)
as a formal coordinate since it starts with \(Q\).
In these charts,
\[
 D_2=Q(Q^2 A_2-2P Q B_2+P^2 C_2)-Q^2 U^2 .
\]
Here the subscript two on \(A_2,B_2,C_2\) denotes the
polynomial remainders before changing to the local charts.
In particular \(D_2/Q\) is a polynomial in \(P,Q\)
and is a scalar \(hP^4\) at \(Q=0\). It follows that
\(C_2|_{Q=0}=hP^2\).

If \(h\ne0\), then \(U|_{Q=0}=0\) on both components.
To see this, generically normalize with
\(q'=F_3^{1/3}\) and \(p'=F_2/(q')^2\), which have
leading terms \(Q,P\). By the bracket with \(F_2,F_3\),
\(D\) has coefficients in \(k(P,Q)\) in these normalized
parameters below order fifteen. The order-one changes
relative to \(q,P\) are zero for \(q'\) and \(U/Q\)
for \(p'\). The coefficient at order three after this
normalization is polynomial: in ordinary coordinates it
is \(D_3-\sigma\partial_Q D_2-(U/Q)\partial_P D_2\),
and every expression here is polynomial because
\(Q\mid D_2\). Its odd degree nine forces it to vanish.
On the \(q=0\) chart the raw order-three term of \(D\)
is \(2PU\cdot hP^2\), since every other product at that
order contains \(q\). On the other hand
\((U/Q)\partial_P D_2=4hP^3U\) there.
Their difference is zero only if \(U\) vanishes on
the component. This proves the claim when \(h\ne0\).

### Two leading functions on a fixed first-output level

Assume now \(h=0\), so \(Q\mid C_2\) as a source
polynomial. On either chart write
\[
 U_* = U|_{Q=0},\quad b_0=B_2|_{Q=0},\quad
 c_1=(C_2/Q)|_{Q=0},\quad
 d_0=C_3|_{Q=0}-\sigma|_{Q=0}c_1 .
\]
In particular \(c_1\) and \(d_0+\sigma c_1\) are restrictions
of the same polynomials on both components.

Fix the actual first component to \(\epsilon S\), where
\(S\) is an independent scalar, retaining \(P\) and the
free coordinate as the chart parameters. This is a regular
formal chart because the leading derivative of the first
component in \(Q\) is \(P^2\), a unit at the generic point.
Then
\[
 q=\epsilon(S/P^2-2U_*/P)+O(\epsilon^2).
\]
Write \(r=P^{-3}\), \(w=P U_*\). The coefficients at
orders two in the second output and three in the third
output are, respectively,
\[
 \mathcal B=r(S^2-3Sw+2w^2)+b_0,\qquad
 \mathcal C=r^2(S-2w)^3+(c_1/P^2)(S-2w)+d_0 .
 \tag{2}
\]
For the third output this uses \(Q=q-\epsilon\sigma\)
inside \(\epsilon^2 C_2\), the source of
\(-\sigma c_1\) in \(d_0\).

At fixed \(S\), the two functions in (2) have zero
two-variable bracket on the chart, by comparing the
order-five determinant coefficient with the Keller order
fifteen. Use the convenient form
\[
       \mathcal B=rS^2-3rwS+b,\qquad
       \mathcal C=r^2S^3-6r^2wS^2+cS+d .
 \tag{3}
\]
Here \(b,c,d\) are rational functions of the chart variables.
Their ordinary degrees are \(4,-2,3\) respectively,
with \(S\) of degree five, \(r\) of degree minus six,
and \(w\) of degree five.

Suppose first \(w\) depends on the free coordinate, so
\(r,w\) are independent. Use \(\partial_r,\partial_w\)
in the function field of the chart, a finite separable
extension of \(k(r,w)\). Let
\(\mathcal D=S\mathcal C-\mathcal B^2\), of \(S\)-degree
at most two by the two leading terms in (3). It too
has zero bracket with \(\mathcal B\). If the two partials
of \(\mathcal B\) are coprime as polynomials in \(S\) over
the chart field, then \(d\mathcal D=h_0\,d\mathcal B\)
with \(h_0\) independent of \(S\), by the degree bound.
Here these differentials hold \(S\) fixed. Exterior
differentiation and the coefficients at \(S^2,S\)
give \(dh_0\wedge dr=dh_0\wedge d(rw)=0\).
Thus \(h_0\) is constant, and homogeneity forces
\(h_0=0\), since the two homogeneous degrees are
eight and four. Now \(\mathcal D\in k[S]\) has degree
eight under scaling with \(S\) of degree five, so
\(\mathcal D=0\). This gives \(b=0\) by setting \(S=0\),
contrary to coprimeness of the partials of \(\mathcal B\),
which would now share \(S\).

Therefore the linear polynomial
\(\partial_w\mathcal B=-3rS+b_w\) has a common root
\(u=b_w/(3r)\) with \(\partial_r\mathcal B\). Equivalently
\[
                  b_r=3wu-u^2,\qquad b_w=3ru .
\]
Euler's formula \(-6r b_r+5w b_w=4b\) gives
\(b=3ru(2u-w)/4\). Differentiating in \(w\) yields
\((4u-w)u_w=5u\), so \(u(w-2u/3)^5\) is constant
in the free coordinate. In fact it belongs to the
scalar field \(k(y)\) on \(x=0\), or \(k(P)\) on \(y=0\).

This invariant is zero. Otherwise the monic equation
for \(u\) obtained from it puts \(u\) in the integrally
closed ring \(k(y)[z]\) on \(x=0\), or \(k(P)[x,x^{-1}]\)
on \(y=0\), because \(w\) belongs to the respective ring.
Then \(u\) and \(w-2u/3\) are both units. In the polynomial
case they are constants in the free coordinate, which is
impossible. In the Laurent case their exponents are
\(5n,-n\) for some integer \(n\), whereas the support
of \(w=P U(x,0,P/x)\) lies in \([-3,3]\).
For \(n\ne0\) the exponent \(5n\) cannot cancel with
\(-n\), and \(n=0\) again makes \(w\) constant.
This proves the claim.

We have either \(b=0\) or \(b=9rw^2/4\).
In the latter case use \(u=3w/2\) and \(\tau=S-u\).
Then
\[
   \mathcal B=r\tau^2,\qquad
   \mathcal C=r^2\tau^3-r^2u\tau^2+\tilde c\,\tau+\tilde d .
\]
The bracket equation is
\(\tau\,\partial_u\mathcal C+2r\,\partial_r\mathcal C=0\)
at fixed \(S\). Its three lower coefficients give
\(\tilde c_u=2r^2u\), \(\tilde d_r=0\), and
\(\tilde d_u=\tilde c-2r\tilde c_r\).
But \(\tilde c=r^2u^2+f\), with \(f_u=0\).
The last equation would give
\(\tilde d_u=-3r^2u^2+f-2r f_r\), impossible
because its \(r\)-derivative must vanish and \(r,u\)
are independent. Thus \(b=0\).

When \(b=0\), divide \(\mathcal B,\mathcal D\) by
\(S\). Now \(\mathcal B/S=r S-3rw\), whose two partials
are coprime. The same bracket and exterior differentiation
argument gives \(d(\mathcal D/S)=h_1\,d(\mathcal B/S)\)
with \(h_1\in k\). Their homogeneous degrees are three
and minus one, so \(h_1=0\) and then \(\mathcal D/S=0\)
by the degree of \(S\). Hence \(\mathcal C=\mathcal B^2/S\).
Comparing (2) with this identity gives
\[
                       c_1=-3(U_*/P)^2,\qquad
                       d_0=2(U_*/P)^3 .
 \tag{4}
\]
Since \(c_1\) is polynomial on its component, \(U_*/P\)
is polynomial there. It is homogeneous linear.

It remains to treat the case \(w\) is constant in the
free coordinate \(f\) of a branch. Then \(U_*\) is
polynomial in \(P\) on \(y=0\), so \(U_*=0\) by odd
degree, or a scalar multiple of \(y^3\) on \(x=0\).
In either case \(U_*/P\) is linear with no \(z\)-term.
In (3) the coefficients of \(S^3,S^2\) in the bracket
give \(c_f=2r b_f\) and \(d_f=-6rw b_f\). Thus
\(c=2rb+c_*\), \(d=-6rwb+d_*\), with constants
for the free derivative. The coefficient of \(S\)
then gives
\[
                b_f\bigl(18rw(rw)_r-2b-(c_*)_r\bigr)=0 .
\]
If \(b_f\ne0\), this makes \(b\) constant in \(f\),
a contradiction. Hence \(b,c,d\) are all constant
in the free coordinate. Formula (2) gives that
\(c_1\) is likewise constant there. Since \(c_1\)
is polynomial of ordinary degree two, it is a
scalar multiple of \(P\) on this component and
vanishes on the intersection line.

Finally \(U_*/P\) is linear on each branch. Their
restrictions to the line \(x=y=0\) have the form
\(\beta z,\delta z\). If both branches have
nonconstant \(w\), comparing \(c_1\) on the line
using (4) gives \(\beta^2=\delta^2\); comparing the
polynomial \(d_0+\sigma c_1\) then gives
\(\beta^3=\delta^3\). Thus \(\beta=\delta\).
If a branch has constant \(w\), its slope on the
intersection line is zero and its \(c_1\) vanishes
there, so (4) gives the same zero slope on the other
branch if needed. The linear forms consequently
come from a single linear \(L(x,y,z)\).
Thus \(U=P L+Q M\) for a linear \(M\), since
\(U-PL\) vanishes on both components of the reduced
divisor. Expanding the displayed \(\widetilde g,\widetilde
p,\widetilde q\) gives exactly the original degree-five
vector, including the shift \(\sigma\).

## Exclusion of the reduced common factor in the pencil

**Theorem SexticKellerTangentConicReducedCommonFactorExclusion.**
Let \(P=xz+y^2\) and \(Q=xy\) in characteristic zero.
There is no three-variable Keller map with ordinary highest
tuple \((P^2Q,PQ^2,Q^3)\).

Work over the algebraic closure. By the first-correction
theorem, a completion differs in ordinary degrees at most
four from
\[
 (\widetilde g\widetilde p^{\,2},
   \widetilde g\widetilde p\widetilde q,
   \widetilde g\widetilde q^{\,2}),
 \qquad
 \widetilde p=P+M,\quad
 \widetilde q=Q+\sigma-L,\quad
 \widetilde g=Q+\sigma+2L ,
\]
where the corrections are linear. Give \(x,y,z\) weights
\(1,1,w\) and increase \(w\) from one to two. For \(1<w<2\)
the reference leading tuple is
\[
                    (x^3 y z^2,\ x^3y^2z,\ x^3y^3).
 \tag{5}
\]
The only possible earlier changes from the degree-at-most-four
errors are a \(z^4\) term in the third component at \(w=3/2\),
and a \(z^4\) term in the second or \(z^3\) times a binary
linear form in the third at \(w=5/3\). These follow by
comparing their weights with \(4+2w,5+w,6\).

At \(w=3/2\) the first two entries of (5) fix \(xy\) and
\(z/y\) up to finitely many choices. The additional \(z^4\)
term in the third would vary on their generic fibres,
contradicting the leading determinant equation.
At \(w=5/3\) the tuple could have form
\[
       (x^3yz^2,\ x^3y^2z+\alpha z^4,\
                  x^3y^3+(\lambda x+\mu y)z^3).
\]
Its determinant is
\[
 2x^2z^5(2\alpha\lambda xz^3-6\alpha\mu yz^3
                    +3\lambda x^4y^2+(6\mu-12\alpha)x^3y^3).
\]
Thus \(\lambda=0\), \(\mu=2\alpha\), and \(\alpha\mu=0\),
so all three coefficients vanish.

### The transition at \(w=2\)

At weights \((1,1,2)\), write \(p=xz,\ q=xy+uz,\ g=xy+vz\)
for the reference leading parameters. The actual weighted
highest tuple is
\[
 \bigl(gp^2+a z^4,\ gpq+z^3(bx+cy),\
       gq^2+z^2(dx^2+e xy+f y^2)+h z^3\bigr)
 \tag{6}
\]
of weights \(8,7,6\). Use the leading zero determinant.
The coefficients at \(x^6y^3z^4,x^6yz^5,x^5y^4z^4,
x^5y^2z^5\), in order, give
\[
        e=2b,\qquad d(u-v)=0,\qquad f=2c-a,\qquad h=2bu .
\]
With these substitutions the coefficient at \(x^5 z^6\)
is \(4bd+2dv(u-v)\), so \(bd=0\). The coefficient at
\(x^2 y^3 z^6\) is \(-24(a-c)^2\), giving \(a=c\).
After that substitution the coefficients at \(x^4yz^6\)
and \(x^3y^2z^6\) are \(-8ad-12b^2\) and \(18au^2\),
respectively. Hence \(b=0,\ ad=0,\ au=0\).

We use the following weighted conic-relation obstruction.
Suppose the highest entries \(A,B,C\) in (6) satisfy
\(AC=B^2\), have relatively closed fraction field
\(k(A,B)\), and their polynomial invariants are precisely
\(k[A,B,C]\). Then no Keller completion exists.
Indeed the actual conic relation \(F_1F_3-F_2^2\)
has weight below fourteen. One can subtract from it
any polynomial in \(F_1,F_2,F_3\) of weight below
fourteen. Its derivative in the third target variable
has weight at most seven, since the only positive-weight
monomials allowed are the three individual components
and \(F_3^2,F_2F_3\). Thus the bracket with \(F_1,F_2\)
of the remainder has weight eight and leading term
\(\kappa A\) for \(\kappa=\det JF\). A highest term
of positive weight \(d<14\) in the remainder has
zero bracket with \(A,B\) at weight \(11+d\).
By the assumed constant ring it can be absorbed
by a polynomial in the actual components. Iterating
would make the remainder constant, a contradiction.

For the applications below the polynomial invariant
assertion follows from relative closedness and surjectivity
of the leading map onto \(AC=B^2\). On \(A\ne0\)
the cone has coordinate ring \(k[A,A^{-1},B/A]\);
a reduced denominator with a zero off the numerator
would lift to a point of the source, an impossibility
for a polynomial pullback. Likewise use the chart
\(C\ne0\). The intersection of the two rings is the
cone ring \(k[A,B,C]\). For example, express them
in \(A=s^2,B=st,C=t^2\). Their intersection belongs
to \(k[s,t]\) and is invariant under simultaneous
change of signs, so lies in the cone ring.

If \(a\ne0\), the equations for (6) give \(u=d=0\).
Then the tuple is \(G(z^2,zy,y^2)\) with
\(G=x^2(xy+vz)+az^2\). Put \(T=y/z\) and \(A=Gz^2\).
On a generic fibre of this pair,
\[
               T x^3+v x^2+a z=A/z^3.
\]
This cubic in \(x\) is irreducible over the rational
function field in \(z\) over the algebraic closure
of \(k(A,T)\). A rational root cannot have a pole at
infinity, since its cubic term would have pole order
at least three and dominate the other terms. A root
finite there cannot cancel \(az\). The generic field
is therefore geometrically integral and \(k(A,T)\)
is relatively closed. Every \(A\ne0\) and finite \(T\)
is realized by taking \(x=0\), \(az^4=A\), \(y=Tz\).
The ray \(A=B=0\) is realized on \(z=0\), where
\(C=x^3y^3\); the origin also occurs. The conic-relation
obstruction applies.

If \(a=0,\ u\ne v\), we have \(d=0\) and (6) is
\(g(p^2,pq,q^2)\). Put \(A=gp^2\) and \(T=q/p\).
On the generic fibre the equations give
\[
             y=(Tx-u)z/x,\qquad
             z^3=\frac{A}{x^2(Tx+v-u)}.
\]
The simple pole at \(x=(u-v)/T\) proves geometric
integrality, so \(k(A,T)\) is relatively closed.
Every \(A\ne0\) and finite \(T\) is realized by choosing
\(x\ne0\) with \(Tx+v-u\ne0\) and solving for \(z,y\).
The opposite ray is again surjective on \(z=0\),
and the origin occurs. The same conic-relation
obstruction applies.

If \(a=0,\ u=v\ne0\), (6) is
\((q p^2,q^2p,q^3+d p^2)\). Choose its last two
entries \(B=q^2p,\ C=q^3+d p^2\). Their parameter
Jacobian in \((p,q)\) is nonzero at generic \(p=0\),
where its value is \(3q^4\). At \(x=0\), \(p=0\)
and \(q=u z\); at \(z=0\), \(p=0\) and \(q=xy\).
The differentials of \(p,q\), and hence those of
\(B,C\), are independent at both generic points.
The rational constant field of the pair is \(k(p,q)\)
since the source field is \(k(p,q)(x)\).
A rational primitive for the leading derivation is,
up to scalar, \(1/(x\det\partial(B,C)/\partial(p,q))\),
using \(\det J(p,q,1/x)=1\).
It has a pole along \(x=0\) and is regular along
\(z=0\). Both divisors dominate \(p=0\) in the
constant-field plane, with order one of \(p\).
No correction by a constant of the derivation can
make the primitive regular at both. This contradicts
the divisor-regular
[leading-pair slice](KELLER_LEADING_PAIR_RATIONAL_SLICE.md).

### Exclusion of the final weighted frame

Only \(a=u=v=0\) remains. In particular
\(\widetilde q,\widetilde g\) are binary, because
the \(z\)-coefficients of both \(\sigma-L\) and
\(\sigma+2L\) vanish. Fix \(y,z\) at weights \(1,2\)
and raise the weight \(r\) of \(x\) starting at two.
The weighted reference leading tuple is (5), of
weights \(3r+5,3r+4,3r+3\).
For \(r\ge2\), the residual terms of degree at most
four and \(x\)-exponent at most three have weight
at most \(3r+2\). Thus the only possible new ties
come from \(x^4\), at \(r=3,4,5\) in the third,
second, first component, respectively. At each
first such tie the unchanged two entries of (5)
fix \(p=xz,q=xy\) up to finitely many choices.
Then \(x\) is still free, and the new \(x^4\)
term contradicts the leading determinant equation.
The map has therefore degree three in \(x\).

Write its highest \(x\)-coefficient as a binary vector
\(K(y,z)\). This has ordinary degree three with
highest part \(y(z^2,zy,y^2)\): it comes from the
product of the three reference coefficients of \(x\),
and a degree-at-most-four error contributes at
most a linear term to \(K\). The \(x^8\) determinant
coefficient gives \(\det(K_y,K_z,K)=0\).
The projective image is thus a curve, containing
the conic given by the highest binary part, since
that part is the limit along the radial direction.
Factor \(K\) on this conic. Its two binary parameters
have degree one and independent linear parts, by
degree three and the highest ratio \(z/y\).
Center them affinely, so \(K=(\ell+c_0)\nu(u_1,u_2)\)
in the two binary variables, with \(\ell\ne0\)
linear and \(\nu(u_1,u_2)=(u_1^2,u_1u_2,u_2^2)\).

Now apply the
[top-coefficient argument for a cubic transverse leading coefficient](SEXTIC_KELLER_NONBINARY_QUARTIC_LINE_FACTOR.md#the-highest-z-coefficient).
That argument from its \(A=(\ell+c)\nu\) frame
onward uses only the degree three in the transverse
variable and the displayed form of its top coefficient.
It first excludes \(c_0\ne0\) by the quadratic
lowest-coefficient descent. For \(c_0=0\), the
determinant coefficients force the constant terms
of the next two transverse coefficients to vanish.
Here the transverse variable is \(x\), so the
centered line \(u_1=u_2=0\) is collapsed, again
contradicting the Keller condition. This finishes
the last weighted case.

## The first correction for a square common factor with a transverse lower member

**Theorem SexticKellerTangentConicSquareCommonFactorTransverseFirstCorrection.**
Let \(P=xz+y^2\) and let \(Q\) be a homogeneous binary
quadratic in \(x,y\) with \(Q(0,y)\ne0\).
For a characteristic-zero three-variable Keller map
with ordinary highest tuple \(x^2(P^2,PQ,Q^2)\), its
degree-five terms can be absorbed by the reference
\[
         \bigl((xP+U_2)^2,\ (xP+U_2)(xQ+V_2),\
                           (xQ+V_2)^2\bigr)
\]
for homogeneous quadratic polynomials \(U_2,V_2\).

Work over the algebraic closure and write \(A=xP\),
\(B=xQ\). The pair \(A,B\) is onto. For generic
nonzero \(x\) write \(r=y/x\), so \(B=x^3Q(1,r)\)
and \(z\) is determined by \(A\). To prescribe nonzero
\(B\) choose \(r\) off the roots and rescale \(x\);
for \(B=0\) choose a finite root of \(Q(1,r)\).
Such a root exists since \(Q(0,y)\ne0\) gives a
nonzero \(r^2\)-coefficient. On the generic fibre
\(x^3=B/Q(1,r)\) is geometrically integral because
a root of the quadratic has multiplicity one or
two. Hence \(k(A,B)\) is relatively closed, and
the polynomial constant ring is \(k[A,B]\) by
surjectivity. In particular all its homogeneous
degrees are multiples of three.

Let \(\Delta=F_1F_3-F_2^2\). As in the conic-relation
argument, its highest coefficient of degree below
twelve belongs to \(k[A,B]\), by the bracket with
\(F_1,F_2\). Thus its terms of degrees eleven and
ten vanish, and its degree-nine term, if nonzero,
belongs to \(k[A,B]\) and is divisible by \(x^3\).
Put \(n=(Q^2,-2PQ,P^2)\). The degree-eleven equation
gives \(n\cdot E=0\) for the degree-five vector \(E\).
By coprimeness of \(P,Q\), write
\[
                        E=(2P U,\ QU+PV,\ 2QV)
\]
with cubics \(U,V\). If \(E_4\) is the degree-four
vector, the degree-ten equation is
\[
                 x^2 n\cdot E_4=(QU-PV)^2 .
\]
Put \(\lambda=Q(0,y)/y^2\ne0\), and
\(S=(QU-PV)/x\), which is a polynomial.
At \(x=0\), \(V=\lambda U\) and \(y^4\) divides
\(S^2\), so \(y^2\) divides \(S\).
Writing \(Q=\lambda y^2+\beta xy+\alpha x^2\)
shows that
\(S|_{x=0}\equiv(\beta y-\lambda z)U|_{x=0}
\bmod y^2\). Indeed the derivatives at \(x=0\)
of \(U,V\) contribute multiples of \(y^2\), and
\(P_x=z\). Thus \(y^2\mid U|_{x=0}\).
There are a binary linear \(L(y,z)\) and quadratics
\(U_2,V_2\) with
\[
                        U=P L+xU_2,\qquad V=Q L+xV_2 .
 \tag{7}
\]

Suppose \(L\ne0\). In degree nine the discriminant
on \(x=0\) is the cross term of \(E,E_4\).
By (7), \(E|_{x=0}=2L(P^2,PQ,Q^2)|_{x=0}\),
so that cross term equals \(2L S^2|_{x=0}\).
It vanishes since the degree-nine invariant is
divisible by \(x^3\). Hence \(S|_{x=0}=0\) and
\(V_2-\lambda U_2\) is divisible by \(x\).

Change target frame by replacing \(Q\) by
\(Q'=Q-\lambda P=x(\alpha x+\beta y-\lambda z)\)
inside the leading Veronese tuple. The determinant
relation for the conic is unchanged by this shear.
Put \(R=\beta y-\lambda z\), and use source
weights \((1/2,1,1)\) on \((x,y,z)\).
The weighted highest tuple in this frame is
\[
                 (x^2+2L)(y^4,\ y^2 xR,\ x^2R^2).
 \tag{8}
\]
Indeed \(U\) starts with \(y^2 L\) in these weights,
and \(V-\lambda U=Q'L+x(V_2-\lambda U_2)\)
starts with \(xRL\). The degree-four error can
compete only in the last entry, where its part
on \(x=0\) vanishes: in this target frame the
degree-ten identity \(n\cdot E_4=S^2\) restricts
there to \(y^4(E_4)_3=0\).

### The generic fibre of the weighted pair

For the first two entries of (8), put \(A\) equal to
the first and \(T=xR/y^2\) equal to their ratio. On a
generic level their relative volume differential must
be exact by the
[weighted leading-pair slice](KELLER_LEADING_PAIR_RATIONAL_SLICE.md).
At fixed \(T\) solve \(z=(\beta y-Ty^2/x)/\lambda\).
If \(s=xy^2\), the equation at fixed \(A\), after
substitution, is
\[
                    s^3+(a y^5-A)s+b y^8=0
 \tag{9}
\]
for constants \(a,b\) over the algebraic closure of the
generic value field. The coefficient \(b\) is a nonzero
scalar times \(T\) times the \(z\)-coefficient of \(L\).
If \(b=0\), then \(a\ne0\) because \(L\ne0\).
The relative volume after fixing \(T\) is proportional
to \(y^2 dx\wedge dy/x=y^2 ds\wedge dy/s\).
Dividing by the differential of \(A\) gives the
relative one-form
\[
                         \frac{y^2\,dy}{3s^2+a y^5-A}.
 \tag{10}
\]

If \(b\ne0\), the locus \(s=0\) in (9) is the point
\(s=y=0\), smooth with denominator \(-A\) in (10).
Where \(s\ne0\) the affine generic level is smooth.
At infinity \(y\) must have a pole since (9) is monic
in \(s\). The two dominant terms are \(s^3\) and
\(b y^8\), so \(s\) has pole order \(8/3\) times
that of \(y\). Indeed any larger order would make
\(s^3\) dominate, and any smaller order would make
\(b y^8\) dominate. The denominator of (10) then
has pole order \(16/3\) times that of \(y\).
Thus (10) is regular at all points of the complete
normalization. A nonzero regular differential cannot
be exact there.

If \(b=0\), the fibre on \(s\ne0\) is
\(s^2+a y^5=A\) and (10) is proportional to
\(y^2dy/s^2\). At a root of \(a y^5=A\) the
parameter \(s\) gives a simple pole with nonzero
residue, since the root is nonzero and simple.
This again prevents exactness. Therefore \(L=0\),
and (7) gives \(U=xU_2,V=xV_2\).
Substitution in the tangent vector is precisely the
degree-five correction of the stated Veronese reference.

## Two correction orders at the double common divisor

**Theorem SexticKellerTangentConicSquareCommonFactorTransverseJetReduction.**
Under the hypotheses of the preceding theorem, write
\[
\lambda=Q(0,y)/y^2,\quad q=Q-\lambda P,\quad
\widehat F=(F_1,F_2-\lambda F_1,F_3-2\lambda F_2+\lambda^2 F_1).
\]
In the Veronese reference for \(\widehat F\) write the
degree-three parameters as \(xP+U_2,xq+V_2'\), where
\(V_2'=V_2-\lambda U_2\). Then \(x\mid V_2'\).
The homogeneous degree-four error relative to this
reference is of the form
\[
                      (2P U,\ qU+P V,\ 2q V),
       \qquad \deg U=\deg V=2,\qquad x\mid V.
\]
The third component of \(\widehat F\) has its degree-six,
degree-five, degree-four and degree-three parts divisible
by \(x^4,x^3,x^2,x\), respectively. In the variables \(y,z\),
the degrees of the three components of \(\widehat F\) are at
most \(4,3,2\), with coefficients in \(k[x]\).

### The first transverse coefficient

Write \(q=x(R+\alpha x)\) with \(R=\beta y-\lambda z\).
The discriminant \(D=\widehat F_1\widehat F_3-\widehat F_2^2\)
has terms of degrees eleven and ten equal to zero, by the
preceding proof. The change of target frame preserves \(D\).
The degree-ten equality contracts the degree-four error of the
Veronese reference against \(x^2(q^2,-2Pq,P^2)\).
Since \(P,q\) are coprime, every vector in this kernel
has the displayed form: the first entry is divisible by \(P\),
the third by \(q\), and the equation gives the middle entry.

Use weights \((1/2,1,1)\). The weighted highest tuple is
\((a^2,ab,b^2)\), with
\[
                     a=xy^2,\qquad b=x^2R+L(y,z),
                 \qquad L=V_2'|_{x=0}.
\]
The leading parameter weights are \(5/2,2\). Terms of the
degree-four error cannot tie the first or second components.
They cannot tie the third either, since its degree-four
error \(2qV\) is divisible by \(x\).
The pair \(ab,b^2\) is independent, since \(a,b\) are
independent. Its parameter Jacobian is \(2b^2\).
Consequently a rational slice for \(ab,b^2\) gives a
nonzero constant multiple of a rational primitive of the
relative differential of \((a,b)\) on each geometric generic
fibre.

Write \(L=u z^2+v yz+\gamma y^2\). On levels \(a=A,b=B\)
solve \(x=A/y^2\). The equation for \(z\) is
\[
 u z^2+(v y-\lambda A^2/y^4)z+
                \gamma y^2+\beta A^2/y^3-B=0 .
\]
The relative differential is, up to sign,
\[
 \frac{dy}{y^2(2u z+v y-\lambda A^2/y^4)}
                    =\frac{y^2\,dy}{W},\qquad
 W=y^4(2u z+v y-\lambda A^2/y^4).
\]
If \(u\ne0\), this is the hyperelliptic model
\[
 W^2=(v^2-4u\gamma)y^{10}+4uB\,y^8
       -(2v\lambda+4u\beta)A^2y^5+\lambda^2 A^4 .
\]
The right side is squarefree over the algebraic closure of
the generic level field. Its constant term is nonzero, and
away from zero a multiple root would make the transcendental
value \(B\) a critical value of the nonconstant rational function
obtained by solving this equation at \(W=0\) for \(B\).
Its degree is eight or ten. The differential is regular at
finite branch points by taking \(W\) as parameter, and at
infinity its possible pole is cancelled by the pole order
of \(W\), at least four. This nonzero regular differential
on a complete smooth curve has no rational primitive.

If \(u=0\) and \(v\ne0\), the curve is rational in \(y\)
and the denominator \(W=v y^5-\lambda A^2\) has simple
nonzero roots. The differential has nonzero residues
there, which also preclude a rational primitive.

It remains to test \(L=\gamma y^2\). If \(\gamma\ne0\),
the differentials of \(a,b\) are independent at the
generic point of \(x=0\), where \(a=0,b=\gamma y^2\):
\(da\) has a nonzero normal coefficient, and \(db\) a
nonzero tangential coefficient. The same is true of
\(ab,b^2\) since \(2b^2\) is a unit there. By the
[divisor-regular leading-pair slice](KELLER_LEADING_PAIR_RATIONAL_SLICE.md)
there must be a rational function \(\phi\), regular at
that divisor, with \(\det J(ab,b^2,\phi)=1\) after
rescaling. On the full field one such primitive is
\[
                            \phi_0=\frac{y^3}{6\lambda a^2 b^2}.
\]
Indeed \(\det J(a,b,y)=\lambda a^2/y^2\).
The source field is \(k(a,b)(y)\), since \(x,z\) are
recovered from \(a,b,y\), so every other primitive differs
from \(\phi_0\) by an element of \(k(a,b)\).
Along \(x=0\), the valuation restricted to \(k(a,b)\)
is the order along \(a=0\): \(a\) has order one and
the residue \(b=\gamma y^2\) is transcendental over
\(k\). Any correction cancelling the order-two pole must
therefore have its product with \(a^2\) reduce into
\(k(\gamma y^2)\). The residue of \(a^2\phi_0\) is
\(y^3/(6\lambda b^2)\), which is outside this field.
No correction makes the slice regular. Thus \(\gamma=0\)
as well, proving \(x\mid V_2'\).

### The next transverse coefficient

Lower the weight of \(x\) from \(1/2\) to \(1/3\) with
\(y,z\) fixed in weight one. Since \(V_2'\) is divisible
by \(x\), the highest reference parameters up to the first
tie are \(a=xy^2,\ b=x^2R\). The only possible ties are
the \(P V\) contribution in the second degree-four error
and \(2qV\) in the third, both at weight \(1/3\).
Indeed the reference weights of the three outputs are
\(2w+4,3w+3,4w+2\) for the moving weight \(w\). The
degrees at most three are below these weights, and the
first degree-four error has weight at most four. In the
second error \(qU\) has weight at most \(w+3\) and
\(PV\) at most four. The third error has weight at most
\(w+3\).

Set \(V_0=V|_{x=0}\). The leading tuple at the tie is
\[
                (a^2,\ a(b+\eta),\ (b+\eta)^2-\eta^2),
                            \qquad \eta=V_0(y,z)/x.
\]
Its determinant is zero since the expected determinant
weight \(8w+7\) is positive. The tuple's field has
transcendence degree at most two, and adjoining \(a\)
and \(\eta\) is algebraic over it. Thus \(\eta\) is
algebraic over \(k(a,b)\). The latter field is relatively
closed, as its source field is \(k(a,b)(y)\), with
\[
                       x=a/y^2,\qquad
                       z=(\beta y-b y^4/a^2)/\lambda .
\]
If \(V_0\ne0\), substitution into \(y^2 V_0(y,z)/a\)
has positive degree in \(y\). The leading degrees
from \(z^2,yz,y^2\) are respectively ten, seven and
four, so cancellation of the highest present one is
impossible. Therefore \(V_0=0\) and \(x\mid V\).

Let \(E_3\) be the degree-three error. The degree-nine
part of \(D\) is
\[
 x^2(q^2,-2Pq,P^2)\cdot E_3
             -2x(qU_2-PV_2')(qU-PV).
\]
The second summand is divisible by \(x^3\). The degree-nine
part of \(D\) is also divisible by \(x^3\), by the
leading constant-ring argument in the preceding theorem.
Reducing after division by \(x^2\) gives
\(y^4(E_3)_3|_{x=0}=0\). This yields the degree-three
divisibility. The degree-four part of the third output is
\((V_2')^2+2qV\), the degree-five part is \(2xqV_2'\),
and the degree-six part is \((xq)^2\), giving the
other stated divisibilities. For the degree bounds in
\(y,z\), the reference parameters have degrees at most
two and one, respectively. In the quartic error the
second entry \(qU+PV\) has degree at most three in
\(y,z\), since \(q,V\) are divisible by \(x\). The third
entry has degree at most two by the divisibilities just
proved, including that of the degree-three part. The
correction pair is unique from the
degree-five vector, so the divisibilities hold over the
original field as well.

## The highest tuple in the transverse variables

**Theorem SexticKellerTangentConicSquareCommonFactorTransverseNormalForm.**
Under the hypotheses of the preceding theorem, make the same target
shear and use \(R=\beta y-\lambda z\) as the third linear source
coordinate, denoted here by \(z\). Let \(A,B,C\) be the highest
parts of the sheared components in the variables \(y,z\) over \(k[x]\).
There is a translation \(X=x+\alpha\) and there are scalars
\(\mu,\delta\) such that
\[
 (A,B,C)=(S^2,ST,T^2+\delta S),\qquad
                  S=Xy^2,\qquad T=X^2z+\mu Xy .
\]
Moreover both quartic tangent corrections \(U,V\) in the preceding
theorem are divisible by \(x\) in the untranslated source frame.
Hence the degree-five and degree-four terms are together absorbed
by a Veronese reference with parameter degrees at most three.

### The transverse discriminant

Write \(R=k[x,y,z]\) for the source ring in this proof.
The weight \((0,1,1)\) gives highest determinant weight
\(4+3+2-2=7\), so \((A,B,C)\) has zero determinant.
The rational slice proof works for this weight: rescale \(y,z\)
by \(\epsilon^{-1}\) and leave \(x\) fixed. The rescaled
polynomials are in \(R[\epsilon]\), and the determinant has
order seven. The same formal rectification of an independent
leading pair therefore yields the rational slice.

The highest powers of \(x\) in \(A,B,C\) are \(x^2y^4,x^3y^2z,x^4z^2\).
The polynomial \(D_\perp=AC-B^2\) is the part of \(y,z\)-degree six
of the ordinary discriminant, whose ordinary degree is at most nine.
Consequently \(\deg_x D_\perp\le3\). On \(y=1\) write \(t=z\) and
use lower-case \(a,b,c,d\) for \(A,B,C,D_\perp\) in this chart.
The two weight-zero fractions \(C^2/A,D_\perp^2/A^3\) are
algebraically dependent: otherwise adjoining the positive-weight
element \(A\) would make their field have transcendence degree
at least three inside \(k(A,B,C)\).
In \(k(t)((1/x))\) take the square root of \(a\) with leading term
\(x\), and put
\[
                         \rho=c^2/a,\qquad h=d/a^{3/2}.
\]
The series \(h\) is algebraic over \(k(\rho)\), regular at \(x=\infty\),
and all its \(1/x\)-coefficients lie in \(k[t]\), since \(a\) is
monic in \(x\). The leading term of \(\rho\) is \(x^6 t^4\).
The residue field of the restricted valuation on \(k(\rho,h)\)
is algebraic over \(k\), by finiteness over \(k(\rho)\).
Work over the algebraic closure for the argument. The residue of
\(h\) is then a scalar \(\delta\).
If \(h-\delta\) has positive order \(n\) and leading coefficient
\(p(t)\ne0\), the element \(\rho^n(h-\delta)^6\) has nonzero
residue \(t^{4n}p(t)^6\), which is not scalar. Thus \(h=\delta\).
After restoring homogeneity, either \(D_\perp=0\), or
\(D_\perp^2=\delta^2 A^3\) with \(\delta\ne0\).

### The nonsquare case

If \(A\) is not a square, the latter identity is impossible by
unique factorization, so \(AC=B^2\). Factor the tuple as
\(g(s^2,st,t^2)\) with \(s,t\) coprime and \(g,s,t\) homogeneous
in \(y,z\). Here the letters \(s,t\) denote polynomial parameters.
The \(x\)-degree of \(s\) is zero or one. If it is one, \(g\)
is independent of \(x\) and its product with the square of the
leading \(x\)-coefficient of \(s\) is \(y^4\). Thus \(g\) is a
square, contrary to the hypothesis. If the \(x\)-degree of \(s\)
is zero, \(s^2\) divides \(y^4\), by the leading coefficient of
\(A\). The leading \(x\)-term of \(t/s=B/A\) is \(xz/y^2\),
so \(y^2\) divides \(s\). After a scalar normalization the tuple
therefore has the form
\[
 g(x)\bigl(y^4,\;y^2(xz+\ell(y,z)),\;(xz+\ell(y,z))^2\bigr)
\]
with \(g\) a squarefree monic quadratic and \(\ell\) linear.
Shift \(x\) by the coefficient of \(z\) in \(\ell\) to a
coordinate \(X\). On fixed generic \(A\) and \(B/A\), eliminate
\(z\) using \(B/A=(Xz+c y)/y^2\). The resulting curve is
\(g_*(X)y^4=A\), and its relative volume differential is
proportional to \(y^3\,dX/X\). Indeed fixing \(B/A\) gives the
two-form \(y^2\,dX\wedge dy/X\), and the derivative in \(y\)
of \(g_*(X)y^4\) is \(4A/y\) on the level.

If \(g_*(0)\ne0\), this differential has a simple pole with
nonzero residue over \(X=0\). Otherwise \(g_*(X)=X(X-d)\)
for \(d\ne0\), and the fourth-root cover is irreducible by
the simple zero at \(d\). The differential has a single pole
of order four over \(X=0\). At \(X=d\), \(dX\) has order
three and \(y^3\) has pole order three. At infinity \(X\)
has pole order two and \(y\) has a simple zero. These orders
follow from the defining fourth-power relation at the two
simple finite roots and at the pole of order two of \(g_*\).

If a primitive exists, project it onto the character
\(y\mapsto\zeta y\) with eigenvalue \(\zeta^3\), where \(\zeta\)
is a primitive fourth root of unity. Its form is \(p(X)y^3\)
with \(p(X)\) rational. It has at most a pole of order three
over \(X=0\) and is regular elsewhere. Thus \(p\) has no finite
pole, must vanish at \(X=d\), and has degree at most one by
the orders at infinity. But
\[
        d((X-d)y^3)=y^3\left(-\tfrac12+\frac{3d}{4X}\right)dX
\]
is not proportional to the required differential. This excludes
the nonsquare case by the leading-pair rational slice.

### The square case and the remaining common factor

Choose the square root \(S\) of \(A\) with leading \(x\)-term
\(xy^2\). The discriminant identity gives \(D_\perp=\delta S^3\)
with the same sign choice as in the formal calculation.
Then \(S^2\mid B^2\), so \(B=ST\) and \(C=T^2+\delta S\)
for a polynomial \(T\). Comparison of \(x\)-coefficients in \(A\)
shows \(S=xy^2+U_2(0,y,z)\). The constant term in this
coefficient comparison gives \(2y^2 U(0,y,z)=0\). Thus \(x\mid U\),
as claimed. Together with \(x\mid V\) this absorbs the quartic
error by adding \(U/x,V/x\) to the two reference parameters.

On the chart \(y\ne0\) write \(s=z/y\), and set
\[
       S=y^2(x+u(s)),\qquad T=y(s L(x)+M(x)),
\]
where \(u(s)=u_2s^2+u_1s+u_0\), \(L\) is monic quadratic and
\(\deg M\le1\). These bounds follow from the weights and the
highest \(x\)-term \(x^2z\) of \(T\).
The pair \(S^2,ST\) is independent and has nonzero parameter
Jacobian generically. On fixed \(S,T\), put \(\xi=T^2/S\)
and eliminate \(y\). With \(b=sL+M\), the level equation is
\(f=b^2-\xi(x+u(s))=0\), and its differential is proportional to
\(dx/f_s\). In fact the volume is \(y\,dx\wedge dy\wedge ds\),
fixing \(T\) gives \(T\,dx\wedge ds/b^2\), and on the
level \(dS=-S\,df/b^2\) at fixed \(T\).

The quadratic discriminant in \(s\) is
\[
 4\xi\bigl(L^2(x+u_0)+u_2M^2-u_1 LM\bigr)
                   +\xi^2\bigl(u_1^2-4u_2(x+u_0)\bigr).
\]
If \(u_1,u_2\) do not both vanish, it is squarefree of degree
five over the geometric generic field. Indeed the two
polynomials of \(x\) multiplying \(\xi,\xi^2\) have degrees
five and at most one, with the second nonzero. Common zeros
are simple for generic \(\xi\); away from them a multiple
root would make \(\xi\) a critical value of their ratio.
Then \(dx/f_s\) is \(dx\) divided by a square root of the
discriminant. It is nonzero and regular at the branch points
and at infinity, so cannot be exact. Therefore \(u\) is constant.

Put \(X=x+u_0\). On the resulting curve take \(r=\sqrt{X}\):
the differential is a nonzero constant times
\(dr/L(r^2-u_0)\). If \(L\) has a simple root different
from \(-u_0\), there is a simple pole with nonzero residue.
If \(L\) has a double root other than \(-u_0\), the denominator
is \((r-\tau)^2(r+\tau)^2\) with \(\tau\ne0\). Its reciprocal
has nonzero residue at \(r=\tau\), since the derivative
of \((r+\tau)^{-2}\) there is nonzero. Thus \(L=X^2\),
leaving \(S=Xy^2\) and \(T=X^2z+(\mu X+\nu)y\).

Finally suppose \(\nu\ne0\). The pair \(S,T\) is coprime and
onto, and its source field is \(k(S,T)(y)\): one solves
\(X=S/y^2\) and then \(z\) on the generic chart. The polynomial
derivation kernel is \(k[S,T]\) by surjectivity and relative
closedness. The three parameter minors of
\((S^2,ST,T^2+\delta S)\) are, up to sign,
\[
                        2T^2-\delta S,\quad 4ST,\quad 2S^2,
\]
and have gcd one. Apply the
[weighted saturation descent](SEXTIC_KELLER_MIXED_DEGREE_LEADING_REDUCTIONS.md#descending-through-the-lower-terms)
with weights \((0,1,1)\). Here the parameter weights are two
and one and the determinant order seven exceeds every output
weight. For clarity the graded-basis step still works with
a zero source weight. Lift a homogeneous \(k\)-basis of
\(R/(S,T)\), including \(1\). The lifts generate \(R\) over
\(k[S,T]\) by induction on weight since both parameter weights
are positive. They are independent because a homogeneous
relation in the \(I=(S,T)\)-adic graded ring has independent
monomials in \(S,T\) over \(R/I\): if
\(\sum_{i=0}^n S^i T^{n-i}r_i\in I^{n+1}\), first absorb
the right side into the coefficients modulo \(I\). In the
resulting zero relation \(S\mid r_0\) by coprimality.
Transfer this multiple to the next coefficient, divide by
\(S\), and induct to obtain \(r_i\in I\) for all \(i\).
Taking minimal total parameter degree in any finite relation
gives independence of the lifts. In the determinant descent
one uses all lifts except \(1\), including any of weight zero,
for the tangent coefficient calculation. This absorbs every
error through loss four, since weights in the source ring are
nonnegative and the determinant order remains larger. The
result would factor through two polynomials, a contradiction.
Thus \(\nu=0\).

The square root, the shift and the remaining scalar coefficients
are determined by the displayed polynomial comparisons. They
therefore belong to the original coefficient field.

## Exclusion of the transverse discriminant scalar

**Theorem SexticKellerTangentConicSquareCommonFactorTransverseVeronese.**
Under the hypotheses of the transverse normal form, its scalar
\(\delta\) is zero. Thus the highest tuple in \(y,z\) over \(k[X]\)
is the Veronese tuple on \(Xy^2,X^2z+\mu Xy\).

Use the translated coordinate, now called \(x\), and put
\[
                 a=xy^2,\qquad b=x^2z+\mu xy.
\]
The source field is \(k(a,b)(y)\) with
\[
                  x=a/y^2,\qquad z=b y^4/a^2-\mu y^3/a .
\]
In the source ring localized at \(a\), the intersection with
\(k(a,b)\) is exactly \(k[a,a^{-1},b]\). This follows either
from these formulas or directly from the resulting ring
identity
\[
              k[x,y,z,a^{-1}]=k[a,a^{-1},b,y,y^{-1}].
\]
The coefficient of \(y^0\) in these Laurent coordinates
also gives a \(k[a,b]\)-linear projection from the
unlocalized source ring to \(k[a,b]\). Indeed the
coefficient of a term with \(y\)-exponent zero coming
from \(x^i y^j z^k\) has power \(i-2k+n=(j+n)/2\) of \(a\),
where \(n\) counts the uses of \(-\mu y^3/a\) in expanding
\(z^k\). This power is nonnegative. The projection preserves
the transverse grading, since \(a,b,y\) have weights \(2,1,1\).
In particular the unlocalized polynomial kernel for
\(\det J(a,b,\cdot)\) is \(k[a,b]\).

Suppose \(\delta\ne0\), and write \(f=(a^2,ab,b^2+\delta a)\)
for the leading map in parameter coordinates. Its normal vector is
\((2b^2-\delta a,-4ab,2a^2)\). The gcd of its entries in the
source ring is \(x\) up to a nonzero scalar: \(y\) cannot divide
the first, and the first has \(x\)-order one. At a first error of
loss \(r\) relative to a two-parameter reference, the linear
determinant equation puts the normal contraction in \(k[a,b]\).
Use the projection to absorb the scalar coefficient vector into
the reference, leaving a tangent error \(f_a U+f_b V\).
Solving by the minors shows \(xU,xV\) are polynomials, with
weights \(2-r,1-r\), respectively. In particular \(V=0\)
when \(r\ge2\). The third error coefficient is \(\delta U+2bV\),
so \(U\) is polynomial even for \(r=1\), since \(b\) is
divisible by \(x\). Every error of loss at least two therefore
admits polynomial absorption. It suffices to handle the
simple pole of \(V\) at loss one.

Use the rescaling in \(\epsilon\) of transverse weights \((0,1,1)\).
In parameter coordinates the scalar correction at loss one is
\(G_\epsilon=f+\epsilon K\), where the three entries of \(K\)
have weights \(3,2,1\) in \(k[a,b]\).
Write \(V=V_++c/x\) with \(V_+\in k[x]\) and \(c\in k\).
After absorbing the polynomial tangent terms, the reference
uses \(\widetilde a=a+\epsilon U,\
\widetilde b=b+\epsilon V_+\). The first remaining error
relative to it is \(f_b c/x\). Write \(E^{(2)}\) for the
polynomial coefficient of the error at order two.

Solve the first two reference output equations exactly in
\(k[x,y,z,a^{-1}][[\epsilon]]\), using the leading parameter
Jacobian \(2a^2\). Denote their parameters by \(a',b'\).
With \(w=c/x\), their expansions in the fraction field are
\[
 a'=\widetilde a+\epsilon^2 P_2+O(\epsilon^3),\qquad
 b'=\widetilde b+\epsilon w+\epsilon^2 Q_2+O(\epsilon^3).
\]
Let \(H=F_{3,\epsilon}-G_{3,\epsilon}(a',b')\). Its order
is at least two, and \(\det J(F_{1,\epsilon},F_{2,\epsilon},H)\)
is the Keller determinant of order seven: the removed
output depends on \(a',b'\), just as the first two do.
The coefficient \(H_2\) therefore belongs to \(k(a,b)\).
It also belongs to the localization at \(a\), and has
transverse weight zero. The ring identity gives
\[
                            H_2\in k[b^2/a].
\]
It is regular at \(x=0\), since \(b^2/a\) has positive order
there.

On the other hand, the first two output comparisons give
\[
\begin{aligned}
 2aP_2&=E^{(2)}_1-K_{1,b}w,\\
 aQ_2&=E^{(2)}_2-bP_2-(U+K_{2,b})w.
\end{aligned}
\]
The first coefficient \(K_1\) is a linear combination of
\(ab,b^3\), so \(K_{1,b}\) is divisible by \(x\). It follows
that \(P_2\) has at most a simple pole and \(Q_2\) at most
a double pole along \(x=0\). The third comparison is
\[
          H_2=E^{(2)}_3-\delta P_2-2bQ_2
                            -(2V_++K_{3,b})w-w^2 .
\]
All terms except \(w^2\) have at most simple poles, since
\(x\mid b\). Regularity of \(H_2\) forces \(c=0\).
This absorbs the loss-one correction as well. The determinant
order seven allows the same tangent descent at each later loss
through four. At that point the tuple factors through two
polynomials, contradicting the Keller determinant. Hence
\(\delta=0\).

## Exclusion when the linear term of the transverse parameter vanishes

**Theorem SexticKellerTangentConicSquareCommonFactorTransverseNonzeroSlope.**
In the transverse Veronese conclusion, \(\mu\ne0\).
Thus a Keller completion in this doubled common-factor branch must
have transverse highest tuple \((S^2,ST,T^2)\) with
\(S=Xy^2,\ T=X^2z+\mu Xy\) and a nonzero linear slope \(\mu\).

Suppose \(\mu=0\), write \(x\) for the translated coordinate, and
let \(a=xy^2,\ b=x^2z\). Work in transverse weights \((0,1,1)\).
The projection and derivation kernel of the preceding proof apply.
In addition,
\[
                   k[x,y,z,b^{-1}]\cap k(a,b)=k[a,b,b^{-1}].
 \tag{11}
\]
Over an algebraic closure the pair maps the chart \(b\ne0\) onto
the parameter chart \(b\ne0\), by taking \(x=1\) and solving
\(y^2=a,z=b\). A reduced denominator with a zero off the numerator
would contradict this surjectivity. This proves (11), and the
coefficient-field assertion also follows by descent.

At each error loss \(r\), project the error onto \(k[a,b]\) as
before and absorb that scalar vector into the reference. The
remaining error is tangent to \(f=(a^2,ab,b^2)\). Its tangent
coefficients \(U,V\) can have poles only along \(x=0\), since the
common divisor of the three source minors is \(x^2\).
In fact they have at most simple poles: \(U\) is the first
error divided by \(2a\), and \(aV\) is the second error minus
\(bU\), which is polynomial since \(x^2\mid b\).
Their weights are \(2-r,1-r\). Only \(U,V\) at loss one
and \(U\) at loss two can therefore have poles.

### Loss one

Write the pole parts as \(u=l(y,z)/x,\ v=c/x\), where \(l\) is
linear and \(c\) scalar. After absorbing the polynomial parts
\(U_1,V_1\), put
\(\widetilde a=a+\epsilon U_1,\ \widetilde b=b+\epsilon V_1\)
and write \(G_\epsilon=f+\epsilon K\) for the reference in
parameters. The three weights of \(K\) are \(3,2,1\).
Solving two output equations identifies the first parameter
perturbation relative to \(\widetilde a,\widetilde b\) as
\((u,v)\). The remaining third-output error starts at order
two, and its order-two coefficient is killed by the leading
derivation. Thus it is in \(k(a,b)\).

First solve \((F_1,F_2)\), whose parameter Jacobian is \(2a^2\).
The normalized coefficient at order two in \(F_3\) is in
\(k[a,a^{-1},b]\) of weight zero, hence regular at \(x=0\)
as it belongs to \(k[b^2/a]\).
Its pole of order two could only come from
\[
                     -\left(v-(b/a)u\right)^2 .
\]
Indeed the linear projection to the normal component is
\((b^2/a^2,-2b/a,1)\), regular at \(x=0\). It annihilates the
order-two tangent parameters. Applying it to the quadratic
term \(-(u^2,uv,v^2)\) gives the displayed term. Products with
the polynomial first corrections and the first derivatives of
\(K\) have at most simple poles. Since \(b/a=xz/y^2\),
the order-two residue is \(-c^2\), so \(c=0\).

Now solve \((F_2,F_3)\), whose parameter Jacobian is \(2b^2\).
By (11) the normalized coefficient at order two in \(F_1\)
has weight two and belongs to \(k[a,b,b^{-1}]\).
The monomials of weight two there are \(a^m b^{2-2m}\)
for \(m\ge0\), of distinct \(x\)-orders \(4-3m\).
The normal projection is now \((1,-2a/b,a^2/b^2)\).
With \(v=0\), the order-two coefficient has pole order at
most two. Its terms of pole order two have residues
\[
 -l^2+2\kappa \frac{y^2}{z}l+e\,\frac{y^4}{z^2}
\]
for scalars \(\kappa,e\). Here the quadratic term contributes
\(-l^2\), the first correction to the second output
contributes \(2(V_1(0)+K_{2,a})y^2l/z\), and the polynomial
order-two error in the third output contributes \(e y^4/z^2\).
There are no other double-pole terms, since \(K_3\) depends
only on \(b\) and the mixed quadratic term in the third
output is zero when \(v=0\). By the listed monomial orders
the residue must be a scalar times \(y^4/z^2\), the residue
from \(m=2\). Consequently \(zl/y^2\) satisfies a monic
quadratic over \(k\). It is a polynomial in \(z/y\) with
zero constant term, so it must vanish. This proves \(l=0\)
and absorbs the first loss.

### Loss two

Project the scalar error at loss two into the reference
and absorb the polynomial part of its first tangent
coefficient. Write the resulting reference as
\(G_\epsilon(\widetilde a,\widetilde b)\) with
\(\widetilde a=a+\epsilon U_1+\epsilon^2 U_2\) and
\(\widetilde b=b+\epsilon V_1\), where \(V_1,U_2\in k[x]\).
The unabsorbed tangent coefficient is \(c/x\) in the
\(a\)-direction. By the parameter weights the reference
polynomials have the forms
\[
\begin{aligned}
 G_1(A,B)&=A^2+\epsilon(k_1AB+k_2 B^3)
                         +\epsilon^2(k_3 A+k_4 B^2),\\
 G_2(A,B)&=AB+\epsilon(\alpha A+k_5 B^2)+\epsilon^2 k_6 B,\\
 G_3(A,B)&=B^2+\epsilon k_7 B+\epsilon^2 k_8 .
\end{aligned}
\]
Here all displayed coefficients are scalars.
Since the third tangent error vanishes and the third output
has weight two, \(F_{3,\epsilon}=G_3(\widetilde a,\widetilde b)\)
exactly. Thus on rectifying \((F_2,F_3)\) the second parameter
stays \(\widetilde b\).
The second-output difference from the reference is
\(\epsilon^2 b c/x+\epsilon^3 s(x)\) for some polynomial \(s\).
Put \(W=V_1+\alpha\) and \(d=s-Wc/x\).
The adjusted first parameter \(a'=\widetilde a+\Delta\)
therefore has
\[
 \Delta=\frac{\epsilon^2 b c/x+\epsilon^3 s}{b+\epsilon W}
        =\epsilon^2 c/x+\epsilon^3 d/b-\epsilon^4 Wd/b^2
          +O(\epsilon^5).
\]
The remainder \(H=F_{1,\epsilon}-G_1(a',\widetilde b)\)
starts in order three. Its coefficients \(H_3,H_4\) are found
from
\[
 G_1(a',\widetilde b)-G_1(\widetilde a,\widetilde b)
 =(2\widetilde a+\epsilon k_1\widetilde b+\epsilon^2 k_3)
                                         \Delta+\Delta^2 .
\]
In particular
\[
 H_3=E^{(3)}_1-2ad/b-(2U_1+k_1 b)c/x
 \tag{12}
\]
with polynomial \(E^{(3)}_1\).
It is a constant of the leading derivation, since the
rescaled determinant has order seven. Its weight is one
and it lies in the ring (11). The monomials of weight
one in that ring have \(x\)-orders \(2-3m\) for \(m\ge0\).
Formula (12) gives pole order at most two, so in fact
\(H_3=h_0 b+h_1 a/b\).

Assume \(c\ne0\). Since \(H_3\) cannot have a double pole,
(12) forces \(d\) to be polynomial and \(W(0)=0\).
Its possible simple-pole residue then gives
\(U_1(0,y,z)=0\): the other terms at this pole are scalar
multiples of \(y^2/z\), and \(U_1(0,y,z)\) is a linear
polynomial. The fourth coefficient is
\[
 H_4= E^{(4)}_1+2a Wd/b^2-(2U_1+k_1 b)d/b
                    -(2U_2+k_1 V_1+k_3)c/x-c^2/x^2 .
\]
It has pole order at most two, with residue at that order
\(2W'(0)d(0)y^2/z^2-c^2\).

For the constant-field test at order four one must remove
the lifted previous coefficient. Subtract
\(\epsilon^3(h_0\widetilde b+h_1 a'/\widetilde b)\) from
\(H\). This subtracts \(h_0 V_1+h_1(U_1/b-aV_1/b^2)\)
from \(H_4\). It depends on the rectified parameters,
so leaves the determinant intact. Its order-four remainder
belongs to \(k[a,b,b^{-1}]\) of weight zero, by the
leading determinant equation and (11). Remove also the
scalar term \(h_1 V_1(0)a/b^2\) at this coefficient.
The resulting weight-zero constant has pole order at most
two. Its residue at order two is a scalar times \(y^2/z^2\)
minus \(c^2\), because \(U_1(0)=0\).
The weight-zero monomials in (11) have orders \(-3m\)
for \(m\ge0\). Thus a weight-zero constant with pole
order at most two must be regular. The displayed residue
cannot vanish for \(c\ne0\). This contradiction absorbs
the second loss, and every later tangent coefficient is
polynomial by the weight bound. The descent factors the
map through two polynomials, contrary to its determinant.

## Exclusion of the transverse Veronese tuple

**Theorem KellerTransverseVeroneseLeadingExclusion.**
Let \(k\) have characteristic zero, let \(m\ge2\) be an integer,
and let \(\mu\in k\). Put \(q=xz+\mu y\), \(a=xy^m\), and \(b=xq\).
There is no three-variable polynomial Keller map with highest
tuple \((a^2,ab,b^2)\) in weights \((0,1,1)\).

**Corollary SexticKellerTangentConicSquareCommonFactorTransverseExclusion.**
Let \(P=xz+y^2\) and let \(Q\) be a homogeneous binary quadratic in
\(x,y\) with \(Q(0,y)\ne0\). The ordinary sextic highest tuple
\(x^2(P^2,PQ,Q^2)\) has no Keller completion.

### The first transverse correction

Write \(R=k[x,y,z]\), \(S=k[a,b]\). In rational coordinates
\(a,b,y\) one has
\[
                    x=\frac{a}{y^m},\qquad
                    z=\frac{b y^{2m}}{a^2}
                                  -\frac{\mu y^{m+1}}{a}.
\]
In particular \(a,b,y\) are algebraically independent,
and \(k(a,b)\) is the rational constant field of
\(\det J(a,b,\cdot)\). The coefficient of \(y^0\) in the
displayed Laurent coordinates defines an \(S\)-linear
projection \(\pi:R\to S\). Indeed a term of \(x^i y^j z^k\)
using the second summand for \(z\) exactly \(n\) times has
\(y\)-exponent \(j-mi+2mk-(m-1)n\).
If this is zero, the exponent of \(a\) is
\(i-2k+n=(j+n)/m\ge0\). The exponent of \(b\) is \(k-n\).
The projection preserves the grading with \(a,b,y\) of
weights \(m,1,1\). It follows also that \(R\cap k(a,b)=S\):
an element in the intersection equals its own coefficient
of \(y^0\).

Let \(E\) be the first transverse correction, with component
weights \(2m-1,m,1\), and write
\(f=(a^2,ab,b^2)\) in parameter coordinates. Put
\(n=f_a\times f_b\). The rescaled determinant has order
\(2m+(m+1)+2-2=3m+1>1\), so the first-order equation is
\[
                            \det J(a,b,n\cdot E)=0 .
\]
This follows by expanding one error row at a time, since
the entries of \(n\) are functions of \(a,b\).
Thus \(n\cdot E\in S\), and \(E-K\), with \(K=\pi(E)\),
is tangent to \(f\). Write
\[
                            E-K=f_a U+f_b V.
\]
The parameter minors have source gcd \(x^2\), since
\(\gcd(a,b)=x\). Solving with those minors permits poles only
at \(x=0\). In fact the poles are at most simple. The first
entry is \(2aU\), and the second is \(bU+aV\), with \(x\mid b\).
The weights of \(U,V\) are \(m-1,0\), respectively. Consequently
\[
              U=U_0(x,y,z)+l(y,z)/x,\qquad
              V=V_0(x)+c/x
\]
for a polynomial \(U_0\) homogeneous of weight \(m-1\), a form
\(l\) of degree \(m-1\), a polynomial \(V_0\), and \(c\in k\).
Since the parameter weights are \(m,1\) and \(m\ge2\),
\[
                          K_2=\alpha a+\beta b^m,\qquad K_3=\eta b
\]
for scalars \(\alpha,\beta,\eta\).

### The highest coefficient along a fixed \(x,q\)

Localize at \(x\) and use \(x,q,y\) as coordinates.
The original determinant in these coordinates is a nonzero
constant divided by \(-x\), in particular independent of \(y\).
The first component has \(y\)-degree \(2m\) and leading
coefficient \(A=x^2\), since the lower transverse terms
have \(y\)-degree at most \(2m-1\).
The second has \(y\)-degree at most \(m\), with coefficient
of \(y^m\) given by
\[
                           B=x^2q+x(V_0+\alpha)+c .
\]
Indeed \(ab\) gives \(x^2q y^m\), and \(aV+\alpha a\)
gives the three remaining contributions to this coefficient.
At fixed \(x,q\), \(bU\) has \(y\)-degree at most \(m-1\),
since \(z=(q-\mu y)/x\), and \(\beta b^m\) is independent of \(y\).
All still lower transverse terms have \(y\)-degree at most
\(m-1\).
Likewise the third component has the exact form
\[
      F_3=b^2+2b W(x)+2c q+g_0(x),\qquad W=V_0+\eta/2,
\]
because its weight-one correction is \(2bV+\eta b\)
and its weight-zero part is a polynomial \(g_0(x)\).

Write \(\{u,v\}_{x,q}=u_xv_q-u_qv_x\).
The coefficient of \(y^{3m-1}\) in the determinant with
rows \(F_3,F_1,F_2\), taken in \(x,q,y\), is
\[
        m B\{F_3,A\}_{x,q}-2m A\{F_3,B\}_{x,q}
                      =-2m x^3\{F_3,B/x\}_{x,q}.
\]
This coefficient must vanish. Set
\[
                H=xq+W(x)+\frac{c}{x},\qquad
                g=g_0-\left(W+\frac{c}{x}\right)^2 .
\]
Then \(B/x=H+\alpha-\eta/2\), \(F_3=H^2+g(x)\), and
\(\{F_3,H\}_{x,q}=x g'(x)\). Hence \(g\) is constant.
But \(H_q=x\ne0\), so the level \(H=0\) has points in the
chart \(x\ne0\), where \(dF_3=2H\,dH=0\).
This contradicts the Keller determinant and proves the theorem.

For the corollary, the preceding transverse normal form and
the exclusion of its discriminant scalar produce a highest
tuple on \(a=Xy^2,\ b=X(Xz+\mu y)\) after the target shear and
source normalizations. These changes preserve the Keller property.
The theorem with \(m=2\) excludes the tuple for every \(\mu\).

## The reduced binary member with a doubled common factor

**Theorem SexticKellerTangentConicSquareCommonFactorReducedBinaryMember.**
Let \(P=xz+y^2\), \(Q=xy\), and suppose a characteristic-zero
Keller map has highest tuple \(x^2(P^2,PQ,Q^2)\). Its degree-five
tuple has the form
\[
                    (2P U,\ Q U+P V,\ 2Q V)
\]
with homogeneous cubics \(U,V\). There are scalars
\((\beta,\gamma)\ne(0,0)\) such that
\[
 U|_{x=0}= y^2(\beta z+\gamma y),\qquad
 x\mid V,\qquad (V/x)|_{x=0}=y(\beta z+\gamma y).
\]
Its highest tuple in weights \((1/2,1,1)\) is
\[
                y^2(x^2+2\beta z+2\gamma y)(y^2,xy,x^2).
\]

Put \(a=xP,\ b=xQ=x^2y\). The pair \(a,b\) is onto
by taking \(x=1\) and solving for \(y,z\); its source
field is \(k(a,b)(x)\) by solving on \(x\ne0\).
Its polynomial constant ring is \(k[a,b]\).
Thus the discriminant \(F_1F_3-F_2^2\) has degree at
most nine, by the first nonzero coefficient argument
used for the transverse member. Its degree-nine part
belongs to \(k[a,b]\) and is divisible by \(x^3\).
Put \(n=(Q^2,-2PQ,P^2)\). The degree-eleven equation
is \(n\cdot E_5=0\), which gives the displayed tangent
form. The degree-ten equation is
\[
                    x^2 n\cdot E_4=(QU-PV)^2.
\]
It gives \(x\mid V\); write \(V=x V_2\). On \(x=0\) it
then becomes
\[
                    y^4(E_4)_3=(yU-y^2 V_2)^2.
\]
Hence \(y\mid U|_{x=0}\); write \(U|_{x=0}=yL_2(y,z)\).

If \(L_2=0\), both cubics \(U,V\) are divisible by \(x\).
The degree-five correction is then absorbed by a
Veronese reference with parameters \(xP+U/x,xQ+V/x\).
Each parameter has degree at most two in \(x\), and the
remaining error has ordinary degree at most four. Thus
the actual map has degree at most four in \(x\). Its
coefficient of \(x^4\) is a translated \((z^2,zy,y^2)\)
plus a constant vector. After translating \(y,z\), the
coefficient of \(x^{11}\) in the determinant forces this
constant vector to vanish, since the normals to the
Veronese conic span all three target directions. This
is the impossible quartic transverse frame: the
[quartic top-coefficient exclusion](SEXTIC_KELLER_NONBINARY_QUARTIC_LINE_FACTOR.md#exclusion-with-a-quadratic-lowest-coefficient)
uses only the transverse degree four, the quadratic
highest coefficient \((z^2,zy,y^2)\), and the Keller
identity from that frame onward. Hence \(L_2\ne0\).

The degree-nine part of the discriminant vanishes on
\(x=0\). There its value is \(2y^3 L_2(E_4)_3\), since
the other entries of the quintic vector vanish on this
plane. Therefore \((E_4)_3|_{x=0}=0\) and the degree-ten
equation gives \(V_2|_{x=0}=L_2\). At weights
\((1/2,1,1)\) the three leading weights are \(5,9/2,4\).
Only the ordinary sextics and the displayed quintic terms
can contribute, since the third quartic has zero restriction
on \(x=0\). The resulting highest tuple is
\[
                       (x^2 y^2+2y L_2)(y^2,xy,x^2).
\]

If the coefficient \(u\) of \(z^2\) in \(L_2\) is nonzero,
complete the square in \(z\) by a linear change preserving
\(y\), to write \(L_2=u Z^2+v y^2\). Fix generic values
\(A\) of the first leading component and \(T=x/y\) of the
ratio of the second to the first. The level equation is
\[
                    A=T^2y^6+2v y^5+2u y^3 Z^2.
\]
After fixing \(T\) the relative area form is a scalar
multiple of \(y\,dy\wedge dZ\). Dividing by the derivative
of the equation in \(Z\), the relative fibre differential
is therefore proportional to \(dy/M\), where \(M=y^2 Z\)
satisfies
\[
                    M^2=\frac{A y-T^2y^7-2v y^6}{2u}.
\]
The polynomial is squarefree of degree seven over the
geometric generic field. The root at zero is simple,
and a repeated nonzero root would make \(A\) a critical
value of \(T^2y^6+2v y^5\). The differential is nonzero
and regular on the smooth projective curve, at simple
branch points and at infinity, where \(y,M\) have pole
orders two and seven. This contradicts the necessary
rational primitive from the leading-pair slice. Thus
\(u=0\) and \(L_2=y(\beta z+\gamma y)\) with at least one
nonzero scalar, as claimed.

## The transverse coefficient in the reduced binary member is nonzero

**Theorem SexticKellerTangentConicSquareCommonFactorReducedBinaryNonzeroTransverse.**
In the preceding reduced binary member theorem, \(\beta\ne0\).
In particular a Keller completion of \(x^2(P^2,PQ,Q^2)\),
with \(P=xz+y^2\), \(Q=xy\), must reach the highest tuple
\[
                       Z y^2(y^2,xy,x^2),\qquad
                       Z=x^2+2\beta z+2\gamma y
\]
in weights \((1,2,2)\), with \(x,y,Z\) weighted source coordinates.

Suppose \(\beta=0\), so \(\gamma\ne0\), and put \(H=x^2+2\gamma y\).
At weights \((1,2,2)\) the leading tuple is \(H y^2(y^2,xy,x^2)\).
Keep the first two weights fixed and increase the weight \(w\)
of \(z\) to its first tie with a leading component.
The original monomial \(x^4z^2\) of the first sextic forces
\(2<w\le3\). At the tie the expected determinant weight
\(24-w\) is positive, so the new leading tuple has rank at
most two. Its specialization at \(z=0\) is the binary leading
tuple, with image dense in the Veronese cone. The new image
closure, being irreducible of dimension at most two and
containing the cone, is therefore that cone.

Write the new tuple as \(g(p^2,pq,q^2)\) with \(p,q\)
coprime weighted polynomials. Restricting to \(z=0\) gives
\(p_0=yh,q_0=xh,g_0=H y^2/h^2\), up to constant
normalizations, since \(x,y\) are coprime. The polynomial
\(H\) is irreducible, so \(h=1\) or \(h=y\).

If \(h=1\), the parameter weights two and one force
\(p=y,q=x\). The common factor \(g\) has weight six.
A \(z\)-term at \(2<w<3\) then requires \(w=5/2\):
the only possibility is \(g=Hy^2+c xz^2\), \(c\ne0\).
At \(w=3\) this \(h=1\) case is impossible, since \(g y^2\)
cannot contain the obligatory \(x^4z^2\).
For \(w=5/2\) fix \(T=x/y\) and the generic first output \(A\).
The relative area after fixing \(T\) is proportional to
\(y\,dy\wedge dz\). Setting \(M=y^2 z\), the generic level
and its differential become
\[
            M^2=\frac{Ay-T^2 y^7-2\gamma y^6}{cT},
                           \qquad \omega\ \propto\ dy/M.
\]
The polynomial on the right is squarefree of degree seven:
zero is a simple root, and a repeated nonzero root would
make \(A\) a critical value of \(T^2y^6+2\gamma y^5\).
The nonzero differential is regular on the complete
normalization and cannot have a rational primitive. This
contradicts the weighted leading-pair slice.

It remains to handle \(h=y\). The parameter weights four
and three and the common-factor weight two permit a tie
only at \(w=3\). Thus
\[
                  g=H,\qquad p=y^2+d xz,\qquad q=xy+e z.
\]
The coefficient of \(x^4 z^2\) in \(g p^2\) gives \(d^2=1\),
so \(d\ne0\).

### The case \(e=0\)

Fix \(T=q/p\) and the first leading output \(A=H p^2\).
Since \(p=xy/T\) and
\(\partial T/\partial z=-T^2 d/y\), the area on the fixed
\(T\)-plane is proportional to \(y\,dx\wedge dy\).
Put \(s=xy\), \(r=xs\), and \(K=A T^2\). The generic
level equation and relative differential, up to nonzero
constants, are
\[
                       s^4=\frac{r(K-r^2)}{2\gamma},
                       \qquad \omega\ \propto\ \frac{dr}{rs}.
\]
Indeed \(dx\wedge dy=dr\wedge ds/r\), \(y=s^2/r\), and
\(K=r^2+2\gamma s^4/r\).

The fourth-root cover is irreducible by a simple root of
\(r(K-r^2)\). The differential has an order-two pole at
the point over \(r=0\), and is regular elsewhere. At the
other two roots, \(dr\) has order three while \(s\) has
order one. At infinity \(r,s\) have pole orders four and
three, so the form is regular there too. If it had a
primitive, projection to the character \(s\mapsto\zeta s\)
with eigenvalue \(\zeta^{-1}\), for a primitive fourth
root \(\zeta\), would give a primitive \(L(r)/s\).
The pole bounds imply \(L\) is polynomial, vanishes at
both nonzero branch values, and has degree at most zero
by the orders at infinity. This is impossible.

### The case \(e\ne0\)

On a generic fixed ratio \(T=q/p\), solve
\[
              z=\frac{T y^2-xy}{e-T d x},\qquad
              p=\frac{y(e y-d x^2)}{e-T d x}.
\]
Pass on the first level \(A=Hp^2\) to \(r^2=H\) with
\(rp=\sqrt A\); equivalently \(r=\sqrt A/p\) there.
Since \(y=(r^2-x^2)/(2\gamma)\), the resulting curve has
equation
\[
 f=r(r^2-x^2)\bigl(e r^2-(e+2\gamma d)x^2\bigr)
                                  -C(e-T d x)=0,
                    \qquad C=4\gamma^2\sqrt A.
\]
The differential is proportional to \(dx/f_r\).
Indeed the fixed-ratio area is \(p\,dx\wedge dy/(e-T d x)\).
On the level its expression in \(x,r\) is proportional to
\(dx\wedge dr/(e-T d x)\), and the derivative of \(rp\)
along that level is the derivative of \(f\) divided by
\(4\gamma^2(e-T d x)\).

This curve is smooth at all affine points for generic \(C,T\).
Off \(e-T d x=0\) it is the generic level of the rational
function obtained by solving \(f=0\) for \(C\): its critical
locus has only finitely many critical values, because in
characteristic zero that function is constant on each
irreducible component of the locus where its differential
vanishes. At \(x=e/(T d)\) the roots in \(r\) are simple
unless \(e+2\gamma d=0\) and \(r=0\). In this exceptional
case \(f_x=C T d\ne0\) there. Thus \(dx/f_r=-dr/f_x\)
extends regularly throughout the affine curve.

When \(e+2\gamma d\ne0\), the homogeneous part of degree five
has five distinct roots on the line at infinity, all in the
chart \(r/x\). At each, \(1/x\) is a parameter and \(f_r\)
has pole order four, while \(dx\) has pole order two.
When \(e+2\gamma d=0\), the two nonzero roots at infinity
are again simple. At the remaining point write
\(u=1/x\), \(v=r/x\); the equation is
\[
                       e v^3(v^2-1)=C(e u^5-T d u^4).
\]
On each normalized branch the orders of \(u,v\) are \(3h,4h\)
for some positive integer \(h\). Hence \(f_r\) has pole order
\(4h\) and \(dx\) has pole order \(3h+1\), so the form is
regular at this point as well. Each component has nonconstant
\(x\), since the leading coefficient in \(r\) is \(e\ne0\).
The resulting differential is nonzero and regular on the
complete normalization, which precludes a rational primitive.
This excludes \(e\ne0\), finishing the proof that \(\beta\ne0\).

## The first lift in the transverse coordinate

**Theorem SexticKellerTangentConicSquareCommonFactorReducedBinaryFirstLift.**
Under the preceding theorem, take \(Z=x^2+2\beta z+2\gamma y\)
as the third source coordinate. Increase its weight from two,
keeping weights one and two on \(x,y\). A first tie must occur
at weight three, and the highest tuple there is
\[
                 Z\bigl(y^4,\ y^2(xy+eZ),\ (xy+eZ)^2\bigr),
                          \qquad e\ne0 .
\]
After constant target shears and an affine linear shift from y to
Y over k[x], the degrees in Z are at most two, two, and three,
with third cubic coefficient e^2. The first two quadratic
coefficients in Z are b_1 Y+b_0 and e Y^2+P(x), respectively,
where deg P is at most three and b_1=0 if P is nonconstant.

The ordinary \(x^4z^2\) term implies the first component has
\(Z\)-degree at least two: replacing \(z\) by its affine
expression in \(Z\) preserves its degree in that variable.
Hence there is a finite first tie. Denote its weight by \(w>2\)
and its highest tuple by \(A,B,C\). Their weights are
\(d_1=8+w,d_2=7+w,d_3=6+w\), and their coefficients of \(Z\)
are \(y^4,xy^3,x^2y^2\), respectively, with no terms of
\(Z\)-degree zero. The expected determinant weight
\(18+2w\) is positive, so the leading field has transcendence
degree at most two.

In fact \(AC=B^2\) at this tie. Choose an integer \(M>0\)
clearing \(d_1,d_2\) and put \(V=B^{M d_1}/A^{M d_2}\).
Both \(V\) and \(R=AC/B^2-1\) are of weight zero, so
\(R\) is algebraic over \(k(V)\), since adjoining \(A\) to
the weight-zero field raises transcendence degree by one.
The order of \(V\) along \(Z=0\) is \(M\), with leading
coefficient \(x^{M d_1}y^{-M(4+w)}\). If \(R\ne0\), its
order there is an integer \(j\ge1\), with coefficient
\(R_0\in(x^2y^6)^{-1}k[x,y]\). The residue of \(R^M/V^j\)
is nonzero and algebraic over \(k\), by finiteness of the
valued extension of \(k(V)\). Its \(y\)-order is therefore
zero. But \(\operatorname{ord}_y R_0\ge-6\) while
\((4+w)j>6\), a contradiction.

Factor the tuple after division by \(Z\) as
\(g(p^2,pq,q^2)\) in coprime weighted parameters \(p,q\).
On \(Z=0\) it is \(y^2(y^2,xy,x^2)\), so
\(p_0=yh,q_0=xh,g_0=y^2/h^2\) after constant
normalization, with \(h=1\) or \(h=y\).
If \(h=1\) then \(p=y,q=x\) and \(g\) has weight four.
If \(h=y\), then \(g=1\) and \(p,q\) have weights four
and three. Only \(Z\)-linear terms can tie in these
parameters or in the common factor, at \(w=3,4\).

### Common-factor lifts

For \(h=1\) the common factor is \(g=y^2+c xZ\) at \(w=3\),
or \(g=y^2+c Z\) at \(w=4\), with \(c\ne0\).
Fix \(T=x/y\) and the first generic level \(A\).
The relative area after fixing \(T\) is proportional to
\(y\,dy\wedge dZ\).
For \(w=3\), the equation is \(A=Z y^4+cT y^3 Z^2\)
and the differential is proportional to \(dy/W\), where
\(W=y^2(2cT Z+y)\) satisfies \(W^2=y^6+4cTAy\).
For \(w=4\), \(A=Zy^4+c y^2 Z^2\), and the same form
has \(W=y(2c Z+y^2)\) with \(W^2=y^6+4c A\).
Both polynomials are squarefree of degree six over the
geometric generic field. The differential is nonzero and
regular on the complete normalization, so has no rational
primitive. This excludes common-factor lifts.

### Parameter lift at weight four

For \(h=y,w=4\), the parameters are \(p=y^2+dZ,q=xy\),
where \(d\ne0\). Fix the ratio \(T=q/p\) and the first
level \(A=Zp^2\). The fixed-ratio area is proportional
to \(p\,dx\wedge dy\), since \(\partial_Z T=-Td/p\).
With \(r=x/y\) the level and differential are
\[
            A=\frac{r^2(r-T)y^6}{dT^3},
                            \qquad \omega\ \propto\ r y^4\,dr.
\]
The sixth-root cover is irreducible by the simple pole
at \(r=T\). A primitive in the character of \(y^4\)
would be \(y^4 L(r)\), with the rational equation, up
to scaling,
\[
                  L'-\tfrac23(2/r+1/(r-T))L=r.
\]
At \(0,T\) the fractional coefficients exclude cancellation
of the leading local power of \(L\); the prescribed
vanishing orders on the right give orders two and one
of \(L\), respectively. There is no pole at another
finite point. At infinity a polynomial of degree \(n>2\)
would give a term of degree \(n-1\) with nonzero
coefficient \(n-2\), so \(\deg L\le2\). This is impossible.

### Parameter lifts at weight three

Now \(p=y^2+d xZ,q=xy+e Z\).
If \(e=0,d\ne0\), fix \(T=q/p,A=Zp^2\). Here the
fixed-ratio area is proportional to \(y\,dx\wedge dy\),
and with \(r=x/y\) the equation and differential are
\[
               A=\frac{r(r-T)y^5}{dT^3},\qquad
                          \omega\ \propto\ y^3\,dr.
\]
The fifth-root cover is irreducible and character projection
gives \(L'-\tfrac35(1/r+1/(r-T))L=1\) for a rational
\(L\). It has no finite poles and must vanish at both
\(0,T\). Comparing degrees at infinity, \(n-6/5\) is
nonzero for every integer degree \(n\), so \(\deg L\le1\),
which is impossible.

Suppose \(d,e\) are both nonzero. Fix \(T=q/p\), put
\(v=p/y\) and \(r=x-Tv\), and set \(\lambda=d/e\),
\(\alpha=\lambda T\). The equations for the level become
\[
             Z=-yr/e,\qquad y=(1+\alpha r)v+\lambda r^2,
                       \qquad r v^2 y^3=K=-eA .
\]
Indeed solving for \(Z\) first gives
\(p=y(e y-d x^2)/(e-T d x)\) and
\(y=(d x^2+(e-T d x)v)/e\).
The relative area is proportional to
\(p\,dx\wedge dy/(e-T d x)\), hence to
\(y v\,dr\wedge dv\) using the displayed formula for \(y\).
The generic level in \(r,v\) is smooth at its finite
points, since \(K\) avoids the finitely many critical
values of the polynomial \(r v^2 y^3\) in characteristic
zero. With \(y=(1+\alpha r)v+\lambda r^2\), its
differential is proportional to
\[
                   \frac{dr}{r y (2y+3(1+\alpha r)v)}.
 \tag{13}
\]
This is the relative area divided by the level differential
and so is regular throughout the smooth affine curve.

Consider a point of the complete normalization where
\(v\) has pole order \(n>0\) and \(r\) is finite. If \(r=0\),
then \(y\sim v\) and \(r\) has zero order \(5n\);
(13) has order \(2n-1\). If \(r\ne0\), then \(y\) has
zero order \(2n/3\) by the equation. Necessarily
\(r=-1/\alpha\) at the point, and \(r+1/\alpha\)
has order \(n\) because \((1+\alpha r)v\) tends to
\(-\lambda r^2\ne0\). The last factor of the denominator
is then a unit. The order of (13) is \(n/3-1\ge0\).

If \(r\) has pole order \(N>0\), put \(s=\operatorname{ord}(v)\).
If \(s<-N\), then \(\operatorname{ord}(y)=s-N\), and the
level equation gives \(5s-4N=0\), impossible. If \(s>-N\),
\(y\) has pole order \(2N\), and the last factor in (13)
also has pole order \(2N\). The form has order \(4N-1\).
If \(s=-N\), cancellation gives \(\operatorname{ord}(y)=N\)
by the level equation, while the last factor still has pole
order \(2N\). The form has order \(N-1\).
Thus (13) is regular everywhere on the complete normalization.
It is nonzero on every component, since \(r\) is not constant
on any component of the nonzero level. It cannot be exact.

The leading-pair rational slice excludes all these cases,
leaving \(d=0,e\ne0\), as claimed.

### The coefficients of the next transverse powers

The weight-three tuple also constrains the ordinary degree in
\(Z\). Write \(A_i\) for the coefficient of \(Z^3\) in output
\(i\) and \(B_i\) for the coefficient of \(Z^2\). Higher powers
are impossible at weights \(11,10,9\). The leading tuple gives
\(A_3=e^2\), \(A_2\) constant, \(A_1=u x+v\), and
\(B_2=e y^2+\) terms of weight at most three in \(x,y\).
Subtract constant multiples of the third output from the first
two to remove \(v,A_2\), preserving their highest forms.
The coefficient of \(Z^7\) in the determinant is
\(3e^2 u\,\partial_y B_2\), since \(A_x=(u,0,0)\),
\(A_y=0\), and the only way to reach that order uses two
highest \(Z^3\) entries and one \(Z^2\) entry.
Therefore \(u=0\). At order \(Z^6\) the determinant equation
then gives \(\det J_{x,y}(B_1,B_2)=0\).

Complete the square in \(B_2\), using \(Y=y\) plus an affine
linear polynomial in \(x\), to write \(B_2=e Y^2+P(x)\),
\(\deg P\le3\). This preserves the weights on \(x,Y,Z\)
and all leading forms. If \(P\) is nonconstant, the generic
fiber of \(B_2\) is geometrically integral: \(s-P(x)\) is
squarefree over the geometric generic field. Thus \(B_1\)
belongs to \(k(B_2)\) and is a polynomial in \(B_2\)
by surjectivity of \(B_2\) over the algebraic closure.
Its weight bound four shows it is affine linear in \(B_2\).
If \(P\) is constant, the determinant equation makes
\(B_1\) a polynomial in \(Y\) of degree at most two.
In either case a constant target shear removing a multiple
of the second output leaves \(B_1=b_1Y+b_0\), where \(b_1=0\)
if \(P\) is nonconstant.

## Exclusion of the reduced binary member with a doubled common factor

**Theorem SexticKellerTangentConicSquareCommonFactorReducedBinaryExclusion.**
Let \(P=xz+y^2\) and \(Q=xy\) over a characteristic-zero field.
The sextic highest tuple \(x^2(P^2,PQ,Q^2)\) has no Keller completion.

In the original coordinates the degree-six terms in the three
outputs are
\[
          x^4z^2+2x^3y^2z+x^2y^4,\qquad
          x^4yz+x^3y^3,\qquad x^4y^2.
\]
The target operations in the first-lift theorem only subtract
multiples of the third output from the first two and a multiple
of the second from the first. Consequently the coefficient of
\(z^2\) in the resulting first output still has ordinary
degree-four part \(x^4\). The second and third outputs have
no \(z^2\) term of ordinary degree six.

In the coordinates of that theorem the first output has
degree at most two in \(Z\) and its coefficient of \(Z^2\) is
\(b_1Y+b_0\). Here \(Z=x^2+2\beta z+2\gamma y\) with
\(\beta\ne0\), and \(Y\) differs from \(y\) by an affine
linear function of \(x\) alone. Returning to the original
source coordinates, its coefficient of \(z^2\) is exactly
\(4\beta^2(b_1Y+b_0)\). This has ordinary degree at most one,
contradicting the degree-four term \(x^4\). This proves the
exclusion.

## The first axial weights in the two square-member frames

**Theorem SexticKellerTangentConicQuadraticOuterSquareMemberFirstTies.**
Put \(P=xz+y^2\). In the frame with ordinary highest tuple
\((x^2P^2,x^4P,x^6)\), the highest tuple in weights \((1,1,w)\)
is \((x^4z^2,x^5z,x^6)\) for \(1<w<4/3\).
At \(w=4/3\) it is either still this tuple or
\[
                           (x^4+c z^3)(z^2,xz,x^2),\qquad c\ne0.
\]
In the frame with ordinary highest tuple \((P^3,P^2x^2,Px^4)\),
the highest tuple for \(1<w<3/2\) is \((x^3z^3,x^4z^2,x^5z)\).
Either it stays this tuple through weight \(5/3\), or there is a
tie at \(w=3/2\) of the form
\[
                        z(x^3+c z^2)(z^2,xz,x^2),\qquad c\ne0 .
\]

### Recovering the two binary variables

We use the normalization argument at the
[first lift of a finite binary leading field](SEXTIC_KELLER_BINARY_REDUCED_LOWER_PENCILS.md#lifting-the-finite-normalization-at-the-first-new-weight).
In the applications below the binary variables \(x,z\) have
weights \(1,w\), \(1<w<2\), and the unused \(y\) has weight one.
Increase the weight of \(y\) to its first tie. The leading
determinant still vanishes: that weight is at most the largest
component weight, and the other two component weights sum to
more than \(1+w\). Thus specialization at \(y=0\) identifies
the two leading image rings. Every rational element of the
binary image field integral over its image ring lifts to a
polynomial of its own weight. Indeed its lift in the new
leading field is integral over the polynomial source ring,
and the source ring is integrally closed. Relations among
these elements and the image coordinates are preserved.

In each exclusion below the binary image field is \(k(x,z)\).
If both \(x,z\) are integral, they therefore lift immediately.
In the other exclusions we recover their polynomial lifts
by divisibility from the lifted integral functions. In either
event the lift of \(x\) equals \(x\), since its weight is one.
The lift of \(z\) can only be \(z\) or \(z+\eta y\), with a
nonzero \(\eta\) possible only at the weight \(w\) of \(y\).
This follows from the restriction to \(y=0\) and the fact
that the moving weight is strictly above one and \(w<2\).
If neither lift depends on \(y\), there is no tie.
Otherwise use the lift of \(z\) as the new \(z\) coordinate
and raise the weight of \(y\) again. The binary image and all
the integral identities are the same in the new coordinates,
and now neither variable can acquire a \(y\)-term.
This also contradicts the existence of a tie for a Keller map.

### The first and second ties

Write \(A,B,C\) for the three leading components.
A monomial \(x^i y^j z^\ell\) of ordinary degree at most
five can tie the reference \(x^{i_0}z^{\ell_0}\) only when
\(\ell>\ell_0\) and \(w=(i_0-i-j)/(\ell-\ell_0)>1\).
Up to \(4/3\) this gives the following possible new terms,
tested in increasing order as their predecessors are removed.
\[
\begin{array}{c|c|l|l|l}
\text{ordinary frame}&w&A&B&C\\ \hline
x^2(P^2,Px^2,x^4)&6/5& & &z^5\\
&5/4& &z^5&yz^4,\ xz^4\\
&4/3&z^5&yz^4,\ xz^4&y^2z^3,\ xyz^3,\ x^2z^3\\ \hline
P(P^2,Px^2,x^4)&5/4& & &z^5\\
&4/3& &z^5&yz^4,\ xz^4
\end{array}
\]
The ordinary degree-six terms with a factor of \(y\) have
smaller weight than their corresponding reference monomials
at every \(w>1\).

At the first tie in either frame, the added term is \(d z^5\)
in \(C\), with \(d\ne0\). The unchanged pair \(A,B\) gives
\(z/x\) and \(x^6\) rationally, and \(C\) then gives \(x^5\).
Thus the image field is \(k(x,z)\).
In the first frame the binary tuple vanishes simultaneously
only at the origin, so \(x,z\) are integral over the image
ring by homogeneous induction. In the second frame let
\(W=xz,\ U=x^2z\). Both are integral by \(W^3=A,\ U^2=B\),
and \(d z^7=C z^2-W U^2\) makes \(z\) integral too.
The preserved relation \(U z=W^2\) forces the lifted \(z\)
to divide the lifted \(W\), since the lifted \(z\) is \(z\)
or \(z+\eta y\), in particular prime. Now \(x=W/z\) lifts
as a polynomial. This excludes both first ties.

At the second tie the unchanged first component \(A\) and
the second component \(B\), possibly with a \(c z^5\) term,
are algebraically independent in \(x,z\): their Jacobian
has the nonzero contribution of their highest \(x\)-terms.
Therefore the coefficient of \(y z^4\) in \(C\) vanishes
by the leading determinant. Both second ties are now
binary. We compute their fields together with the third
tie of the first frame.

### The binary field at the remaining ties

At weight \(4/3\) in the first frame, write the possible
corrections as \(a z^5\), \(z^4(b y+c x)\), and
\(z^3(d y^2+e xy+f x^2)\), respectively. The leading
determinant coefficients of \(x^9 z^5,x^8y z^5,x^4 y z^8\)
are
\[
          -12b+6e,\qquad 12d,\qquad 50ad-6be-28cd.
\]
Hence \(b=d=e=0\). The surviving tuple is binary.

In these three binary situations use weights \(k,k+1\)
on \(x,z\), with \(k=4\) at \(5/4\), and \(k=3\) at
\(4/3\). Put \(t=z^k/x^{k+1}\), \(\rho=z/x\).
The components have the form
\[
 A=\frac{\rho^D}{t^6}(1+a t),\qquad
 B=\frac{\rho^{D-1}}{t^6}(1+c t),\qquad
 C=\frac{\rho^{D-2}}{t^6}(1+f t).
\]
Here \(D=26\) or \(21\) at the second ties in the first
and second frames, respectively, and \(D=20\) at the
third tie in the first frame. At either second tie \(a=0\).
The two weight-zero fractions
\[
              R_1=AC/B^2=\frac{(1+a t)(1+f t)}{(1+c t)^2},
       \qquad R_2=B^D/A^{D-1}
                   =\frac{(1+c t)^D}{(1+a t)^{D-1}t^6}
\]
recover \(t\) whenever \(R_1\) is nonconstant. To see this,
a proper field they generate has index two, since
\(R_1\) has degree at most two. Its involution must preserve
\(R_2\). If \(a,c\) are distinct and nonzero, the two
poles at \(0,-1/a\) have distinct orders \(6,D-1\),
and the zeros at \(-1/c,\infty\) have distinct orders
\(D,5\). If exactly one of \(a,c\) vanishes, the two
poles have distinct orders and there is a unique zero
at \(-1/c\) or \(\infty\), respectively. If \(a=c\ne0\),
there is a unique pole at zero and zeros of distinct
orders one and five. Each case forces at least three
fixed points and hence the identity. If \(a=c=0\) and
\(R_1\) is nonconstant, it has degree one.
Knowing \(t\), one recovers \(\rho\) from \(A/B\), so the
binary image field is \(k(x,z)\). If \(R_1\) is constant,
then \(a=c=f\), since \(R_1(0)=1\). This gives no tie
at the second weights and exactly the displayed exceptional
binomial tuple at the third weight of the first frame.

We first rule out the birational cases in the first frame
and the case with a nonzero second-component correction
in the second frame.
At \(5/4\) in the first frame, \(A=x^4z^2\),
\(B=x^5z+c z^5\), \(C=x^6+f xz^4\).
If \(c\ne0\), the tuple vanishes simultaneously only
at the origin; its binary variables are integral.
If \(c=0,f\ne0\), \(W=x^2z\) is integral since
\(W^2=A\), and \(x\) is integral by
\[
                         x^{13}+f W^4-C x^7=0.
\]
In the lifted ring \(x\) stays the linear variable \(x\).
The relation gives \(x^7\mid W^4\), hence \(x^2\mid W\),
and \(z=W/x^2\) lifts as a polynomial. This excludes the
second tie of the first frame.

At \(4/3\) in the second frame with \(c\ne0\), \(W=xz\)
is integral since \(W^3=A\), and \(z\) is integral
by \(c z^7=B z^2-W^4\). Then \(U=x^2z\) is integral
by \(U^2=B-c z^5\). The relation \(U z=W^2\) again
recovers \(x\) from the lifted \(z\).
If \(c=0,f\ne0\), the second tie still permits
\((x^3z^3,x^4z^2,x^5z+f xz^4)\).

At \(4/3\) in the first frame, if \(a\ne0\) and \(R_1\)
is nonconstant, the tuple has no simultaneous zero
outside the origin. On the axes this follows from
\(a\ne0\) and the \(x^6\) in \(C\). Off the axes a common
zero would require \(a=c=f\). Thus the binary variables
are integral. If \(a=0,f\ne0\), set \(W=x^2 z\). The
relations \(A=W^2\) and \(x^{10}+f W^3-C x^4=0\)
make \(W,x\) integral, and the lifted relation forces
\(x^2\mid W\). If \(a=f=0,c\ne0\), use \(C=x^6\),
\(A=W^2\) and \(x^7 B=x^{10}W+cW^4\) for the same
conclusion. This finishes the early exclusions in the first frame.

### Completing the early ties in the smooth-common-factor frame

In the residual at \(4/3\) in the smooth-common-factor
frame, \(W=xz\) and \(U=x^2z\) are integral, and the
image field is \(k(x,z)\), as computed above. The identity
\[
                    C W^3 U^3=U^7+dW^{10}
\]
implies that the lifted \(W,U\) have exactly the same
prime support: a prime dividing either one must divide
the other. The lifted \(W\) cannot be irreducible, since
then \(U\) would be a power of it up to a scalar, contrary
to their weights \(7/3,10/3\). Its irreducible factors
can only specialize to the two factors of \(xz\), with
their respective weights, since each has positive weight
and nonzero specialization. Thus \(W=xZ\) with
\(Z=z\) or \(Z=z+\eta y\). The specialization and equal
support then give \(U=x^2 Z\). The leading outputs now
depend only on \(x,Z\), and the same argument applies
at a putative next tie after using \(Z\) as a coordinate.
This excludes the last case at \(4/3\).

The next possible tie in this frame is at \(3/2\), with
corrections \(a z^5,\ z^4(b y+c x),\
z^3(d y^2+e xy+f x^2)\).
The determinant coefficients of \(x^7z^7,x^6y z^7,
x^3 y z^9\) are
\(-12b+6e,\ 12d,\ 40ad-6be-18cd\).
Thus \(b=d=e=0\). The preceding binary field calculation
applies with \(k=2,D=15\), and shows birationality
unless \(a=c=f\). In the exceptional case the tuple
is \(z(x^3+c z^2)(z^2,xz,x^2)\) as stated.

In the birational case, if \(a\ne0\), use \(W=xz\)
and the identities
\[
 W^3=A-a z^5,\qquad
 B z^2=W(A+(c-a)z^5),\qquad
 C z^4=W^2(A+(f-a)z^5).
\]
At least one of \(c-a,f-a\) is nonzero.
Cubing the corresponding identity gives an equation for
\(z\) with nonzero scalar leading coefficient, so \(z\)
and then \(W\) are integral. In the lifted ring the middle
identity forces the lifted prime \(Z\) to divide \(W\):
at \(Z=0\) its right side is \(W^4\).
Thus \(x=W/Z\) lifts too.

If \(a=0\), then \(W^3=A\), and \(U=x^2z\) is integral
by \(U^5-BU^3+c W^7=0\). The remaining identity
\[
                           C W^3 U=U^5+f W^7
\]
shows that any prime dividing the lifted \(W\) also divides
\(U\). The converse follows from the monic equation
for \(U\) if \(c\ne0\) and from the last identity if
\(f\ne0\). At least one of \(c,f\) is nonzero. Again the
weights, now \(5/2,7/2\), and the specialization \(W=xz\)
force \(W=xZ,\ U=x^2 Z\). This excludes the birational
cases at \(3/2\).

Finally, if there is no tie up to \(3/2\), the only
possible tie at \(5/3\) is \(C=x^5z+d z^4\) with
\(d\ne0\). From \(A,B\) one has \(t=z/x,x^6\), and
\(C\) also gives \(x^4\). Thus the image field
is \(k(x^2,t)\). Its elements \(W=xz,\ S=z^2\)
are integral: \(W^3=A\), and
\(C S=W B+d S^3\). In the lifted polynomial ring
\(B S=W^4\), so every prime of \(S\) divides \(W\).
The lifted \(S\) cannot be irreducible since its weight
\(10/3\) exceeds the weight \(8/3\) of \(W\).
It factors into two lifts of \(z\), each of weight
\(5/3\) and equal to \(z\) or \(z+\eta y\).
They cannot be distinct, since their product would divide
\(W\) and exceed its weight. Thus \(S=Z^2\) and \(Z\mid W\),
so \(W=x Z\) by the weight and specialization of \(W/Z\).
These identities recover the original binary tuple in \(x,Z\),
and preclude a first \(y\)-tie, or a second one after the
coordinate change. This excludes \(5/3\) as asserted.

## Exclusion of the binomial axial ties

**Theorem SexticKellerTangentConicQuadraticOuterSquareMemberBinomialExclusion.**
Neither exceptional binomial in the preceding axial tie
theorem has a Keller completion in its stated ordinary
frame. Thus \((x^2P^2,x^4P,x^6)\) has no axial tie
through \(w=4/3\), and \((P^3,P^2x^2,Px^4)\) has none
through \(w=5/3\).

Scale the respective weights to \((3,3,4)\) or \((2,2,3)\).
Use \(A,B,C\) for the highest tuple. At the proposed
binomial tie it is \(g(z^2,xz,x^2)\), with
\[
        g=x^4+c z^3\quad\hbox{or}\quad g=z(x^3+c z^2),
                 \qquad c\ne0 .
\]
In each case the binary common factor \(g\) is squarefree:
\(x^4+c z^3\) or \(x^3+c z^2\) is irreducible over
the algebraic closure by the coprime exponents, and
in the second case it is coprime to \(z\).

Raise the weight \(W\) of \(y\) alone. In the first frame
the ordinary monomial \(x^2y^4\) of the first component
must tie by \(W=7/2\), since its weight is \(6+4W\)
and the binary first weight is twenty. In the second
frame \(y^6\) must tie by \(W=5/2\), since the binary
first weight is fifteen. At the first tie the three
output weights still have sum larger than the source
weights, so the highest map has rank at most two.
Its specialization at \(y=0\) has dense image in the
Veronese cone. The new irreducible image closure,
of dimension at most two, is therefore that cone.
Write the new tuple as \(G(p^2,pq,q^2)\) with
coprime weighted parameters. At \(y=0\) they give
\(p_0=hz,\ q_0=hx\), where \(h^2\) divides the
squarefree \(g\), hence \(h\) is scalar.

In both frames the first tie is strictly below the
weight of \(z\). The parameter lifts must then be
\(p=z,q=x\). In fact the weight of \(p\) is four or
three, respectively, and a positive power of \(y\)
can neither have that weight nor attain it with an
additional factor of \(x\) or \(z\) at the allowed
moving weights. The smaller-weight \(q\) cannot involve
\(y\) at all.

In the first frame \(G\) has weight twelve. In the
range \(3<W\le7/2\), no new monomial
\(x^i z^j y^\ell\) can have weight twelve: for
\(\ell=1\) the weight \(W\) would have to be integral,
for \(\ell=2\) the only possible \(2W\) is seven
leaving five, and for \(\ell=3\) the only possible
\(3W\) is ten leaving two. Neither five nor two
is a nonnegative combination of three and four.
Higher powers of \(y\) already exceed the weight.
Thus the required tie cannot occur.

In the second frame \(G\) has weight nine.
For \(2<W\le5/2\) the only possible added terms are
\[
             y^4\text{ at }9/4,\qquad
             x y^3\text{ at }7/3,\qquad
             x^2 y^2\text{ at }5/2 .
\]
This follows by solving \(2i+3j+W\ell=9\) with
\(\ell>0\). Write the added term as \(d x^{4-n}y^n\),
\(d\ne0\), where \(n=4,3,2\). Fix the generic first
output \(A=G z^2\) and ratio \(T=x/z\). On that level,
\[
            A=T^3 z^6+c z^5+d T^{4-n}z^{6-n} y^n .
\]
The fixed-ratio area is proportional to \(z\,dy\wedge dz\),
so the relative fibre differential is proportional to
\(dz/(z^{5-n}y^{n-1})\).

Let \(D(z)=A-T^3 z^6-c z^5\), which has nonzero constant
term and six simple nonzero roots for generic \(A,T\).
The \(n\)-th root cover is irreducible by a simple
root of \(D\). If \(n=4\), then \(y^4\) is a nonzero
constant times \(D/z^2\). Over zero, if \(z\) has
order \(2h\), then \(y\) has order \(-h\) and
\(dz/(z y^3)\) has order \(3h-1\). At a root
of \(D\), \(y\) has order one and \(dz\) has
order three. At infinity \(y,z\) have equal
pole orders, and \(dz/(z y^3)\) is regular.
Thus it is a nonzero regular differential on the
complete normalization.
For \(n=3\), set \(M=zy\). The model is \(M^3\)
equal to a nonzero scalar times \(D(z)\), with
relative differential \(dz/M^2\), again nonzero
and regular at the finite branch points and
infinity. For \(n=2\), put \(M=z^2 y\). Now
\(M^2\) is a nonzero scalar times \(D(z)\) and
the differential is \(dz/(z M)\), with simple
poles of nonzero residue at \(z=0\) since \(D(0)\ne0\).
None of these differentials has a rational primitive.
The weighted leading-pair slice gives a contradiction
in every case.

## The next axial weight in the doubled square-member frame

**Theorem SexticKellerTangentConicDoubledSquareMemberWeightThreeHalves.**
For an ordinary highest tuple \((x^2P^2,x^4P,x^6)\),
\(P=xz+y^2\), there is no axial tie before \(w=3/2\).
At weight \(3/2\), the highest tuple in \((1,1,w)\)
is either binary of the form
\[
           (x^4 z^2+b xz^4,\ x^5z+e x^2z^3,\
                       x^6+j x^3z^2+f z^4)
\]
or has a nonzero \(y\)-term and is
\[
                (x^4+z^2(i y+e x))(z^2,xz,x^2),\qquad i\ne0.
\]
In the latter case affine normalizations in \(x,y\)
give the complete highest tuple in the \(z\) variable
as \((y z^4,xy z^3,x^2 y z^2)\), of \(z\)-degrees
\(4,3,2\).

After the preceding exclusions there are no possible
ties before \(3/2\), by testing
\(w=(i_0-i-j)/(\ell-\ell_0)\) as in the first tie
table. At \(3/2\) the three possible corrections are
\[
 \bigl(z^4(a y+b x),\
       z^3(c y^2+d xy+e x^2),\
       f z^4+ z^2(g y^3+h xy^2+i x^2y+j x^3)\bigr).
\]
The determinant coefficients of \(x^{10}z^4\),
\(x^9y z^4\), and \(x^8 y^2z^4\) give \(a=2d-i\),
\(h=2c\), and \(g=0\), respectively. The coefficient
of \(x^4y^3z^6\) then gives \(-12ch=0\), hence
\(c=h=0\). The coefficient of \(x^6y z^6\) becomes
\(-12(d-i)^2\), so \(a=d=i\).
If \(i=0\) this is the binary tuple in the statement.
If \(i\ne0\), the coefficients of \(y z^{10}\)
and \(x^3y z^8\) give \(f=0\), \(b=2e-j\).
The coefficient of \(x^4 z^8\) is then
\(10i(e-j)^2\). Therefore \(b=e=j\), giving
the stated common factor.

In the nonbinary case take \(Y=i y+e x\) as a source
coordinate. For \(w>3/2\) the new reference tuple
is \((Y z^4,xY z^3,x^2Y z^2)\).
There is no higher \(z\)-degree in the first component;
in the second, only a constant multiple of \(z^4\)
can tie next; and in the third, only affine-linear
multiples of \(z^3\) can tie. These are the only possible
terms of higher \(z\)-degree, by the previously excluded
ties and the coefficients at \(3/2\).
At \(w=2\), with \(y\) now denoting \(Y\), the most general
tying tuple is
\[
          (y z^4,\ xy z^3+c z^4,\
                            x^2y z^2+(d x+e y)z^3).
\]
Its determinant is \(z^9(4cdz+(8c-4d)xy+e y^2)\);
hence \(c=d=e=0\). The constant \(z^3\) term in the
third component can tie only at \(w=3\). At fixed
first and second outputs of the reference, \(z\)
is free and \(x/z\) is fixed. The third reference
output is constant on this level, so adding a
nonzero multiple of \(z^3\) contradicts the leading
determinant. Thus the \(z\)-degrees are \(4,3,2\).

It remains to normalize the three highest \(z\)
coefficients, allowing terms that never tied at
positive weights. After translating \(Y\), these
coefficients have the form
\[
       A=y,\qquad B=x(y+r)+s y+t,\qquad
       C=x^2(y+u)+\text{terms of \(x\)-degree at most one}.
\]
The bounds follow from the ordinary degree limits
and the listed coefficients at \(3/2\).
The zero-weight leading tuple
\((z^4 A,z^3 B,z^2 C)\), of weights \((0,0,1)\),
has rank at most two by the leading determinant.
Thus its weight-zero fractions \(R=AC/B^2\) and
\(V=B^4/A^3\) are algebraically dependent. At the
valuation \(x=\infty\), \(V\) has pole order four.
Consequently the residue of \(R\), which is
\(y(y+u)/(y+r)^2\), is scalar: the residue
extension of the finite field extension of \(k(V)\)
is algebraic. Comparing coefficients gives
\(r=u=0\).

If \(R\ne1\), then \(R-1\) has order \(n=1\) or
\(2\) at \(x=\infty\), since its numerator has
\(x\)-degree less than two. With leading coefficient
\(r_0(y)\), the element \((R-1)^4 V^n\) has nonzero
residue \(r_0(y)^4 y^n\), which must be scalar.
The order at \(y=0\) makes this impossible for
\(n=1,2\). Thus \(y C=B^2\). In particular \(y\)
divides \(B\), and translating \(x\) by the scalar
\(s\) gives the asserted coefficients.

## Weight two in the smooth-common-factor frame

**Theorem SexticKellerTangentConicSmoothSquareMemberWeightTwo.**
For the ordinary highest tuple \((P^3,P^2x^2,Px^4)\) with
\(P=xz+y^2\), there is no axial tie before \(w=2\).
At weight two the highest tuple in \((1,1,w)\) is either
binary of the form
\[
       (x^3z^3+b xz^4,\ x^4z^2+e x^2z^3+f z^4,\
                         x^5z+j x^3z^2+\ell xz^3)
\]
or of the form \(z(x^3+z(i y+e x))(z^2,xz,x^2)\)
with \(i\ne0\). In the latter case its complete axial
highest tuple is \((y z^4,xy z^3,x^2y z^2)\) after affine
normalizations, of \(z\)-degrees \(4,3,2\).

The axial monomial test gives no tie between \(5/3\)
and two. At weight two, after the previous exclusions,
the possible corrections are
\[
 \bigl(z^4(a y+b x),\
        z^3(c y^2+d xy+e x^2)+f z^4,\
        z^2(g y^3+h xy^2+i x^2y+j x^3)+z^3(k y+\ell x)\bigr).
\]
In increasing elimination order the determinant
coefficients of \(x^8z^6,x^7y z^6,x^6y^2 z^6\)
give \(a=2d-i\), \(h=2c\), and \(g=0\).
The coefficient of \(x^3y^3 z^7\) then gives
\(-12ch=0\), and that of \(x^5y z^7\) reduces to
\(-12(d-i)^2\). Thus \(c=h=0\), \(a=d=i\).
Now the coefficients of \(x^6z^7\) and \(x^3y z^8\)
are, respectively,
\[
                  5i(2e-j-b)+6k,\qquad
                  5i^2(2e-j-b)+7ik .
\]
Subtracting \(i\) times the first from the second
gives \(ik=0\), and the first then gives \(k=0\)
also if \(i=0\). For \(i=0\) this leaves the
stated binary tuple. For \(i\ne0\) it gives
\(b=2e-j\). The coefficient of \(xy z^9\)
becomes \(4i^2(2f-\ell)\), and that of \(x^4z^8\)
then reduces to \(10i(e-j)^2\).
Thus \(\ell=2f\) and \(b=e=j\). At \(z^{10}\)
the coefficient is now \(8if^2\), so \(f=\ell=0\).
This is the nonbinary tuple in the statement.

In that case take \(Y=i y+e x\). Above weight two
the new reference is \((Y z^4,xY z^3,x^2Y z^2)\).
The only possible further tie is a constant
multiple of \(z^3\) in the third component, at
\(w=3\): the ordinary support and all the previous
ties exclude any higher \(z\) terms in the other
components and any nonconstant \(z^3\) coefficient
in the third. As before the first two leading
outputs fix \(x/z\) and leave \(z\) free, excluding
this tie. The highest \(z\)-coefficients have the
same bounds and forms as in the doubled square
member just proved: the first is affine linear
with nonzero \(Y\)-coefficient, the second has
quadratic part \(xY\), and the third has cubic
part \(x^2Y\). The weight-zero argument there
therefore gives the same axial normal form.

## The binary tie in the doubled square-member frame

**Theorem SexticKellerTangentConicDoubledSquareMemberBinaryRefinement.**
Let \(P=xz+y^2\). For the ordinary highest tuple
\((x^2P^2,x^4P,x^6)\), suppose a Keller map has a nontrivial
binary tie at weights \((1,1,3/2)\). Its weighted tuple
must be
\[
                 (x^4z^2,x^5z,x^6+jx^3z^2),\qquad j\ne0 .
\]
The next tie obtained by increasing the weight of \(y\)
occurs at \(5/4\). There the highest tuple is
\[
                         (x^2P^2,x^4P,x^6+jxP^2).
\]

Use the binary form from the weight-three-halves theorem,
with coefficients \(b,e,j,f\) in the order displayed there.
Write \(t=z^2/x^3\), \(\rho=z/x\); with the binary weights
scaled to two and three the leading outputs are
\[
 \frac{\rho^{14}}{t^6}(1+bt),\quad
 \frac{\rho^{13}}{t^6}(1+et),\quad
 \frac{\rho^{12}}{t^6}(1+jt+ft^2).
\]
Put \(R=AC/B^2\) and \(S=B^{14}/A^{13}\), where
\(A,B,C\) are these outputs. Then
\[
 R=\frac{(1+bt)(1+jt+ft^2)}{(1+et)^2},\qquad
 S=\frac{(1+et)^{14}}{(1+bt)^{13}t^6}.
\]
The image field is \(k(R,S,A/B)\), since the three
outputs are recovered from those elements. The degree
of the source field over it is \([k(t):k(R,S)]\):
\(A/B\) is a nonzero scalar multiple of \(\rho\)
over \(k(t)\) and hence is transcendental over \(k(t)\).

### Nonconstant projective ratio

If \(R\) is nonconstant, the binary image field is
\(k(x,z)\) except when \(b=e=j=0,\ f\ne0\), where it
is \(k(x^2,xz,z^2)\).
Indeed the index divides \(\deg R\le3\) and \(\deg S\).
If \(b\ne0\) and \(b\ne e\), the two poles of \(S\)
have orders six and thirteen, so its degree is
nineteen. If \(b=0,e\ne0\), its poles at zero and
infinity have distinct orders six and eight. The
only possible proper index is two, and its involution
preserving \(S\) must fix both poles and the unique
zero \(-1/e\), which is impossible. If \(b=e\ne0\),
then \(S=(1+et)/t^6\), and \(R\) has degree at most
two by cancellation. The unique pole and the two
zeros of distinct orders again preclude an involution.
Finally, if \(b=e=0\), an index-two involution
preserving \(t^{-6}\) must send \(t\) to \(-t\).
This requires \(R=1+ft^2\) with \(f\ne0\), and \(A/B=\rho\)
then gives the stated fixed field.

At any first lift in \(y\), its moving weight is in
\((1,5/4]\), since the ordinary \(x^2 y^4\) term in
the first output ties by weight \(5/4\). The highest
image ring specializes isomorphically to the binary
one, by the first-tie argument used above. In particular
integral elements of the binary image field lift to
polynomials and all rational identities persist.
A polynomial lift of \(x\) or \(z\), of weight one
or \(3/2\), cannot contain \(y\) in this weight range.
Therefore a tie is impossible whenever both lifts
are polynomial.

In the even-field exception the binary tuple vanishes
simultaneously only at the origin. Thus \(X=x^2\),
\(W=xz\) and \(Z=z^2\) are integral and lift. By
weight \(X\) lifts unchanged, since no monomial in
\(y\) has the required weight two, even after
multiplication by a power of \(x\) or \(z\).
The lifted equation \(W^2=XZ\) then makes \(W/x\)
a polynomial lift of \(z\). The expressions
\(A=X^2 Z,\ B=X^2 W,\ C=X^3+f Z^2\) show that
no tie occurs, a contradiction.

Consider the birational cases with \(f\ne0\).
If the leading tuple has no common zero outside the
origin, both variables lift by integrality.
A common zero must lie off the axes and satisfy
\(b=e\ne0,\ C=(x^3+b z^2)(x^3+h z^2)\), with
\(h=f/b\ne0\). Put \(W=x^3+b z^2,\ U=W x^3,\ V=W z^2\).
The identities
\[
                C=U+hV,\quad W^2=U+bV,\quad A^3W=UV^3
\]
make \(W,U,V\) integral. If \(b\ne h\), solve for
\(U,V\) as affine polynomials in \(W^2,C\); the last
identity then gives a monic equation for \(W\) after
division by a nonzero scalar. If \(b=h\), use \(W^2=C\)
and substitute \(U=C-bV\) to make \(V\) integral,
then \(U\). The lifted \(W\) is prime since it is
homogeneous of positive weight and specializes to the
irreducible \(x^3+b z^2\). The third and second identities
force it to divide both \(U\) and \(V\). Now \(x=A/V\)
has cube \(U/W\), and \(z=B/(W x^2)\) has square \(V/W\).
Thus the rational lifts \(x,z\) are polynomial, excluding
this case too.

If \(f=0\) and \(R\) is nonconstant, \(x\) is integral
over the image ring. For \(j=0\) use \(C=x^6\). For \(j\ne0\), put
\(M=C-x^6=jx^3 z^2\). Then
\[
 j^2 x^5 A=M\bigl(b C+(j-b)x^6\bigr),\qquad
 j^3 x^5 B^2=M\bigl(e C+(j-e)x^6\bigr)^2 .
\]
At least one of \(b,e\) differs from \(j\), since \(R\) is
nonconstant. The corresponding equation has a nonzero
scalar leading coefficient in \(x\), so makes \(x\) integral.

The lift of \(x\) is unchanged. If \(b\ne0\), a pole
of the rational lift of \(z\) is impossible, since
\(b xz^4\) would dominate \(x^4z^2\) in \(A\) and
give a pole. This holds also at the prime \(x\), where
the order of \(x\) is one. If \(b=0,e\ne0\), use
\(e x^2z^3\) and \(x^5z\) in \(B\) in the same way.
Thus there is no tie in either case. With \(b=e=0\),
necessarily \(j\ne0\). The identity
\(C-x^6=jx^3 z^2\) shows that the lift of \(z\) can
have only a simple pole at \(x=0\).
Thus \(xz\) lifts to a polynomial of weight \(5/2\),
and the only possible \(y\)-term in it is \(d y^2\),
at weight \(5/4\). Substitution into \(B=x^4(xz)\)
gives \(d=1\) from its ordinary \(x^4y^2\) term, and
then \(A=x^2P^2,B=x^4P,C=x^6+jxP^2\), as asserted.

### Constant projective ratio

If \(R\) is constant, its value is one. The factor identity
gives either \(b=e=j,\ f=0\) or
\(b=0,\ j=2e,\ f=e^2\).
In the first case \(b\ne0\) for a nontrivial tie, and the
binary leading tuple is \(x(x^3+b z^2)(z^2,xz,x^2)\).
The first \(y\)-lift still has rank two with conic image:
it contains the binary conic image on \(y=0\). The
specialized common factor is squarefree, so in the
conic factorization the parameters restrict to \(z,x\)
and cannot acquire a \(y\)-term at moving weights
in \((1,5/4]\). Their common factor has weight four.
The only possible new term in it is \(d zy^2\), at
weight \(5/4\); this follows by solving
\(2i+3k+\ell W=8\) for \(\ell>0,\ 2<W\le5/2\)
in weights scaled by two. The resulting first output
cannot contain \(x^2y^4\), which is obligatory at
that weight. Thus this case is impossible.

In the second case \(e\ne0\), the tuple is the
Veronese on the coprime parameters \(x^2z,x^3+e z^2\).
At the first \(y\)-lift its image is again conic, with
constant common factor by specialization. The second
parameter cannot change at weights in \((1,5/4]\).
The first parameter can acquire \(d y^3\) at \(7/6\)
or \(dxy^2\) at \(5/4\), with no other possibility by
the weights \(7,6\) in weights scaled by two.
The former produces \(y^6\) in the first component,
contrary to the original sextic. In the latter, the
ordinary coefficient of \(x^4 y^2\) in the second
component gives \(d=1\). The lifted parameters are
therefore \(a=xP,\ b=x^3+e z^2\).

On generic fixed levels \(a=a_0,\ b=b_0\), the relative
differential is proportional to \(dx/(xyz)\), by
the determinant in \(y,z\) of the two level equations
\[
                    z^2=(b_0-x^3)/e,\qquad
                    y^2=a_0/x-xz .
\]
This differential is regular and nonzero on the complete
normalization. Over \(x=0\), \(z\) is a unit, \(x\)
has order two, \(y\) has pole order one, and \(dx\)
has order one. At the simple roots of \(b_0-x^3\),
\(y\) is generically nonzero, and the orders of \(dx\)
and \(z\) agree. The zeros of \(a_0/x-xz\) are simple
on the first double cover, since \(a_0\) avoids the
critical values of the nonconstant function \(x^2 z\).
Ramification there likewise cancels the factor \(y\)
in the denominator. At infinity the pole orders
of \(x,z,y\) are \(4,6,5\), while \(dx\) has order
of pole five. Thus the form cannot have a rational
primitive, contradicting the leading-pair slice for
the Veronese outputs. This excludes the constant-ratio
case and completes the proof.

## The lower \(x\)-weights of the mixed doubled frame

**Theorem SexticKellerDoubledSquareMixedLowXWeights.**
Under the mixed binary refinement, suppose \(j\ne0\) and
the highest tuple at weights \((4,5,6)\) is
\((x^2P^2,x^4P,x^6+jxP^2)\). After a constant target shear
of the second output by the third, the highest tuple in
weights \((s,5,6)\), for every \(0<s<4\), is
\[
                        (x^2y^4,x^4y^2,jxy^4).
\]
In particular the second output has the form
\[
                             q(x)y^2+u(x)y+v(x)z+w(x)
\]
with \(q,v\) monic of degrees four and five.

For \(s<4\) just below four the displayed tuple is highest.
Write \(A,B,C\) for its three entries, and retain the notation
for the entries at any tie. The other monomials of ordinary
degree six are below the reference weights for \(0<s<4\).
For the terms of degree at most five, use the strict weight
bounds from the binary \((1,1,3/2)\) tuple and the mixed
\((4,5,6)\) tuple, except for the explicitly retained terms
of those tuples. A monomial \(x^i y^h z^k\) can tie a
reference \(x^p y^q\) from below only if \(i<p\), at
\(s=(5h+6k-5q)/(p-i)\). The positive candidates are
\[
\begin{array}{c|l|l|l}
s&A&B&C\\ \hline
15/4&&y^5&\\
11/3&&xy^3z&\\
7/2&y^3z^2&z^4,\ x^2 y z^2&\\
10/3&&xy^4&\\
13/4&&yz^3&\\
3&y^4z,\ xyz^3&y^2z^2,\ x^2y^2z&yz^3\\
11/4&&y^3z&\\
8/3&&xz^3&\\
5/2&y^5&y^4,\ x^2y^3&\\
7/3&&xyz^2&\\
2&z^4,\ xy^2z^2&z^3,\ xy^2z,\ x^3z^2&y^2z^2\\
7/4&&yz^2&\\
5/3&&xy^3&\\
3/2&yz^3&y^2z&\\
5/4&&y^3&\\
1&y^2z^2,\ xy^3z&x^2z^2,\ x^3yz&y^3z\\
2/3&&xz^2&\\
1/2&y^3z&z^2,\ x^2yz&\\
1/3&&xyz&\\
1/4&&yz&
\end{array}
\]
The \(j xP^2\) terms with more than one power of \(x\)
are below \(jxy^4\) at all these weights.
The expected determinant weight is \(6s+39>0\) at
each tie, so each leading tuple has rank at most two.

An isolated \(z\)-dependent term in the second component
cannot tie: the first and third unchanged entries fix
\(x,y\) algebraically, leaving \(z\) free. At \(10/3\)
the second component can tie by a multiple of \(xy^4\);
remove it by subtracting the corresponding multiple of
the third output. This does not change the tuple at any
greater weight, since the third weight is then strictly
smaller than the second. The prior strict weight bounds
for the second component's unlisted terms still hold,
since the third component has smaller weights in both
prior frames. The introduced ordinary \(x^6\) term in
the second component is below its reference for \(s<4\).
In the rest of the proof we use this sheared second
component.

### Simultaneous ties involving \(z\)

At \(7/2\) and \(1/2\), fix the generic value \(D=C/j\)
and set \(x=D/y^4\). At \(7/2\) the first two outputs are
\[
        A=D^2/y^4+a y^3 z^2,\qquad
        B=D^4/y^{14}+b z^4+dD^2 z^2/y^7.
\]
At \(1/2\) they have the same expressions with \(z^2\)
replaced by \(z\). If \(a=0\), existing independence
in \(y\) forces \(b=d=0\). Otherwise solve for \(z^2\)
or \(z\), and substitute into \(B\):
\[
 B=\frac{b A^2}{a^2 y^6}
      +\left(\frac{d}{a}-\frac{2b}{a^2}\right)
                  \frac{A D^2}{y^{10}}
      +\left(1-\frac{d}{a}+\frac{b}{a^2}\right)
                  \frac{D^4}{y^{14}}.
\]
This cannot be independent of \(y\), as comparison of the
three distinct powers forces \(b=d=0\) and then leaves
a nonzero last term.

At \(s=3\) write the corrections as
\((a y^4z+b xyz^3,\ c y^2z^2+d x^2 y^2z,\ e yz^3)\).
In order, the coefficients of
\(x^5y^6z^2,x^4y^3z^5,x^4y^9,x^2y^9z,y^9z^2\)
in the determinant are
\[
            42bj-36e,\quad 6be,\quad
            14aj-4dj,\quad 10adj-8cj,\quad 6acj.
\]
Since \(j\ne0\) these force all five coefficients to vanish.

At \(s=2\) write the corrections as
\((a z^4+b x y^2 z^2,\ c z^3+dxy^2z+e x^3z^2,\
f y^2z^2)\).
The coefficients of \(x^3y^9\) and \(x^2y^7z^2\)
give \(d=c=0\). Put \(F=f/j\). Four more coefficients
give
\[
       7b-2e-6F=0,\quad
       56a+20be+4bF-32eF=0,\quad
       48ae+32aF-4beF=0,\quad aeF=0 .
\]
If \(a\ne0\), the last two give \(e=F=0\), a contradiction
from the first two. If \(a=0\), then \(beF=0\);
substituting each possibility in the first two again gives
\(b=e=F=0\).

At \(s=3/2\), the first two outputs on \(D=C/j\) are
\(D^2/y^4+a y z^3\), \(D^4/y^{14}+b y^2z\).
For \(b\ne0\) elimination of \(z\) makes the first
vary with \(y\) at fixed \(B\): the coefficient of
\(B^3\) is \(a/(b^3 y^5)\), forcing \(a=0\) and then
leaving \(D^2/y^4\). For \(b=0\) the second fixes \(y\),
so \(a=0\) directly.

At \(s=1\) put \(u=xy,\ t=z/u\).
The three outputs can be written
\(y^2u^2 M(t),\ u^4 N(t)/y^2,\ y^3 u L(t)\)
where \(M,N\) are quadratics or lower with constant
term one, and \(L\) is affine linear with constant term
\(j\). At fixed \(C\), the ratio \(C^6B^2/A^7\)
equals \(N^2 L^6/M^7\). It must be constant in \(t\),
since otherwise it recovers \(t\) algebraically, and
\(A\) then recovers \(y\) algebraically, giving rank
three. But this ratio being constant forces
\(\deg M\le1\) by degree comparison, and then \(M=1\)
by parity of root multiplicities. Thus \(N=1,L=j\).

### Binary ties

For each binary tie in this proof, the first further tie
on raising the unused \(z\)-weight has leading image
specializing isomorphically to that of the binary tuple.
Indeed the expected determinant weight remains positive:
the two fixed source weights sum to less than the sum of
any two output weights, and the first moving weight is
at most the largest output weight. Thus one can use the
same normalized integral lifts as in the preceding
axial weight arguments. For these binary ties the
initial \(z\)-weight is larger than either fixed-variable
weight, so a recovered polynomial lift of \(x\) or \(y\)
cannot involve \(z\).

At \(s=15/4\) or \(5/4\), the binary second entry would
be \(x^4 y^2+b y^m\) with \(b\ne0\) and \(m=5\) or \(3\),
respectively. The first and third entries give \(x,y^4\)
rationally. Over this field the second recovers \(y\),
since in the fourth-root extension with \(y^4\) fixed
the expression \(x^4 y^2+b y^m\) has trivial stabilizer.
Write \(D=C/j\). The equation \(B y^{14}=D^4+b y^{m+14}\)
makes \(y\) integral, so it lifts unchanged. The rational
lift of \(x=D/y^4\) cannot have a pole on \(y=0\), since
\(x^4y^2=B-b y^m\) is polynomial. Thus both lifts are
polynomial, contradicting the further tie.

At \(s=5/2\) the binary tuple is
\[
              (x^2 y^4+a y^5,\ x^4 y^2+d x^2 y^3+b y^4,\
                 j x y^4).
\]
Put \(t=x^2/y,\ \rho=y/x\). The output weights in
binary weights \(1,2\) are \(10,8,9\), and the
zero-weight ratios \(AB/C^2,\ C^{10}/A^9\) are
\[
       \frac{(t+a)(t^2+dt+b)}{j^2 t},
                      \qquad \frac{j^{10}t^5}{(t+a)^9}.
\]
For \(a\ne0\) they recover \(t\). A proper index
would be three, by the two degrees, hence the second
ratio would be a function of the first. The only pole
of the second ratio is at \(-a\), where the first
ratio is zero. The fibre of the first ratio there
would have to be a triple point, forcing
\(t^2+dt+b=(t+a)^2\). But then the poles of the
first ratio at zero and infinity have orders one
and two, whereas the second ratio has zeros there
of orders five and four. These orders cannot
come from a function of the first ratio.
If \(a=0\), the two ratios become a quadratic
and a nonzero constant times \(t^{-4}\). Their
field can have a proper index only when \(d=0\),
by the involution \(t\mapsto-t\).
Once \(t\) is recovered, \(A/C\) recovers \(\rho\),
so the binary field is full except possibly for
\(a=d=0\).

For \(a\ne0\) or \(b\ne0\) in the full-field case, \(y\)
is integral. Indeed with \(D=C/j\) the corresponding
identities are \(A y^4=D^2+a y^9\) if \(a\ne0\), or
\(B y^6=A^2+dA y^5+b y^{10}\) if \(a=0,b\ne0\).
Then the lift of \(y\) stays \(y\).
From \(A\) the rational lift of \(x\) can have poles
only at \(y=0\), and \(B\) precludes them: if \(x\)
has a pole, the \(x^4y^2\) term gives a strictly
lower, negative order than either correction.
Thus there is no tie.

The full-field case with \(a=b=0,d\ne0\) and the
tie at \(s=5/3\) can both use \(W=xy^2,V=xy\).
In both cases the field is full: when \(a=0\) the
ratio \(A/C\) gives \(x\), then \(C\) gives \(y^4\),
and \(B\) has a nonzero odd power of \(y\) in addition
to its even term, with trivial stabilizer under
the fourth roots of unity. In both cases \(W^2=A\)
and \(C V^2=jW^3\). For \(s=5/2\) the identity is
\(V^6+d W^3 V-BW^2=0\); for \(s=5/3\), where
\(B=x^4 y^2+dxy^3\), it is
\(V^6+dW^4-BW^2=0\). Thus \(V,W\) lift integrally.
They have equal prime support by these identities.
The lift of \(V\) cannot be prime, since the weights
of \(W,V\) are \(5,3\) or \(7,4\) after scaling the
\(x\)-weight to one. Its factors restrict to \(x,y\),
each of which has weight below the moving \(z\)-weight.
So \(V=xy\) and \(W=xy^2\) still, precluding a tie.

The only remaining case is \(a=d=0,\ b\ne0\) at
\(5/2\). Its image field is \(k(x,y^2)\) by the
involution above, or directly from \(A,C,B\).
Put \(H=y^2\), so \(A=x^2H^2,\ B=x^4H+bH^2,\ C=jxH^2\).
The identity \(B H^3=A^2+bH^5\) makes \(H\) integral.
In a lift the second identity for \(B\) forces the
rational \(x=A/(C/j)\) to be polynomial: the possible
lift of \(H\) is \(y^2\), \(y^2+\eta z\) at \(z\)-weight
four, or \(y^2+\eta xz\) at \(z\)-weight three.
At a root of \(H\) it has multiplicity at most two
and \(x^4 H\) cannot have a pole. Elsewhere the same
identity precludes a pole as well. The weights force
\(x\) to lift as \(x\). If \(H=y^2\), there is no tie.
If \(H=y^2+\eta z\), it is a coordinate and one can
raise the now unused \(y\)-weight. Each possible lift
of \(H\), of weight four, adds a term linear in \(y\),
either \(xy\) at weight three or \(y\) at weight four.
Take it again as the coordinate and continue. Beyond
four no lift can depend on \(y\), a contradiction.

Finally suppose \(H=y^2+\eta xz\), \(\eta\ne0\).
In parameters \(x,H\),
\[
               \det\frac{\partial(B,C)}{\partial(x,H)}
                          =j H^2(7x^4-2bH)=\sigma.
\]
At generic \(x=0\) the two differentials of \(B,C\)
are independent, since \(\sigma\) is a unit there
and \(x,H\) have independent differentials there.
A rational primitive for their Jacobian derivation
is \(-y/(\eta x\sigma)\). Every correction belongs
to \(k(x,H)\), since the source field is \(k(x,H)(y)\).
At \(x=0\) the pole order is one and its residue
contains \(y\) over the residue field \(k(H)=k(y^2)\).
It cannot be cancelled by an invariant. This
contradicts the divisor-regular leading-pair slice.

These exclusions leave the reference tuple for every
\(0<s<4\). Taking \(s=0\) shows the support of the
second output has \(5\deg_y+6\deg_z\le10\):
for every monomial, pass to the limit of the
weight bound \(4s+10\). Thus it is quadratic in \(y\)
and affine in \(z\), with no \(yz\)-term. Its
original ordinary coefficients of \(y^2,z\) were
\(x^4,x^5\) in degree six, respectively, and the
target shear changed the ordinary sextic only by
a multiple of \(x^6\). This proves the assertion.

## Residues and coordinate completion for the affine transverse output

**Lemma KellerQuadraticAffineOutputRootRestriction.**
Let a three-variable Keller map in characteristic zero have
an output \(S=v(x)z+q(x)y^2+u(x)y+w(x)\) with \(v,q\)
nonzero. Over the algebraic closure every root of \(v\)
is a root of \(q\). If \(v,u\) have a common root \(r\),
then \(v'(r)=q'(r)=u'(r)=0\) and \(w'(r)\ne0\).

**Corollary SexticKellerDoubledMixedCoprimeAffineOutputInversion.**
In the lower \(x\)-weight theorem the map is invertible
if \(\gcd(v,u)=1\). Consequently any noninvertible map
in that branch has a root \(r\) of \(v,u,q\) at which
\(v,q,u\) all vanish to order at least two and \(w'(r)\ne0\).

### The two branches over a root of the transverse coefficient

Work on a generic level \(S=t\), extending constants
to the algebraic closure of \(k(t)\) as needed.
The Keller identity gives a regular one-form on this level
whose exterior derivative equals the relative area form.
Indeed the other two outputs give a pair \(G_1,G_2\)
with \(dG_1\wedge dG_2\) a nonzero constant multiple of
that area form; take the same multiple of \(G_1\,dG_2\)
with the inverse normalization. This is the residue
argument used in
[Wang–Yan, Lemma 3.14](https://arxiv.org/pdf/2610.36495)
for an affine-variable component.

On \(v\ne0\) the level ring is \(K[x,v^{-1},y]\), with
\(K\) the extended constant field. Up to the choice of
sign, its area form is \(dx\wedge dy/v(x)\), with primitive
\(-y\,dx/v(x)\). Every closed regular one-form on this
open set is an exact form plus \(h(x)dx\), with
\(h\in K[x,v^{-1}]\). This follows by integrating the
\(dy\) coefficient in the polynomial variable \(y\),
and using closedness to see that the remaining
\(dx\) coefficient is independent of \(y\).

Suppose \(v(r)=0\) but \(q(r)\ne0\), and let \(k\) be
the multiplicity of \(r\) in \(v\). With \(z\) fixed
constant the equation of the generic level gives
two regular formal branches \(y_\pm(x)\) over \(x=r\),
since the generic quadratic has distinct roots there.
Regularity of the global primitive implies that the
residues of \(y_\pm(x)dx/v\) on these branches are equal:
the exact term has zero residue, and \(h(x)dx\) has
the same residue on both. Write
\(G=w-u^2/(4q)\). The difference \(y_+-y_-\) is
\[
             2q^{-1/2}\bigl(t-G(x)-v(x)z\bigr)^{1/2}.
\]
The residue uses only its terms through order \(k-1\)
in \(x-r\), so the \(vz\) term does not contribute.
For \(k=1\) the residue difference is nonzero.
For \(k\ge2\), if \(G'(r)\ne0\), the coefficient
of the lowest power of \(t-G(r)\) in the residue
difference, namely its power \(1/2-(k-1)\), is nonzero.
It is contributed only by taking the linear term
\(-G'(r)(x-r)\) a total of \(k-1\) times in the
binomial expansion and multiplying by the nonzero
constant term of \(2q^{-1/2}(x-r)^k/v\).
All higher powers of \(t-G(r)\) cannot cancel this
one. Therefore equal residues force \(G'(r)=0\).
But then \(x=r,\ y=-u(r)/(2q(r))\) is a critical point
of \(S\), since \(v'(r)=0\) and \(S_x=G'(r)\) there.
Both cases are impossible in a Keller map.

If \(r\) is also a root of \(u\), then \(S_y=S_z=0\)
on all of \(x=r\). The polynomial
\(S_x=q'(r)y^2+u'(r)y+v'(r)z+w'(r)\) must have no
zero there over the algebraic closure. Its nonconstant
coefficients must therefore vanish and its constant
coefficient must be nonzero, as asserted.

### The coprime coordinate completion

If \(\gcd(v,u)=1\), the substitution \(y\mapsto
q y^2+u y+w\) is invertible over \(k[x]/(v)\).
The coefficient \(u\) is a unit there and \(q^5=0\)
since every root of the degree-five polynomial \(v\)
is a root of \(q\). The inverse is polynomial of
degree at most five in its argument \(t\). For example
solve \(y=u^{-1}(t-w)-u^{-1}q y^2\) as a formal
series in the nilpotent \(q\); the coefficient of
\(q^i\) has degree at most \(i+1\) in \(t\), and
the series stops after \(i=4\). Uniqueness follows
by subtracting two solutions, since the multiplier
\(u+q(y_1+y_2)\) is a unit.
Lift the coefficients of this inverse to degree at
most four in \(x\), obtaining \(L(x,t)\) of total
degree at most nine.

Then
\[
              (x,y,z)\longmapsto
                    \left(x,S,\frac{y-L(x,S)}{v(x)}\right)
\]
is a coordinate completion. The last entry is polynomial
by the inverse modulo \(v\). With output coordinates
\((x,t,h)\) it has inverse
\[
           y=L(x,t)+v(x)h,\qquad
           z=\frac{t-q(x)y^2-u(x)y-w(x)}{v(x)}.
\]
The numerator is divisible by \(v\) by the same inverse
identity. The inverse degrees in \(y,z\) are at most
nine and seventeen, respectively: before dividing
by the degree-five \(v\), the bounds are
\(4+2\cdot9,\ 5+9,\ 6,\ 1\) in the numerator of the
third coordinate.

By the lower \(x\)-weight theorem every monomial
\(x^i y^h z^k\) of either remaining output has
\(5h+6k\le20\) and \(i+h+k\le6\). Under the completion
its degree is at most \(i+9h+17k\le54\).
The new map fixes the coordinate \(t\) and is still
Keller. On each fixed \(t\)-plane the remaining two
outputs are invertible by the retained planar
degree-100 theorem of Moh. The
[planar specialization and descent argument](LOW_DEGREE_COORDINATE_KELLER_MAPS_IN_THREE_VARIABLES.md)
then proves invertibility over the original field
using this actual composed degree bound.

## Exclusion of the mixed doubled frame

**Theorem SexticKellerDoubledSquareMixedExclusion.**
In characteristic zero, no three-variable Keller map of
ordinary degree at most six in the mixed doubled square-member
frame of the lower \(x\)-weight theorem has the highest tuple
\[
             (x^2P^2,\ x^4P,\ x^6+jxP^2),
             \qquad P=xz+y^2,\quad j\ne0
\]
in weights \((4,5,6)\).

Work over the algebraic closure. After the shear of the lower
\(x\)-weight theorem, write \(F_2=qy^2+uy+vz+w\), with \(q,v\)
monic of degrees four and five. The same theorem, by taking
the limit in weights \((s,5,6)\) as \(s\) decreases to zero,
bounds each of \(F_1,F_3\) by weight 20 in \((0,5,6)\).
Their parts of weight 20 are \(p(x)y^4,r(x)y^4\), with \(p\)
monic quadratic and \(r\) linear with leading coefficient \(j\).

### Increasing the transverse weight

Increase the weight of \(z\) from six to ten, keeping the other
two weights fixed. The only possible ties for \(F_1,F_3\)
before ten are \(z^3\) at \(20/3\) and \(yz^2\) at \(15/2\),
with coefficients in \(k[x]\). This follows by solving
\(5h+Wl=20\) for \(6<W<10\) with \(5h+6l\le20\) and \(l>0\).
At each of these weights the expected determinant weight is
\(45-W>0\), so the highest tuple has rank at most two.
On the generic level \(q y^2=s\) the two other highest
outputs are
\[
             s^2p/q^2+a(x)\zeta,\qquad
             s^2r/q^2+b(x)\zeta ,
\]
with \(\zeta=z^3\) or \(yz^2\), respectively. Here \(x,\zeta\)
are algebraically independent over \(k(s)\), and the
Jacobian in them must vanish. Its \(\zeta\) coefficient
would force \(a/b\) constant if both were nonzero. Subtracting
that constant multiple would then make a nontrivial constant
linear combination of \(p/q^2,r/q^2\) constant in \(x\).
If just one of \(a,b\) is nonzero, the same Jacobian forces
the corresponding other rational function constant. Both
possibilities contradict their distinct positive orders at
infinity, six and seven. This eliminates both ties.

At weight ten put \(U=y^2,V=z\). The highest first and third
outputs are binary quadratics in \(U,V\), and the second is
\(qU+vV\). On \(qU+vV=s\) put \(H=U/s\) and substitute
\(V=s(1-qH)/v\). Dividing the other two outputs by \(s^2\)
gives a pair of polynomials in \(H\) of degree at most two
over \(k(x)\). Their joint image over \(k(s)\) has dimension
at most one by the rank condition. If either depends on
\(H\), the image is a curve, defined by a single irreducible
polynomial over \(k(s)\). That polynomial remains irreducible
over the purely transcendental extension \(k(s,x)\) and
defines the image there as well. In fact the image there
still has dimension one, and the extended relation ideal
is a prime of height one contained in its relation ideal.
The six monomials of total degree at most two in the two
outputs give polynomials in \(H\) of degree at most four.
They are linearly dependent over \(k(s,x)\), so the defining
equation has degree at most two.

At \(H=1/q\) the two outputs are \(p/q^2,r/q^2\).
Their six monomials of degree at most two have distinct
orders \(0,6,7,12,13,14\) at \(x=\infty\). They are
linearly independent over \(k(s)\), a contradiction.
Thus both polynomials are independent of \(H\). Undoing
the substitution gives highest polynomials
\[
              p(y^2+(v/q)z)^2,\qquad
              r(y^2+(v/q)z)^2 .
\]
In particular \(r(v/q)^2\in k[x]\). A pole of \(v/q\) at
a finite point would require \(r\) to vanish there to order
at least two, which it cannot. Thus \(v=qT\) with \(T\)
monic linear. Write \(P_T=y^2+Tz\). The terms of \(F_1-pP_T^2\)
and \(F_3-rP_T^2\) have weight at most three in \((0,1,2)\).

### The bracket on \(F_2=s\)

On \(v\ne0\) the level \(F_2=s\) has coordinates \(x,y\)
with \(z=(s-w-uy-qy^2)/v\). Write \(A,C\) for the
restrictions of \(F_1,F_3\). The chain rule gives
\[
                A_x C_y-A_y C_x=\kappa/v(x)
\]
for some \(\kappa\in k^\times\), with derivatives at fixed
\(s\). Since \(P_T=(s-w-uy)/q\), each restriction is of the
form
\[
       a_3(x)y^3+a_2(x)y^2+
       (b_0(x)+b_1(x)s)y+d_0(x)+d_1(x)s+d_2(x)s^2.
\]
Indeed the possible remainder monomials in \(y,z\) are
\(y^3,yz,y^2,z,y,1\). The cubic coefficient can have a
finite pole only at \(T=0\), and at most a simple pole.
Only \(yz\) contributes such a pole, through
\(z=(s-w-uy)/v-y^2/T\).
The coefficients \(d_2\) of \(A,C\) are \(p/q^2,r/q^2\),
respectively. Each \(b_1\) tends to zero at \(x=\infty\):
for \(A\) it is \(-2pu/q^2\) plus the coefficient of
\(yz\) divided by \(v\), and for \(C\) use \(r\) in place
of \(p\). The degree bounds \(\deg u\le5\) and
\(\deg [yz]F_i\le4\) give the assertion.

These properties persist under constant invertible
linear combinations of the pair. For each nonzero
combination \(d_2\) is nonzero and tends to zero at
infinity. If both outputs have the same \(y\)-degree
\(m\ge2\), the coefficient of \(y^{2m-1}\) in the bracket
makes their leading coefficients proportional over \(k\).
They belong to \(k(x)\), and the vanishing coefficient is
\(m\) times their derivative determinant. A constant shear
therefore lowers one of the degrees. If one degree is zero,
the bracket forces the other to be one. It remains to
consider degrees \((3,2)\), degrees \((m,1)\) for
\(m=2,3\), and the case of both degrees at most one.

### Degrees three and two

Call the cubic output \(G\) and the quadratic \(D\), with
leading coefficients \(a,b\in k(x)^\times\). The coefficient
of \(y^4\) in their bracket gives \(2a'b-3ab'=0\), so
\(a^2/b^3\in k^\times\). Since \(k\) is algebraically closed,
write \(b=h^2,a=\eta h^3\) with \(h\in k(x)^\times\)
and \(\eta\in k^\times\). The pole bound on \(a\) forces
\(h\in k[x]\), since any finite pole of \(h^3\) would have
order at least three. Set
\[
        Y=hy+(b_0+b_1s)/(2h),\qquad D=Y^2+d(x,s).
\]
In \(x,Y\) coordinates write
\[
                 G=\eta Y^3+c(x,s)Y^2+e(x,s)Y+f(x,s).
\]
All four coefficients \(c,d,e,f\) are polynomial in \(s\)
over \(k(x)\). The new bracket is a nonzero constant
divided by \(vh\). Its coefficients at \(Y^3,Y^2,Y\)
give \(c_x=0,\ 2e_x=3\eta d_x,\ f_x=c d_x\).
Consequently \(c\in k[s]\), \(e-(3\eta/2)d\in k[s]\)
and \(f-cd\in k[s]\). The constant term of the bracket,
\(-e d_x\), is independent of \(s\) and nonzero. Thus
\((e^2)_x=3\eta e d_x\) is independent of \(s\) and nonzero.
If \(e\) had positive degree \(n\) in \(s\), comparing the
coefficients of \(e^2\) from \(s^{2n}\) down through \(s^n\)
would force all coefficients of \(e\) into \(k\). Indeed
the leading coefficient has constant square, and each
successive coefficient occurs linearly times twice the
leading coefficient, with all other summands already
constant. This would give \((e^2)_x=0\). Hence \(e\)
is independent of \(s\), and the \(s^2\) coefficient of
\(d\) is constant in \(x\).

Put \(L=b_1/(2h)\) for the quadratic output. This tends to
zero at \(x=\infty\). The coefficient of \(sy^2\) in
\[
                  G=\eta Y^3+c(s)D+e(x)Y+(f-cd)(s)
\]
is \(h^2(3\eta L+[s]c)\). It vanishes by the form of \(G\)
before the change of variable. Therefore \(L\) is constant
and is zero. But the \(s^2\) coefficient of \(d\) is
\(d_2-L^2\), which would now make the nonconstant \(d_2\)
of \(D\) constant. This is impossible.

### One degree at least two and one degree one

Write the linear output as \(D=b(x,s)y+d(x,s)\) and the
other as \(G\) of degree \(m=2\) or \(3\). Use \(x,D\) as
coordinates in the function field over \(k(s)\). The
bracket gives \(b\,\partial_x G|_D=\text{constant}/v\).
In particular the constant coefficient of \(G\) as a
polynomial in \(D\) is a rational primitive in \(x\)
of a nonzero constant divided by \(vb\).

Write \(b=b_0+b_1s\). If \(b_1\ne0\) and \(b_0/b_1\)
is nonconstant, the generic equation \(b_0/b_1=-s\)
has simple finite solutions away from all fixed zeros
and poles of \(v,b_0,b_1\). This follows from separability
in characteristic zero and the transcendence of \(s\).
The differential \(dx/(vb)\) has a nonzero residue at
each such point in the algebraic closure of \(k(s)\), so
it cannot have a rational primitive. If \(b_1\ne0\) and
\(b_0/b_1\) is constant, specialize \(s\) to its negative.
On this level \(D\) is independent of \(y\), while \(G\)
still has degree \(m\), since its leading coefficient
is independent of \(s\). The bracket cannot be a nonzero
constant in \(y\) on that level. The specialized identity
is valid because the original restrictions only use
denominators in \(x\). Thus \(b=b_0(x)\ne0\).

In \(x,D\) coordinates the positive-degree coefficients of
\(G\) are now independent of \(x\), by the bracket. They
are polynomial in \(s\), since the substitution for \(y\)
only divides by \(b(x)\). The coefficient \(H_m\) of
\(D^m\) is a nonzero constant in \(k\) by the form of \(G\)
before substitution. The coefficient of \(y^{m-1}\) in
\(G\) is
\[
                    b^{m-1}\bigl(mH_m d+H_{m-1}(s)\bigr).
\]
It has degree at most one in \(s\) in the original form.
Comparing the \(s^2\) coefficient again forces the
nonconstant \(d_2\) of \(D\) to be constant.

### Both degrees at most one

Use the original \(A,C\), since the constant linear changes
do not alter the condition that both degrees are at most one.
Let their \(y\) coefficients be \(\rho_0+\rho_1s\) and
\(\sigma_0+\sigma_1s\). The coefficient of \(ys^2\) in
the bracket gives \(\rho_1'\sigma_1-\rho_1\sigma_1'=0\).
A nonzero pair \((\rho_1,\sigma_1)\) would therefore be
a rational multiple of a nonzero constant vector. But the
coefficient of \(s^3\) independent of \(y\) gives
\(\sigma_1(p/q^2)'-\rho_1(r/q^2)'=0\). The two derivatives
are linearly independent over \(k\) by their distinct
orders at infinity. Thus \(\rho_1=\sigma_1=0\).
Repeating the same comparison at \(y\) and \(s^2\) gives
\(\rho_0=\sigma_0=0\), contradicting the nonzero bracket.
This completes the exclusion.

## Exclusion of the binary tie in the smooth-common-factor frame

**Theorem SexticKellerTangentConicSmoothSquareMemberBinaryExclusion.**
Suppose a three-variable characteristic-zero Keller map of
ordinary degree at most six has ordinary highest tuple
\((P^3,P^2x^2,Px^4)\), where \(P=xz+y^2\). If the highest
tuple in weights \((1,1,2)\) is binary in \(x,z\), it is
\[
                         (x^3z^3,\ x^4z^2,\ x^5z).
\]
In this case the highest tuple at weights \((2,3,4)\)
is \((P^3,P^2x^2,Px^4)\).

Work over the algebraic closure. With \(x,z\) held at
weights one and two, the binary tuple from the weight-two
theorem is
\[
  A_0=x^3z^3+b_0xz^4,\quad
  B_0=x^4z^2+e_0x^2z^3+f_0z^4,\quad
  C_0=x^5z+j_0x^3z^2+\ell_0xz^3.
\]
Increase the weight \(w\) of \(y\) from one to \(3/2\).
The possible ties strictly before \(3/2\) are
\[
 \begin{array}{c|c|c}
 w&B&C\\ \hline
 5/4&& y^4z\\
 4/3&y^3z^2&xy^3z\\
 7/5&&y^5 .
 \end{array}
\]
For the table, test \(i+w h+2l=d\) at \(h>0\)
for \(d=9,8,7\), using the ordinary degree bound and
the specified terms in degree six. At each tie the
leading determinant vanishes since \(9+8+7-(3+w)>0\).

The three binary entries are pairwise algebraically
independent. In the variables \(x,t=z/x^2\) their weights
are \(9,8,7\), and the orders at \(t=0\) of their coefficient
functions are \(3,2,1\). For each pair the ratio with
the powers of \(x\) cancelled thus has nonzero \(t\)-order.
The isolated ties in \(C\) are impossible by independence of
\(A_0,B_0\). At \(4/3\), fix \(A_0\) at a generic value
and use \(y^3\) as a transcendental variable over the
function field of the binary level. The two other entries
are affine in \(y^3\). If both coefficients are nonzero,
the coefficient of \(y^3\) in their relative determinant
forces the ratio of the two coefficients to be constant
along the binary level. This ratio is proportional to
\(z/x\), which is algebraically independent of \(A_0\).
Indeed at fixed \(z/x\) the polynomial \(A_0\) has
nonzero terms of distinct \(x\)-degrees six and possibly
five. If only one coefficient is nonzero, the unchanged
entry would be constant along the level of \(A_0\), again
impossible. This eliminates every tie before \(3/2\).

### The highest tuple at weights \((2,3,4)\)

At \(w=3/2\) all the tying terms have even \(y\)-degree.
Use \(x,P,z\) as polynomial generators of \(k[x,y^2,z]\).
In these variables the general tuple is
\[
 \begin{aligned}
 A&=P^3+aPz^3+b xz^4,\\
 B&=x^2P^2+cP^2z+d xPz^2+e x^2z^3+f z^4,\\
 C&=x^4P+i xP^2+h x^2Pz+g Pz^2+j x^3z^2+\ell xz^3 .
 \end{aligned}
\]
This support follows by the same exponent test at
\(w=3/2\) and the substitution \(y^2=P-xz\). The
determinant in \(x,P,z\) vanishes because the determinant
in \(x,y,z\) is its product with \(2y\).
The following coefficients determine the tuple:
\[
\begin{array}{c|r}
 x^3P^5 & 6(2c-h)\\
 P^6 &3ci\\
 x^4P^4 z &12(2d-j)\\
 xP^5 z &6(ch+di-2g)\\
 x^5P^3 z^2 &-18(a-2e)\\
 x^6P^2 z^3 &-24b\\
 x^5P z^5 &2(9ae-7bd-10bj).
\end{array}
\]
In particular \(b=0\), \(a=2e\), and \(ae=0\), giving
\(a=b=e=0\). At these values the coefficient of \(x^3P^3 z^3\)
is \(12(dj+4f)\) and that of \(x^2P^2 z^5\) is \(36fj\).
Together with \(j=2d\) these give \(d=f=j=0\).
The coefficient of \(x^2P^4 z^2\) is now \(-18\ell\),
giving \(\ell=0\). The remaining relations give
\(h=2c,\ g=c^2,\ ci=0\). Thus, with \(Q=x^2+cz\),
the tuple is
\[
                       (P^3,\ P^2 Q,\ P Q^2+i xP^2),
                       \qquad ci=0.
\]

### The two relative differentials

If \(c\ne0\), apply the weighted leading-pair slice to
\(P^3,P^2Q\). Their differential wedge is \(3P^4\)
times \(dP\wedge dQ\), so on generic \(P=S,Q=T\) this
would give a rational primitive of the relative form
for \(P,Q\). Eliminating \(z\) gives
\[
                          y^2=S+(x^3-Tx)/c .
\]
The cubic is squarefree for generic \(S,T\).
As \(dP\wedge dQ\wedge dx=2cy\,dx\wedge dy\wedge dz\),
the relative form is proportional to \(dx/y\).
It is nonzero and regular on the complete normalization.
At a finite ramification point \(y\) is a parameter;
at infinity \(x,y\) have pole orders two and three.
Hence this form has no rational primitive, a contradiction.

Suppose \(c=0\) and \(i\ne0\). Use the first and third
highest entries \(P^3,G=x^4P+i xP^2\). At the generic
point of \(x=0\) their differentials are independent,
as \(P=y^2\) there and \(G_x\) at fixed \(P\) is
\(iP^2\) there. The weighted leading-pair slice now
gives a rational slice regular at that divisor.
But in the source field \(k(x,P)(y)\),
\[
                    J(P^3,G,y)=3P^3x(4x^3+iP).
\]
The derivation \(J(P^3,G,\cdot)\) annihilates \(k(x,P)\)
and is a nonzero constant multiple of \(\partial_y\)
over that field. Any slice with nonzero scalar value
\(\kappa\) is consequently of the form
\[
               \frac{\kappa y}{3P^3x(4x^3+iP)}+H(x,P),
                         \qquad H\in k(x,P).
\]
The first term has a simple pole at \(x=0\), with
coefficient \(\kappa y/(3iP^4)\) relative to the
parameter \(x\). This coefficient is outside \(k(P)\)
in the residue field, since \(P=y^2\) there.
The divisor valuation on \(k(x,P)\) is the \(x\)-adic
valuation with residue field \(k(P)\). Thus \(H\) cannot
cancel the pole. This contradicts the regular slice.

It follows that \(c=i=0\). Specializing the resulting
highest tuple to \(y=0\) gives the three uncorrected binary
entries at \(w=1\), proving the assertion.

## The axial weights five-halves and three in the smooth frame

**Theorem SexticKellerTangentConicSmoothSquareMemberWeightThreeReduction.**
In the binary case of the smooth-common-factor weight-two
theorem, the highest tuple in \((1,1,w)\) stays
\((x^3z^3,x^4z^2,x^5z)\) for \(2<w<3\).
At \(w=3\) it is
\[
                          (x^3z^3,\ x^4z^2,\ x^5z+i x^2z^2).
\]
In weights \((1,2,3)\) the highest tuple is
\((P^3,P^2x^2,Px^4+iP^2)\), with the same \(i\).

By the binary exclusion at weight two, all lower ordinary
terms have weights strictly below \(18,16,14\), respectively,
in \((2,3,4)\). As \(w\) increases above two, a lower ordinary
monomial \(x^r y^s z^l\) ties the reference monomial
\(x^{r_0}z^{l_0}\) only for \(l>l_0\) and
\(w=(r_0-r-s)/(l-l_0)\). The bounds \(r+s+l\le5\) and
\(2<w\le3\) give a single tie before three, namely \(z^3\)
in the third entry at \(5/2\).

### Excluding the tie at five-halves

Hold \(x,z\) at weights \(1,5/2\) and increase the weight of
\(y\) from one to \(7/4\). Test \(r+w_y s+5l/2=d\) for
\(d=21/2,9,15/2\), using the strict lower-term bounds at
\((2,3,4)\). The possible terms before \(7/4\) are \(yz^3\)
in the second entry and \(xyz^2\) in the third at \(3/2\),
and \(y^3z\) in the third at \(5/3\). At \(7/4\) the only
tying \(y\)-terms are the ordinary terms completing \(P^3\),
\(P^2x^2\) and \(Px^4\). Both lists are checked by the
[retained script](verify_sextic_smooth_square_weight_three.py).

At the first tie, on a geometric generic level of
\(A_0=(xz)^3\) the function \(xz\) is a nonzero scalar.
The two binary entries \(x^4z^2\) and \(x^5z+\delta z^3\)
are both nonconstant in \(x\) there. The ratio of the two
potential \(y\)-coefficients is proportional to \(z/x\),
also nonconstant. The coefficients of the relative
determinant in \(y\) therefore exclude both terms, by the
same affine calculation as for the weight-two binary tie.
Independence of the first two entries eliminates the
isolated tie at \(5/3\). The highest tuple at \(7/4\)
is consequently \((P^3,P^2x^2,Px^4+\delta z^3)\).
Its determinant is nonzero unless \(\delta=0\), since the
first two entries are independent in \(x,P\) and \(z\)
is independent of \(x,P\). This excludes the tie at \(5/2\).

### The highest tuple at axial weight three

The axial exponent test gives the possible corrections to
\((x^3z^3,x^4z^2,x^5z)\) at \(w=3\) as
\[
                    (a z^4,\ (d x+u y)z^3,\
                                 (i x^2+v xy+\lambda y^2)z^2).
\]
In their determinant the coefficients at \(x^7z^6\),
\(x^6yz^6\), and \(x^3yz^7\) are
\(6(v-2u),12\lambda,2(16a\lambda-6d\lambda-3uv)\).
Thus \(u=v=\lambda=0\).

Now increase the weight of \(y\) to two, keeping \(x,z\)
at weights one and three. The same test gives only one
intermediate tie, \(y^3z\) in the third entry at \(5/3\).
The first two binary entries are algebraically independent:
with \(t=z/x^3\) the ratio cancelling their powers of \(x\)
has nonzero order at \(t=0\), since \(3/12\ne2/10\).
The isolated tie is therefore impossible. At \(w_y=2\)
the allowed polynomials are exactly
\[
 \begin{aligned}
 A&=P^3+a z^4+b y^3z^2+c xy z^3,\\
 B&=x^2P^2+d xz^3+e y^2z^2+f y^5+g xy^3z+h x^2y z^2,\\
 C&=x^4P+i x^2z^2+j y z^2+k y^4+\ell xy^2z
                  +m x^2y^3+n x^3yz.
 \end{aligned}
\]
Their determinant must vanish as its expected weight is
\(12+10+8-6=24>0\).

### The coefficient calculation over the rationals

Let \(I\subset\mathbb Q[a,b,c,d,e,f,g,h,i,j,k,\ell,m,n]\)
be the ideal of all coefficients of this determinant in
\(x,y,z\). The following entries or powers are in \(I\):
\[
 \begin{gathered}
 b-2h+n,\quad c-2h+n,\quad 2f-m,\quad 2g-m-n,\\
 (m-n)^4,\quad j^3,\quad (i+k-\ell)^3,\quad (2h-n)^2,\quad
 e^3,\quad d^3,\quad a^2,\quad (n^2-8k+4\ell)^3 .
 \end{gathered}
\]
The retained script constructs the determinant and its
coefficient ideal over \(\mathbb Q\), computes polynomial
witnesses for these memberships, and verifies the matrix
product recovering all twelve expressions from the ideal
generators. Consequently each displayed expression
vanishes at every solution over a field of characteristic
zero, and the same is true of the bases of the powers.
The script requires Python 3 and Singular on the path.
With Singular 4.3.2, run
\[
 \texttt{python3 research/jacobian/monodromy-theory/verify\_sextic\_smooth\_square\_weight\_three.py}
\]
from the repository root. The success line records both
tie tables and the twelve coefficient identities.

Taking \(Q=x^2+h y\), the relations give
\[
                              (A,B,C)=(P^3,P^2Q,PQ^2+iP^2).
\]
For example \(m=n=2h\), \(k=h^2+i\) and
\(\ell=h^2+2i\); all of \(a,b,c,d,e,j\) vanish.
This also verifies that the resulting family satisfies
the highest determinant equation.

If \(h\ne0\), the first two outputs again give the
relative differential for \(P,Q\) up to a nonzero scalar
on each geometric generic component. On the level
\(P=S,Q=T\), solve \(y=(T-x^2)/h\) and \(z=(S-y^2)/x\).
Since \(dP\wedge dQ\wedge dx=-hx\,dx\wedge dy\wedge dz\),
the relative differential is proportional to \(dx/x\)
on this rational curve. Its nonzero residue at \(x=0\)
contradicts the leading-pair slice. Thus \(h=0\), giving
the asserted tuple at \((1,2,3)\). Specialization at \(y=0\)
gives the asserted weight-three axial tuple.

## Complete axial degrees when the remaining coefficient is nonzero

**Corollary SexticKellerTangentConicSmoothSquareMemberNonzeroAxialCoefficient.**
Under the weight-three reduction, if \(i\ne0\), the
highest tuple in weights \((1,1,w)\) for every \(w>3\)
is \((x^3z^3,x^4z^2,i x^2z^2)\). In particular the
three degrees in \(z\) are \(3,2,2\).

The weight-three conclusion bounds the degree in \(z\) of
the first and third entries by three and two, respectively.
Indeed any term of higher \(z\)-degree already exceeded
their bounds at weight three or had a coefficient forced
to zero there. In the second entry the only possible term
of higher \(z\)-degree is \(\delta z^3\). It ties at \(w=4\).

To exclude this coefficient, hold \(x,z\) at weights one
and four, and raise the weight of \(y\) to \(5/2\). The
reference weights are \(15,12,10\). In the first two
entries, a monomial of ordinary degree at most five has
weight in \((1,2,3)\) strictly below \(12,10\), respectively.
In the third entry the bound is strictly below eight
apart from \(iP^2\). If a lower ordinary monomial
\(x^r y^s z^l\) with \(s>0\) ties for \(1<w_y\le5/2\),
its weight at \((1,5/2,4)\) must reach the corresponding
reference weight. No such monomial does. For example,
writing \(b=r+2s+3l\), its weight at \((1,5/2,4)\) is
\(b+s/2+l\). The integer bounds \(b\le11,9,7\) give \(s/2+l\le b/3\) from \(2s+3l\le b\).
Thus \(b+b/3\) is strictly below \(15,12,10\) for
\(b\le11,8,7\), respectively. For
the second entry at \(b=9\), equality requires \(s=0\),
which accounts only for the already included \(z^3\).
The ordinary terms completing \(P^3,P^2x^2\) and the
correction \(iP^2\) tie exactly at \(w_y=5/2\). Thus
the highest tuple in weights \((2,5,8)\) is
\((P^3,P^2x^2+\delta z^3,iP^2)\).

Suppose \(\delta\ne0\). Fix generic \(P=S\) and
\(P^2x^2+\delta z^3=T\). On this level put \(U=Sx\);
the curve tower is
\[
                       U^2=T-\delta z^3,\qquad
                       y^2=(S^2-Uz)/S .
\]
The second component and \(P^3\) are independent.
Their relative differential on the generic curve is
a nonzero scalar multiple of \(dz/(xy)\), by the minor in
\(x,y\) of \(P\) and \(P^2x^2+\delta z^3\).
The first double cover is the smooth model of
\(U^2=T-\delta z^3\). The function \(Uz\) has a single pole
of order five there, so the second cover is branched at
infinity. Its finite branch points are simple since
\(S^2\) is a generic value of \(Uz\) on the first curve.
At these finite branch points \(x\) is nonzero, and the
ramification cancels the factor \(y\) in the denominator.
At \(x=0\), the first ramification cancels \(x\), and
\(y\) is a unit. At infinity on the final curve the pole
orders of \(z,x,y\) are \(4,6,5\), respectively. Thus the
form is nonzero and regular on the complete normalization.
It cannot have a rational primitive, contrary to the
weighted leading-pair slice. Hence \(\delta=0\), and the
three complete highest \(z\)-terms follow.

## Relative differentials at the next transverse weight

**Lemma KellerConicQuadraticFibreDifferential.**
Let \(k\) have characteristic zero, let \(p,q,r\in k[x]\)
have degrees \(2,\le2,\le3\), respectively, with zero
polynomials permitted for \(q,r\), and suppose \(p\) is monic.
Put \(U=y^2+xz\) and \(B=pU^2+qUz+rz^2\).
If \(q,r\) are not both zero, the relative Jacobian
differential of \(U,B\) has no rational primitive on at
least one geometric component of its generic fibre.

Work over the algebraic closure of the function field of
the generic values \(U=S,B=T\). Write \(D=pS^2-T\).
Eliminating \(U\) leaves
\[
                         rz^2+qSz+D=0,\qquad y^2=S-xz .
\]
The relative form, up to sign, is
\(\omega=dx/(2y(2rz+qS))\). This follows by taking the
minor in \(y,z\) at fixed \(x\). The expressions below
are in the function fields of the components and the
orders refer to their smooth projective models. The
parameters \(S,T\) are independent, so they avoid all
fixed critical values of the rational functions used
in the argument.

### The case \(r=0\)

Here \(q\ne0\). With \(W=qy\) the model is
\[
                  W^2=qN/S,\qquad N=S^2(q+xp)-Tx ,
                  \qquad \omega=dx/(2S W).
\]
The cubic \(N\) has simple roots: outside zero they
come from a generic level of the nonconstant rational
function \((q+xp)/x\), and at zero any root is simple
by the coefficient of \(T\). A nonzero root of \(q\)
is never a root of \(N\). If \(q\) has a double root
outside zero, or a simple root at zero, then \(qN\)
has a double root there and \(\omega\) has a simple
pole of nonzero residue. If \(q\) is squarefree and
prime to \(N\), the squarefree polynomial \(qN\)
has degree at least three, and \(dx/W\) is nonzero
and regular on the complete normalization. The only
case left is \(q=\alpha x^2\) with \(\alpha\ne0\).
Then \(N=xN_1\) and \(V=W/x\) satisfies a squarefree
cubic equation \(V^2=\alpha xN_1/S\). The form
\(dx/(xV)\) has a double pole at the single point over
zero and is regular elsewhere. A primitive would
have just one simple pole, which would give a degree-one
map of this elliptic curve to the projective line.
This is impossible.

### Finite points when \(r\ne0\)

Put \(W=2rz+qS\), so
\[
             W^2=\Delta=q^2S^2-4rD,\qquad
             z=(W-qS)/(2r),\qquad \omega=dx/(2yW).
\]
Away from the zeros of \(r\), every zero of \(\Delta\)
is simple because it comes from a generic level of
\(p-q^2/(4r)\). A constant such function contributes
no zeros. At a simple zero the ramification cancels
\(W\) in the relative differential. If \(y^2\) also
vanishes there, \(x/(2r)\) is a unit and the \(W\)-term
in \(S-xz\) makes that zero simple on the first cover.
The second ramification then cancels \(y\). At other
finite points with \(rW\ne0\) a zero of \(y^2\) is
simple since it comes from a generic level of
\(p+q/x+r/x^2\). Indeed the condition \(z=S/x\)
in the quadratic equation gives precisely that level,
and the quadratic has nonzero derivative in \(z\)
there. A constant rational function again gives
no such zeros.

At a root of \(r\) where \(q\) is nonzero, \(W\) is a
unit. The two branches of \(z\) are one finite branch
and one with a pole of order \(\operatorname{ord}r\).
Outside zero these give a nonzero finite value or a
pole of \(y^2\), respectively, for generic \(T\).
At zero the finite branch gives \(y^2=S\) in the
residue field. On the pole branch a zero of \(y^2\)
can occur only when \(r\) has order one and
\(q(0)=-r'(0)\). In that case write the pole branch
as \(-qS/r-z_{\mathrm{fin}}\). The coefficient of
\(x\) in \(S+xqS/r+xz_{\mathrm{fin}}\) depends
nontrivially on \(T\) through
\(z_{\mathrm{fin}}(0)=(T-p(0)S^2)/(q(0)S)\). The
zero is simple and the differential is regular.

At a common root \(a\) of \(r,q\), put
\(d=\operatorname{ord}_a r\) and \(e=\operatorname{ord}_a q\),
allowing \(e=\infty\). Measure orders relative to \(x-a\),
so that fractional orders are allowed. The order of
\(W\) is \(\min(d,2e)/2\), as \(D(a)\ne0\) generically.
If \(d\le2e\), both roots \(z\) have order \(-d/2\).
For \(d=2e\), the same holds since \(D(a)\ne0\)
precludes cancellation in either leading root.
If \(d>2e\), necessarily \(d=3,e=1\) and the two
orders of \(z\) are \(-1,-2\), by balancing the linear
and constant terms and then the quadratic and linear
terms. For \(a\ne0\), \(y\) has half the order of \(z\).
In every one of these cases the sum of the orders of
\(y,W\) is strictly below one. If the ramification
index of \(x-a\) is \(E\), the order of \(\omega\) is
\(E(1-\operatorname{ord}_{x-a}(yW))-1\), hence is
nonnegative.

For \(a=0\) and \(d=1\), \(y\) is a unit and the
differential is regular. For \(d=2\), \(xz\) is finite
and its value cannot be \(S\) for generic \(T\), as
the leading equation at \(z=S/x\) contains \(T\) with
nonzero coefficient. Thus \(y\) is a unit, \(W\) has
order one and there is a simple pole. For \(d=3,e=1\)
the same reasoning gives a simple pole on the branch
of order \(-1\) for \(z\) and regularity on the other
branch. It remains to treat \(d=3,e\ge2\) at zero.
Then \(r=\rho x^3\) and \(q=\gamma x^2\) with
\(\rho\ne0\). The curve of \(W/x\) satisfies
\[
                   (W/x)^2=x\bigl(\gamma^2S^2x
                                    -4\rho p S^2+4\rho T\bigr).
\]
This is a squarefree cubic, and \(x=0\) is a branch
point. On that elliptic curve \(y^2=S-xz\) has a pole
of order one there, so there is a single point over
zero on the final curve, with ramification index four
over \(x\). The orders of \(x,W,y\) at this point are
\(4,6,-1\), respectively. Thus \(\omega\) has a
double pole there. This final curve has positive genus
because it covers the elliptic curve.

### Infinity when \(r\ne0\)

Write \(d=\deg r\), \(e=\deg q\), taking \(e=-\infty\)
if \(q=0\), and \(D_\infty=\deg\Delta\).
Let \(v\) be the growth exponent of a root \(z\)
relative to \(x\) at infinity. If \(2e<d+2\), the
quadratic and constant terms give \(v=(2-d)/2\) and
\(D_\infty=d+2\). If \(2e>d+2\), the two possible
values are \(2-e,e-d\) and \(D_\infty=2e\).
If \(2e=d+2\), both roots have growth \((2-d)/2\)
and \(D_\infty\ge d\) by the independent coefficient
\(T\) of \(r\). These balances also follow directly
from the product of the roots and the expression
for \(W\). All values of \(v\) are strictly greater
than \(-1\) since \(d\le3,e\le2\). Hence \(y\) grows
as \(x^{(1+v)/2}\). The denominator \(yW\) grows
as \(x^{(1+v+D_\infty)/2}\), and the displayed bounds
give \(v+D_\infty\ge1\). Thus at infinity the form has
at most simple poles, which are themselves obstructions
if present. In the case \(d=3,e\ge2\) of the preceding
paragraph with a double pole at zero, \(D_\infty=5\)
and \(v=-1/2\), so there are no poles at infinity.

These calculations leave either a simple pole of nonzero
residue, a nonzero regular differential on a complete
curve, or just one double pole on a positive-genus
component. Each possibility rules out a rational
primitive. This proves the lemma.

## The highest tuple at transverse weights \((0,1,2)\)

**Theorem SexticKellerTangentConicSmoothSquareMemberTransverseReduction.**
Under the smooth-common-factor weight-three reduction with
\(i\ne0\), divide the third output by \(i\). After a
translation of \(x\), the highest tuple at weights
\((0,1,2)\) is one of
\[
       \bigl(U^3,(x^2+c)U^2,U^2\bigr),
       \qquad
       \bigl(U^3,m y z^2,U^2\bigr)
       \quad(m\ne0),\qquad U=y^2+xz.
\]
The second case requires a nonzero \(y\)-coefficient
in the coefficient of \(z^2\) of the second output.

We use the leading-pair slice also at nonnegative
weights. Its formal proof applies unchanged: the
rescaling by powers of the formal parameter is defined
for weight zero and each rescaled output still has
nonnegative parameter exponents. The same recursive
rectification and determinant identity therefore
produce the slice and its divisorial regularity.

### The coefficient of highest \(z\)-degree

By the complete axial degrees and the weight-\((1,2,3)\)
tuple, the highest outputs in weights \((0,0,1)\) are
\[
               (a z^3,b z^2,c z^2),\qquad
      a=a_0(x)+\alpha y,\quad b=b_0(x)+L(x)y,\quad c=c(x),
\]
where \(a_0,b_0,c\) are monic of degrees \(3,4,2\),
respectively, and \(\deg L\le1\). The leading determinant
vanishes, so \(a^2/c^3\) and \(b/c\) are algebraically
dependent. Indeed the two ratios are independent of
\(z\), and \(c z^2\) is transcendental over \(k(x,y)\).

If \(\alpha\ne0\), use \(x,a\) as plane coordinates
and write \(b=\lambda(x)a+D(x)\) with \(\lambda=L/\alpha\).
The vanishing bracket of the two ratios forces
\((D/c)'=0\) and \(2\lambda'c+\lambda c'=0\), by
the two coefficients in \(a\). The latter forces
\(\lambda=0\), since \(\lambda^2c\) is a polynomial
and \(c\) is quadratic. Then \(D=b_0\) has degree
four and \(D/c\) cannot be constant. Thus \(\alpha=0\).

If \(L\ne0\), the ratio \(b/c\) depends on \(y\), so
\(a_0^2/c^3\) must be constant. Monicity gives
\(a_0^2=c^3\), whence \(a_0=\tau^3,c=\tau^2\) for a
monic affine \(\tau\). Translate \(x\) so that \(\tau=x\).
On the generic level of \(xz\) and \(b z^2\) one solves
for \(z\) and then \(y\), and the relative form is
proportional to \(x\,dx/L(x)\). This uses the minor
of the two level equations in \(y,z\), which is
\(-x L(x)z^2\). The same form, up to a nonzero scalar,
comes from the leading pair \(((xz)^3,b z^2)\).
If \(L\) has a nonzero root, this differential has a
nonzero residue there. The leading slice excludes
that possibility. Therefore \(L=m\ne0\) or
\(L=\ell x\) with \(\ell\ne0\).

### Raising the weight of \(y\) when \(L=0\)

Give \(z\) weight two and increase the weight of \(y\)
from zero to one, keeping \(x\) of weight zero.
For \(L=0\), the only possible intermediate tie is
\(\delta y^3z\) in the second entry at \(2/3\).
This follows by testing \(ws+2l=6,4,4\) for the
three entries, using their strict lower-term
bounds in \((1,2,3)\). If \(a_0,c\) are independent
in the sense \(a_0^2/c^3\) is nonconstant, the
unchanged entries \(a_0 z^3,c z^2\) are independent
in \(x,z\); hence \(\delta=0\). Otherwise translate
\(x\) so that \(a_0=x^3,c=x^2\). On the generic level
\(xz=S,\ b_0(x)z^2+\delta y^3z=T\), if \(\delta\ne0\)
one obtains
\[
                    y^3=\frac{Tx^2-S^2b_0(x)}{\delta Sx}.
\]
The relative form is a nonzero scalar multiple of
\(dx/y^2\), from the minor in \(y,z\). The rational
function on the right has only simple finite zeros,
since away \(x=0\) they are a generic level of
\(b_0/x^2\). At zero it has a simple pole, is a unit,
or has a simple zero, according to the order of
\(b_0\). At infinity it has a pole of order three.
Consequently the form is regular at the finite
ramification points and at infinity, and is nonzero.
This excludes \(\delta\ne0\).

At weight one put \(t=y^2/z\). The three highest
outputs are \(z^3 A(x,t),z^2 B(x,t),z^2 C(x,t)\),
where \(A\) is monic cubic in \(t\), \(C\) is monic
quadratic in \(t\), and \(B\) is quadratic in \(t\)
with monic quadratic coefficient in \(x\) at \(t^2\).
The respective highest ordinary forms of these
polynomials in \(x,t\) are
\((x+t)^3,\ x^2(x+t)^2,\ (x+t)^2\). These bounds again
follow from the support at \((1,2,3)\).
The two ratios \(A^2/C^3,\ B/C\) are dependent.
At \(t=\infty\), the first tends to one and the second
has nonconstant residue, the coefficient of \(t^2\)
in \(B\). If \(A^2/C^3-1\) were nonzero, it would have
positive valuation. A valuation on a function field
of transcendence degree one that is nontrivial on a
nonconstant function has residue field algebraic over
the base field. Here this follows by restricting
first to the rational function field of that function
and using finiteness of the resulting extension.
The nonconstant residue of \(B/C\) contradicts it.
Thus \(A^2=C^3\). Monicity in \(t\) gives
\(A=(t+\tau)^3,C=(t+\tau)^2\) with \(\tau=x+\text{constant}\).
Translate \(x\) to absorb the constant.

The resulting highest tuple is \((U^3,B_*,U^2)\) with
\[
                       B_*=p(x)U^2+q(x)Uz+r(x)z^2,
\]
where \(p\) is monic quadratic, \(\deg q\le2\) and
\(\deg r\le3\). The two latter bounds follow by
cancellation of the degree-four leading polynomial
\(x^2(x+t)^2\). If \(q,r\) are not both zero, the
preceding relative differential lemma contradicts
the slice for \(U^3,B_*\).
With \(q=r=0\), if \(p'(0)\ne0\), the latter leading
pair has independent differentials at the generic point
of \(x=0\), so has a regular slice there. In
\(k(x,U)(y)\) its Jacobian derivation takes \(y\)
to \(3U^4p'(x)x\). A primitive for any nonzero scalar
therefore has the form of a nonzero scalar times
\(y/(3U^4p'(x)x)\) plus a function of \(x,U\).
Its simple pole at \(x=0\) cannot be cancelled by
that function, since its coefficient contains \(y\)
over the residue field \(k(U)=k(y^2)\).
Thus \(p'(0)=0\), giving the first tuple in the theorem.

### Raising the weight of \(y\) when \(L\ne0\)

Here the second highest weight as \(y\) increases
is \(4+w\), from \(L y z^2\). The same support test
shows there are no ties for \(0<w<1\). At weight one,
after the translation already made, put \(t=y^2/z\)
and write the highest outputs as
\[
             z^3 A,\qquad y z^2(\ell x+m+n t),
                        \qquad z^2 C,
\]
where
\[
 \begin{aligned}
 A&=(x+t)^3+a t^2+bxt+ct,\\
 C&=(x+t)^2+d t .
 \end{aligned}
\]
Exactly one of \(\ell,m\) is nonzero.
Put \(L_1=\ell x+m+n t\). The fractions
\(A^2/C^3\) and \(t^2 L_1^4/C^5\) are dependent by
the rank bound. After cancelling the common factors
in their bracket, it is equivalent to the equation
\[
\begin{split}
0={}&(2C A_x-3A C_x)
          \bigl(2L_1C+t(4C(L_1)_t-5 L_1 C_t)\bigr)\\
 &-(2C A_t-3A C_t)\,t(4C(L_1)_x-5 L_1 C_x).
\end{split}
\]
The six coefficients needed below are
\[
\begin{array}{c|l}
 x^5t &2\ell(2b-3d)\\
 x^4t &6m(2b-3d)\\
 x^4t^2 &2(6a\ell+6b\ell-2bn-18d\ell+3dn)\\
 x^3t &8cm\\
 x^3t^2 &2(14am+bd\ell+10bm+8c\ell-8cn-3d^2\ell-36dm)\\
 x^2t^2 &2(3bdm+2cd\ell+12cm-9d^2m).
\end{array}
\]
If \(m\ne0\), these give \(2b=3d,c=0,a=b,d=0\),
hence \(a=b=c=d=0\). If \(\ell\ne0\), they give
\(2b=3d,a=b,cd=0\).
When \(c=0\) the coefficient of \(x^3t^2\) then gives
\(d=0\). When \(c\ne0\), \(d=0\) and the same
coefficient gives \(n=\ell\). This last tuple
has the form \((U^3+cV^2,\ell UV,U^2)\) with \(V=yz\).
The pair \(U^2,\ell UV\) has independent differentials
at the generic point of \(z=0\). The source field is
\(k(U,V)(y)\), since \(z=V/y\), \(x=(U-y^2)/z\).
Its relative slice is a nonzero scalar times
\(y/(U^2 V)\) plus a rational function of \(U,V\),
since \(J(U,V,y)=-V\).
At \(z=0\) the residue of that function after multiplication
by \(V\) lies in \(k(U)\) whenever it has at most a simple
pole. The coefficient from the first term contains \(y\)
with \(y^2=U\), so this pole cannot be cancelled. The
regular slice excludes \(c\ne0\).

The remaining polynomials are \(A=U^3,C=U^2\) in source
coordinates, with second entry \(B=y z(L z+n y^2)\).
For \(L=\ell x\), on \(U=S\) this entry equals
\[
                         y(S-y^2)(\ell(S-y^2)+ny^2)/x.
\]
The source field is \(k(U,B)(y)\), and
\(J(U^3,B,y)=-3U^2B\). The primitive proportional to
\(y/(U^2B)\) has a simple pole at \(z=0\) if \(n\ne0\),
or at \(x=0\) if \(n=0\), on which \(B\) has order
one. In either case \(U^3,B\) have independent
differentials there, and \(U=y^2\) in the residue field.
The pole cannot be cancelled by an invariant, again
contradicting regularity.

Finally let \(L=m\ne0\). The generic curve on \(U=S,B=T\)
is \(y z(mz+n y^2)=T\), with \(x=(S-y^2)/z\).
If \(n\ne0\), put \(W=2mz+n y^2\) and \(V=yW\).
Apart from a nonzero scalar, the relative form is
\[
 \frac{dy}{y z W}
       =\frac{1}{2T}\left(dy+n\,\frac{y^3\,dy}{V}\right),
                    \qquad V^2=n^2 y^6+4mTy .
\]
The polynomial defining the double cover is squarefree.
The odd differential \(y^3\,dy/V\) is regular on the
affine curve. Any rational primitive is regular there,
and its odd part is \(Vh(y)\) with polynomial \(h\).
Its differential has numerator
\(f h'+f'h/2\) over \(V\), with \(f=n^2 y^6+4mTy\).
For \(h\ne0\) the numerator has degree at least five,
so it cannot equal \(y^3\). This excludes \(n\ne0\),
leaving \(B=m y z^2\), as asserted.

## Exclusion of the transverse cusp with a linear coefficient in \(y\)

**Corollary SexticKellerTangentConicSmoothSquareMemberTransverseCuspExclusion.**
A characteristic-zero Keller map has no highest tuple
\((U^3,m y z^2,U^2)\) of weights \((0,1,2)\) with
\(U=y^2+xz\) and \(m\ne0\). Thus only the first tuple
of the transverse reduction survives.

Apply the positive-weight correction part of the
[squarefree cusp argument](SEXTIC_KELLER_ALIGNED_CUSP_EXCLUSION.md#divisibility-of-each-square-correction)
with \(W=U\), \(H=m y z^2\), cusp outputs denoted
\(P,Q\) of weights four and six, and third output
denoted \(T\) of weight five. Here \(W\) is prime
and coprime to \(H\), the pair \((W,H)\) is onto
by taking \(z=1\), and \(k(W,H)\) is relatively
closed because \(k(x,y,z)=k(W,H)(z)\).
The rescaled determinant has order
\(4+6+5-(0+1+2)=12>9\).

For clarity, that part of the cusp argument only
uses the weight bounds through loss \(r\le3\) and
regularity of the chart along the generic point
of \(W=0\) at fixed \(H\). They hold at these
nonnegative weights. The formal rescalings still have
nonnegative exponents. The divisor \(W=0\) is reduced
and the restriction of \(H\) to it is onto by \(z=1\).
In particular the regular fixed-level chart exists
at its generic point. The rational and polynomial
kernels are \(k(W,H)\) and \(k[W,H]\), respectively,
by the field and surjectivity conditions. Although
the weight-zero source ring is \(k[x]\), the weight-zero
part of this polynomial kernel is \(k\). At any
loss \(r\le3\) the kernel comparison gives the same
leading tangent error as in that argument, and the
square and cube pole comparisons are at orders \(2r\)
and \(3r\), both strictly below twelve. Their scalar
conclusion on \(W=0\) makes the square correction
divisible by \(W\) because its weight is positive
and less than five. Thus the correction of the root
is polynomial, and iteration removes all
positive-weight square errors.

Consequently there is a polynomial \(V\) of weight
two with highest form \(U\) such that \(P=V^2+A(x)\).
The reference for \(Q\) at this point is a polynomial
in \(V\) of degree three plus a scalar multiple of
\(T\), with remaining weighted degree at most two.
At loss four the linear contraction of the
determinant gives
\[
                  2 E_{Q,2} - 3U A(x)\in k[U,H].
\]
Its weight is two, so it is a scalar multiple of \(U\).
Subtract a scalar multiple of \(P\) and of \(T\) from
\(Q\) to remove the quadratic reference term and
the reference term in \(T\). The subtraction of
the \(A(x)\) part of \(P\) is absorbed below weight
two. We therefore have
\[
                 Q=V^3+\tfrac32 V A(x)+cV+E,
                            \qquad \deg_{(0,1,2)} E\le1 .
\]
The exact differential wedge is
\[
 dP\wedge dQ=
    A'(x)\bigl(\tfrac32 A(x)+c\bigr)\,dx\wedge dV
      +2V\,dV\wedge dE+A'(x)\,dx\wedge dE .
\]
If the leading part of \(E\) has weight \(j\in\{0,1\}\),
the only contribution at weight \(4+j+5\) of the
three-form obtained by wedging with \(dT\) is
\(2U\,dU\wedge dE_j\wedge dH\). The volume form has
weight three, so this contribution is zero. Hence
\(E_j\in k[U,H]\). For weight one this forces zero,
and for weight zero it forces a constant. Removing
that constant makes \(E=0\).

The Keller determinant is now, up to sign,
\[
                 A'(x)\bigl(\tfrac32 A(x)+c\bigr)\,J(x,V,T).
\]
It is zero if \(A\) is constant. Otherwise the middle
factor is nonconstant, so the product cannot be a
nonzero scalar. Both cases are impossible.

## Exclusion of the nonzero weight-three coefficient

**Theorem SexticKellerTangentConicSmoothSquareMemberNonzeroCoefficientExclusion.**
In the binary case of the smooth-common-factor weight-two
theorem, the coefficient \(i\) in the weight-three reduction
is zero.

Suppose \(i\ne0\). Use the transverse reduction and its
cusp exclusion, divide the third output by \(i\), and
subtract the scalar multiple of the third output from
the second that removes \(c\) in the transverse tuple.
Denote the outputs in the order third, second, first by
\(P,B,Q\). With \(U=y^2+xz\), the highest tuples in the
two weight systems are then
\[
\begin{array}{c|ccc}
 &P&B&Q\\ \hline
 (0,1,2)&U^2&x^2U^2&U^3\\
 (1,2,3)&U^2+x^4U/i&x^2U^2&U^3 .
\end{array}
\]
The translation of \(x\) in the transverse reduction
preserves the second row. The subtraction from \(B\)
also preserves it, since the weight of \(P\) is lower
than that of \(B\) in that row.

### Completing the square in \(U\)

Over \(k(x)\), replace \(z\) by \((U-y^2)/x\) and
complete the square in \(U\). Put \(V=U+l(x)y+m(x)\)
so that
\[
   P=V^2+A,\qquad B=x^2P+bV+c,\qquad
   Q=V^3+dV^2+eV+f .
\]
Here all coefficients belong to \(k(x)[y]\). The
bounds in \((0,1,2)\) give
\(\deg_y A,\deg_y c\le3\), \(\deg_y b,\deg_y d\le1\).
In the variables \(x,y,V\), the determinant is
\(\kappa/x\) for a nonzero scalar \(\kappa\).

The weights of \(x,y,U\) induced by the second row
are \(1,2,4\) in the Laurent polynomial ring allowing
\(x^{-1}\). Put \(M=x^4/(2i)\). The highest part of
\(l y+m\) in these weights is \(M\) and that of \(A\)
is \(-M^2\). The second row also gives highest
coefficient \(-2x^2M\) for \(b\), by forming
\(B-x^2P\), and highest coefficient \(-3M\) for \(d\)
by expanding \((V-M)^3\). In particular at \(y=0\)
the asymptotic highest powers of \(A,b,d\) in \(x\)
are \(-M^2,-2x^2M,-3M\), respectively.

Compare the coefficients of \(V^5\) and \(V^4\) in
\(J_{x,y,V}(P,B,Q)\). The first is \(4x d_y\).
With \(d_y=0\), the second is
\(-2(3x A_y-2x e_y+b_y d')\), where a prime denotes
\(d/dx\) for coefficients in \(k(x)\).
It follows that \(d\in k(x)\) and
\(e=3A/2+h b+E\) with \(h=d'/(2x)\) and
\(E\in k(x)\). Thus
\[
         Q=dP+h(B-x^2P)+V(P+A/2+E)+G
\]
for some \(G\in k(x)[y]\). The highest power of \(h\)
at \(x=\infty\) is \(-6M/x^2\). In particular \(h'\ne0\).

In this form the coefficients of \(V^3,V^2\) in
the determinant give
\[
\begin{split}
  2xG_y&=(h'b+E')b_y,\\
  0&=-x(3A+2E)A_y-2\partial_y((h'b+E')c)
                            +2(b_xG_y-b_yG_x).
\end{split}
\tag{*}
\]
Both identities follow by differentiating the displayed
expressions in \(x,y,V\), with \(V\) held fixed in the
first two partial derivatives. In particular
\(G=(h'b^2+2E'b)/(4x)+g_0(x)\). In the second
identity, the terms other than \(-x(3A+2E)A_y\) have
degree at most three in \(y\). Hence \(\deg_y A\le2\).

### Excluding quadratic \(A\)

If \(\deg_y A=2\), its leading coefficient \(A_2\)
is nonzero. The second identity in \((*)\), at \(y^3\),
forces both \(b_y\ne0\) and a nonzero cubic coefficient
\(c_3\) of \(c\). Fix \(P=S\) at a generic scalar. On
\(V^2=S-A\) there are two expansions at \(y=\infty\)
with \(V\sim\pm\eta y\), where \(\eta^2=-A_2\).
Expand over \(k(x,S,\eta)\), and differentiate in \(x,y\)
at fixed \(S\), extending the derivations to \(V\).
The relative bracket of \(B,Q\) is, up to sign,
\(\kappa/(2xV)\), by the change from \(V\) to \(P\).
Its coefficient at \(y^5\) is zero on each branch.
The leading \(y^3\) coefficients of \(B,Q\) on the
two branches are \(c_3,hc_3\pm A_2\eta/2\).
The vanishing of the bracket coefficient says that
both ratios \(h\pm A_2\eta/(2c_3)\) have zero
\(x\)-derivative. Their sum contradicts \(h'\ne0\).

### Excluding affine nonconstant \(A\)

Assume \(A\) has degree one in \(y\). On the level
\(P=S\), set \(u=A=S-V^2\). The change first from
\(y\) to \(u\) and then from \(u\) to \(S\) gives a
relative bracket in \(x,V\) equal up to sign to
\(\kappa/(x A_y)\), independent of \(V,S\).

If \(b_y=0\), then \(G_y=0\). The second identity
in \((*)\) forces \(h'b+E'\ne0\) and makes \(c\)
quadratic in \(u\) with nonzero quadratic coefficient
\(-3x/(4(h'b+E'))\). On the level \(P=S\), \(B\)
has degree four in \(V\) and \(Q-hB\) has degree
three. In the relative bracket the term
\(-h' B B_V\) has degree seven, while the bracket of
\(B,Q-hB\) has degree at most six. This is impossible.

Thus write \(b=\beta(x)u+\gamma(x)\) with \(\beta\ne0\).
By integrating the second identity in \((*)\) in
\(y\), its terms other than \(\partial_y((h'b+E')c)\)
give a polynomial of degree at most three. Since
\(h'b+E'\) is affine and nonconstant in \(y\), \(c\)
has degree at most two. Use
\(c=c_2 u^2+c_1 u+c_0\) and \(G=g_2u^2+g_1u+g_0\)
with coefficients in \(k(x)\). The first identity
in \((*)\) gives \(g_2=h'\beta^2/(4x)\ne0\).
On \(P=S\) the two outputs are
\[
\begin{aligned}
 B&=x^2S+(\beta u+\gamma)V+c,\\
 Q&=dS+h(B-x^2S)+(S+u/2+E)V+G .
\end{aligned}
\]
For their bracket assign weights two and one to
\(S,V\), so \(u=S-V^2\) has weight two. The
\(V\)-derivatives here are at fixed \(S\).

If \(c_2\ne0\), the leading bracket of the weight-four
parts \(c_2u^2,(h c_2+g_2)u^2\) forces the ratio of
their coefficients to be constant, say \(\lambda\in k\).
Subtract \(\lambda B\) from \(Q\) and put \(H=h-\lambda\).
Then \(g_2=-H c_2\) and the weight-three part of the
sheared \(Q\) is \((S+K u)V\), where \(K=1/2+H\beta\).
The highest bracket with \(c_2u^2\), after dividing
by \(u^2\), is
\[
                  c_2'(1+K)S+(4c_2K'-3c_2'K)V^2.
\]
If \(c_2'\ne0\) the first coefficient gives \(K=-1\),
making the second nonzero. Thus \(c_2'=K'=0\).
As \(H\beta\) is a nonzero scalar and \(g_2=H'\beta^2/(4x)\),
the relation \(g_2=-H c_2\) gives \(H'/(xH^3)\)
a nonzero constant value. This is impossible at
infinity because \(H\) has growth \(x^2\).

If \(c_2=0\), the highest bracket is between
\(\beta uV\) and \(g_2u^2\). Dividing it by \(u^2\)
gives \(-4\beta' g_2 V^2-\beta g_2'(S-3V^2)\).
Thus \(\beta\) and \(g_2\) are nonzero constants and
\(h'=4xg_2/\beta^2\). Set \(K=1/2+h\beta\).
The next bracket, in weight five, between the parts
of weights \((3,3)\) and \((2,4)\), after division by
\(uV\), is
\[
                (2\beta K'-8xg_2)S
                          +(-3\beta K'-4g_2c_1')u .
\]
It gives \(c_1'=-3x\). For weight four, write
\(F=h c_1+g_1\). The relevant parts are
\[
\begin{array}{c|cc}
 \text{weight}& B&Q\\ \hline
 4&0&g_2u^2\\
 3&\beta uV&(S+Ku)V\\
 2&x^2S+c_1u&dS+F u\\
 1&\gamma V&(h\gamma+E)V.
\end{array}
\]
After taking the bracket and substituting \(S=u+V^2\),
the coefficient of \(u^2\) in weight four is
\(-3\beta xh-\beta F'-3x/2\).
The coefficient of \(uV^2\), using the preceding
one, is \(2c_1 K'-4g_2\gamma'\).
Thus \(\gamma'=2xc_1/\beta\); since \(c_1'=-3x\),
\(\gamma\in k[x]\) has degree at most four.
But at \(y=0\), \(\gamma=b-\beta A\) has highest
power \(\beta M^2\) in \(x\), of degree eight.
This contradiction finishes the affine case.

### Excluding \(A\) independent of \(y\)

It remains to consider \(A\in k(x)\). On the generic
level \(P=S\), one has \(V^2=S-A(x)\), with \(V\)
constant for the \(y\)-derivative. If \(b_y=0\),
then \(Q=hB+L\) with \(L\) independent of \(y\).
The relative bracket in \(x,y\) is
\(-B_y(h'B+L_x)\). If \(B_y=0\) this vanishes;
otherwise it has positive \(y\)-degree since \(h'\ne0\).
Both contradict the bracket \(\pm\kappa/(2xV)\).

If \(b_y\ne0\), use \(b\) itself as the affine coordinate
in place of \(y\) at fixed \(x\). The relative bracket
is then \(\pm\kappa/(2xV b_y)\), independent of \(b\).
As in the affine case, \((*)\) gives \(\deg_b c\le2\).
Writing \(c=c_2 b^2+c_1 b+c_0\), the quadratic
coefficient of \(G\) in \(b\) is \(g_2=h'/(4x)\ne0\).
If \(c_2\ne0\), the leading bracket gives
\((h c_2+g_2)/c_2=\lambda\in k\).
After subtracting \(\lambda B\), the linear coefficient
in \(b\) of \(Q\) is \(W=(h-\lambda)(V+c_1)+E'/(2x)\).
Since \(h-\lambda\ne0\), this is nonzero as a function
of the generic \(S\). Leading coefficients in the
bracket now give \(2c_2 W_x-c_2'W=0\).
Using \(V_x=-A'/(2V)\), comparison of \(V^{-1}\)
forces \(A'=0\). If \(c_2=0\), the leading
coefficients of the linear \(B\) and quadratic \(Q\)
are \(V+c_1,g_2\), and their relation
\(2g_2(V+c_1)_x-g_2'(V+c_1)=0\) gives the same
conclusion. These comparisons are valid since \(V\)
is transcendental over \(k(x)\) before specializing
\(S\). Finally, \(A'=0\) makes \(P=V^2+A\) critical
at \(V=0\) on the open set where the \(x\)-denominators
are invertible. Such points exist since the coefficient
of \(z\) in \(V\) is \(x\). This contradicts the
Keller condition and proves \(i=0\).

## The affine output in the binary smooth-common-factor frame

**Theorem SexticKellerTangentConicSmoothSquareMemberAffineOutput.**
In the binary case of the smooth-common-factor weight-two
theorem, the highest tuple in \((1,1,w)\) is
\((x^3z^3,x^4z^2,x^5z)\) for all \(w>3\).
The three \(z\)-degrees are \((3,2,1)\), and the third
output is
\[
                         S=v(x)z+q(x)y^2+u(x)y+w(x),
\]
where \(q,v\) are monic of degrees four and five.

By the nonzero-coefficient exclusion, the highest tuple
in \((1,2,3)\) is \((U^3,x^2U^2,x^4U)\), with
\(U=y^2+xz\). For increasing axial weight the monomial
test used above gives exactly the possible corrections
\[
              (0,\delta z^3,(a x+b y)z^2)
\]
at \(w=4\). The first two highest outputs are independent
in \(x,z\), as is seen from the nonconstant ratio
cancelling their \(x\)-powers in the parameters \(x,z/x^4\).
Their leading determinant with the third entry forces
\(b=0\).

Hold \(x,z\) at weights one and four and raise the weight
of \(y\) to \(5/2\). Testing \(r+ws+4l=15,12,9\) for
\(s>0\) with the strict lower-term bounds of \((1,2,3)\)
gives no earlier tie. At \(5/2\), the ordinary terms
complete \(U^3,x^2U^2,x^4U\), and the only new \(y\)-term
is \(j y^2z\) in the third output.
In the variables \(x,U,z\) the tuple is
\[
         (U^3,\ x^2U^2+\delta z^3,\
                         x^4 U+jUz+(a-j)xz^2).
\]
Its leading determinant is zero. After dividing the
determinant in \(x,U,z\) by \(-3U^2\), the coefficient
polynomial is
\[
              2j xU^3+4(a-j)x^2 U^2 z
                   -12\delta x^3 U z^2-3\delta(a-j)z^4.
\]
Hence \(j=a=\delta=0\). This excludes the axial tie at
four, and the ordinary tuple persists in \((2,5,8)\).

Above four the only possible axial tie is \(c z^2\) in
the third output at weight five. Hold \(x,z\) at weights
one and five and raise the weight of \(y\) to three.
The same test, now with \((2,5,8)\) and the bounds
\(18,14,10\), gives no intermediate tie. At three the
general tuple is
\[
 \begin{aligned}
 A&=U^3+a y z^3,\\
 B&=x^2U^2+b y^3z+d x y z^2,\\
 C&=x^4U+c z^2+e x y^3+f x^2 y z.
 \end{aligned}
\]
Its determinant is zero. The coefficients at
\(x^8z^5,\ x^7y^2z^4,\ x^6y^4z^3,\ x^5y^6z^2,\
x^4y^8z,\ x^3y^{10}\), up to nonzero rational factors,
are the linear expressions
\[
 \begin{gathered}
 a-2d+f,\quad 4a+6b-2d-3e-2f,\quad
 11a+14b-18d-12e+2f,\\
 3a+3b-11d-9e+4f,\quad
 6b+8d+12e-7f,\quad 4b+3e-2f.
 \end{gathered}
\]
Their vanishing gives \(a=e=0,\ b=d=f/2\).
Now the coefficients at \(x^5 y z^5\) and \(x^3z^6\)
are \(-6(4c+2d^2)\) and \(6cd\), respectively.
Thus \(c=d=0\), ruling out the last axial tie and
recovering the ordinary tuple in weights \((1,3,5)\).

The complete axial degrees and the bounds in
\((1,3,5)\) give the highest tuple in \((0,0,1)\) as
\[
                    (a(x)z^3,\ b(x,y)z^2,\ c(x,y)z),
\]
where \(a,b_0,c_0\) are monic of degrees \(3,4,5\),
and \(b=b_0(x)+\beta y,\ c=c_0(x)+L(x)y\) with
\(\deg L\le1\). The two ratios \(a/c^3,b/c^2\) are
algebraically dependent. If \(L\ne0\), the first
has order three at \(y=\infty\). If \(\beta\ne0\),
the second has order one, and the quotient of
the first by the cube of the second has residue
\(aL^3/\beta^3\). If \(\beta=0\), the second has
order two, and the quotient of the square of the
first by the cube of the second has residue
\(a^2/b_0^3\). These are nonconstant in \(x\), which
is impossible for a nontrivial valuation of a field
of transcendence degree one over \(k\). Thus \(L=0\).
Now \(a/c_0^3\) is a nonconstant function of \(x\),
so dependence of the ratios forces \(\beta=0\).

The third output is consequently \(v(x)z+q(x)y^2+
\alpha y^3+u(x)y+w(x)\), by the weight-\((1,3,5)\)
bounds. Keep \(x,z\) of weights zero and two and
increase the \(y\)-weight from zero to one. The
only possible tie before one is \(\alpha y^3\) in
the third output at \(2/3\). The first two highest
outputs \(a(x)z^3,b_0(x)z^2\) are independent since
\(a^2/b_0^3\) is nonconstant. Hence \(\alpha=0\).
The ordinary highest tuple gives the asserted degrees
and monicity of \(q,v\).

### The two highest quadratic parameters

**Theorem SexticKellerTangentConicSmoothSquareMemberQuadraticParameters.**
In the affine-output setting, \(q\mid v\).
Write \(v=qT\) with \(T\) monic affine and put
\(U_T=y^2+Tz\). In weights \((0,1,2)\) the
highest tuple is of the form
\[
       (U_T U_D^2,\ p(x)U_T U_D,\ q(x)U_T),
                      \qquad U_D=y^2+Dz ,
\]
where \(p\) is monic quadratic and \(D\) is monic
affine. If \(D\ne T\), then \(q=p^2\).
The map is invertible if \(\gcd(v,u)=1\). Any
noninvertible completion in this branch requires
a common root of \(v,u,q\), all vanishing to order
at least two there, with \(w'\ne0\) at the root.

By the preceding weight bounds, the highest outputs
\(A,B,C\) in \((0,1,2)\) are homogeneous polynomials
in \(Y=y^2,z\) over \(k[x]\) of degrees three, two,
and one, respectively. At \(z=0\) they are
\(Y^3,p(x)Y^2,q(x)Y\), with \(p\) monic quadratic by
the ordinary highest tuple. Also \(C=qY+vz\).
The highest determinant vanishes. On the symbolic
level \(C=s\), use \(H=Y/s\) and \(z=s(1-qH)/v\);
the polynomials \(A/s^3,B/s^2\) then have degrees at
most three and two in \(H\), over \(k(x)\).
Their image over \(k(s)\) has dimension at most one.

If both polynomials are independent of \(H\), evaluate
them at \(H=1/q\). This gives \(A=(Y+(v/q)z)^3\)
and \(B=p(Y+(v/q)z)^2\). Polynomiality of the
coefficient of \(z^3\) in \(A\) gives \(q\mid v\),
and this is the case \(D=T\).

Otherwise the image is a curve with an irreducible
equation of total degree at most three over \(k(s)\).
Indeed the equation stays irreducible over \(k(s,x)\)
and still defines the image, since at least one
polynomial depends on \(H\). Over this latter field,
the polynomial parametrization has degree at most
three in \(H\). A generic affine line therefore
meets the curve at most three times. This follows by
pulling it back to a polynomial of degree at most
three in \(H\): the parametrization is finite onto
the curve since \(H\) is integral over the algebra
of the output polynomials, and a generic line
has as many distinct intersections as the degree
of the curve in characteristic zero.

At \(H=1/q\) the coordinates of that image are
\(a=1/q^3,\ b=p/q^2\). Their orders at \(x=\infty\)
are twelve and six. If \(q\ne p^2\), then
\(a-b^2=(q-p^2)/q^4\), of order \(12+n\) with
\(1\le n\le4\). From the ten monomials in \(a,b\)
of total degree at most three, replace \(a,ab,a^2\)
by \(a-b^2,\ ab-b^3,\ a^2-ab^2\). The resulting
orders are all distinct: the three replacements
have orders \(12+n,18+n,24+n\), and the other seven
have orders \(0,6,12,18,24,30,36\). This contradicts
existence of the curve equation. Thus \(q=p^2\).
Now \(a=b^2\) on \(H=1/q\). As \(b\) is nonconstant
in \(x\), the image curve contains the parabola
and is therefore that parabola. In source variables
this says \(AC=B^2\).

Over \(k(x)[Y,z]\), the linear form \(C\) divides
\(B\). Writing \(B=CL\) gives \(A=CL^2\), where
\(L=Y/p+\lambda z\) by the restriction at \(z=0\).
Put \(T=v/q\) and \(D=p\lambda\). Then
\[
             C=q(Y+Tz),\quad B=p(Y+Tz)(Y+Dz),
                         \quad A=(Y+Tz)(Y+Dz)^2.
\]
Let \(\alpha,\gamma\in k[x]\) be the coefficients of
\(Y^2z,Yz^2\) in \(A\).
They are \(\alpha=T+2D\) and \(\gamma=D^2+2TD\).
Thus \(4\gamma=\alpha^2+2\alpha T-3T^2\).
A finite pole of \(T\) would give a pole of this
polynomial, since the highest pole order comes
from \(-3T^2\). Consequently \(q\mid v\) and \(D\)
is polynomial. Both are monic affine: \(T=v/q\) is
monic affine, and \(\alpha\) has degree at most one
with leading coefficient three by the ordinary
coefficient of \(y^4z\) in \(U^3\). This proves the
displayed form.

The quadratic affine-output root restriction applies
to \(S\). If \(\gcd(v,u)=1\), its coordinate completion
argument gives the same bounds nine and seventeen
for the inverse \(y,z\) coordinates. Indeed \(q,v\)
have the same degrees and root restrictions as in
that argument. For a monomial \(x^r y^s z^l\) of
either other output, \(s+2l\le6\) and
\(r+s+l\le6\). Its composed degree is at most
\(r+9s+17l\le54\). The retained planar degree-100
theorem and the same specialization and descent
then give inversion. The remaining root assertion
is the other conclusion of the quadratic affine-output
root restriction.

**Corollary SexticKellerTangentConicSmoothSquareMemberDistinctParameterInversion.**
In the quadratic-parameter reduction, if \(D\ne T\) then
the Keller map is invertible. Thus a noninvertible
completion has \(D=T\) and \(u\) shares a root with \(v\).

Suppose \(D\ne T\) and the map is noninvertible.
By the affine-output root restriction and the coprime
inversion argument, there is a common root \(r\) of
\(q,u,v\) where each vanishes to order at least two
and \(w'(r)\ne0\), over the algebraic closure.
On \(x=r\), the third output is constant and has
\(x\)-derivative \(w'(r)\), so the first two restricted
outputs form a planar Keller pair. Since \(q=p^2\)
in this case, \(p(r)=0\). In weights \((1,2)\) on
\(y,z\) the first restricted output has weight-six
highest form
\[
                     (y^2+T(r)z)(y^2+D(r)z)^2
\]
and the second restricted output has weight at most
three. The difference \(D-T\) is a nonzero scalar,
so \(D(r)\ne T(r)\).

The second output cannot be constant. If its weight
is one, it is affine in \(y\) and the planar Keller
determinant forces the derivative of the first
output in \(z\) to be constant. This is incompatible
with its weight-six highest form, which depends on
\(z\). For weight two or three, the highest bracket
must vanish. Two independent variables and the
weighted Euler relation then make the first highest
form a scalar times the cube or square, respectively,
of the second. In detail, contraction of the wedge
of the two logarithmic differentials with the Euler
derivation gives the zero differential of their ratio
with the two weights cancelled.

The displayed form cannot be a cube, since its two
distinct factors have multiplicities one and two
when neither parameter vanishes. If one parameter
vanishes, its factor is \(y^2\), and the multiplicity
of the other, irreducible factor is still one or two.
The displayed form cannot be a square unless
\(T(r)=0\). In that case \(D(r)\ne0\), and the
second highest form would be a nonzero scalar times
\(y(y^2+D(r)z)\). The second restricted output
therefore has \(z\)-degree one with nonconstant
affine coefficient in \(y\), since its weight is
at most three. Such a polynomial has a critical
point: take the root of its coefficient of \(z\)
and solve for \(z\) to annul its \(y\)-derivative.
This contradicts the planar Keller condition in
the only remaining case.

## Exclusion of the coincident quadratic parameters

**Theorem SexticKellerTangentConicSmoothSquareMemberCoincidentParameterExclusion.**
Let \(k\) have characteristic zero. Let \(T\in k[x]\)
be monic affine, \(p,q\in k[x]\) monic of degrees
two and four, and \(u,w\in k[x]\). Put
\(U=y^2+Tz\). There is no Keller map
\[
             (U^3+E_A,\ pU^2+E_B,\ S),\qquad
                  S=qU+uy+w,
\]
with \(E_A,E_B\) of weights at most five and three,
respectively, in \((0,1,2)\).

Work over the algebraic closure. Fix the third
output \(S=s\) and use \(x,y\) on the open set
\(v=qT\ne0\), where
\[
     U=(s-uy-w)/q,\qquad z=(s-qy^2-uy-w)/v .
\]
Denote the first and second restricted outputs by
\(A,B\), in \(k(x)[s,y]\). Their bracket at fixed
\(s\) is
\[
                       A_x B_y-A_y B_x=\kappa/v
\tag{1}
\]
for a nonzero scalar \(\kappa\), up to its sign.
The bounds on \(E_A,E_B\) give the following list
of possible powers, with coefficients in \(k(x)\):
\[
 \begin{array}{c|l}
 A& y^5,\ y^4,\ y^3,\ s y^3,\ y^2,\ s y^2,\
                    y,\ sy,\ s^2y,\ 1,\ s,\ s^2,\ s^3\\
 B& y^3,\ y^2,\ y,\ sy,\ 1,\ s,\ s^2 .
 \end{array}
\]
For example the terms of weight five in \(E_A\)
are multiples over \(k[x]\) of \(y^5,y^3z,yz^2\),
and the terms of weight three in \(E_B\) are
multiples of \(y^3,yz\). The coefficients of \(s^3\)
and \(s^2\) in the respective constant-in-\(y\)
terms are \(a=1/q^3,\ b=p/q^2\).

### Total degree in \(s,y\)

Give \(y\) weight one and decrease the weight of
\(s\) from above \(5/3\). At \(5/3\), only \(y^5\)
can tie \(a s^3\) in the first output, and the
second highest output is \(b s^2\). The expected
highest bracket is zero by (1). Since \(b\) is
nonconstant, this eliminates \(y^5\).
At \(3/2\), the only ties are \(s y^3\) in the first
output and \(y^3\) in the second. Put \(t=y^3/s^2\).
Dividing by \(s^3,s^2\) gives highest expressions
\(a+r(x)t,\ b+d(x)t\). Their bracket in \(x,t\)
must vanish. If both \(r,d\) are nonzero, the
coefficient of \(t\) forces \(r/d\) constant
and the constant coefficient would make a
nontrivial scalar combination of \(a',b'\) zero.
Those derivatives have distinct orders at infinity,
thirteen and seven. If only one coefficient is
nonzero the other derivative would vanish.
Hence \(r=d=0\). The isolated \(y^4\) tie in the
first output at \(4/3\) is likewise impossible
because \(b'\ne0\).

Consequently \(A,B\) have total \((s,y)\)-degree at
most three and two. Their homogeneous parts of these
degrees, after putting \(\theta=y/s\) and dividing
by \(s^3,s^2\), are polynomials \(f(x,\theta),g(x,\theta)\)
of degrees at most three and two in \(\theta\), with
\(f(x,0)=a,\ g(x,0)=b\). Their \(x,\theta\) Jacobian
vanishes by (1).

If one depends on \(\theta\), the image curve has
an irreducible equation of total degree at most
three, by the same polynomial parametrization
argument as in the quadratic-parameter reduction.
Its specialization at \(\theta=0\) contains
\((a(x),b(x))\). If \(q\ne p^2\), the ten monomials
of total degree at most three in \(a,b\) are
independent by the valuation comparison there:
\(a-b^2\) has order \(12+n\) with \(1\le n\le4\)
and replaces \(a\) in the order-twelve pair,
with the corresponding replacements at orders
eighteen and twenty-four. Thus \(q=p^2\), and
then \(a=b^2\) makes the image curve the parabola.
So \(f=g^2\). Since \(\deg_\theta f\le3\), \(g\)
has degree one, say \(g=b+n(x)\theta\) with \(n\ne0\).

### Excluding a nonconstant top on the symbolic level

In this case \(B\) is linear in \(y\), with
coefficient \(ns+n_0(x)\). Its square determines
the top of \(A\), which is quadratic in \(y\)
with \(y^2\)-coefficient \(n^2s+r_0(x)\).
Using \(x,B\) as coordinates over \(k(s)\)
makes \(A\) polynomial in \(B\). By (1) its
\(x\)-derivative at fixed \(B\) is
\(\kappa/(v(ns+n_0))\). In particular the constant
coefficient gives a rational primitive in \(x\)
for this rational function.

If \(n_0/n\) is nonconstant, its generic level
\(-s\) has simple finite points away from the
fixed zeros and poles of \(v,n\), by separability.
At each such point \(dx/(v(ns+n_0))\) has a
simple pole of nonzero residue. This is impossible
for the differential of a rational primitive.
Hence \(s_0=-n_0/n\in k\). On specializing to
\(s=s_0\), (1) remains valid since the original
denominators are polynomials in \(x\) alone.
Then \(B\) is independent of \(y\) and (1)
forces \(A\) to be linear in \(y\). Hence
\(r_0=-n^2s_0\).

Write \(t=s-s_0\). There are rational functions
\(d_j,m_j\in k(x)\) such that
\[
 \begin{aligned}
 B&=tny+bt^2+d_1t+d_0,\\
 A&=tn^2y^2+(2bn t^2+m_1t+m_0)y
                              +\text{a polynomial in }t .
 \end{aligned}
\]
In the second line the polynomial in \(t\) has
coefficients in \(k(x)\) and is independent of \(y\).
The coefficient of \(B^2\) in \(A\) written in
\(x,B\) is \(1/t\). Comparing the linear coefficient
gives \(A=B^2/t+(\alpha(x)+\beta(x)/t)B+L(x,t)\).
Its \(x\)-derivative at fixed \(B,t\) is
\(\kappa/(vnt)\), so \(\alpha,\beta\in k\).
The function \(L\) is a Laurent polynomial in \(t\)
with at most a simple pole at \(t=0\). Its only
coefficient possibly nonconstant in \(x\) is the
one at \(t^{-1}\). Computing it at \(y=0\) gives
\(-d_0^2-\beta d_0\). Thus
\[
         tA=B^2+(\alpha t+\beta)B+tK(t)-G(x),\qquad
        K\in k[t],\quad G=d_0^2+\beta d_0,
\tag{2}
\]
and \(G\) is nonconstant by (1).

Restore the polynomial map in (2), with \(t=S-s_0\).
The right-hand expression for \(G(x)\) is now
polynomial in the source, so \(G(x)\in k[x]\).
Also \(d_0\in k[x]\), since it is rational and
satisfies a monic quadratic over \(k[x]\).
Choose a root \(r\) of \(2d_0+\beta\); it exists
since \(G\), and hence \(d_0\), is nonconstant.
Then \(G'(r)=0\). At each point on the source
plane \(x=r\), differentiate (2) and use the
invertible differential of the Keller map.
The derivative of the target polynomial
\(B^2+(\alpha t+\beta)B+tK(t)-tA\) must vanish.
Its partials in \(A,B\) force \(t=0,B=-\beta/2\)
throughout the plane. The two outputs \(S,B\)
would be constant there, contradicting the
Keller condition.

### Excluding a top independent of \(y/s\)

If \(f,g\) are independent of \(\theta\), the actual
\(B\) has degree at most one in \(y\), with
coefficient \(N(x)\), and \(A\) has degree at
most two in \(y\). Its \(y^2\) coefficient is
independent of \(s\) and its \(y\) coefficient
is at most affine in \(s\). If \(N=0\), (1)
makes \(A\) linear in \(y\) and forces
\(-A_y B_x=\kappa/v\). But \(B_x\) has \(s\)-degree
two with coefficient \(b'\ne0\), which is impossible.

If \(N\ne0\), express \(A\) as a polynomial in
\(B\) at fixed \(x,s\). By (1) its positive-degree
coefficients are independent of \(x\), and the
coefficient \(\lambda\) of \(B^2\) is in \(k\)
since the \(y^2\) coefficient is independent of \(s\).
Its coefficient \(C_1(s)\) of \(B\) is polynomial in
\(s\), as seen by matching \(y\)-coefficients.
The \(s^2y\) coefficient in \(A\) vanishes, hence
\(2\lambda b+[s^2]C_1=0\) upon division by \(N\).
The function \(b\) is nonconstant, so \(\lambda=0\).
Thus \(A_y=a_0(x)+a_1(x)s\). The coefficient of
\(sy\) in (1) makes \(a_1/N\) constant in \(x\),
and the coefficient of \(s^3\) at \(y=0\) gives
\(Na'-a_1b'=0\). The distinct orders of \(a',b'\)
at infinity again give a contradiction.
This proves the exclusion.

**Corollary SexticKellerTangentConicSmoothSquareMemberBinaryInversion.**
Every characteristic-zero sextic Keller map in the
binary case of the smooth-common-factor weight-two
theorem is invertible.

Indeed the affine-output and quadratic-parameter
reductions apply. Coincident parameters are
excluded by the theorem just proved, while
distinct parameters give inversion by the
distinct-parameter corollary.

## The first corrections in the nonbinary axial tuple

**Theorem SexticKellerTangentConicNonbinaryAxialFirstCorrection.**
Let \(F\) be a characteristic-zero Keller map of
ordinary degree at most six with complete highest
\(z\)-terms \((y z^4,xy z^3,x^2y z^2)\) in an affine
source frame. Write \(a,b,c\in k[x,y]\) for the
coefficients of \(z^3,z^2,z\) in its three outputs.
Then
\[
                 x^2a-2xb+c=\lambda x^3y
                       \qquad(\lambda\in k)
\]
and \(b(x,0)=x a(x,0),\ c(x,0)=x^2 a(x,0)\).
In weights \((0,-1,1)\) its highest tuple is therefore
\[
                (Wz^3,xWz^2,x^2Wz),
                            \qquad W=yz+a(x,0).
\]

Put \(H=x^2a-2xb+c\) and
\(\mathcal D=x\partial_x-4y\partial_y\).
The coefficient of \(z^7\) in the full determinant
is \(-y(\mathcal D+1)H\). To obtain it, use the
three highest \(z\)-terms except for replacing one
by its term of degree one less in \(z\), and sum
the three determinants. Thus \((\mathcal D+1)H=0\).
By the ordinary degree bound, \(H\) has degree at
most five. A monomial \(x^i y^j\) in its kernel
requires \(i-4j+1=0\). Among nonnegative \(i,j\)
with \(i+j\le5\), only \((3,1)\) qualifies.
This proves the claimed formula for \(H\).

Write \(a_0,b_0,c_0\) for the restrictions at \(y=0\).
Then \(c_0=2xb_0-x^2a_0\). At \(y=0\) the
coefficient of \(z^6\) in the determinant is
\[
  3x^2 a_0 b_0'-2x^2 b_0 a_0'-3x a_0 c_0'
       +x c_0 a_0'+2 b_0 c_0'-c_0 b_0'.
\]
Indeed every term involving a replacement of degree
two less in \(z\) contains a remaining factor of
\(y\). The two-replacement terms surviving on the
plane give the displayed expression. Substitution
of \(c_0=2xb_0-x^2a_0\) factors it as
\(2\Delta(x\Delta'+2\Delta)\), where
\(\Delta=xa_0-b_0\). Since \(\Delta\) is polynomial
and the field has characteristic zero, its
vanishing forces \(\Delta=0\).
Thus \(b_0=xa_0,\ c_0=x^2a_0\). The claimed
highest tuple follows by comparing \(z\)-degree
minus \(y\)-degree in each entry.

## The doubled square-member frame through axial weight two

**Theorem SexticKellerTangentConicDoubledSquareMemberAxialThroughTwo.**
For a degree-at-most-six characteristic-zero Keller map
with ordinary highest tuple \((x^2P^2,x^4P,x^6)\),
where \(P=xz+y^2\), the highest tuple in weights
\((1,1,w)\) is \((x^4z^2,x^5z,x^6)\) for
\(1<w<2\). At \(w=2\) the highest tuple is binary
in \(x,z\).

### The remaining tie at three-halves

The nontrivial binary tie at \(3/2\) was already
excluded. In the nonbinary case, replace \(iy+ex\)
by \(y\). The leading tuple is then
\[
                     (Wz^2,Wxz,Wx^2),\qquad W=x^4+y z^2 .
\]
Take the leading pair \(B=Wxz,C=Wx^2\).
Its differentials are independent at the generic
point of \(z=0\): their transverse and tangential
terms are \(x^5\,dz, 6x^5\,dx\).
Thus there must be a rational slice regular along
that divisor. In the source field
\(k(B,C)(x)\), put \(t=B/C=z/x\). Directly
\[
                         J(B,C,x)=-C xz^2 .
\]
A primitive for unit value of this derivation is
\(1/(2Cz^2)\), and all such primitives differ
by a function of \(B,C\). Multiplying the displayed
primitive by \(B^2\), its residue at \(z=0\) is
\(x^4/2\). The restricted valuation on \(k(B,C)\)
is \(B\)-adic, with residue field \(k(C)=k(x^6)\),
so no invariant can cancel the pole. This
excludes the nonbinary tie.

### The binary tie at five-thirds

With the tuple unchanged at \(3/2\), the axial
monomial test for ordinary correction terms of
degree at most five gives only one possible tie
before \(w=2\), at \(5/3\). Its form is
\[
             (x^4z^2,\ x^5z+a z^4,\ x^6+(b x+c y)z^3).
\]
The first two entries are independent in \(x,z\),
so the leading determinant gives \(c=0\).

Suppose \((a,b)\ne(0,0)\). Put \(t=z^3/x^5\) and
\(\rho=x^2/z\). Then \(x=t\rho^3,z=t^2\rho^5\).
For the leading outputs \(A,B,C\), the ratios are
\[
        AC/B^2=(1+b t)/(1+a t)^2,\qquad
        B^{11}/A^{10}=(1+a t)^{11}/t^3 .
\]
When \(a\ne0\), their degrees as rational functions
are at most two and exactly eleven, and the first
is nonconstant. Thus they generate \(k(t)\).
When \(a=0,b\ne0\), the first alone is affine
nonconstant. The powers of \(\rho\) in the
outputs are \(22,20,18\), so the image field is
\(k(t,\rho^2)\), containing \(X=x^2,J=xz,Z=z^2\).

Raise the \(y\)-weight from one to the first tie,
at most \(4/3\) by the ordinary term \(x^2y^4\).
Specialization at \(y=0\) identifies the new leading
image ring with the binary one. Indeed the leading
determinant still vanishes, and the specialized ring
has the same dimension two. Integral elements of
the binary image field therefore lift to polynomials,
with their original weights, by the normalization
argument in the first axial tie analysis.

When \(a\ne0\), the binary tuple vanishes
simultaneously only at the origin, so \(X,J,Z\)
are integral. When \(a=0,b\ne0\), the identities
\[
          X^5(C-X^3)^2=b^2 A^3,\qquad
          J^3=X(C-X^3)/b
\]
first make \(X\) integral and then make the rational
lift of \(J\) polynomial. In either case the
weight-two \(X\) must lift as \(x^2\). The
weight-\(8/3\) element \(J\) can only acquire
a multiple of \(y^2\), at weight \(4/3\).
But \(J^2=XZ\) forces \(x\) to divide the
lift of \(J\). In the second case where \(Z\)
is initially rational, this divisibility follows
instead from \(J^3=X(C-X^3)/b\). Hence \(J=xz\)
and \(Z=z^2\) in the lift. As the three entries
are polynomial in \(X,J,Z\), this contradicts
the existence of the tie.

### Nonbinary tuples at weight two

At weight two the monomial test, after the
preceding exclusions, gives the highest tuple
on \(x\ne0\) as
\[
              (x^8t^2 F,\ x^7t G,\ x^6 H),
                    \qquad s=y/x,\quad t=z/x^2,
\]
where
\[
 \begin{aligned}
 F&=1+t a(s)+t^2 A_0,\\
 G&=1+t b(s)+t^2 B(s),\\
 H&=1+t c(s)+t^2 D(s)+t^3 C_0.
 \end{aligned}
\]
Here \(\deg a,\deg D\le2,\ \deg b\le3,\ 
\deg B\le1,\ \deg c\le4\), and \(A_0,C_0\)
are scalars. Put \(R=FH/G^2,\ V=G^8/(t^6 F^7)\).
The leading determinant implies their dependence.
In particular a rational function of \(R,V\)
of order zero at \(s=\infty\) has scalar residue
whenever the valuation of their field is nontrivial.

For the preliminary coefficient identities, the
matrix with rows
\[
       (8F,2F+tF_t,F_s),\quad
       (7G,G+tG_t,G_s),\quad
       (6H,tH_t,H_s)
\]
has determinant zero. Use subscripts for the
coefficients of powers of \(s\). Its coefficient
of \(t\) gives \(c_4=0,\ c_3=2b_3,\ c_2=2b_2-a_2\)
and \(c_1=2b_1-a_1\). At \(t^2s^5\) the remaining
term is \(36b_3^2\); at \(t^2s^3\), after \(b_3=0\),
it is \(24(a_2-b_2)^2\). Thus \(b_3=c_3=0\)
and \(a_2=b_2=c_2\). If \(a_1=a_2=0\),
the coefficient of \(t^2s\) also gives
\(D_2=b_1^2\).

Suppose first \(a_2\ne0\). At \(s=\infty\),
\(V\) has pole order two with leading coefficient
\(a_2/t^5\), and the residue of \(R\) is
\(1+tD_2/a_2\). Hence \(D_2=0\).
Unless \(D_1=2B_1\), the function \(R-1\) has order
one and \((R-1)^2V\) has residue proportional to
\(t^{-3}\), a contradiction. Thus write
\(G=F+L,\ H=F+2L+K\), where \(L\) is of degree
at most one in \(s\) and \(K\) is independent of
\(s\). Now \(FH-G^2=FK-L^2\) has degree at most
two in \(s\). If its leading coefficient in \(s\)
is \(N(t)\) in degrees two, one, or zero,
respectively, the order-zero test functions are
\((R-1)V,\ (R-1)^2 V^3,\ (R-1)V^2\).
Their residues, up to nonzero scalars, are
\(N/t^7,\ N^2/t^{19},\ N/t^{12}\).
The degrees of \(N\) are at most four, five, five,
respectively. None of the residues can be constant.
Thus \(F K=L^2\). If \(K\ne0\), the discriminant
of \(F\) in \(s\) must vanish. But its linear term
in \(t\) is \(-4a_2 t\). Therefore \(K=L=0\) and
\(F=G=H\).

Next let \(a_2=0,\ a_1\ne0\).
If \(G\) is independent of \(s\), then
\(c_1=-a_1\) and \(R\) has pole order \(m=2\) or
\(3\), according to whether \(D_2=0\).
The residue of \(R^7 V^m\) has \(t\)-order
\(-12\) or \(-18\) at zero. Indeed the leading
coefficient of \(V\) at \(s=\infty\) has
\(t\)-order \(-13\), and that of \(R\) has
\(t\)-order \(m\). This is impossible.
If \(G\) depends on \(s\), then \(V\) has pole
order one. If \(D_2\ne0\), the residue of \(R/V\)
is proportional to \(t^6/(b_1+t B_1)^{10}\),
nonconstant. Thus \(D_2=0\), and the residue
of \(R\) is
\[
                         R_0=a_1(c_1+t D_1)/(b_1+t B_1)^2.
\]
If it is zero, then \(c_1=D_1=0\) and \(b_1=a_1/2\).
The residue of \(RV\) has \(t\)-order \(-6\),
again impossible. A nonzero constant \(R_0\)
forces \(B_1=D_1=0\) by polynomial degrees in
\(t\). The numerator of \(R-R_0\) has \(s\)-degree
at most one. As above, a nonzero leading
coefficient \(N(t)\) of \(s\)-degree one or zero
gives residue proportional to \(N/t^7\) for
\((R-R_0)V\), or \(N/t^{12}\) for
\((R-R_0)V^2\), respectively. Its degree in
\(t\) is at most five. Thus \(R=R_0=1\),
where the last equality follows at \(t=0\).
Then \(a_1=c_1=b_1\), and \(F=G=H\) since
\(F(H-2G+F)=(G-F)^2\) and the latter two
differences are independent of \(s\).

In either of these cases with \(F=G=H\) the
leading tuple is \((Wz^2,Wxz,Wx^2)\) with
\[
                   W=x^4+(a_0x^2+a_1xy+a_2y^2)z+A_0z^2 .
\]
For \(a_2\ne0\), fix generic \(C=Wx^2\) and
\(r=z/x=B/C\). Complete the square in \(y\) with
\(\widetilde y=y+a_1 x/(2a_2)\) and put
\(M=x^2\widetilde y\). Then
\[
 M^2=\frac{x}{a_2 r}
         \left(C-x^6-(a_0-a_1^2/(4a_2))r x^5-A_0r^2 x^4\right).
\]
This polynomial is squarefree of degree seven
for generic values. The relative form of the
leading pair \(B,C\) is proportional to \(dx/M\),
a nonzero regular differential on the complete
curve. If \(a_2=0,a_1\ne0\), the source field
is \(k(B,C)(x)\). As in the weight-three-halves
obstruction, the pair \(B,C\) is regular at the
generic point of \(z=0\). Here
\(J(B,C,x)=-a_1 C x^2 z\), and the primitive
\(1/(2a_1 C xz)\) has a simple pole along \(z=0\).
Multiplying by \(B\) gives residue
\(1/(2a_1x^2)\) outside \(k(C)=k(x^6)\).
Both cases contradict the slice.

It remains to take \(a_1=a_2=0\). If \(G\) is
independent of \(s\), then \(V\) is a nonconstant
function of \(t\) alone and dependence forces \(H\)
to be independent of \(s\) as well. Otherwise
\(V\) has pole order eight at \(s=\infty\), and
\(D_2=b_1^2\). If \(b_1=0\), then \(G_1=B_1t^2\)
and \(H_1=D_1t^2\), with no quadratic term.
If \(H_1\ne0\), the residue of \(R^8V\) has
\(t\)-order \(-6\); if \(H_1=0\), that of \(R^4V\)
has the same order. This is impossible.
Thus \(b_1\ne0\), and constancy of the residue
of \(R\) gives \(F=(1+kt)^2\), \(k=B_1/b_1\),
by substituting \(t=0\) in that residue.
In particular the residue is one. If \(R\ne1\)
and the numerator of \(R-1\) has leading coefficient
\(N(t)\) in \(s\)-degree one or zero, respectively,
the residues of \((R-1)^8V\) or \((R-1)^4V\)
have respective \(t\)-orders
\(8\operatorname{ord}_0 N-14\) or
\(4\operatorname{ord}_0 N-6\), which cannot
be zero. Hence \(FH=G^2\), giving
\(G=(1+kt)(1+(j+b_1s)t)\) and
\(H=(1+(j+b_1s)t)^2\).
For \(k\ne0\) this uses divisibility of \(G\)
by \(1+kt\); for \(k=0\) its degree in \(t\)
is at most one because \(H\) has degree at most
three.

In source variables the tuple is
\(((zL)^2,zLM,M^2)\) with \(L=x^2+kz\),
\(M=x^3+(j x+b_1y)z\).
If \(k\ne0\), on generic \(zL=r,M=s\) set
\(N=xz\). Then \(N^2=rz-kz^3\) is elliptic,
and the relative form is proportional to
\(dz/(zN)\). This form has one double pole,
at \(z=0\), and no other poles. A primitive
would give a degree-one map from the elliptic
curve to the projective line. If \(k=0\), use
\(B=zLM,C=M^2\), which are independent at the
generic point of \(z=0\). With \(P_0=x^2z\),
\[
                         J(B,C,x)=-2b_1 M^2 P_0 .
\]
The source field is \(k(P_0,M)(x)\).
A primitive \(-x/(2b_1M^2P_0)\) has a pole
whose coefficient at \(P_0=0\) contains \(x\)
over \(k(M)=k(x^3)\); no invariant cancels it.
The contradiction finishes the weight-two
nonbinary exclusion.

### Returning to the conic in the binary doubled frame

**Theorem SexticKellerTangentConicDoubledSquareMemberWeight234Reduction.**
Under the doubled square-member hypotheses, the
highest tuple in weights \((2,3,4)\) is
\[
             (x^2P^2,\ x^4P+kxP^2,\ x^6+jx^3P+eP^2),
                                 \qquad ke=0 .
\]

At axial weight two the binary tuple has form
\[
       (x^4 z^2+a x^2 z^3+A z^4,\
        x^5z+b x^3z^2+B xz^3,\
        x^6+c x^4z+C x^2z^2+D z^3).
\]
Raise the weight of \(y\) from one to \(3/2\),
holding the weights of \(x,z\) at one and two.
The allowed monomials which could tie, using the
previous weight bound and ordinary support, are
\[
\begin{array}{c|lll}
 w_y &F_1&F_2&F_3\\ \hline
 6/5&&&y^5\\
 5/4&&y^4z&xy^4\\
 4/3&y^3z^2&xy^3z&x^2y^3,\ y^3z\\
 7/5&&y^5&\\
 3/2&y^4z,\ xy^2z^2,\ x^2y^4,\ x^3y^2z&
     y^2z^2,\ xy^4,\ x^2y^2z,\ x^4y^2&
     y^4,\ xy^2z,\ x^3y^2 .
\end{array}
\]
The isolated ties are impossible, since each pair
of unchanged entries is independent. This is seen
by comparing their weights \(8,7,6\) and orders
\(2,1,0\) in \(z/x^2\) at zero. At \(5/4\), on
the generic level of the first output, the
relative bracket of the two affine functions of
\(y^4\) would make the ratio of their \(y^4\)
coefficients constant along that binary level
if both were nonzero. Their ratio is proportional
to \(z/x\), which is independent of the first
output. If just one coefficient were nonzero,
the other binary output would be dependent on
the first one. This rules out the tie at \(5/4\).

At \(4/3\), write \(Y=y^3\). Denote the four
potential corrections by \(qYz^2,r xYz,u x^2Y+
vYz\). The determinant in \(x,Y,z\) vanishes.
Taking the coefficients at \(x^{10}z^2\) and
\(x^6Yz^2\) gives \(q-2r+u=0,\ qr-2qu+ru=0\),
hence \(q=r=u\). The coefficient at \(Y^2z^3\)
gives \(qrv=0\). If \(q=0\), the coefficient at
\(x^8z^3\) now gives \(v=0\). Otherwise \(v=0\),
and the coefficients at \(x^4Yz^3,x^2Yz^4,Yz^5\)
give \(c=2b-a,\ C=2B-A,\ D=0\).
The coefficients at \(x^6z^4,x^2z^6\) then reduce
to \(10q(a-b)^2,\ 8q(A-B)^2\). Thus the
candidate tuple is \((Wz^2,Wxz,Wx^2)\) with
\(W=x^4+a x^2 z+A z^2+q y^3\).
At generic levels \(z/x=t,\ Wx^2=s\),
the variable \(M=xy\) satisfies
\[
                          M^3=x(s-x^6-a t x^5-A t^2 x^4)/q .
\]
The relative form of the leading pair of last
outputs is proportional to \(x\,dx/M^2\).
The polynomial is squarefree for generic \(s,t\).
At the simple finite branch points this form
is regular; at infinity \(x,M\) have pole orders
three and seven. Hence the form is nonzero and
regular on the complete curve. This excludes the
tie at \(4/3\).

At \(3/2\) use coordinates \(x,P,z\) in the
subring \(k[x,y^2,z]\). The allowed tuple is
\[
 \begin{aligned}
 A_1&=x^2P^2+p_1 P^2z+p_2xPz^2+p_3x^2z^3+p_4z^4,\\
 B_1&=x^4P+p_5xP^2+p_6x^2Pz+p_7Pz^2+p_8x^3z^2+p_9xz^3,\\
 C_1&=x^6+p_{10}x^3P+p_{11}x^4z+p_{12}P^2+
                          p_{13}xPz+p_{14}x^2z^2+D z^3 .
 \end{aligned}
\]
Its determinant in \(x,P,z\) is zero. Write
\(a=p_1,d=p_5,b=p_6,h=p_7,g=p_8,j=p_{10},e=p_{12}\).
Taking the coefficients at \(x^{11}z^2,x^{10}Pz,
x^9z^3,x^9P^2,x^8Pz^2,x^6P^3,x^7P^2z,x^5P^2z^2\),
in that order, gives
\[
\begin{gathered}
 p_3=0,\quad p_2=2g,\quad p_4=g^2,\quad
 p_{11}=2b-a,\quad p_9=g(b-a),\\
 p_{13}=(a j-14ad+4jb+4bd)/6,\qquad
 p_{14}=(a-b)^2+(j-4d)g+2h,\\
 D=4adg-4dbg-a p_{14}+2bh+p_{13}g.
\end{gathered}
\]
Let \(I\) be the rational ideal of coefficients of the
determinant after these substitutions. The following
polynomials to the indicated powers belong to \(I\).
\[
\begin{array}{c|r}
 \text{polynomial}&\text{power}\\ \hline
 aj&6\\
 ag&2\\
 ad&4\\
 a^2-ab-dg+h&3\\
 e(ab+gj-2h)&3\\
 b(j-2d)&6\\
 bj^2+8ae-4be&6\\
 g(dj-2e)&3\\
 g(d^2-e)&3
\end{array}
\]
The [retained checker](verify_doubled_square_weight_234.py)
constructs the rational coefficient ideal and
verifies polynomial identities witnessing all nine
memberships. It also checks both tie lists and the
eight eliminations. Run it with Python 3 and Singular
4.3.2 from the repository root.

These equations give five families. For \(a\ne0\)
one has \(j=g=d=0\) and \(h=ab-a^2\).
If \(e\ne0\), then \(b=2a\). If \(a=0\),
then \(h=dg\); with \(g\ne0\) one gets \(e=d^2\)
and either \(j=2d\) or \(b=d=e=0\).
With \(a=g=0\), either \(b=0\), or \(j=2d\)
and \(e=d^2\). Substitution yields the following
list, allowing overlaps:
\[
\begin{array}{l|l}
 \text{highest tuple}&\text{parameters}\\ \hline
 (P^2L,PL^2,L^3+eP^2)&L=x^2+a z\\
 (P^2L,PLN,LN^2)&L=x^2+a z,\ N=x^2+(b-a)z\\
 (x^2P^2,x^4P+dxP^2,x^6+jx^3P+eP^2)&\\
 (W^2,x^3W,x^6+jx^2W)&W=xP+g z^2\\
 (W^2,WV,V^2)&W=xP+g z^2,\ V=x^3+dP+b xz .
\end{array}
\]

### Generic-fibre differentials of the five families

All generic levels in the following calculations may
be taken over the algebraic closure of their parameter
field. A nonzero scalar multiple on a fixed level
does not change exactness. We also use that for a
smooth double cover \(M^2=f(x)\) with \(f\) squarefree,
the cover branches at the simple roots and at
infinity when \(\deg f\) is odd. In the residue-zero
cases with a double pole at just one point of a
positive-genus curve, a primitive would have a
single simple pole and hence degree one.

**The first family.** If \(a\ne0\), the first two
outputs give the same relative differential up to
scalar as \(P,L\), since their bracket in \(P,L\)
is \(3P^2L^2\). On generic \(P=s,L=t\) one gets
\(y^2=s+(x^3-tx)/a\) and relative form \(dx/y\).
This is regular and nonzero on the elliptic model.
If \(a=0\) the family is already the third one.

**The second family.** Write
\(L=x^2+\alpha z,\ N=x^2+\beta z\).
When \(\alpha=\beta\), use the first family.
Otherwise fix \(P/N=s,\ LN^2=t\), furnished by
the leading pair of the two last outputs.
Use \(r=N\), \(M=xr\). Then
\[
 M^2=\frac{\alpha r^3-\beta t}{\alpha-\beta},
 \qquad
 y^2=sr+\frac{M}{\alpha-\beta}
                 -\frac{Mt}{(\alpha-\beta)r^3},\qquad
 \omega\ \propto\ \frac{dr}{M y}.
\]
The form follows from the Jacobian of
\((P/N,LN^2,N)\) in \((P,L,N)\), which is \(N\),
and \(J(P,L,N)\), which is proportional to \(xy\).
If \(\alpha\beta\ne0\), the first cover is elliptic.
At its finite ramification points \(dr/M\) is regular;
at \(r=0\) the pole of \(y^2\) only increases the
order of the differential. Other finite zeros of
\(y^2\) are simple for generic \(s\): they come
from a generic level of the rational function
obtained by solving the displayed equation for
\(s\). At infinity of the final curve, the pole
orders of \(r,M,y\) are \(4,6,3\), respectively.
Hence the form is everywhere regular and nonzero.

If \(\alpha=0,\ \beta\ne0\), \(M\) is a nonzero
scalar on each first component. The double cover
in \(y\) of the \(r\)-line has four simple finite
zeros in addition to its odd poles at zero and
infinity, so has genus two. The form has a double
pole at infinity and is regular elsewhere.
If \(\beta=0,\ \alpha\ne0\), normalize \(r=\ell^2\),
\(M=\ell^3\). The differential becomes a scalar
multiple of \(d\ell/(\ell^2 y)\). The poles of
\(y^2=s\ell^2+\ell^3/\alpha-t/(\alpha\ell^3)\)
are of order three at \(\ell=0,\infty\).
The form is regular there (order zero at zero),
and the generic finite zeros are simple. Thus
this case also gives a nonzero regular form.
All nontrivial cases contradict the slice.

**The third family.** If \(de\ne0\), the last
two outputs as functions of \(x,P\) have independent
differentials at \(x=0\), with Jacobian
\(2deP^3\) there. Their derivation annihilates
\(k(x,P)\), and this is the full invariant field
since \(k(x,y,z)=k(x,P)(y)\). The derivative of
\(y\) equals the Jacobian in \(x,P\) times
\(J(x,P,y)=-x\). A slice would therefore be a
nonzero scalar times \(y\) divided by that
derivative, plus an invariant. Its simple pole
at \(x=0\) cannot be removed since its coefficient
involves \(y\) outside \(k(P)=k(y^2)\).
Thus \(de=0\).

**The fourth family.** For \(g\ne0\), the first
two outputs fix \(W\) and \(x\) on each geometric
component. The conic
\(xy^2+x^2z+g z^2=W\) has relative form
proportional to \(dz/y\). At each point above
\(z=\infty\), \(y\) has pole order one and the
form has a nonzero residue. If \(g=0\), the
tuple is already in the third family.

**The fifth family.** Fix generic \(W=r,V=s\)
using the first two outputs, whose bracket in
\(W,V\) is \(2W^2\). Write
\(\alpha=d,\ \beta=b,\ \gamma=g\).
If \(\alpha\gamma\ne0\), elimination yields
\[
 \begin{aligned}
 M&=2\alpha\gamma z-\beta x^2,\\
 M^2&=(\beta^2+4\alpha\gamma)x^4-4\alpha\gamma s x+
                                                 4\alpha^2\gamma r,\\
 y^2&=(s-x^3-(\alpha+\beta)xz)/\alpha,\qquad
                      \omega\ \propto\ dx/(y M).
 \end{aligned}
\]
The first cover is a squarefree quartic or linear
cover. Its finite ramification cancels the factor
\(M\) in the differential. Finite zeros of \(y^2\)
are simple for generic \(r,s\). Indeed on \(y^2=0\)
the two fixed-level polynomials restrict to a map
of the \(x,z\) plane, and any rank-deficient locus
there has nondominant image in characteristic zero.
If that plane map is nondominant, there are no
zeros on the generic fibre at all.

For the quartic case, at infinity
\(z=\zeta x^2+O(x^{-1})\) on either first branch.
If the cubic term of \(y^2\) cancels, then
\(\gamma=\alpha+\beta\), \(\zeta=-1/(\alpha+\beta)\),
and the constant term of \(y^2\) is
\(s/(2\alpha+\beta)\ne0\).
Thus \(y\) is a unit or has a pole, and \(dx/M\)
is regular. For the linear case the first cover
has \(x,M\) of pole orders two and one. If the
cubic term in \(y^2\) cancels, then
\(\beta=-2\alpha,\ \gamma=-\alpha\); the next
term is a nonzero multiple of \(xM\), of pole
order three. Otherwise the pole order of \(y^2\)
is six. In both cases the pole of \(y\) cancels
the possible pole of \(dx/M\) after pullback.
Hence the form is nonzero and regular.

If \(\alpha\beta\ne0\) but \(\gamma=0\),
the equation is
\[
               y^2=(x^3-s)/\beta+
                        (\alpha+\beta)r/(\beta x),
               \qquad \omega\ \propto\ dx/(x^2 y).
\]
For \(\alpha+\beta\ne0\) it is a double cover
of genus two, with exactly one double pole of
\(\omega\) at \(x=0\). This can be read off
from \((xy)^2=x(x^4-sx+(\alpha+\beta)r)/\beta\).
For \(\alpha+\beta=0\), the identity
\[
                d(y/x)=\frac{x^3/2+s}{\beta x^2y}\,dx
\]
reduces \(\omega\) modulo exact differentials
to a nonzero scalar multiple of \(x\,dx/y\).
The latter is nonzero in the cyclic-cover
remainder on the squarefree cubic.

If \(\alpha=0,\ \beta\ne0\), the equation is
\[
 y^2=r/x-(s-x^3)/\beta-\gamma(s-x^3)^2/(\beta^2x^3),
                   \qquad \omega\ \propto\ dx/(x^2y).
\]
When \(\gamma=0\), the same genus-two analysis
gives a unique double pole. For \(\gamma\ne0\),
\(y^2\) has a pole of order three at zero.
At infinity it grows as \(x^3\) or is a unit,
with limiting value \(s/\beta\) in the latter
case \(\gamma=\beta\). The form is regular at
both points and at the generic simple zeros
of \(y^2\). If \(\alpha=\beta=0,\gamma\ne0\),
use the conic of the fourth family.
Finally if \(\alpha\ne0,\beta=\gamma=0\), the
tuple is in the third family with
\(d=\alpha,e=\alpha^2\), already excluded.
When all three parameters vanish the tuple is
the uncorrected one. This concludes the reduction.

## An affine output below two quartic transverse parts

**Lemma KellerTwoQuarticTransversePartsAffineOutputExclusion.**
In characteristic zero, give \(x,y,z\) weights
\((0,3,4)\). Suppose the first two outputs of a
polynomial map have weights at most twelve, with
\(y^4\)-coefficients \(p(x),r(x)\) of degrees two
and at most one, respectively, with \(r\ne0\).
Suppose the third output is
\(q(x)y^2+u(x)y+v(x)z+w(x)\), where \(v\ne0\)
and \(\deg q\ge3\). Then the determinant cannot
be a nonzero constant.

Increase the \(z\)-weight from four to six. At
four the only extra highest terms in the first
two entries are multiples of \(z^3\). The only
possible intervening tie is \(yz^2\) at \(9/2\).
At either weight fix \(q y^2=S\). The first two
highest entries are \(S^2p/q^2+a(x)\zeta\),
\(S^2r/q^2+b(x)\zeta\), where \(\zeta\) is the
tying monomial. The leading Jacobian in \(x,\zeta\)
vanishes because the expected determinant weight
is positive. If both coefficients are nonzero, its
\(\zeta\) coefficient makes their ratio constant,
and its remaining coefficient would make a
nonzero scalar combination of \(p/q^2,r/q^2\)
constant. If only one coefficient is nonzero,
the other rational function would be constant.
Neither is possible by their distinct positive
orders at \(x=\infty\). Both ties disappear.

At six normalize by fixing \(q y^2+v z=S\)
and setting \(H=y^2/S\). The first two highest
entries, divided by \(S^2\), are quadratic
polynomials in \(H\) over \(k(x)\). Their image
over \(k(S)\) has dimension at most one. If
one depends on \(H\), the image is a curve of
degree at most two by the polynomial-parametrization
argument for the [mixed doubled frame](#increasing-the-transverse-weight).
At \(H=1/q\) the two values are \(p/q^2,r/q^2\).
Writing their orders at infinity as \(m,n\),
one has \(0<m<n<2m\), since \(\deg q\ge3\).
Thus the orders \(0,m,n,2m,m+n,2n\) of the six
monomials of degree at most two are distinct.
No curve equation of degree at most two can hold.
Both polynomials are therefore independent of \(H\).
It follows that the weight-\((0,1,2)\) highest
tuple is
\[
           (pU^2,rU^2,qU),\qquad U=y^2+(v/q)z .
\]
From \(r(v/q)^2\in k[x]\) and \(\deg r\le1\),
a finite pole of \(v/q\) is impossible. Put
\(T=v/q\in k[x]\).

The bracket on the level of the actual third
output now has the form from the mixed
doubled-frame proof. In detail, after solving
for \(z\) and substituting \(U=(S-u y-w)/q\),
both outputs have form
\[
 a_3(x)y^3+a_2(x)y^2+(b_0(x)+b_1(x)S)y
                     +d_0(x)+d_1(x)S+d_2(x)S^2 .
\]
The respective \(d_2\)'s are \(p/q^2,r/q^2\),
and the bracket is a nonzero constant divided by
\(v(x)\), as a polynomial identity over \(k(x)[S]\).
The degree reductions and the arguments for
[a linear output](#one-degree-at-least-two-and-one-degree-one)
and [both degrees at most one](#both-degrees-at-most-one)
use only this form, the distinct orders of the
\(d_2\)'s, and the bracket. The remaining
cubic-quadratic case also applies without
requiring \(T\) affine or any bound on the
degree in \(x\) of \(b_1\).
Indeed use the rational scale \(h\) in the
[cubic-quadratic argument](#degrees-three-and-two),
write \(D=Y^2+d(x,S)\) with \(Y=hy+(b_0+b_1S)/(2h)\),
and let \(G=\eta Y^3+c(S)D+e(x)Y+f(S)\).
The comparison of the \(S y^2\) coefficient gives
\(b_1/(2h)\in k\). The \(S^2\) coefficient
of \(d\) is constant in \(x\) by the same bracket
calculation. Adding \((b_1/(2h))^2\) shows that
the original \(d_2\) of \(D\) is constant,
a contradiction. This proves the lemma.

## Removing the mixed lower corrections in small \(x\)-weights

**Theorem SexticKellerDoubledSquareNonzeroMixedCoefficientsExclusion.**
Under the doubled square-member reduction at weights
\((2,3,4)\), one cannot have \(k j\ne0\) when \(e=0\).

Divide the second and third entries by \(k,j\),
respectively, in the leading computations. Lower
the \(x\)-weight \(s\) from two while keeping the
\(y,z\)-weights three and four. The reference tuple
is \((x^2y^4,xy^4,x^3y^2)\). Testing
\(si+3h+4l=2s+12,s+12,3s+6\) using the support
at \((2,3,4)\) gives the following list. It is
tested from larger to smaller \(s\), removing
the terms ruled out at each step.
\[
\begin{array}{c|lll}
 s&A&B&C\\ \hline
 5/3&&&yz^2\\
 3/2&yz^3,\ y^5&&xy^3\\
 4/3&&&y^2z\\
 1&y^2z^2,\ xy^3z&y^3z&y^3,\ xz^2,\ x^2yz\\
 2/3&&&z^2\\
 1/2&y^3z&&xyz\\
 1/3&&&yz .
\end{array}
\]
The isolated transverse terms in the third entry
are excluded by fixing the independent first two
highest entries, which fix \(x,y\) algebraically.

### The tie at \(s=3/2\)

Fix the leading value \(S=xy^4\) of the middle
entry. The first and third on this level are
\[
            S^2/y^4+a y z^3+b y^5,\qquad
                            S^3/y^{10}+cS/y .
\]
The latter is nonconstant in \(y\), so vanishing
of the leading determinant forces \(a=0\).

If \(b\) or \(c\) is nonzero, the resulting binary
leading image field is \(k(x,y)\). For \(b\ne0\)
the two rational functions on the fixed middle level
have degrees nine and ten in \(y\). For \(b=0,c\ne0\)
the first determines \(y^4\) and the second has
degree ten. Any common proper field then has
index two and must be fixed by \(y\mapsto-y\),
which fails for the term \(cS/y\).

Raise the \(z\)-weight from eight with \(x,y\)
of weights three and six, a common scaling of
the current frame. The first new tie is at most
nine, since the ordinary \(x^4z^2\) in the first
entry ties then. Specialization at \(z=0\)
identifies the image rings, as in binary
normalization lifting, since both are domains
of transcendence degree two. If \(b\ne0\), the
identity \(b y^9-Ay^4+B^2=0\) makes \(y\)
integral, so it lifts unchanged (its weight is
below that of \(z\)). Any pole of the lifted
\(x=B/y^4\) must therefore be along \(y=0\).
A pole of order \(m\ge1\) there would give a
pole of order \(3m-2\) in \(x^3y^2+cxy^3\),
since its first term has strictly lower valuation
than the second. Hence \(x\) also lifts as a
polynomial and is unchanged.

If \(b=0,c\ne0\), use the rational lifts of
\(x,y\). At a prime divisor where \(x\) has pole
order \(m\ge1\), polynomiality of \(x^2y^4\)
requires that \(y\) have positive order.
Polynomiality of \(x^3y^2+cxy^3\) then forces
\(\operatorname{ord} y\ge3m/2\), since its first
term has strictly lower valuation than the
second. It follows that \(B=xy^4\) vanishes to
order at least \(5m\) there. The divisor can be
taken weighted homogeneous, as the rational
lifts are homogeneous. Its defining prime
specializes to a nonzero nonconstant polynomial
at \(z=0\), since \(B|_{z=0}=xy^4\ne0\).
This would give a factor of multiplicity at
least five in \(xy^4\), impossible.
Thus \(x\) is polynomial and unchanged by its
weight. A pole of \(y\) is now incompatible
with \(B=xy^4\), since \(x\) is prime of
multiplicity one. Hence \(y\) is also unchanged.
In either case the unchanged rational parameters
preclude the first \(z\)-tie. Therefore \(b=c=0\).

### The tie at \(s=1\)

The potential highest tuple is
\[
 \left(
 x^2y^4+a y^2z^2+bxy^3z,\
 xy^4+c y^3z,\
 x^3y^2+d y^3+f xz^2+g x^2yz
 \right).
\]
The coefficients of its determinant at
\(x^4y^8,\ x^3y^7z,\ x^2y^6z^2,\ xy^5z^3,\ y^4z^4\)
give, after division by nonzero scalars,
\[
 \begin{gathered}
 5b-4c-2g,\quad 10a+bc+3bg-6cg-4f,\quad
 6ac+8ag-bcg+bf-8cf,\\
 4acg+6af-3bcf,\qquad acf .
 \end{gathered}
\]
If \(c=0\), the first two identities give
\(g=5b/2,\ f=5a/2+15b^2/8\).
The others give \(b(20a+f)=0\) and \(af=0\),
so \(a=b=f=g=0\).
If \(c\ne0,f=0\), the fourth identity gives
\(ag=0\). If \(g=0\), the first and third
give \(b=4c/5,a=0\), contradicting the second.
If \(a=0\) and \(g\ne0\), the third gives \(b=0\)
and the first \(g=-2c\), again contradicting
the second. Finally if \(acf=0\) with \(cf\ne0\),
then \(a=0\), and the fourth gives \(b=0\);
the third now contradicts \(cf\ne0\).
Thus \(a,b,c,f,g=0\).

If \(d\ne0\), the binary tuple gives \(x\) and
\(y^4\) by the first two entries. The last entry
\(x^3y^2+d y^3\) has degree three in \(y\) over
\(k(x)\), so the degrees four and three give
\(y\) in the image field. It is integral by the
last entry. At the first new \(z\)-tie, of weight
at most five by the ordinary first output,
\(y\) therefore lifts unchanged. The rational
lift \(x=B/y^4\) cannot have a pole along \(y=0\),
since \(x^3y^2+d y^3\) is polynomial. Its
weight then forces it unchanged, a contradiction.
Thus \(d=0\).

### The last ties and the affine obstruction

At \(s=1/2\), on the level \(S=xy^4\) the two
other entries have the form
\(S^2/y^4+a y^3z,\ S^3/y^{10}+b S z/y^3\).
If \(a\ne0\), eliminating \(z\) with the
first entry gives the second as a combination
of \(A y^{-6}\) and \(y^{-10}\), with coefficients
\(bS/a\) and \((1-b/a)S^3\). They cannot both
vanish, contrary to the leading rank condition.
If \(a=0\) the unchanged first entry forces
\(b=0\). This eliminates all positive \(s\)
ties. At \(s=0\) the first two outputs have
weight at most twelve, with \(y^4\)-coefficients
of degrees two and one. The third output
has weight at most six and therefore has form
\(q(x)y^2+u(x)y+v(x)z+w(x)\).
Here \(q,v\) are nonzero of degrees three and
four from the \(j x^3P\) term in the initial
\((2,3,4)\) highest tuple. These are precisely
the hypotheses of the two-quartic transverse
parts exclusion, giving the contradiction.

**Theorem SexticKellerDoubledSquareQuarticThirdCorrectionExclusion.**
Under the doubled square-member reduction at weights
\((2,3,4)\), one cannot have \(e\ne0\).

Here \(k=0\) by the preceding reduction. Lower the
\(x\)-weight \(s\) from two with the \(y,z\)-weights
three and four. The reference tuple is
\((x^2 y^4,x^4y^2,e y^4)\).
The possible ties, tested in order using the
strict bounds at \((2,3,4)\), are as follows.
\[
\begin{array}{c|ll}
 s&A&B\\ \hline
 7/4&&y^3z\\
 5/3&&xyz^2\\
 3/2&yz^3,\ y^5&z^3,\ y^4,\ x^2y^3\\
 4/3&&xy^2z\\
 5/4&&yz^2\\
 1&y^2z^2,\ xy^3z&y^2z,\ xy^3,\ x^2z^2,\ x^3yz\\
 3/4&&y^3\\
 2/3&&xz^2\\
 1/2&y^3z&z^2,\ x^2yz\\
 1/3&&xyz\\
 1/4&&yz .
\end{array}
\]
The third entry stays \(e y^4\) throughout the
descent, as every other term has positive
\(x\)-weight at two or strictly smaller
\((0,3,4)\)-weight. The isolated transverse
ties are impossible because the unchanged
entries fix \(x,y\) algebraically.

### The tie at \(s=3/2\)

Fix \(y\) generically using the third entry.
The first two are
\[
        x^2 y^4+a y z^3+b y^5,\qquad
             x^4 y^2+c z^3+d y^4+g x^2 y^3.
\]
If \(a\ne0\), eliminating \(z^3\) using the
first leaves an uncancelled \(x^4y^2\) in the
second, contradicting the leading rank.
Then \(a=0\) forces \(c=0\).
Subtract \(d/e\) times the third entry from
the second by a constant target shear.
The only introduced degree-six term is an
\(x^6\) term, below the reference at all
\(s<2\).

Put \(L=x^2+b y\). On the binary level,
\(L=A/(C/e)\) and \(y^4=C/e\). The second
entry is
\[
           y^2 L^2+(g-2b)L y^3+(b^2-bg)y^4.
\]
If \(b\ne0,\ g\ne2b\), then \(y\) belongs
to the image field: in the extension of
\(k(L,y^4)\) obtained by adjoining \(y\),
no nontrivial fourth-root scaling preserves
the cubic term. If \(b=0,g\ne0\), the same
conclusion follows from degrees four and
three in \(y\) over \(k(L)\).

In either case \(y\) is integral, as is
\(K=x^2y\) from \(B=K^2+g y^2 K\).
Raise the \(z\)-weight from eight with
\(x,y\) of weights three and six. The first
tie is at most nine by the ordinary first
entry. The integral lifts are \(y\) unchanged
and possibly \(K=x^2y+\delta xz\), the only
new term allowed by its weight twelve.
But the identity \(A=y^3K+b y^5\) cannot
supply the ordinary \(x^4z^2\) required
at nine. If the tie were earlier, neither
generator could lift nontrivially. Both
possibilities are impossible.

It remains at this tie to check \(b\ne0\),
\(g=2b\). Here \(Y=y^2\) belongs to the image
field since \(B=Y L^2-b^2Y^2\), with \(L\)
and \(Y^2\) already in the field.
The integral \(Y\) lifts either unchanged
or as \(y^2+\delta xz\) at nine. If it stays
unchanged, \(yL\) in the source rational
field is polynomial by \(B=Y L^2-b^2Y^2\).
The same weight-twelve bound on \(yL\) as
on \(K\) prevents the first entry \(A=y^3(yL)\)
from acquiring the required ordinary term.
If the lift of \(Y\) has \(\delta\ne0\), it
is prime, and \(Y L^2\) polynomial implies
\(L\) is polynomial. Its weight makes it
unchanged. Then \(A=Y^2 L\) and the ordinary
terms force \(\delta=1\). But the leading
pair \(C,A\) on the generic level of \(Y,L\)
gives a differential proportional to \(dx/x\):
here \(Y=xz+y^2,\ L=x^2+b y\), and one
solves first for \(y\) and then for \(z\).
The nonzero residue contradicts the slice.
Thus both \(b,g\) vanish.

### The ties at \(s=1\) and \(s=3/4\)

On the generic \(y\)-level at \(s=1\), write
\(t=z/(xy)\). The first two entries become
\[
                 x^2y^4 L(t),\qquad
                 x^4y^2 M(t)+xy^3 N(t),
\]
with \(L(0)=M(0)=1\) and degrees at most two,
two, and one, respectively. The coefficients
of the relative bracket in \(x,t\) give
\(M=L^2\) and \(N^2/L\) constant when \(N\)
is nonzero. In particular \(L\) has degree
at most one. If \(L=1+b t\) with \(b\ne0\),
then \(N=0\) and the leading tuple is
\((y^3 R,R^2,e y^4)\) with \(R=x^2 y+b xz\).
At fixed \(y,R\) the relative form for the
pair \(C,A\) is proportional to \(dx/x\),
by solving for \(z\). This excludes \(b\ne0\).
If \(L=1\), then \(M=1\), \(N\) is a scalar
\(d\), and the tuple is
\((x^2y^4,x^4y^2+dxy^3,e y^4)\).

Suppose \(d\ne0\) and put \(H=xy^2\).
The field \(k(H,y)\) over \(k(H^2,y^4)\)
has the sign action on \(H\) and the fourth-root
action on \(y\). The entry \(B=H^4/y^6+d H y\)
is invariant exactly under the simultaneous sign.
Thus the image field contains \(Y=y^2\) and
\(J=yH\). They are integral, using the third
entry for \(Y\) and
\(J^4+dY^5 J-Y^5B=0\) for \(J\).
Raise the \(z\)-weight from four to the first
tie, at most five. The lift of \(Y\) is \(y^2\)
or \(y^2+\delta xz\) at five. In the latter
case with \(\delta\ne0\), \(B\) is no longer
needed to see divisibility: \(A=J^2/Y\)
polynomial gives \(Y\mid J\).
Then \(J/Y\) has weight four and is unchanged,
equal to \(xy\), giving \(A=Yx^2y^2\).
This cannot supply the ordinary \(x^4z^2\).
If \(Y=y^2\), then \(J/y\) is polynomial of
weight seven. Its only possible \(z\)-correction
is a multiple of \(x^2z\) at five. The
ordinary first entry forces \(J/y=xP\) there,
but \(B=J^4/Y^5+dJ\) would then have a pole
at \(y=0\). This excludes \(d\ne0\).

At \(s=3/4\) the only new tie is \(d y^3\)
in the second entry. If \(d\ne0\), the binary
tuple gives \(x^2=A/(C/e)\); over this parameter
\(C\) gives \(y^4\) and the second entry has
degree three in \(y\), so \(y\) is in the
image field. With \(K=x^2 y\) the two integral
parameters are \(y,K\), since
\[
                  (A,B,C)=(y^3K,K^2+d y^3,e y^4).
\]
Use weights three and twelve for \(x,y\),
and raise the \(z\)-weight from sixteen to
the first tie, at most twenty-one. Then \(y\)
is unchanged, and the only possible new
term in \(K\) is \(\delta z\) at weight
eighteen. Without it the tie is impossible.
With \(\delta\ne0\), take \(x,y,K\) as source
coordinates, with weights \(3,12,18\).
Raise the now unused \(x\)-weight. At every
first tie the image ring again lifts by its
integral parameters \(y,K\). They can acquire
only terms in \(x\) alone for \(y\), and terms
in \(x\) and linear in \(y\) for \(K\).
Take the lifted pair as new coordinates,
preserving the binary tuple on \(x=0\).
The possible weights for these changes are
finite, and the \(x\)-weight strictly increases
at each step. A tie must recur because the
Keller map depends on \(x\); its first tie
has \(x\)-weight at most the largest component
weight, so the leading determinant still
has positive expected weight at every step.
This gives a contradiction and excludes the tie.

### The tie at \(s=1/2\)

Fix \(y\) generically and write \(t=z/x^2\).
The first two highest entries are
\(x^2y^4+a y^3z,\ x^4y^2+b z^2+c x^2 y z\).
The relative bracket forces the second to be
the square of \(R=x^2y+a z\): the ratio of
its \(t\)-polynomial to the square of the
\(t\)-polynomial of \(R/x^2\) is constant
and equals one at \(t=0\). If \(a=0\) there
is no tie. Otherwise take \(x,y,R\) as source
coordinates of weights \(1,6,8\), scaling the
current weights by two. The highest tuple is
\((y^3R,R^2,e y^4)\).
In its image field the elements \(Y=y^2,J=yR\)
are integral, since \(Y^2=C/e,\ J^2=YB\);
also \(A=YJ\). They are rational there by
these formulas, for instance \(Y=A^2/(BC/e)\)
and \(J=A/Y\).

Raise the unused \(x\)-weight to the first
tie. The lifted \(Y\) is monic quadratic in
\(y\), as its weight is twelve and the other
weights are positive. If it is squarefree
then \(Y\mid J\) by \(J^2=YB\). Specializing
at \(x=0\) would make \(y^2\mid yR\), impossible.
Thus \(Y\) is the square of \(y+\ell\) for a
polynomial \(\ell\) of weight six in \(x\)
alone, since \(R\) has weight eight. Use
\(y+\ell\) as coordinate. The quotient
\(J/(y+\ell)\) is polynomial of weight
eight and specializes to \(R\). It differs
from \(R\) only by terms in \(x\) and at most
linear in the new \(y\). Taking it as new
coordinate preserves the same binary tuple.
Again only finitely many increasing \(x\)
weights permit these changes, whereas a
first tie must recur, with positive expected
determinant weight by the same bound.
This excludes \(a\ne0\).

All positive \(s\) ties have been eliminated.
At \(s=0\), the first and third entries
have weight at most twelve, with \(y^4\)
coefficients of degrees two and zero.
The second has weight at most six, so is
\(q(x)y^2+u(x)y+v(x)z+w(x)\).
The ordinary terms give \(\deg q=4,\deg v=5\),
as the target shear introduced only the
ordinary pure power \(x^6\).
The two-quartic transverse parts exclusion
now applies after putting this affine output
last. This proves \(e=0\).


## The third-output mixed correction

**Theorem SexticKellerDoubledSquareThirdMixedCorrectionExclusion.**
In the doubled square-member weight-\((2,3,4)\) reduction, the
coefficient \(j\) is zero.

Suppose \(j\ne0\). The preceding exclusions give \(k=e=0\).
Work over the algebraic closure and scale the third output by \(1/j\).
Its weight-\((2,3,4)\) highest part is then \(x^6/j+x^3P\), with
\(P=xz+y^2\). Denote the highest parts at subsequent weights by
\(A,B,C\). At each weight used below the expected determinant
weight, the sum of the three output weights minus the sum of
the source weights, is positive. Thus the highest triple has
rank at most two. For a first tie obtained by increasing the
weight of a source variable absent from a binary triple, the
image field identifies with that of the binary specialization.
Indeed, both images have dimension two and the latter is a
closed subset of the former after taking closures. Rational
parameters of the binary image therefore lift as homogeneous
rational functions of the same weights.

### Axial weights up to three

At weights \((1,1,w)\) with \(w>2\), the reference triple is
\((x^4z^2,x^5z,x^4z)\). Comparison with the weight-\((2,3,4)\)
support gives the following possible ties through \(w=4\).
The ordinary degree bound is six.
\[
\begin{array}{c|lll}
w&A&B&C\\ \hline
5/2&&z^3&\\
3&xz^3,\ yz^3&x^2z^2,\ xyz^2&xz^2,\ yz^2\\
4&z^3&xz^2,\ yz^2&z^2
\end{array}
\]
The support comparisons and determinant identities in this proof
are reproduced by the [exact coefficient checker](check_sextic_doubled_square_third_mixed.py).
Run it from the repository root with
`.venv/bin/python research/jacobian/monodromy-theory/check_sextic_doubled_square_third_mixed.py`.
The internal replay used SymPy 1.14.0 and Singular 4.3.2 over the rationals.

At \(w=5/2\), suppose \(B=x^5z+\lambda z^3\) with
\(\lambda\ne0\). The image field is \(k(x,z)\): use \(z=A/C\),
\(x^4=C/z\), and \(x^5=(B-\lambda z^3)/z\).
Increase the weight of \(y\) from one to its first tie, which
occurs at or before \(7/4\) because \(x^2y^4\) occurs in the
first output. Use \(x,z\) also for their rational lifts.
At a pole of \(z\), the identity \(C=x^4z\) makes \(x\) regular,
and \(\lambda z^3\) then dominates \(B\), impossible.
At a pole of \(x\) of order \(m\ge1\), regularity of \(B\)
forces \(\operatorname{ord}z\ge5m\) and thus
\(\operatorname{ord}A\ge6m\). A homogeneous prime dividing \(A\)
specializes at \(y=0\) to a nonzero polynomial of positive
weight, and its sixth power cannot divide \(x^4z^2\).
The lifts are polynomial. The lift of \(x\), of weight one,
is unchanged. The only possible changes to the lift of \(z\)
are adding a multiple of \(y^2\) at \(w_y=5/4\) or of \(xy\)
at \(w_y=3/2\). Either change produces a term of ordinary degree
eight in \(x^4z^2\) after substitution. Hence there is no first
tie, contradicting the presence of \(x^2y^4\).
This excludes \(\lambda\ne0\).

At \(w=3\) write the possible highest triple as
\[
 \bigl(x^4z^2+(a x+\beta y)z^3,\quad
       x^5z+(c x^2+\delta xy)z^2,\quad
       x^4z+(e x+\varphi y)z^2\bigr).
\]
Three determinant coefficients give
\[
 -\beta-4\delta+6\varphi=0,\qquad
 -5\beta\delta+5\beta\varphi+2\delta\varphi=0,\qquad
 \beta\delta\varphi=0 .
\]
The last equation sets one variable to zero, the second sets
another to zero, and the first sets the third to zero.

The resulting binary image field is \(k(x,z)\) for all
\(a,c,e\). In fact, put \(\theta=z/x^3\). The ratios
\[
 R=B^7/C^8=\frac{(1+c\theta)^7}{\theta(1+e\theta)^8},
 \qquad S=A^4/B^5=\frac{\theta^3(1+a\theta)^4}{(1+c\theta)^5}
\]
have orders \((-1,3)\) at \(\theta=0\). No other place
restricts to the same place of \(k(R,S)\). At a finite
nonzero pole of \(R\), if \(c\ne e\) its pole order is
eight and the order of \(S\) is zero or four. If \(c=e\)
its pole order is one and \(S\) has a pole. At infinity
\(R\) can have a pole only if \(e=0,c\ne0\). Its pole
order is then six and the order of \(S\) is two or minus
two. Since the pole at zero is simple and all residue
fields are \(k\), the degree formula for places yields
\(k(R,S)=k(\theta)\). Now \(B/C=x(1+c\theta)/(1+e\theta)\)
recovers \(x\), and then \(z\).

Increase \(w_y\) from one to the first tie, at most two
because \(x^2y^4\) occurs in the first output. For a
homogeneous prime write \(u=\operatorname{ord}x\) and
\(v=\operatorname{ord}z\) for the rational lifts.
If \(u<0\) and \(v>3u\), regularity of \(B\) gives
\(v\ge-5u\). Then \(A\) has order at least \(-6u\),
impossible since \(x z^2(x^3+a z)\) has no factor of
multiplicity six. If \(u<0\) and \(v<3u\), the term of
\(A\) of smallest valuation is unique and has a pole.
If \(u<0\) and \(v=3u\), cancellation in all three
outputs requires \(a=c=e\ne0\). Regularity of \(A\)
forces \(\operatorname{ord}(x^3+a z)\ge-7u\), so \(C\)
has order at least \(-3u\). Its binary specialization
\(xz(x^3+a z)\) is squarefree, which is again impossible.
These specialization arguments use that the relevant
output restricts to a nonzero polynomial. Thus the
weight-one lift of \(x\) is polynomial and is \(x\)
itself.

If \(a\ne0\) or \(e\ne0\), the lift of \(z\) is
polynomial, since at any pole of \(z\) the \(a xz^3\)
or \(e xz^2\) term has a strictly smaller and negative
valuation than the other term in its output. Here
\(\operatorname{ord}x\) is zero or one. A change to \(z\)
of weight three can only add a multiple of \(y^2\) at
\(w_y=3/2\) or of \(xy\) at \(w_y=2\).
It gives a term of degree eight in \(x^4z^2\), while
the substituted \(a xz^3\) has degree at most seven.
No first tie is possible, so \(a=e=0\).

If \(c\ne0\), the only possible pole of \(z\) is along
\(x=0\) by \(C=x^4z\). The identity for \(B\) bounds it
by a simple pole, so \(U=xz\) lifts polynomially.
Its weight is four. Its possible corrections containing
\(y\) are \(y^3\) at \(4/3\), \(xy^2\) at \(3/2\),
and \(x^2y,y^2\) at two. Any ordinary cubic correction
would give a term of degree eight in \(A=x^2U^2\).
The only first tie is therefore at two and the lift is
\(U=xz+y^2=P\), with coefficient fixed by the ordinary
term \(2x^3y^2z\) of \(A\).
At \((1,2,3)\) we would have highest triple
\((x^2P^2,x^4P+cP^2,x^3P)\).

For this putative triple, lower the weight \(s\) of \(x\)
from one, with \(y,z\)-weights two and three, comparing
to the reference \((x^2y^4,c y^4,x^3y^2)\).
The only ties for \(0<s<1\) are
\[
\begin{array}{c|lll}
s&A&B&C\\ \hline
2/3&&&z^2,\ y^3\\
1/2&z^3,\ y^3z&&xyz\\
1/3&&&yz .
\end{array}
\]
At \(s=2/3\) rank removes \(z^2\). If \(C=x^3y^2+d y^3\)
with \(d\ne0\), the binary image field is \(k(x,y)\).
The first two entries recover \(x^2,y^4\), and no
nonidentity element of the group \(\mu_2\times\mu_4\)
fixes \(C\), by comparing its two monomials.
Use scaled weights \((2,6,9)\) and raise \(w_z\) to
the first tie, at most ten by the term \(x^4z^2\).
The lift of \(y\) is polynomial by \(y^4=B/c\) and
cannot change by weight. The identity \(A=x^2y^4\)
allows the lift of \(x\) poles only along \(y=0\),
where \(x^3y^2+d y^3\) would then have a pole.
Thus \(x\) also lifts polynomially and cannot change.
This contradicts the tie. At \(s=1/2\) fix \(y\)
and use the bracket of the resulting first and third
outputs in \(x,z\). If the tying coefficients are
\(a',b',c'\), respectively, the \(z^2,z^1,z^0\)
coefficients, after nonzero monomial factors are removed,
give \(a'=0,\ b'c'=0,\ 2c'-3b'=0\).
At \(s=1/3\) the two unchanged outputs already have
rank two in \(x,y\), excluding \(yz\).

At \(s=0\) the first two outputs have weight at most eight
with \(y^4\)-coefficients of degrees two and zero, the latter
equal to \(c\). The third has form \(q(x)y^2+u(x)y+v(x)z+w(x)\)
with \(\deg q=3\) and \(\deg v=4\). Apply the proof of the
[two quartic transverse parts lemma](#an-affine-output-below-two-quartic-transverse-parts)
with initial weights \((0,2,3)\) instead of \((0,3,4)\).
The only initial extra highest term is \(yz^2\),
and there is no intervening tie before \(w_z=4\).
The comparison at four is again the weight-\((0,1,2)\)
comparison, so the rest of that proof applies unchanged.
This excludes \(c\ne0\) and completes the axial analysis
through \(w=3\).


### Axial weights four and five

At \(w=4\) the first and third outputs are already independent
in \(x,z\). Thus the possible \(yz^2\) term of \(B\) vanishes.
The binary triple is
\[
      (x^4z^2+a z^3,\quad x^5z+b xz^2,\quad x^4z+c z^2).
\]
Increase the \(y\)-weight \(t\) through \(5/2\). The possible
ties from the starting support and the axial exclusions are
\[
\begin{array}{c|lll}
t&A&B&C\\ \hline
5/3&&y^3z&\\
2&y^2z^2,\ x^2yz^2&xy^2z,\ x^3yz&y^2z,\ x^2yz\\
9/4&&y^4&\\
7/3&xy^3z&x^2y^3&xy^3\\
12/5&y^5&&\\
5/2&x^2y^4,\ x^3y^2z&y^2z,\ x^4y^2&x^3y^2 .
\end{array}
\]
At each row with only one monomial the other two binary
entries are independent, as seen from their two-variable
bracket at the smallest \(z\)-power. This removes that tie.

At \(t=2\), write
\[
\begin{aligned}
 A&=z^2(x^4+d y^2+f x^2y+a z),\\
 B&=xz(x^4+g y^2+h x^2y+b z),\\
 C&=z(x^4+i y^2+\ell x^2y+c z).
\end{aligned}
\]
The determinant equations imply that all \(y\)-coefficients
are zero, or that \(a=b=c=g=i=0,\ \ell=h,\ f=2h,\ d=h^2\).
For the algebraic verification, let \(I\) be the ideal of the
determinant coefficients over \(\mathbb Q\). The checker
reduces each product of a generator of \((d,f,g,h,i,\ell)\)
and a generator of \((a,b,c,g,i,\ell-h,f-2h,d-h^2)\)
to zero modulo a Gröbner basis of \(\sqrt I\). Over the
algebraic closure this forces one of the two stated cases.
In the second case, put \(Q=x^2+h y\). If \(h\ne0\), the
last two entries are \(B=x^3 z Q\), \(C=x^2 z Q\).
They fix \(x\) and \(zQ\) on the generic fibre and hence
the function field of the fibre is \(k(B,C)(z)\).
The relative Jacobian differential for that pair is
\[
     \frac{dz}{J(B,C,z)}=\frac{1}{C x^2 h}\frac{dz}{z}.
\]
Its residue at \(z=0\) is nonzero, contradicting the
[leading-pair rational slice](KELLER_LEADING_PAIR_RATIONAL_SLICE.md).
Thus the tie at \(t=2\) disappears.

At \(t=7/3\) use \(Y=y^3\) and denote the three new
coefficients by \(d,e,f\). The determinant in \(x,Y,z\)
has coefficients \(-d-4e+6f,\ -5de+5df+2ef,\ def\)
at \(x^9z^2,\ x^6Yz,\ x^3Y^2\), respectively.
Their vanishing forces \(d=e=f=0\).

At \(t=5/2\) the ordinary leading terms fix the
coefficients that recombine into \(P\). The tuple in
\(x,P,z\) is
\[
 (x^2P^2+a z^3,\quad x^4P+b' xz^2+g Pz,\quad x^3P+c z^2)
\]
where \(b'=b-g\). The change to \(x,P,z\) has determinant
\(2y\), so the determinant in \(x,P,z\) vanishes. Its
coefficients at \(x^6Pz^2,\ x^5P^2z,\ x^4P^3,\ x^2Pz^3\)
give \(3a=0,\ 4(2b'-3c)=0,\ 4g=0,\ -9ag-4b'c=0\).
Consequently \(a=b'=c=g=0\), proving that at
\((1,5/2,4)\) the highest tuple is
\((x^2P^2,x^4P,x^3P)\).

The latter weights bound the transverse degrees by
\((2,2,1)\), with the coefficient of \(z^2\) in the
second output at most constant in \(x,y\). If this
coefficient is \(\lambda\ne0\), the highest triple at
\((1,1,5)\) is
\((x^4z^2,x^5z+\lambda z^2,x^4z)\). Its image field is
again \(k(x,z)\), by \(z=A/C\) and recovery of \(x^4,x^5\).
Raise the \(y\)-weight to the first tie, at most three
by \(x^2y^4\). At an \(x\)-pole of order \(m\), \(C\)
first forces \(\operatorname{ord}z\ge4m\), and \(B\)
then forces \(\operatorname{ord}z\ge5m\). The latter
gives order at least \(6m\) for \(A\) and contradicts
binary specialization. Thus the lift of \(x\) is
polynomial and unchanged. At any \(z\)-pole, the term
\(\lambda z^2\) now dominates \(B\), so \(z\) also
lifts polynomially. A change in that lift must have
ordinary degree at least two, since no linear monomial
in \(x,y\) has weight five. Such a change gives ordinary
degree at least eight in \(A=x^4z^2\). There can be no
first tie. This proves that the transverse degrees are
\((2,1,1)\).

### The exact bracket on the third-output fibres

Let \(a_2,b_1,c_1\in k[x,y]\) be the respective highest
\(z\)-coefficients. The weight-\((1,5/2,4)\) bound makes
each at most linear in \(y\). On \(y=0\) their respective
\(x\)-degrees are \(4,5,4\).
The highest \(z\)-coefficient of the determinant of the
actual map vanishes, so the triple
\((a_2 z^2,b_1 z,c_1 z)\) has rank at most two.
Fixing \(c_1 z\) gives algebraic dependence of
\(a_2/c_1^2\) and \(b_1/c_1\) in \(x,y\).

If \(T=b_1/c_1\) depends on \(y\), express \(a_2/c_1^2\)
as a polynomial over \(k(x)\) in \(T\) of degree at
most two. Indeed, for \(c_1\) nonconstant in \(y\),
both \(1/c_1\) and \(a_2/c_1\) are affine in \(T\)
by linear fractional inversion. For \(c_1\) independent
of \(y\), \(T\) and \(a_2/c_1^2\) are affine in \(y\).
The Jacobian in \(x,T\) then forces every coefficient
of this polynomial to be constant. But on \(y=0\)
its value has order four at \(x=\infty\), while
\(T\) has a pole there. This is impossible. Thus
\(T\in k(x)\), and \(T\) is nonconstant by its pole
at infinity. The same Jacobian forces \(a_2/c_1^2\)
to be independent of \(y\). Since \(a_2,c_1\) have
\(y\)-degree at most one and \(a_2\ne0\), it follows
that \(c_1\) and hence both \(a_2,b_1\) belong to
\(k[x]\).

On a level \(S\) of the actual third output solve
for \(z\) over \(k(x)\). The other outputs become
\[
 G=p(x)S^2+r(x,y)S+u(x,y),\qquad
 H=q(x)S+v(x,y),
 \qquad p=a_2/c_1^2,\quad q=b_1/c_1,
\]
where \(r,u,v\in k(x)[y]\). Their \(x,y\) bracket at
fixed \(S\) equals \(\det JF/c_1\), a nonzero function
of \(x\). The \(S^2\) term gives \(r=h v+i(x)\) with
\(h=p'/q'\). The \(S\) term gives
\[
             u=\frac{h'}{2q'}v^2+\frac{i'}{q'}v+j_0(x).
\]
The constant term in \(S\) is then
\[
      \left[\left(\frac{h'}{2q'}\right)'v^2+
             \left(\frac{i'}{q'}\right)'v+j_0'\right]v_y .
\]
It is nonzero, so \(v\) has positive \(y\)-degree, and
it has \(y\)-degree zero, forcing
\((h'/q')'=0\). In characteristic zero this gives
\(h=\alpha q+\beta\) and \(p=\alpha q^2/2+\beta q+\gamma\)
for constants \(\alpha,\beta,\gamma\).
The order of \(p\) at \(x=\infty\) is four and \(q\)
has a pole there, a contradiction. This proves \(j=0\).


## The second-output mixed correction

**Theorem SexticKellerDoubledSquareSecondMixedCorrectionInversion.**
In the doubled square-member weight-\((2,3,4)\) reduction,
if \(k\ne0\), the map is invertible.

By the preceding exclusions, \(j=e=0\).
Work over the algebraic closure and scale the second output by \(1/k\).
Write \(A,B,C\) for the highest entries at the weight under
consideration. At weights \((s,3,4)\) just below \(s=2\) they are
\((x^2y^4,xy^4,x^6)\).
The following are all possible ties as \(s\) decreases strictly
between two and zero. In this table the second entry has been
scaled as above.
\[
\begin{array}{c|lll}
s&A&B&C\\ \hline
11/6&&&yz^2\\
9/5&&&xy^3\\
7/4&&&x^2yz\\
5/3&&&y^2z\\
8/5&&&xz^2\\
3/2&yz^3,\ y^5&&y^3,\ x^2y^2,\ x^4y\\
7/5&&&xyz\\
4/3&&&z^2,\ x^3z\\
6/5&&&xy^2\\
7/6&&&yz\\
1&y^2z^2,\ xy^3z&y^3z&y^2,\ x^2z,\ x^3y\\
4/5&&&xz\\
3/4&&&x^2y\\
2/3&&&z\\
3/5&&&xy\\
1/2&y^3z&&y
\end{array}
\]
The [exact checker](check_sextic_doubled_square_second_mixed.py)
enumerates the support using ordinary degree at most six and
the stated weight-\((2,3,4)\) highest tuple. It also checks the
subsequent tie lists and the determinant coefficients used below.
Run it from the repository root with
`.venv/bin/python research/jacobian/monodromy-theory/check_sextic_doubled_square_second_mixed.py`.
The internal replay used SymPy 1.14.0.

At each tie here and below the expected determinant weight
is positive, so the highest entries have rank at most two.
When only \(C\) can tie and the tying monomials involve \(z\),
the first two entries fix \(x,y\) algebraically and all those
coefficients vanish.

### Binary ties

When a binary tie below changes the leading tuple, increase
the weight of the unused variable \(z\) to the first further tie.
In a frame \((s,3,4)\) with \(A\) containing \(x^2y^4\), this
occurs at or before \(w_z=6-s\), by the ordinary term \(2x^3y^2z\).
The expected determinant weight remains positive to this bound.
The leading image field at the tie identifies with its binary
specialization at \(z=0\) by the rank argument of the preceding
section. Use this identification for the rational lifts below.

First suppose \(A=x^2y^4,\ B=xy^4\) and
\(C=x^6+\lambda x^a y^b\) with \(\lambda\ne0\) and \(b\) odd.
This covers the isolated binary ties at \(9/5,3/4,3/5,1/2\),
once the possible term of \(A\) at \(1/2\) is removed. At that
last weight \(B,C\) are independent in \(x,y\), by the contribution
of \(x^6\) to their bracket, and thus the term \(y^3z\) of \(A\)
vanishes. In each odd case \(A/B=x\) and \(B/x=y^4\), and the
last entry recovers \(y^b\), so the binary image field is \(k(x,y)\).
If the lift of \(x\) has a pole, \(B\) forces the lift of \(y\)
to have a zero there. The \(x^6\) term then gives a pole of \(C\),
since \(a<6\). Thus the polynomial lift of \(x\) is unchanged by
weight. A pole of \(y\) could occur only at \(x=0\), where
\(\operatorname{ord}x=1\) and \(xy^4\) would have a pole.
Thus \(y\) also lifts polynomially and is unchanged by weight.
This contradicts the further tie.

At \(s=3/2\), \(B,C\) are independent in \(x,y\), since the
term \(x^6\) of \(C\) gives a nonzero monomial of their bracket
with \(B=xy^4\). This removes \(yz^3\) from \(A\).
The remaining possible binary tie is
\[
   (x^2y^4+a y^5,\quad xy^4,\quad
                 x^6+c y^3+d x^2y^2+e x^4y).
\]
If \(a\ne0\), put \(\theta=y/x^2\) and \(T=A/B\).
Two weight-zero ratios are
\[
     R=C/T^6=\frac{1+e\theta+d\theta^2+c\theta^3}{(1+a\theta)^6},
       \qquad S=B/T^9=\frac{\theta^4}{(1+a\theta)^9}.
\]
The point \(\theta=0\) is the only point above its place in
\(k(R,S)\), since the only other zero of \(S\) is infinity,
where \(R=0\) instead of one. The degree of the extension
therefore divides four. The unique pole of \(S\) has order
nine, so the degree also divides nine. Consequently
\(k(R,S)=k(\theta)\), and \(T=x(1+a\theta)\) recovers \(x\)
and then \(y\). For their lifts an \(x\)-pole forces positive
order of \(y\) by \(B=xy^4\), and then \(x^6\) gives \(C\)
a pole. A \(y\)-pole of order \(m\) forces
\(\operatorname{ord}x\ge4m\) by \(B\), so \(a y^5\)
gives \(A\) a pole. Thus both lifts are polynomial and
unchanged, contradicting the further tie. This excludes \(a\ne0\).
With \(a=0\), if \(c\ne0\) or \(e\ne0\), the binary field is
again full: the expression \(C\) has a trivial stabilizer
in the fourth-root extension obtained by adjoining \(y\) to
\(k(x,y^4)\). The same \(x\)-pole argument applies. Once \(x\)
has lifted unchanged, \(B\) excludes a pole of \(y\).
Thus only \(d x^2y^2\) in \(C\) can survive this tie.

At \(s=1\) write the tuple as
\[
\begin{aligned}
 A&=x^2y^4+a y^2z^2+b xy^3z,\qquad
 B=xy^4+c y^3z,\\
 C&=x^6+d y^2+e x^2z+f x^3y .
\end{aligned}
\]
The determinant coefficients at \(x^7 y^6,x^6y^5z,x^4y^7,
x^2y^5z^2,xy^8\) respectively give \(b=c,\ a=0,\ e=cf/4,\
c^2 e=0,\ cd=0\), substituting the earlier equations as needed.
If \(c=0\), this is a binary tie with \(e=0\).
When \(f\ne0\) the field is full by the fourth-root argument,
and regularity of the lifts follows as at \(3/2\) with \(a=0\).
Thus for \(c=0\) only \(d y^2\) can survive.

The quadratic binary corrections in \(C\) remaining at
\(s=3/2,6/5,1\) are of the common form
\[
               C=x^6+\lambda x^m y^2,\qquad m=2,1,0
\]
respectively. Assume \(\lambda\ne0\).
The binary field is \(k(x,Y)\) with \(Y=y^2\).
At an \(x\)-pole, \(B=xY^2\) again forces \(Y\) to vanish,
leaving a pole in \(C\). Thus \(x\) lifts as the unchanged
polynomial \(x\), and \(B=xY^2\) excludes any pole of \(Y\).
The lift of \(Y\), of weight six, can acquire a term containing
\(z\) at a weight \(4<w_z\le6-s\) only in the form of a
multiple of \(xz\), at \(w_z=6-s\). The ordinary coefficient
of \(x^3 y^2 z\) in \(A=x^2Y^2\) fixes that multiple to one.
The new highest tuple is
\[
                      (x^2P^2,xP^2,x^6+\lambda x^mP).
\]

### The three quadratic corrections

For \(m=2\) scale the weights of the last tuple to \((1,2,3)\)
and the third output by \(1/\lambda\). Lower \(w_x=s\)
from one. Relative to \((x^2y^4,xy^4,x^2y^2)\) the only
positive intermediate tie is at \(s=1/2\):
the first entry can add \(a z^3+b y^3z\), the third \(e yz\).
The determinant coefficients at \(x^2 y^5 z^2,x^2y^8,y^7z\)
force \(a=0,\ 3b=2e,\ be=0\). Thus the tie disappears.
At \(s=0\) the first two entries have weight at most eight,
with \(y^4\)-coefficients \(p,r\) of degrees two and one.
The third output is
\[
                         q(x)y^2+u(x)y+v(x)z+w(x)
\]
with \(\deg q=2,\ \deg v=3\).
Use the proof of the
[two quartic transverse parts lemma](#an-affine-output-below-two-quartic-transverse-parts)
starting at weights \((0,2,3)\).
Here the orders of \(p/q^2,r/q^2\) at infinity are two and
three, so the same inequalities \(0<2<3<4\) distinguish all
six orders of the monomials of degree at most two.
The only extra initial term is \(yz^2\), with no intervening
tie before the weight-\((0,1,2)\) comparison. The remainder
of the lemma's proof uses the stated orders and \(\deg r\le1\).
This excludes \(m=2\).

For \(m=0\), the weights \((1,3,5)\) of the last tuple already
put the third output in the displayed affine form, with \(q\)
constant nonzero and \(v\) linear nonconstant.
The [affine-output root restriction](#residues-and-coordinate-completion-for-the-affine-transverse-output)
requires each root of \(v\) to be a root of \(q\), which is
impossible.

For \(m=1\), keep \(y,z\)-weights \(3,24/5\) and lower
\(w_x=s\) from \(6/5\), scaling the third output by \(1/\lambda\).
Relative to \((x^2y^4,xy^4,xy^2)\) the possible ties are
\[
\begin{array}{c|lll}
s&A&B&C\\ \hline
9/10&y^3z&&\\
3/5&xy z^2&y z^2&\\
3/10&y z^2&& .
\end{array}
\]
At the isolated ties the unchanged \(B,C\) are independent
in \(x,y\), so the tying coefficient is zero.
At \(s=3/5\) write the corrections as \(a xy z^2,b y z^2\).
On a fixed generic level \(C=t\), the first two entries
are \(A=t^2+a t z^2/y\), \(B=t y^2+b y z^2\).
If \(a\ne0\), eliminating \(z^2\) makes
\(B=y^2(t+b(A-t^2)/(at))\) nonconstant at a generic fixed
value of \(A\), contradicting rank. Thus \(a=0\).

If \(b\ne0\), the leading pair \(B,C\) is independent.
Over its generic values the fibre has equation
\(z^2=\lambda_0/y+\mu_0 y\), with
\(\lambda_0=B/b,\ \mu_0=-t/b\).
The relative Jacobian differential is \(dy/(2b y^3 z)\).
Modulo exact differentials, \(dy/(y^3 z)\) is the nonzero
constant multiple \(-\mu_0/(3\lambda_0)\) of \(dy/(y z)\),
since
\[
       d(z/y)=-\frac{3\lambda_0 y^{-3}+\mu_0 y^{-1}}{2z}\,dy .
\]
The latter form is regular on the complete normalization.
The double cover is branched at zero, infinity and the two
simple zeros of \(\lambda_0/y+\mu_0 y\). At zero the orders
of \(y,z,dy\) are \(2,-1,1\), at infinity they are
\(-2,-1,-3\), and at either simple zero they are \(0,1,1\).
Every other point is regular for the displayed form.
It is nonzero, so it cannot have a rational primitive on
the complete curve: the derivative of a function with a
pole has a pole, and a regular function on the complete
curve is constant. This contradicts the
[leading-pair rational slice](KELLER_LEADING_PAIR_RATIONAL_SLICE.md)
and excludes \(b\ne0\).

At \(s=0\) the third output now has the displayed affine form
with \(q\) linear and \(v\) quadratic. Using the already
established frame at \((1,5/2,4)\), the other coefficients
have \(\deg u\le3,\ \deg w\le6\).
The affine-output root restriction gives \(v\) proportional
to \(q^2\), and \(u\) is a unit modulo \(v\), since a common
root would force a multiple root of \(q\).
Modulo \(v\) the substitution \(y\mapsto qy^2+uy+w\) has a
polynomial inverse in the argument \(t\) of degree at most
two by \(q^2=0\). In fact it is
\((t-w)/u-(q/u)((t-w)/u)^2\) in this quotient.
Lift its coefficients to degree at most one in \(x\),
giving \(L(x,t)\) of total degree at most three.
As in the coprime coordinate completion above,
\[
                 (x,C,(y-L(x,C))/v)
\]
is a polynomial coordinate system.
Its inverse components in \(y,z\) have degrees at most
three and five: use \(y=L+vh\) and
\(z=(t-qy^2-uy-w)/v\) in coordinates \((x,t,h)\).
The numerator in the last quotient has degree at most seven
and is divisible by \(v\).
Consequently the composed Keller map fixes \(t\) and the
two other outputs have degree at most thirty. The
[planar specialization and descent argument](LOW_DEGREE_COORDINATE_KELLER_MAPS_IN_THREE_VARIABLES.md)
with the retained planar bound proves invertibility over the
original field. This settles the quadratic correction for \(m=1\).

### Completion of the lower-weight analysis

If \(c\ne0\) at \(s=1\), the coefficient identities give
\(a=d=e=f=0,\ b=c\). For \(s<1\) just below that tie the
reference tuple is \((cxy^3z,cy^3z,x^6)\).
The only further positive ties are in the third output:
\(xz,\ x^2y,\ z,\ xy,\ y\) at \(4/5,3/4,2/3,3/5,1/2\),
respectively. At each, the first two entries fix \(x\)
and \(y^3z\). The tying monomial is nonconstant in \(y\)
at fixed \(x,y^3z\), so its coefficient is zero.
Thus in this case \(C\) comes from an actual third output
depending only on \(x\). If \(c=0\) at \(s=1\) and no
quadratic correction survives anywhere, the binary
exclusions likewise leave \(x^6\) as the highest third
output down to arbitrarily small positive \(s\). This
again forces the actual third output to depend only on \(x\).
Its \(x\)-degree is six, so its derivative divides the
Jacobian determinant and cannot be a nonzero constant.
This contradiction exhausts the case \(k\ne0\).


## Axial weights of the uncorrected doubled frame

**Theorem SexticKellerDoubledSquareUncorrectedAxialThroughThree.**
Suppose the highest tuple at \((2,3,4)\) is
\((x^2P^2,x^4P,x^6)\).
Then for \((1,1,w)\) with \(2<w<3\) the highest tuple is
\((x^4z^2,x^5z,x^6)\), and at \(w=3\) it has the binary form
\[
                 (x^4z^2+a xz^3,\ x^5z+c x^2z^2,\
                           x^6+f x^3z+e z^2).
\]

The support bound at \((2,3,4)\), together with ordinary
degree at most six, gives the following alternatives.
\[
\begin{array}{c|lll}
w&A&B&C\\ \hline
5/2&&z^3&xz^2,\ yz^2\\
3&xz^3,\ yz^3&x^2z^2,\ xyz^2&z^2,\ x^3z,\ x^2yz
\end{array}
\]
Here \(A,B,C\) denote the highest entries and the reference
tuple for the table is \((x^4z^2,x^5z,x^6)\).
The support comparisons and determinant coefficients are
reproduced by the [exact checker](check_sextic_doubled_square_uncorrected_axial.py),
run from the repository root with
`.venv/bin/python research/jacobian/monodromy-theory/check_sextic_doubled_square_uncorrected_axial.py`.
The internal replay used SymPy 1.14.0.

At \(w=5/2\), independence of the first two binary entries
forces the \(yz^2\) coefficient of \(C\) to vanish.
Suppose the remaining coefficients \(b,c\) in
\[
             (x^4z^2,\ x^5z+bz^3,\ x^6+c xz^2)
\]
are not both zero. Put \(\theta=z^2/x^5\).
The zero-weight ratios are
\[
          B^2/(AC)=\frac{(1+b\theta)^2}{1+c\theta},
          \qquad C^3/A^2=\frac{(1+c\theta)^3}{\theta^2}.
\]
They recover \(\theta\). For \(b,c\ne0\) they have coprime
degrees two and three if \(b\ne c\), and the first is linear
if \(b=c\). For \(b=0,c\ne0\) the first has degree one.
For \(b\ne0,c=0\) the second generates \(k(\theta^2)\)
and the first is not even in \(\theta\).
Consequently the binary image field recovers \(X=x^3\) from
\(A/C\) and \(\theta\), and \(V=x^2z\) from \(B=X V(1+b\theta)\).
The three entries belong to \(k(X,V)\), with equations
\[
            A=V^2,\qquad B=XV+b V^3/X^2,\qquad C=X^2+c V^2/X .
\]

Raise the weight of \(y\) from one to its first tie, at most
\(7/4\) by \(x^2y^4\) in \(A\).
The expected determinant weight stays positive, so the
binary image field lifts as in the preceding sections.
The lift of \(V\) is polynomial of ordinary degree at most
three by \(A=V^2\). The lift of \(X\) is polynomial by the
equation for \(C\): at an \(X\)-pole the \(X^2\) term has a
strictly smaller valuation than \(c V^2/X\).
The only possible changes to \(X\), of weight three, add
a multiple of \(y^2\) when \(w_y=3/2\).
The possible changes to \(V\), of weight \(9/2\) and
degree at most three, add a multiple of \(y^3\) at \(3/2\)
or \(xy^2\) at \(7/4\).
The \(y^3\) term would introduce a forbidden ordinary \(y^6\)
in \(A\). At \(3/2\) the changed \(X=x^3+\delta y^2\), with
\(\delta\ne0\), is coprime to the unchanged \(V=x^2z\),
contradicting \(X^2\mid V^3\) if \(b\ne0\) and \(X\mid V^2\)
if \(c\ne0\). At \(7/4\) the ordinary term \(2x^3y^2z\) in
\(A=V^2\) forces \(V=xP\) while \(X=x^3\). This again
contradicts the same divisibility for either nonzero
coefficient, by the order along \(x=0\). Thus no first tie
is possible for \((b,c)\ne(0,0)\).

At \(w=3\) write
\[
\begin{aligned}
 A&=x^4z^2+(a x+b y)z^3,\qquad B=x^5z+(c x^2+d xy)z^2,\\
 C&=x^6+e z^2+(f x^3+g x^2y)z .
\end{aligned}
\]
The determinant coefficients at \(x^{10}z^3\) and \(x^6 yz^4\)
give \(b-2d+g=0,\ bd-2bg+dg=0\). Thus \(b=d=g\).
If their common value \(\beta\) is nonzero, the coefficients
at \(yz^6,\ x^7z^4,\ x^4 z^5\) force \(e=0,\ a+f=2c\)
and \(a c-2a f+c f=0\), hence \(a=c=f\).
In that case the tuple is
\[
              (A,B,C)=W(z^2,xz,x^2),\qquad
                     W=x^4+a xz+\beta yz .
\]

To exclude this case put \(t=x/z\). The field generated by
\(A,B\) is \(k(A,t)\), relatively algebraically closed in
the full source field \(k(A,t)(z)\), since
\[
                  y=(A-t^4 z^6-a t z^4)/(\beta z^3).
\]
Its intersection with \(k[x,y,z]\) is precisely \(k[A,B,C]\).
In fact an element of the intersection is polynomial in
\(A\) over \(k(t)\): in \(k(t,z)[y]\) the variable \(A\)
is a linear coordinate, and a reduced denominator
depending on \(A\) would persist under this extension.
Write the resulting polynomial as \(\sum_j a_j(t)A^j\).
No \(a_j\) has a pole at a finite \(t=t_0\). Otherwise,
after multiplying by the maximum power of \(t-t_0\)
in those poles and reducing along \(x=t_0 z\), the terms
of maximal pole order give a nonzero polynomial in \(A\)
with constant coefficients. The reduction of \(A\) is
nonconstant, because its coefficient of \(y\) is
\(\beta z^3\), so this contradicts regularity.
Over the algebraic closure all the \(a_j\) are therefore
polynomials. Along \(z=0\) the term \(t^i A^j\) has
order \(2j-i\), with leading coefficient \(x^{4j+i}\).
Terms with a given order have distinct leading powers
of \(x\), so any negative order would give a pole.
Thus \(0\le i\le2j\) for each term, and each such monomial
is a product of \(A,tA=B,t^2A=C\).

Use the determinant identity for
\(D=F_1F_3-F_2^2-Q(F)\), where \(Q\) has output-weight
less than sixteen for weights \(10,8,6\).
The weighted degree of \(D\) is below sixteen since
the highest relation is \(AC=B^2\). We have
\[
                  J(F_1,F_2,D)=(F_1-\partial_3 Q(F))\det JF .
\]
The right-hand side is nonzero of weight ten, because
\(\partial_3 Q(F)\) has weight less than ten.
If \(D\) has positive weight \(d\), the expected determinant
weight on the left is \(13+d>10\).
Its highest part is then algebraic over \(k(A,B)\), hence
belongs to \(k[A,B,C]\) by the intersection just proved.
It can be removed by adding to \(Q\) a homogeneous output
polynomial of weight \(d\), and iterating decreases the
weight of \(D\) to zero. A constant \(D\) contradicts the
nonzero determinant identity. Thus \(\beta=0\), proving the
claimed binary form.


## Lifting the binary tuple at axial weight three

**Theorem SexticKellerDoubledSquareUncorrectedAxialThreeTransverseReduction.**
Under the uncorrected weight-\((2,3,4)\) hypotheses, write
the binary weight-\((1,1,3)\) highest tuple as in the
preceding theorem. Then \(a=e=0\). If the axial tuple is
corrected, so \((c,f)\ne(0,0)\), the highest tuple at
\((1,2,3)\) is
\[
               (x^2P^2,\ x^4P+cP^2,\ x^6+f x^2P).
\]

The [axial checker](check_sextic_doubled_square_uncorrected_axial.py)
verifies the possible terms in the polynomial lifts and the
Veronese fibre identities below. Work over the algebraic closure. For
\(\theta=z/x^3\) set \(M=1+f\theta+e\theta^2\). The two
zero-weight ratios are
\[
 R=\frac{B^2}{AC}
   =\frac{(1+c\theta)^2}{(1+a\theta)M},
 \qquad
 S=\frac{A^2}{BC^2}
   =\frac{\theta^3(1+a\theta)^2}{(1+c\theta)M^2}.
\]
Raise the \(y\)-weight from one to the first tie, at most
two because of the ordinary \(x^2y^4\) term in the first
entry. The expected determinant weight stays positive,
and the image of the first tie is the same as the binary
image by specialization at \(y=0\), as in the previous
lifting argument. Identities in the binary image field
lift to the same identities among rational functions
at the first tie.

### Nonconstant projective ratio

Suppose \(R\) is nonconstant. The place \(\theta=0\)
has \(R=1\) and \(S=0\) of order three. It is the only
place over its restriction to \(k(R,S)\). In fact the
only other possible finite zero of \(S\) is at
\(1+a\theta=0\). For \(R\) to be a unit there, both
\(1+c\theta\) and \(M\) must vanish to order one, making
\(S\) have a pole there. Infinity can have \(R=1,S=0\)
only if \(a=0\) and \(e=c^2\ne0\); then the order of
\(S\) there is two. If these places at zero and infinity
restrict to the same place, the order of \(S\) downstairs
divides three and two, so the two ramification indices
are three and two. Their sum exceeds \(\deg R\le2\).
Thus the field index divides three. If it is three, then
\(\deg R=3\), \(a,e\ne0\), there is no cancellation in
\(R\), and \(S\in k(R)\). When \(c\ne0\), the numerator
root of \(R\) and infinity are both zeros of \(R\).
At the former \(S\) has a pole and at the latter \(S\)
is finite. When \(c=0\), \(R\) has its sole zero at
infinity of order three, but \(S\) there has a pole of
order one. Both are impossible for \(S\in k(R)\).
We conclude that \(\theta\) is recovered.

The image field is therefore \(k(X,V)\) with \(X=x^2\)
and \(V=xz=\theta X^2\): indeed \(B/C\) recovers \(X\)
once \(\theta\) is known. The lifted identities are
\[
       A=XV^2+aV^3/X,\qquad B=X^2V+cV^2,\qquad
                  C=X^3+fXV+eV^2/X .
\]
Let \(u=\operatorname{ord}X<0\) at a putative pole of
\(X\), and write \(v=\operatorname{ord}V\) there.
If \(v>2u\), the \(X^3\) term of \(C\) dominates;
if \(v<2u\), the lowest-order nonzero term of \(A\)
has negative order and does not cancel. If \(v=2u\),
the expressions
\[
      A=X^5\theta^2(1+a\theta),\quad B=X^4\theta(1+c\theta),
                        \quad C=X^3 M(\theta)
\]
require a common zero at \(\theta=-1/a=-1/c\) of the
three parenthetical factors. Here \(e\ne0\), since
\(e=0\) would give \(a=c=f\), hence \(R=1\).
Regularity of \(A\) requires
\(\operatorname{ord}(1+a\theta)\ge -5u\).
If the root of \(M\) is simple, \(C\) has order at
least \(-2u\), contradicting its squarefree binary
specialization. If the root of \(M\) is double, \(C\)
has order at least \(-7u\), while its binary
specialization has maximal factor multiplicity two.
In either case any prime dividing \(C\) specializes
to a nonconstant factor at \(y=0\), since it is
homogeneous of positive weight and \(C|_{y=0}\ne0\).
This excludes an \(X\)-pole.

The polynomial lift of \(X\) has weight two, so it is
\(x^2\) or \(x^2+\delta y\) at \(y\)-weight two.
If at least one of \(a,c,e\) is nonzero, \(V\) is
polynomial: at a \(V\)-pole the corresponding term
\(aV^3/X,cV^2\) or \(eV^2/X\) is the unique term
of smallest, negative order in its output.
Otherwise \(A=XV^2\) allows a pole of \(V\) only
when \(X=x^2\), and then at most a simple pole along
\(x=0\). Put \(H=xV\) in that case, so that
\[
                  (A,B,C)=(H^2,x^3H,x^6+f xH).
\]
The polynomial lift \(H\), of weight five and
ordinary degree at most three, can gain \(y^3\) at
\(5/3\) or \(xy^2,yz\) at two. The former gives
a forbidden \(y^6\) in \(A\). The latter must be
\(H=xP+\mu yz\), with the coefficient of \(xy^2\)
fixed by the ordinary \(2x^3y^2z\) term in \(A\).
If \(\mu\ne0\), the pair \(A,B\) fixes \(H,x\)
algebraically and has relative differential proportional
to \(dy/(x^2+\mu y)\): on the fixed levels,
\(H=(x^2+\mu y)z+xy^2\).
The simple pole at \(x^2+\mu y=0\) has nonzero
residue, contradicting the leading-pair slice.
Thus \(H=xP\) and \(V\) is polynomial also in this
case.

A polynomial lift of \(V\), of weight four, differs
from \(xz\) only by possible \(y^3\) at \(4/3\),
\(xy^2\) at \(3/2\), or \(x^2y,y^2\) at two.
A cubic addition makes \(XV^2\) of ordinary degree
eight, while \(aV^3/X\), when nonzero, has degree
seven. This contradicts the degree bound on \(A\).
Thus the first tie is at two and \(V=P\), with the
coefficient of \(y^2\) fixed by the ordinary leading
part of \(A\). If \(a\) or \(e\) is nonzero, the
identities require \(X\mid V^3\) or \(X\mid V^2\).
This is impossible: \(P\) has no factor \(x\), and
does not vanish identically on \(x^2+\delta y=0\)
if \(\delta\ne0\).
We have \(a=e=0\). If \(X=x^2+\delta y\) with
\(\delta\ne0\), on the generic pair \(X,V\) one
can solve for \(y=(X-x^2)/\delta\) and then for \(z\)
from \(V=P\). The relative differential is proportional
to \(dx/x\), since \(J(X,V,x)=\delta x\).
For the actual leading pair \(A=XV^2,B=X^2 V+cV^2\)
the additional Jacobian in \(X,V\) is nonzero and
constant on each generic fibre. Again the residue is
nonzero. Thus \(X=x^2\), proving the displayed
reduction for \(R\) nonconstant.

### Constant projective ratio with a common factor

If \(R\) is constant, it equals one. The identity
\((1+c\theta)^2=(1+a\theta)M\) gives either
\(a=c=f\ne0,\ e=0\), or \(a=0,\ f=2c,\ e=c^2\).
In the first case the binary tuple is
\[
                        W(z^2,xz,x^2),\qquad W=x(x^3+a z).
\]
At the first tie its lift factors as
\(\widetilde W(Z^2,XZ,X^2)\), by the conic relation and
unique factorization. On \(y=0\), the parameters
specialize to a common multiple of \(z,x\),
whose square divides \(W\). Since \(W\) is squarefree,
normalize so that they specialize to \(z,x\) and
\(\widetilde W\) specializes to \(W\). The weights
of \(X,Z,\widetilde W\) are \(1,3,4\).
Thus \(X=x\), and any change of \(Z\) has ordinary
degree at least two, giving degree at least eight
in \(\widetilde W Z^2\). So \(Z=z\). A tie at two is
impossible because \(\widetilde W z^2\) cannot contain
the ordinary \(x^2y^4\) term. The only earlier
changes to \(\widetilde W\) are a multiple of \(y^3\)
at \(4/3\), or of \(xy^2\) at \(3/2\).

Put \(t=z/x=B/C\) on the generic fibre. In the
first case with coefficient \(\lambda\ne0\) the
curve and relative differential, up to a nonzero
fibre constant, are
\[
        y^3=\frac{C-x^6-a t x^4}{\lambda x^2},
                       \qquad \frac{dx}{x y^2}.
\]
In the second case they are
\[
        y^2=\frac{C-x^6-a t x^4}{\lambda x^3},
                       \qquad \frac{dx}{x^2 y}.
\]
Indeed \(C=x^2\widetilde W\), and replacing \(z\)
by \(tx\) makes the volume acquire a factor \(x\).
The generic numerator has six simple nonzero roots,
where ramification cancels the denominator involving
\(y\). On the triple cover the orders of \(x,y,dx\)
at zero are \(3,-2,2\), and at infinity they are
\(-3,-4,-4\). On the double cover those orders are
\(2,-3,1\) and \(-2,-3,-3\), respectively.
In both cases the indicated differential is regular
and nonzero on the complete normalization, so it
has no rational primitive. Since \(B=tC\), the same
obstruction applies to the leading pair \(B,C\).
This excludes the squarefree common-factor case.

### Constant projective ratio with coprime parameters

In the remaining constant-ratio case, \(c=0\)
gives the uncorrected axial tuple. Assume \(c\ne0\).
The binary tuple is the Veronese tuple on
\((U,V)=(x^2 z,x^3+c z)\). The parameters are
coprime, so the same factorization and specialization
argument lifts them to polynomial \(U,V\) of weights
five and three and ordinary degrees at most three.
An earlier tie in \(U\) could only add \(y^3\) at
\(5/3\), giving forbidden \(y^6\) in \(U^2\).
An earlier change to \(V\) is \(\delta y^2\) at
\(3/2\). On generic levels \(U=u,V=v\) this gives
\(z=u/x^2\), \(y^2=(v-x^3-cu/x^2)/\delta\), and the
relative form is proportional to \(dx/(x^2 y)\).
It has a simple pole with nonzero residue above
\(x=0\), since \(y\) has a simple pole there and the
double cover is unramified. It follows that the
first tie would have to be at two, with
\[
                  U=xP+\mu yz,\qquad V=x^3+c z+\nu xy.
\]
The ordinary \(2x^3y^2z\) term in \(U^2\) again
fixes the \(xy^2\) coefficient of \(U\).

On \(V=v\) eliminate \(z\). On a level \(U=u\)
the resulting equation and relative form are
\[
\begin{aligned}
 c u &=x^2(v-x^3)+K(x)y+Dxy^2,\\
 K(x)&=\mu v-(\mu+\nu)x^3,\qquad D=c-\mu\nu,
 \qquad \omega\ \propto\ \frac{dx}{K(x)+2Dxy}.
\end{aligned}
\]
For \(D\ne0\), the curve is a double cover with square
root of
\[
                 Q(x)=K(x)^2+4D(x^6-vx^3+c u x).
\]
Its degree is six or three: simultaneous vanishing of
the coefficients at six and three gives
\((\mu-\nu)^2+4c=0,\ \mu^2-\mu\nu+2c=0\).
Because \(c\ne0\) these imply \(\nu=-\mu,\ c=-\mu^2\),
contrary to \(D\ne0\).
The polynomial is generically squarefree. At \(x=0\)
it has no root unless \(\mu=0\), when the root is simple
by the coefficient \(4Dcu\). At nonzero \(x\), a multiple
root would require \(u\) to be a critical value of a
fixed nonconstant rational function of \(x\) over \(k(v)\),
by dividing the equation \(Q=0\) by \(x\) and using
the term linear in \(u\). This is impossible for
transcendental \(u\).
The form \(dx/\sqrt{Q}\) is regular and nonzero on the
complete double cover, including at infinity in either
degree, excluding \(D\ne0\).

If \(D=0\), then \(\mu\nu=c\ne0\). For
\(\mu+\nu\ne0\), the fibre is rational and \(dx/K(x)\)
has nonzero residues at the simple roots of \(K\).
For \(\nu=-\mu\) and \(c=-\mu^2\) put \(L=x^2+\mu y\).
The variables \(x,L,V\) are polynomial coordinates,
and \(cU=LV-x^5\). Now \(J(U,V,x)=\mu V\).
For the leading pair \(A=U^2,B=UV\), the rational
slice with derivative \(\kappa\ne0\) would be
\[
                       \frac{\kappa x}{2\mu U^2 V}
                                  +h(U,V)
\]
for some rational invariant \(h\), since the source
field is \(k(U,V)(x)\). The differentials of \(A,B\)
are independent at the generic point of \(V=0\),
where \(cU=-x^5\). A regular slice is therefore
required there by the divisor-regular version of the
[leading-pair slice](KELLER_LEADING_PAIR_RATIONAL_SLICE.md).
The pole coefficient of the displayed primitive,
\(\kappa x/(2\mu U^2)\) at \(V=0\), is outside
\(k(U|_{V=0})=k(x^5)\). Any pole coefficient of an
invariant lies in that field, by taking its Laurent
expansion in \(V\) and using that \(U|_{V=0}\) is
transcendental over \(k\). Thus no pole cancellation
is possible. This excludes the last constant-ratio
case with \(c\ne0\) and completes the proof.


## Excluding the isolated third-output correction at \((1,2,3)\)

**Theorem SexticKellerDoubledSquareUncorrectedAxialThreeIsolatedThirdExclusion.**
In the corrected \((1,2,3)\) tuple of the preceding theorem,
\(f\ne0\) implies \(c\ne0\).

Suppose \(c=0,f\ne0\), and scale the third output by
\(1/f\). Lower the \(x\)-weight \(s\) from one, with
\(y,z\)-weights two and three. Just below one the leading
tuple is \((x^2y^4,x^4y^2,x^2y^2)\). The possible
intermediate ties are
\[
\begin{array}{c|lll}
s&A&B&C\\ \hline
3/4&&y^2z&\\
2/3&&xz^2,\ xy^3&\\
1/2&z^3,\ y^3z&z^2,\ y^3,\ x^2yz&yz\\
1/3&&xyz&\\
1/4&&yz&
\end{array}
\]
These use the weight-\((2,3,4)\) and \((1,2,3)\)
supports, as verified along with the determinant
coefficients below by the
[axial checker](check_sextic_doubled_square_uncorrected_axial.py).
The expected determinant weight is positive throughout.
An isolated change to the second entry involving \(z\)
vanishes, since the unchanged first and third entries
are independent in \(x,y\). This also removes \(xz^2\)
at \(2/3\), even if \(xy^3\) occurs.

Suppose \(B=x^4y^2+\lambda xy^3\) at \(2/3\), with
\(\lambda\ne0\). The binary image field is \(k(Y,W)\)
with \(Y=y^2,\ W=xy\). In fact \(Y=A/C\),
\(W^2=C\), and \(B=W^4/Y+\lambda W Y\) recovers \(W\).
Raise the \(z\)-weight to the first tie, at most
\(10/3\) by the ordinary \(x^3y^2z\) term of \(A\).
At this tie the image identifies with its binary
specialization. The lift of \(W\) is polynomial by
\(W^2=C\), and has weight \(8/3\), so it is still \(xy\).
At any prime other than \(x,y\), the identity \(A=W^2Y\)
excludes a pole of \(Y\), and \(B=W^4/Y+\lambda W Y\)
excludes a zero. The rational lift of \(Y\) is
therefore a constant times a monomial in \(x,y\),
and the specialization at \(z=0\) makes it \(y^2\).
There is no first tie, a contradiction.

At \(s=1/2\) write
\[
\begin{aligned}
 A&=x^2y^4+a z^3+b y^3z,\qquad
 B=x^4y^2+p z^2+d y^3+e x^2yz,\\
 C&=x^2y^2+\ell yz.
\end{aligned}
\]
The determinant coefficients at \(x^5y^3z^2,x^5y^6,
x^3 y^5 z,x y^4 z^2\) give \(a=0,\ e=3\ell-b,\ 
p=3\ell^2-3b\ell+b^2\) and then
\(-3(b-\ell)^3=0\) after removing a nonzero scalar factor.
Thus \(b=\ell,e=2\ell,p=\ell^2\).

If \(\ell\ne0\), set \(H=x^2y+\ell z\) as a
polynomial coordinate in place of \(z\), of weight
three. The tuple is \((y^3H,H^2+d y^3,yH)\),
independent of \(x\). With \(y,H\)-weights two and
three, raise the \(x\)-weight to the first tie.
Because the map depends on \(x\), the tie occurs
at or before weight nine, the largest output weight.
The expected determinant weight stays positive
to this bound.

If \(d\ne0\), the binary image field is \(k(y,H)\),
by \(y^2=A/C,\ y^3=(B-C^2/y^2)/d\) and \(H=C/y\).
The lift of \(y\) is polynomial by
\(B y^2=C^2+d y^5\), and then \(H\) is polynomial by
\(H^2=B-d y^3\).
If \(d=0\), the lifted relation is \(AB=C^3\).
Each irreducible factor of \(B\) divides \(C\) and
specializes at \(x=0\) to a nonconstant power of \(H\)
up to a scalar. Since \(C|_{x=0}=yH\), there is
exactly one such factor, specializing to \(H\),
and its exponent in \(B\) is two by \(B|_{x=0}=H^2\).
Thus the lifted \(B\) is the square of a polynomial
specializing to \(H\), and that polynomial divides \(C\).
Taking this square root as the lifted \(H\), and the
quotient \(C/H\) as the lifted \(y\), gives the same tuple.

In either case the polynomial lifts have weights two
and three, and specialize to \(y,H\) at \(x=0\).
The changes allowed by weight are a polynomial in \(x\)
for \(y\) and a polynomial in \(x,y\) of \(y\)-degree
at most one for \(H\), with every new term involving
\(x\). The lifted \(y,H\) can therefore be taken as
new polynomial coordinates, preserving the same binary
tuple at the new \(x\)-weight. The possible \(x\)-weights
of all nontrivial changes above \(1/2\) form a finite
set: they are \(2/m,3/m,1/m\) with \(m\) a positive
integer. But a first tie must still occur at a larger
weight at each step by the Keller condition.
This contradiction excludes \(\ell\ne0\).

If \(\ell=0,d\ne0\), the binary image field at \(s=1/2\)
is \(k(x^2,y)\), by \(y^2=A/C\),
\(y^3=(B-C^2/y^2)/d\), and \(x^2=C/y^2\).
Raise the unused \(z\)-weight to the first tie,
at most \(7/2\) by the ordinary \(x^3y^2z\) term.
The lift of \(y\) is polynomial by \(B y^2=C^2+d y^5\).
The lift of \(H=x^2y\) is polynomial by
\(H^2=B-d y^3\). Their weights two and three force
both lifts to be unchanged, so \(x^2\) is also
unchanged. This contradicts the tie.
Thus the reference tuple persists through all
positive \(s<1\).

At \(s=0\), write the actual second and third outputs as
\[
\begin{aligned}
 F_2&=b(x)z+p(x)y^2+u(x)y+w(x),\\
 F_3&=d(x)z+q(x)y^2+r(x)y+h(x).
\end{aligned}
\]
Indeed each of their weights at \((0,2,3)\) is at most
four. Their transverse coefficients have degrees
\(\deg b=5,\ \deg d=3\) from the ordinary \(x^5z\)
term and the \((1,2,3)\) term \(x^3z\), respectively.
The first output has weight at most eight, with
\(y^4\)-coefficient of \(x\)-degree two.
In particular its \(z^2\)-coefficient has degree at
most four in \(x,y\).

On \(F_3=S\) solve for \(z\) over \(k(x)\) and write
\[
            G=p_2(x,y)S^2+p_1(x,y)S+p_0(x,y),\qquad
            H=\lambda(x)S+v(x,y),
                 \qquad \lambda=b/d
\]
for the restrictions of \(F_1,F_2\).
Their \(x,y\) bracket at fixed \(S\) is a nonzero constant
divided by \(d\).
The \(S^3\) coefficient gives
\(\lambda'\partial_y p_2=0\), so \(p_2\) depends only
on \(x\). The \(S^2\) and \(S\) coefficients give
\[
         p_1=g v+i(x),\qquad
         p_0=\frac{g'}{2\lambda'}v^2+
                   \frac{i'}{\lambda'}v+j(x),\qquad
         g=p_2'/\lambda'.
\]
The constant term is
\[
  \left[\left(\frac{g'}{2\lambda'}\right)'v^2
                +\left(\frac{i'}{\lambda'}\right)'v+j'\right]v_y .
\]
Since this is nonzero of \(y\)-degree zero,
\(v\) has positive \(y\)-degree and both displayed
derivatives of coefficients vanish. In particular
\(p_2\) is a constant-coefficient polynomial of degree
at most two in \(\lambda\).
But \(p_2\) is the \(z^2\)-coefficient of \(F_1\) divided
by \(d^2\), so it vanishes at infinity, while
\(\lambda\) has a pole there. Thus \(p_2=0\).
The equations now give \(g=0,\ i=\alpha\lambda+\beta\)
and \(G=\alpha H+\beta S+j(x)\) for constants
\(\alpha,\beta\). This is the actual identity
\(F_1-\alpha F_2-\beta F_3=j(x)\), impossible because
\(F_1\) has a nonzero \(y^4\)-coefficient and the
other two are at most quadratic in \(y\).
This excludes \(c=0,f\ne0\).
