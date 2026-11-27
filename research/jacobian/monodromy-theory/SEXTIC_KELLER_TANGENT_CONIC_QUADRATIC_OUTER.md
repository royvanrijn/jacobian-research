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
