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
