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
