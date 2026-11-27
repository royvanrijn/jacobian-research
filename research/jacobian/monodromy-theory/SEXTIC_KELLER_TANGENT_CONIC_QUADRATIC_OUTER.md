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
