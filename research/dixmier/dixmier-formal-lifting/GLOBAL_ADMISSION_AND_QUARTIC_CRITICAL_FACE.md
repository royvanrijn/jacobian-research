# Global admission and a quartic critical-face obstruction

Fix characteristic-zero constants, the quantum parameter h, and the exact
positional Long image R=X(2-3XQ). The global centralizer is
[QEP1's finite Ore algebra](EXACT_POSITION_CENTRALIZER_AND_COMPLETION.md)
K[[h]][X,Q][B;h Delta], where

\[
 B=3X^2P+2(1-3XQ)Z,\qquad
 \Delta=3X^2\partial_X+2(1-3XQ)\partial_Q.
\]

Global polynomial admission allows every negative Laurent support point
on the generic relative fiber. In particular it allows the first
centered quartic term on QCL1's critical boundary. The family constructed
below nevertheless has no finite scalar mate: its quartic critical face
has polynomial centralizer and a separate all-order mate obstruction.
This is a restricted extension of QCL1, not a nonexistence theorem for
all globally admitted negative Laurent corrections.

## Exact global coordinates and support realization

In [QPH1's chart](QUINTIC_PHASE_AND_BRANCH_OBSTRUCTION.md), set
r=R, v=X^{-1}. Its operator a satisfies [a,v]=h/v and commutes
with r. The exact identities, with coefficients on the left, are

\[
 \begin{split}
 F&=X^{-1}+Q+3XQ^2,\\
 A_0&=Q+3XQ^2-XB/3=a-v,\\
 B&=-3va+9v^2-5rv^3+r^2v^4.                 \tag{1}
 \end{split}
\]

Thus A_0 and C=R A_0-2/3 are global finite operators commuting
with R. No localization is needed to define them.
The exact globally admitted quadratic seed is

\[
 S_b=-3R A_0^2+(18XQ-8)A_0+9Q-3hRX
                       =4a-3ra^2+2v.                    \tag{2}
\]

Its reduction is Long's second classical image. Formula (2) follows by
expanding a=v+A_0 and using [A_0,v]=h/v, rv=2-3XQ and
6v-3rv^2=9Q. It asserts an exact operator identity; the displayed
h correction accounts for the ordering.

Over the generic constant field containing r, write

\[
 z=a-\frac2{3r},\qquad C=r(z-v),\qquad
 S_b-\frac4{3r}=-3rz^2+2v.                               \tag{3}
\]

The centering in (3) is precisely QCL1's critical normalization.
For integers n>=1, j>=0 and m<=0, the global monomial

\[
                         h^n X^{-m} C^j                 \tag{4}
\]

becomes h^n r^j v^m(z-v)^j. Under the finite rescaling
h=s^5, v=s^2V, z=sA it equals

\[
 s^{5n+2m+j} r^j V^m(A-sV)^j,\qquad [A,V]=1/V.           \tag{5}
\]

Every coefficient of the finite Ore expansion on the right is
polynomial in s. Its s^0 term is V^m A^j. Thus the least
QCL1 weight after dividing the position by s^2 is exactly
5n+j+2m-2, realized by the single normal monomial
h^n r^j v^m z^j. All other normal monomials have greater weight.
This proves the support assertion for all n,m,j; a finite scan is
not being used. Multiplication by the generic nonzero scalar r^j
does not restrict that support. Consequently global admission alone
cannot force every correction into QCL1's quadratic-face cone.

## A globally admitted quartic family

For any nonzero constant epsilon, consider the actual global operator

\[
                   S_\epsilon=S_b+\epsilon h^2X^6C^4.    \tag{6}
\]

It commutes with R, has the prescribed Long classical image, and
has generic momentum order four. After (3) and the rescaling in (5),
its finite integral rescaling has critical reduction

\[
 P=gA^2+\mu V+\rho V^{-6}A^4,\quad
       [A,V]=\eta/V,\qquad g\mu\rho\eta\ne0,             \tag{7}
\]

with g=-3r, mu=2, rho=epsilon*r^4 and eta=1.
Only the quartic monomial from the correction lies on that face;
all its ordering companions have strictly greater weight.
The first nonquadratic critical term has even parameter order h^2.
The exact chart expansion of the finite operator can also contain
higher odd powers of h from normal ordering; no global even-series
claim is inferred from the leading critical term.

The following theorem excludes (7) with arbitrary nonzero constants.
It also allows g=0, but that extra case is not needed here.

## A local divisor criterion for all finite partner orders

Use canonical Euler coordinates [H,V]=eta V, eta!=0, on the finite
rational differential algebra k(V)[b;eta partial_V], H=Vb.
Expand coefficients at V=0. Suppose a finite P has least Laurent
weight -N<0, with leading term V^{-N}C(H), where C is a nonzero
polynomial of degree j>=2. On some root-translation orbit
alpha+eta*Z, record the finite divisor of C by the Laurent polynomial
c(t). Assume c(xi)!=0 for a primitive Nth root xi. Then

\[
                         \operatorname{Cent}(P)=k[P],    \tag{8}
\]

and P has no finite rational-coefficient mate with nonzero scalar
commutator. The mate exclusion also holds over every finite
algebraic coefficient cover. Here the centralizer in (8) is in the
finite rational-coefficient differential algebra, not its division ring.

A nonzero finite operator T has least Laurent exponent m in Z and
a leading term V^m D(H), where D is a nonzero polynomial. Its degree
is at most ord_b T. If the coefficient at V^{m-N} in [P,T] must
vanish, exact multiplication gives

\[
              C(H+m\eta)D(H)=D(H-N\eta)C(H).              \tag{9}
\]

For a commuting T this condition always holds; for a nonzero
constant mate it holds whenever m<N. On the selected root orbit,
(9) gives the finite-divisor identity

\[
                  (t^{-m}-1)c(t)=(t^N-1)d(t).             \tag{10}
\]

Evaluating at xi forces N to divide m. Write e=deg D.
The next coefficient at H=infinity in (9) is
eta*(jm+Ne) times the nonzero product of leading coefficients.
Thus e=-jm/N, so m<=0. If m=-Nk, k>=0, the leading coefficient
of P^k satisfies the same difference identity. The rational ratio
with D is Neta-periodic, hence constant: finite root and pole
sets cannot have nonzero additive periods in characteristic zero,
and a nonconstant polynomial cannot be periodic either.

Subtract a scalar multiple of P^k to remove T's least Laurent term.
It strictly raises m. Repeating removes all negative m and then
the scalar term at m=0 in finitely many steps, because the integer
least exponents increase. No differential-order bound on P or T
is needed. Individual subtractions remain finite; their orders
need not stay below the original order of T. For a commuting
operator a nonzero remainder with positive m contradicts e=-jm/N.
This proves (8). If constants were extended to analyze roots, the
subtraction scalar descends by comparing a nonzero leading
coefficient, so (8) holds over the original constants.

For a nonzero constant mate, a residual 0<m<N contradicts (9)--(10).
If m>N, the entire commutator has positive V-valuation. At m=N,
its V^0 coefficient is

\[
                    C(H+N\eta)D(H)-D(H-N\eta)C(H).       \tag{11}
\]

This has degree j+e-1>=1, with nonzero leading coefficient
Neta*(j+e) times the product of leading coefficients.
It cannot be a scalar. This excludes every finite mate over k(V).
For a finite algebraic coefficient cover, apply the normalized
coefficient trace first, as in
[QPC1](UNIVERSAL_P_POWER_OBSTRUCTION.md): trace commutes with
the uniquely extended derivation, and returns a rational-coefficient
mate with the same scalar commutator.

There is a simple useful corollary in the chart [a,v]=eta/v.
If a finite Laurent polynomial has a unique least canonical weight
m-2j, realized by lambda*v^m*a^j with m<0,j>=2, it has polynomial
centralizer and no finite mate. Indeed H=v^2a and

\[
 v^m a^j=v^{-N}\prod_{\ell=0}^{j-1}(H-2\ell\eta),
                   \qquad N=2j-m>2j.                    \tag{12}
\]

The root character is sum_{ell=0}^{j-1}xi^{2ell}, which is nonzero:
the order of xi^2 is N/gcd(N,2)>j. All other terms must have
strictly greater canonical weight m'-2j'; a least coefficient-v
exponent alone is not sufficient.

The hypothesis matters. In ordinary coordinates [b,V]=eta,
P=b^j+V has bottom weight -j and character
1+t+...+t^{j-1}, which vanishes at a primitive jth root.
It has the actual mate -b, with commutator eta.

## The quartic face and all isolated global pole corrections

For (7), H=V^2A gives exactly

\[
 P=\rho V^{-14}C_4(H)+gV^{-4}H(H-2\eta)+\mu V,
 \quad C_4(H)=\prod_{i=0}^3(H-2i\eta).                    \tag{13}
\]

The bottom degree is four and its character is 1+t^2+t^4+t^6.
At a primitive fourteenth root this is nonzero: xi^2 is primitive
of order seven and 1+x+x^2+x^3 has degree less than Phi_7.
The local criterion proves (8) and no mate for this face.

More generally, every isolated global correction (4) with m<0,j>=2
added to S_b has a unique least canonical weight m-2j<-4.
The all-a term of (z-v)^j has that weight; replacing one a by
v increases weight by three, and quantum contractions preserve
canonical weight because [a,v]=eta/v. The quadratic seed has
least weight -4. Thus (12) excludes all such single-monomial
families directly, for every n and without a momentum-order cutoff.
In particular h*X^3*C^3 is also excluded. Global admission permits
their support; the finite-mate obstruction is an additional argument.

## Parameter-pole removal excludes the admitted family

The same conclusion holds for every finite integral rescaling
whose critical reduction satisfies the local divisor criterion,
even with higher generic momentum order. In particular it holds
when the critical reduction is exactly (7). A hypothetical finite mate may first be traced
to rational position coefficients. Its coefficients have finitely
many s-Gauss poles. At a highest pole s^{-N}, the residue commutes
with P, so (8) identifies it with a finite polynomial q(P).
Subtract s^{-N}q of the integral rescaled position. This retains
the commutator and strictly lowers the pole order. After finitely
many steps an integral mate would reduce to a nonzero constant
mate of P, just excluded. This is exactly QCL1's finite
parameter-pole argument, now with the quartic face (7).

Consequently S_epsilon in (6) has no finite rational-position
or algebraic-cover scalar mate. In particular it has no globally
admitted finite mate. The more general assertion applies to all
actual finite exact-R Long images whose critical normalization
and rescaling have precisely face (7).

For comparison, generic orders alone are weaker. A regular
prescribed mate of (6) would need generic order at least five:
order three is excluded by QSO34; at order four subtract its
constant-leading multiple c*S_epsilon. The remainder cannot
have order three, hence has order at most two. Since the
quartic coefficient of S_epsilon has h-valuation two and the
mate's quartic coefficient has valuation at least one, c has
valuation at least minus one. Every cubic coefficient of
S_epsilon has valuation at least two, so the subtraction cannot
leave a classical cubic mate. The quartic-face theorem excludes
all partner orders instead of only these first possibilities.

## The square ordering is also closed

A small change in ordering genuinely leaves the local criterion.
Set

\[
 W=X^3C^2,\qquad S_{\rm sq}=S_b+\epsilon h^2 W^2.         \tag{14}
\]

These are actual global finite operators commuting with R.
The classical reduction is still Long's prescribed second image,
and the generic order is four. Its critical face is

\[
        P_{\rm sq}=-3rA^2+2V+
                     \epsilon r^4(V^{-3}A^2)^2.           \tag{15}
\]

Writing C_2(H)=H(H-2eta), its lowest term is

\[
 \epsilon r^4 V^{-14}C_2(H-7\eta)C_2(H).                 \tag{16}
\]

The roots are 0,2eta,7eta,9eta, on one translation orbit.
Its character (1+t^2)(1+t^7) vanishes at every primitive
fourteenth root, so the local criterion does not apply.
This is the lowest coefficient of an actual operator square;
replacing it by the noncommutative fourth power in (6) changes
that divisor character. The distinction cannot be dropped.

The vanishing character only limits the local divisor criterion.
The [all-order square-face residue proof](SQUARE_ORDERED_QUARTIC_FACE.md)
now proves that (15) has rational polynomial centralizer and no finite
scalar mate over any finite algebraic coefficient cover.
[QOI1](ORDERED_QUARTIC_INFINITY_OBSTRUCTION.md) proves the same
conclusion independently and permits more scalar orderings.
Finite parameter-pole removal therefore closes S_sq at every finite
partner order. It is no longer a remaining candidate.

## A globally admitted mixed cubic-quartic face

Keep W=X^3C^2 and let lambda be a nonzero h-independent scalar.
The exact identities

\[
 \begin{split}
 [C,v]&=rh/v,\\
 vW&=X^2C^2,\\
 Wv&=X^2C^2+2rhX^4C-r^2h^2X^6
 \end{split}                                             \tag{17}
\]

show that both products vW and Wv are globally admitted.
Indeed [C^2,v]=rh(Cv^{-1}+v^{-1}C), while
Cv^{-1}=v^{-1}C-rh v^{-3}; multiplying by v^{-3}
gives the stated ordering correction. Define the fully polynomial
operator E=2X^2C^2+2RhX^4C-R^2h^2X^6=vW+Wv and put

\[
 S_\lambda=S_b-3\lambda h\{RE+CW+WC\}
                             -3R\lambda^2h^2W^2.         \tag{18}
\]

Every factor in (18) belongs to the actual global centralizer of R.
Consequently S_lambda is a finite globally admitted operator commuting
with R. It has Long's prescribed classical reduction and generic
momentum order four. On the generic fiber its leading coefficient
is -3lambda^2*h^2*r^5*v^{-6}, which is nonzero.

Since z=v+C/r and g=-3r, equations (3) and (17) give exactly

\[
 S_\lambda-\frac4{3r}
      =-3r\bigl(z+\lambda hW\bigr)^2+2v.                 \tag{19}
\]

This identity retains the ordered cross terms zW+Wz; replacing them
by twice one product would change the operator. Under (5),
W=s^{-4}r^2V^{-3}(A-sV)^2, so the integral critical reduction is

\[
 P_{\rm mix}=-3r\bigl(A+\lambda r^2V^{-3}A^2\bigr)^2+2V.
                                                                  \tag{20}
\]

With weights wt(V)=1 and wt(A)=3, Ore contractions lower weight
by five. The top commutative weighted symbol of (20) is

\[
             -3r\bigl(\xi+\lambda r^2V^{-3}\xi^2\bigr)^2.
\]

In particular its cubic coefficient has the nonzero top term
-6lambda*r^3*V^{-3}, of total weight six. QOI1's strict quartic
cone has top symbol s^2V^{-6}xi^4+g xi^2 and permits no cubic
term of weight six; its scalar ordering terms have lower weight.
Thus (20) lies outside that proved cone. Its lowest canonical
quartic coefficient still has the roots in (16), so the earlier
divisor character vanishes as well.

The admission and symbol calculations do not themselves exclude a mate.
A separate [twisted-adjoint theorem](TWISTED_ADJOINT_SQUARE_OBSTRUCTION.md)
now excludes every finite scalar mate for both (18) and (20), in every
partner order and on finite algebraic coefficient covers. The rational
quadratic operator in (19) has leading coefficient lambda*h*r^2*v^-3.
Its odd v-divisor forces every finite rational mate to have even order,
whereas its coefficient-fixed adjoint would turn one into a nonzero
skew-adjoint mate of odd order. This contradiction closes the mixed
family without treating the quadratic operator as a canonical coordinate.
No polynomial-centralizer conclusion for this mixed square is needed.

The globally admitted nonsymmetric perturbation
S_(lambda,nu)=S_lambda+nu*h*X^2*C, with nu an h-independent scalar,
is also closed. For nu nonzero its
critical face adds nu*r*V^-2*A to (20), breaking that adjoint.
The separate universal residue argument in the same theorem proves
polynomial centralizer and no finite mate for the critical face;
QCL1's finite parameter-pole descent then excludes every finite
generic mate of this exact global perturbation. These closures do
not classify arbitrary mixed negative-Laurent corrections.

## Scope and evidence

The global monomials (4) disprove a universal favorable Laurent
support restriction from admission alone. The divisor criterion
excludes isolated negative pole monomials of momentum order at
least two and their same-face parameter deformations. The
square-ordered quartic case is closed by its all-order residue
proof and QOI1. The coupled mixed face (20), although outside those
criteria, is closed by the independent adjoint theorem, and its specified
nonsymmetric first-order perturbation by a universal residue certificate.
Other mixed negative-Laurent faces, deformed first images and nonscalar
carriers remain outside these results.

[The exact controls](check_global_quartic_face.py) verify the global
chart identities and commutators, support realization, the isolated
quartic divisor test, polynomial-power identities, the positive
Airy-type control, the square ordering's vanishing character, and
the exact mixed-product, completed-square and critical-symbol
identities. Universal centralizer and mate exclusions use the
written proofs. The mixed admission controls establish the identities;
its mate exclusion uses the separately linked adjoint and residue proofs.
No actual Long mate, bounded lift, DC1/DC2 resolution, formal
verification, independent whole-proof replay or external review
is claimed.
