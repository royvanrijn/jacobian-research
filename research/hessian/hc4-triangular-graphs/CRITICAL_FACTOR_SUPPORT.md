# Where critical-factor units can vanish

15 November 2026. The localized-unit rank budget can be made sensitive
to each critical factor's actual permitted leading levels. An empty
support excludes a rational root independently of special-fibre
irreducibility. The resulting cubic coefficient test closes further
S3 outputs with arbitrary leading multiplicities and parameter degrees.

## 1. Critical-factor support and all subset rank inequalities (HC4MYID0)

Let K=(f,g) be complex plane Keller, B=C[s,t], L=C(s,t), and
P(T,z)=sum_(i=0)^m b_i(T)z^i, m>=2, with a=b_m!=0. Suppose
P(f,h)+g=0 for a rational h in L. The algebraic part of
[HC4MYIC1](LEADING_VALUATIONS_AND_CRITICAL_CURVES.md#2-rational-labels-force-affine-line-critical-curves-hc4myic1)
gives the polynomial u=a(f)*h, the monic equation

\[
 Q(T,U,V)=a(T)^{m-1}\big(P(T,U/a(T))+V\big),\qquad
 J(f,u)Q_U(f,u)=-a(f)^{m-1}J(f,g).                     \tag{1}
\]

Here Q_U is independent of V, with constant leading U coefficient m.
Factor it over C(T), keeping distinct monic irreducible factors:

\[
                   Q_U=m\prod_{i=1}^{ell} H_i(T,U)^{n_i}.
                                                                    \tag{2}
\]

Each H_i belongs to C[T,U]. One can use Gauss's lemma: a
factorization over C[T,U] has leading coefficients whose product is
the unit m, so each can be normalized to one; its factors remain
irreducible over C(T). This counts irreducible critical components
over C(T), not all roots after a parameter extension.

Define the allowed finite support of the i-th factor by

\[
                   S_i=\{c\in C:a(c)=0, H_i(c,0)=0\}.   \tag{3}
\]

Then the polynomial v_i=H_i(f,u) has every irreducible zero divisor
over one of those levels, and no others. For any nonempty subset I
of the distinct factors, the full f-fibre component counts r_c satisfy

\[
                 |I|\ \leq\sum_{c\in\bigcup_{i\in I}S_i}(r_c-1).
                                                                    \tag{4}
\]

All degrees and factor multiplicities are unrestricted. In particular

\[
       \gcd(a(T),H_i(T,0))=1\text{ for some }i
             \quad\Longrightarrow\quad
                  P(f,Z)+g\text{ has no rational root in }L.         \tag{5}
\]

No irreducibility hypothesis on any special Keller fibre is needed
for (5). Support in (3) is a necessary allowance; the theorem does
not assert that every permitted level or every fibre component is used.

### Specialization determines the support

By (1), every v_i is a unit of B[1/a(f)]. Hence every source prime
divisor of v_i lies over some leading zero c. If D is such a
divisor, reducing modulo D puts both Q(c,u,0) and Q_U(c,u) equal
to zero. At a leading zero their equations simplify to

\[
 Q(c,U,0)=U^{m-1}(U+b_{m-1}(c)),\qquad
 Q_U(c,U)=\frac{d}{dU}Q(c,U,0).                         \tag{6}
\]

The only possible common root is U=0. If b_(m-1)(c)!=0, the other
root -b_(m-1)(c) of Q is simple; if that coefficient is zero, every
root is already zero. Thus u=0 modulo D, and H_i(c,0)=0. This
proves (3). Equivalently one can use
U Q_U-(m-1)Q=U^m modulo a, which makes the common-root conclusion
immediate. This step retains the full normalized root equation, not
only its derivative.

If S_i is empty, v_i has no prime divisor in B and is a nonzero
constant. But J(f,u)!=0 in (1), so f,u are algebraically independent;
the positive-U-degree polynomial H_i(f,u) cannot be constant.
This contradiction proves (5).

### Independence gives the subset inequalities

The classes of the v_i are independent even modulo C(f)^*. Indeed
u is transcendental over C(f). A relation
product_i H_i(f,u)^(e_i) in C(f)^* would be a rational-function
identity in that variable, and its order at each distinct irreducible
factor H_i forces e_i=0.

For any subset I, all these units lie in the localization over
S_I=union_(i in I) S_i. The exact unit quotient of
[HC4MYHZ0](LEADING_FIBRE_UNIT_GATE.md#1-the-exact-localized-unit-defect-hc4myhz0)
is the direct sum over c in S_I of Z^(r_c)/Z(1,...,1), of rank
sum_(c in S_I)(r_c-1). Its embedding into the larger localization
does not change these blocks. Independence of the |I| classes
therefore proves (4). The total rank bound alone can be weaker:
excess defect on a level excluded from a subset's support cannot
help that subset. No matching theorem or S-unit height theorem is
needed for this direct lattice argument.

### Two all-degree no-root consequences

For every m>=2, the polynomial

\[
                  P=a(T)z^m+b(T)z^{m-1}+d(T)
\]

has no rational source root whenever gcd(a,b)=1. Its normalized
derivative is U^(m-2)*(mU+(m-1)b); the factor
H=U+(m-1)b/m has constant term coprime to a. This is a no-root
statement, not an irreducibility claim in arbitrary higher degree.

More generally suppose a has only one distinct zero c, b_(m-1)(c)
is nonzero, and Q_U splits into linear factors over C(T). The
specialization in (6) has one nonzero simple critical root
-(m-1)b_(m-1)(c)/m. One polynomial critical branch must take that
value at c, so its support is empty. Again (5) excludes every
rational source root, with arbitrary leading multiplicity.

## 2. Cubic closure by a critical-factor support gcd (HC4MYID1)

Let P=a z^3+b z^2+c z+d have coefficients in C[T], a!=0, and
suppose its nonzero critical discriminant is a polynomial square:

\[
              b^2-3ac=H^2,\qquad H\ne0.                 \tag{7}
\]

If either

\[
                    \gcd(a,b-H)=1\quad\text{or}\quad
                    \gcd(a,b+H)=1,                      \tag{8}
\]

then every closed stable embedding A^(r+2)->A^(r+3), r>=1,
carrying a plane Keller pair through (T,xy-P(T,z)) forces that
pair invertible. All parameter, leading-coefficient and embedding
degrees are unrestricted. No special K1-fibre or nilpotent root-group
hypothesis is used.

In fact the normalized derivative has distinct polynomial roots
beta_+=(-b+H)/3 and beta_-=(-b-H)/3. Condition (8) makes one of
their supports in (3) empty. A reducible source cubic would have a
rational root, contrary to (5). The retained
[root-transitivity and global filtration proof](LEADING_VALUATIONS_AND_CRITICAL_CURVES.md)
then gives the polynomial inverse. This does not infer arbitrary
splitting-group preservation from transitivity.

This closes two useful unrestricted families:

- If a=alpha*(T-c0)^q, alpha!=0, q>=1, b(c0)!=0 and (7)
  holds, one of b+H or b-H is nonzero at c0, so (8) holds.
  The leading multiplicity q has no bound.
- Every sparse cubic a(T)z^3+b(T)z^2+d(T) with gcd(a,b)=1
  closes, with no restriction on the number or multiplicities of
  leading zeros. When b!=0, take H=b and use gcd(a,2b)=1.
  If b=0, coprimality makes a a nonzero constant and the retained
  monic output theorem already closes the case.

For example, for arbitrary q>=1 and arbitrary R,d in C[T],

\[
 P=T^q z^3+z^2-\frac{2R+T^q R^2}{3}\,z+d,
                 H=1+T^q R                              \tag{9}
\]

satisfies (7)-(8). The branch beta_-=-(2+T^q R)/3 has no
permitted leading zero, even when a has high multiplicity. These
are genuinely S3 root covers: P+V is irreducible, and the
V-discriminant of its z-discriminant is 16*H^6!=0, so the
z-discriminant has simple zeros and is not a square. This family
does not meet the preceding critical-genus gate, since its critical
discriminant is a square; its quadratic coefficient is not divisible
by a, and its leading zero drops the degree only to two.

The usual genuine closed embedding j=(s,1,t+P(s,w),w), with an
unused v, still has its literal inverse. The independent quadratic
T*v+y*z+x^2 pulls back to 1+s*v+t*w+w*P(s,w), with actual unit
Hessian and the retained explicit gradient inverse. These controls
have invertible gradients and do not construct nonautomorphic K.

### The remaining square-discriminant branch needs both supports

When (8) fails, define

\[
 S_+=\{a=0,b-H=0\},\qquad S_-=\{a=0,b+H=0\}.
\]

A hypothetical rational source label must have defect at least one
over each S_+ and S_-, and at least two over their union, by (4).
If gcd(a,b)=1 these two sets are disjoint. Thus two distinct
reducible full f-fibres are necessary, one in each support; a single
fibre with many components cannot replace them.

This is a real residual condition rather than an empty one at the
coefficient level. For a=T*(T-1), b=1, c=-4/3 and H=1-2T,
the critical supports are {0} and {1}, both nonempty. The new gcd
test does not close that coefficient family. No Keller realization
or actual reducible source fibre is asserted for this control.

## 3. Two reduced singleton supports force a constant critical gap (HC4MYID2)

There is a further global restriction beyond the rank test. Suppose
two distinct rational critical branches beta_i,beta_j have source
units v_i=u-beta_i(f) and v_j=u-beta_j(f) which are reduced
polynomials. Suppose each has just one allowed leading level, c_i
and c_j respectively. Then

\[
                         beta_i(T)-beta_j(T)\in C^*.     \tag{10}
\]

The coefficient and source degrees are unrestricted. This condition
is on actual reduced polynomial units, not merely their zero sets.

Indeed v_i divides f-c_i, because all its distinct prime divisors
divide that reduced fibre. At a source zero of v_i, the proof of
(3) gives f=c_i and u=0. Consequently v_j is the same constant
kappa_ij=-beta_j(c_i) on the whole zero set of v_i. Since v_i
is reduced, v_i divides v_j-kappa_ij. Interchanging i,j gives
v_j dividing v_i-kappa_ji. The v_i and v_j are nonconstant;
ordinary total degrees therefore force them to have equal degree
and both quotients to be scalars. Write v_j=lambda*v_i+kappa_ij.
Algebraic independence of f,u forces lambda=1, since otherwise
(1-lambda)u would be a rational function of f. It follows that
beta_i-beta_j is a constant, nonzero because the two branches
are distinct. This proves (10), without an embedding-rectification
or coordinate theorem.

### A cubic two-level corollary

For a cubic with D=H^2!=0, assume a is square-free with at most
two distinct zeros and gcd(a,b)=1. If H is nonconstant, every
carried output embedding forces K invertible. Constant a is already
monic; one distinct leading zero closes by Section 2. With two
zeros, empty support again closes by Section 2. Otherwise the two
critical supports are disjoint singletons.

The required reducedness of both critical units follows from the
Keller identity, rather than being added as a hypothesis here. Let
D_source divide one critical unit over its allowed simple leading
level. Then u=0 on D_source, so the tangent derivation J(f,.) kills
that restriction, making v_D(J(f,u))>=1. The other critical unit
is nonzero at that level: b(c)!=0 and H(c)^2=b(c)^2 give distinct
critical values zero and -2b(c)/3. In the identity

\[
                 3J(f,u)v_+v_-=-a(f)^2 J(f,g),
\]

the right side has valuation two, since f-c is reduced and a has
a simple zero. The selected critical unit already has valuation
at least one, so its valuation is exactly one. This holds at each
of its prime divisors; both units are reduced. Apply (10): their
branch difference is 2H/3, contradicting nonconstant H.

In particular the coefficient example above,
P=T*(T-1)*z^3+z^2-4*z/3 with H=1-2T, is closed by this
additional polynomial-pencil argument, although it passed the support
gcd test. This closes it without asserting any actual reducible
Keller fibres.

The constant-gap boundary is still nonempty at the coefficient level.
For a=T*(T-1), b=1-2T and c=4/3, the critical discriminant is
H^2=1. Both supports are again singletons, but the critical gap is
constant, so (10) alone yields no contradiction. Section 4 supplies
the missing exclusion for this fixed model. Repeated leading
zeros also require a separate reducedness argument; the valuation-two
step above is not silently reused for higher multiplicity.

### The remaining two-level constant-gap case has one fixed output model

There is no need to search arbitrarily high coefficient degrees in that
particular residual. Suppose a is square-free with two zeros c0,c1,
gcd(a,b)=1 and D=H0^2 with H0 a nonzero constant. An empty
support is already closed. Otherwise order the zeros and choose the
sign of H0 so that b(c0)=H0 and b(c1)=-H0. Put

\[
 T=c0+Delta*tau,\quad Delta=c1-c0,\qquad
 a(c0+Delta*tau)=alpha*tau*(tau-1),\quad alpha\ne0.
\]

The transformed quadratic coefficient has the unique form
b=H0*(1-2tau)+a*B(tau), with B polynomial. Set rho=B/3.
The polynomial shift z=Z-rho removes that entire extra term.
The new linear coefficient is forced by D=H0^2 to equal
4H0^2/(3alpha), and the new constant coefficient is a polynomial
d_tilde(tau). Rescaling Z=(H0/alpha)*zeta gives

\[
 P(c0+Delta*tau,(H0/alpha)zeta-rho(tau))
      =kappa P_*(tau,zeta)+d_{tilde}(tau),\qquad
 kappa=H0^3/alpha^2,                                   \tag{11}
\]

where the fixed polynomial is

\[
 P_*(tau,zeta)=tau(tau-1)zeta^3+(1-2tau)zeta^2+4zeta/3.
                                                                    \tag{12}
\]

For clarity the actual conjugating maps on output data are

\[
 Psi(tau,X,Y,zeta)=(c0+Delta*tau,\ kappa X,\ Y,
                          (H0/alpha)zeta-rho(tau)),
\]
\[
 Theta(tau,V)=(c0+Delta*tau,\ kappa V-d_{tilde}(tau)).
\]

Both are polynomial automorphisms with written inverses obtained by
solving these affine triangular formulas. Their target Jacobian is
Delta*kappa!=0, and

\[
       (T,xy-P)\circ Psi=Theta\circ(tau,XY-P_*).          \tag{13}
\]

Thus a closed embedding j carrying K becomes Psi^(-1) composed
with j, carrying Theta^(-1) composed with K through the fixed pair
(tau,XY-P_*). Closedness and the Keller condition are preserved,
and K is invertible exactly when the transformed pair is. Extra
ambient variables are left fixed. Conversely the same maps recover
the original output data. This is an exact all-degree equivalence
for this constant-gap two-simple-leading-level branch.

This normalizes output pairs and their embedding data. It does not
assert invariance of a fixed quadratic potential or of the actual
source Hessian under a nonlinear change; fixed-potential compatibility
must still be tracked literally. The next argument closes this fixed
model at the level of the output obstruction.

## 4. The fixed cubic is closed by a second polynomial chart (HC4MYID3)

For the fixed P_* in (12), P_*(f,Z)+g is irreducible over L for
every complex plane Keller pair (f,g). Consequently every closed
stable embedding carrying that pair through (T,xy-P_*) forces it
invertible. Combined with Section 3, this closes the entire cubic
class with square-free a having at most two zeros, gcd(a,b)=1
and nonzero square critical discriminant, in arbitrary coefficient
and embedding degrees. The constant-gap normal form is a reduction
to this proof, rather than a remaining endpoint.

Suppose instead that h in L were a rational root. Rename the original
source coordinates s1,s2 in this section. Set u=f(f-1)h and
t=3u+1-2f. The normalized critical units are

\[
             v_+=u-2f/3=(t-1)/3,\qquad
             v_-=u-2(f-1)/3=(t+1)/3.
\]

Their supports are respectively {0} and {1}. The simple-level
valuation argument of Section 3 makes both reduced. Thus v_+
divides f and v_- divides f-1. Define the polynomials

\[
                  p=f/v_+,\quad q=(f-1)/v_-,\quad A=p-q.
\]

Subtracting f-1 from f gives t(p-q)-(p+q)=3. Solving this
linear identity yields

\[
 p=\frac{A(t+1)-3}{2},\qquad q=\frac{A(t-1)-3}{2},
\]
\[
 f=\frac{A(t^2-1)-3t+3}{6},\qquad
 u=\frac{A(t^2-1)}{9},\qquad h=\frac{A}{pq}.            \tag{14}
\]

All formulas are identities in the source field; p,q are nonzero
polynomials. Since f,u are algebraically independent by (1), the
polynomial map Phi=(A,t):A2->A2 is dominant and generically finite.
Substitution of (14) in the rational root equation, followed by
clearing its denominators, gives the exact polynomial identity

\[
 9p^2q^2 g+N(A,t)=0,\qquad
 N(A,t)=A\big(A^2(t^2-1)-9At+27\big).                  \tag{15}
\]

Consider the irreducible target hyperbola

\[
                          Gamma:\ A(t+1)=3.
\]

At a source point mapping to Gamma, the displayed formula for p
makes p=0, so (15) makes N=0. On Gamma, A is nonzero and

\[
                         N=3A(A+3).
\]

Hence the only possible image point on Gamma is (A,t)=(-3,-2).
The rest of this hyperbola is omitted by the generically finite
polynomial map Phi.

The retained explicit
[Jelonek--Lason Theorem 3.2](https://arxiv.org/html/1411.5011v2)
now rules this out, exactly as in HC4MYIC1. Clear algebraic-equation
denominators to put the nonproperness set inside one target algebraic
curve. A generic omitted point of Gamma is nonproper and belongs
to no other component of that curve, so the theorem supplies a
nonconstant polynomial map A1->Gamma. Such a map cannot exist:
its polynomial coordinate functions would satisfy A(s)*(t(s)+1)=3,
making both factors units of C[s] and therefore constants. This
contradiction excludes the rational root.

A reducible cubic over L has a rational root, so irreducibility is
proved. The retained root-transitivity/global-filtration argument
then closes the fixed output. The normal-form equivalence (13)
closes every constant-gap two-simple-leading-level case, while
Section 3 already closes the nonconstant square-gap cases. No
properness of the original Keller map or of the auxiliary Phi was
assumed, and no Hessian invariance was used.

The nonsquare critical-discriminant case of square-free degree one
also remains separate: its quadratic critical factor must not be
counted as two rational branches over C(T). General higher-degree
intransitive factors, other output geometries, both plane endpoints
and unrestricted HC4/JC2 remain open.

The [source context](sources/critical-support-context.json) records
the literature lookup and retained inputs. Section 4 reuses the
explicit qualitative nonproperness theorem retained with HC4MYIC1;
no new external theorem or function-field S-unit height bound is
added. The [checker](verify_critical_support.py)
and [receipt](receipt-critical-support.json) give exact finite controls
separate from the unrestricted proofs. No formal verification,
independent replay, external review or novelty priority is asserted.
