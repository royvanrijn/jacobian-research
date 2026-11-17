# Complete weight (1,-1,-3) degree-seven exclusion (KCUB18)

Over every characteristic-zero field, every polynomial Keller map in three
variables of source/output weights (1,-1,-3) and ordinary degree at most seven
is a polynomial automorphism. There is no support, collision, or coefficient
nonvanishing assumption. This closes a weighted class, not a global minimum.
It supplies no new counterexample.

The proof combines written polynomial degree/divisibility arguments with seven
exact rational identities. The independent checker reconstructs the complete
20-slot source determinant, every boundary equation and all six complete linear
parameterizations, then multiplies every identity. The total is 5,728 multiplier
monomials, maximum multiplier degree 13. The written case coverage and the
imported plane theorem are not proof-assistant formalizations.

## Complete model and first-coefficient split

The three distinct weights force a diagonal linear part; constant terms are
impossible. Its diagonal is invertible because the Jacobian is a nonzero
constant. Output scaling normalizes the linear part to the identity and the
Jacobian to one, preserving degree and support. Work over an algebraic closure;
the eventual unique polynomial inverse descends to the original field.

Set u=ab, w=a^3c. On a!=0 write F=(af,a^-1g,a^-3h), where the complete model is

```
f=1+c1*u+c2*u^2+c3*w+c4*u^3+c5*u*w,
g=u+c6*u^2+c7*w+c8*u^3+c9*u*w+c10*u^4+c11*u^2*w+c12*w^2,
h=w+c13*u^3+c14*u*w+c15*u^4+c16*u^2*w+c17*u^5
    +c18*w^2+c19*u^3*w+c20*u*w^2.
```

The quotient Jacobian identity is

```
g*(fu*hw-fw*hu)+3*h*(gu*fw-gw*fu)+f*(gu*hw-gw*hu)=1.       (1)
```

The independent replay compares (1) to the full original-source determinant,
and enumerates all ordinary-degree-two-through-seven atoms from the weights.
It imposes no support cap or explicit collision.

The cleared identity in `weighted-113/full/certificate.json` proves
c5^13 is in the ideal of all coefficients of (1)-1. Consequently c5=0.
Write f=1+mu*u+P*u^2+nu*w+R*u^3.
If nu!=0, independent diagonal conjugacies in b,c normalize nu=-1 and,
when mu!=0, mu=-1; mu=0 is its other stratum. They preserve the complete class.
If nu=0 use the univariate argument below, including mu=0.

## Polynomial graph boundaries when nu!=0

Now f=N(u)-w, N=1+mu*u+P*u^2+R*u^3. Put v=f and
G(u,v)=g(u,N(u)-v), H(u,v)=h(u,N(u)-v). The chain rule transforms (1) into

```
(G+v*Gv)*Hu-Gu*(3*H+v*Hv)=1.                              (2)
```

At v=0, G0 H0'-3G0' H0=1. For two nonzero polynomials of degrees m,n,
the leading coefficient of this Wronskian is proportional to n-3m.
If n!=3m its degree is m+n-1, so m+n=1. Otherwise n=3m.
Neither polynomial can be zero. For R!=0 the source bounds are m<=6,n<=7,
so m<=2,n<=6. For R=0 (including P=0) the bounds are m<=4,n<=5,
so m<=1,n<=3. In particular the cubic graph must retain the possible
balanced degrees (2,6); the stricter quadratic bound cannot be imported there.

For both mu=-1 and mu=0 the six branches are R!=0, R=0/P!=0, R=P=0.
The saved input for each branch appends all necessary higher boundary
coefficients to all native determinant equations. It then reduces their
complete affine-in-tail subsystem, of rank six, leaving nine free tails.
The replay performs independent exact row operations; only a monomial in the
explicitly nonzero branch parameter can be inverted. The denominators are
R^2, P, 1 and saturations are R^3, P^2, 1 respectively. Thus no denominator
root is lost: it belongs to a separately listed lower-degree branch.

Five branches have full-ring unit identities:

| mu | graph | multiplier monomials | max degree |
| --- | --- | ---: | ---: |
| -1 | R!=0 | 3032 | 9 |
| 0 | R!=0 | 579 | 11 |
| 0 | R=0, P!=0 | 187 | 6 |
| -1 | R=P=0 | 1578 | 4 |
| 0 | R=P=0 | 115 | 4 |

These exclude those branches. The remaining mu=-1 quadratic graph is excluded
next by a written degree split and a small exact unit identity. Its original
20-second reduced, parameter-field and factor-cleaned queries remain recorded
as timeouts; none is used as an exclusion.

## Remaining quadratic graph

Here N=1-u+P*u^2, P!=0. Source support and the necessary boundary equations give

```
G=G0(u)+v*G1(u)+K*v^2,      G0=A*u+B,
H=H0(u)+v*H1(u)+v^2*H2(u),
deg G1<=2, deg H0<=3, deg H1<=3, deg H2<=1,
```

with constant K. At source (u,w)=(0,0), equivalently (u,v)=(0,1),

```
G=H=0, Gu-Gv=1, -Hv=1, Hu-Hv=0,
[u^2] H(u,N(u))=0.                                        (3)
```

The last condition uses the missing pure-u^2 term in the original h.
The v^4,v^3,v^2,v^1,v^0 coefficients of (2) are respectively

```
3K H2'=0,
3K H1'+2G1 H2'-5G1' H2=0,
3K H0'+2G1 H1'+G0 H2'-4G1' H1-5A H2=0,
2G1 H0'+G0 H1'-3G1' H0-4A H1=0,
G0 H0'-3A H0=1.                                          (4)
```

The replay checks these identities and all degree bounds directly.

If K!=0, the first equation makes H2=C constant.
If C!=0, the second integrates to H1=lambda*G1+d,
lambda=5C/(3K)!=0. In the third equation a quadratic G1 would give a
nonzero cubic term, so deg G1<=1. That same equation gives deg H0<=2.
The last equation forces deg H0<=1: if A=0 it is immediate; if A!=0
a polynomial degree two term cannot balance degree three. The degree-one
coefficient in the third equation now forces G1 constant, hence H1 constant.
Conditions (3) give H0'=H1+2C=-1. The last equation then forces A=0,
and the third gives -3K=0, a contradiction.

If C=0, the second equation gives H1=d constant, and (3) gives d=-1.
The third gives 3K H0'=4d G1', so deg H0<=2, and the boundary equation
again gives deg H0<=1. Conditions (3) give H0'=d=-1, H0(0)=1;
the boundary equation gives A=0. The fourth equation, with H0=1-u and
G1 affine, forces G1'=0. The third then gives K=0, a contradiction.
Thus K=0.

If H2!=0, the v^3 equation is 2G1 H2'-5G1' H2=0. Unless G1=0,
leading degrees give 2 deg H2=5 deg G1; the bounds force both constant.
Thus G is affine, also when G1=0. Conditions (3) put it in the form

```
G=(1+k)*u-k+k*v.
```

Here is a degree-free exclusion of every such affine G, for arbitrary H.
If 1+k=0, equation (2) becomes (1-2v)Hu=1, impossible for polynomials.
Otherwise set q=G and A=1+k. The vector field on the left of (2) sends
q to A*q and v to -A*v. Thus

```
A*(q*Hq-v*Hv-3H)=1.
```

All polynomial solutions are -1/(3A)+sum_j Cj*q^(j+3)*v^j.
At the source q=0,v=1 their value is -1/(3A), contradicting H=0.
Consequently H2=0 as well.

With K=H2=0, the v^2 equation becomes G1 H1'-2G1' H1=0.
If G1=0 apply the affine exclusion. Otherwise H1=c*G1^2.
Condition -Hv=1 at the source gives beta=G1(0)!=0 and c=-1/beta^2.
Since deg H1<=3, G1=alpha*u+beta is affine. From (3) and the constant
coefficient of the boundary equation obtain

```
A=-3-beta, alpha=2beta+4,
G0=A*u-beta, G1=alpha*u+beta, H1=-(G1/beta)^2,
H0=1+(2alpha/beta-1)*u+(alpha^2/beta^2-2alpha/beta+P)*u^2+T*u^3.
```

Only P,T,beta remain. Clear the beta denominators of the u coefficients of
E0=G0 H0'-3A H0-1 and
E1=2G1 H0'-3alpha H0+G0 H1'-4A H1. Up to rational scalar factors they are

```
P*beta^3+3P*beta^2-3T*beta^3+8beta^2+40beta+48,
P*beta^2-3beta^2-9beta-8,
T*(beta+2),
P*beta^3+2P*beta^2+3T*beta^3-4beta^3-20beta^2-32beta-16,
2P*beta^2-5beta^2-12beta-4.
```

`weighted-113/core/certificate.json` expresses 1 in their ideal, with eight
multiplier monomials and maximum degree two. The replay reconstructs these
five polynomials from E0,E1 and verifies both the source conditions and this
identity. Its selective substitution changes only u,v and retains all Laurent
parameter powers, including beta denominators. The earlier source-condition
audit incorrectly passed Laurent coefficients through a polynomial-only
substitution helper; that use was corrected and the entire checker replayed.
The preceding receipt is retained as `validation-before-laurent-fix.json`. This excludes the final quadratic branch.

## Univariate first quotient when nu=0

Let f=phi(u)=1+mu*u+P*u^2+R*u^3 be nonconstant.
Equation (1) is

```
phi'*(g*hw-3*gw*h)+phi*(gu*hw-gw*hu)=1.                    (5)
```

A multiple root of phi would make its left side zero, so phi is squarefree.
At each root the nonzero constant Wronskian in w has source degrees at most
two in both g,h. Leading degrees force both at most one. In the original
source g's w^2 coefficient is constant and h's is affine in u.

If deg phi>=2, its at least two distinct roots force both w^2 coefficients
to vanish globally. Write g=A(u)w+B(u), h=D(u)w+E(u), with deg A<=2,
D(0)=1 and B(0)=E(0)=0. If A=0, (5) reads (phi*B)'*hw=1.
Both factors are constant, and integration gives phi*B=k*u with k!=0,
impossible because nonconstant phi with phi(0)=1 cannot divide u.
For A!=0, the w coefficient in (5) is

```
-2phi'*A*D+phi*(A'*D-A*D')=0,
```

so phi^2*D=k*A for a nonzero constant k. Its left degree is at least four,
but deg A<=2, impossible.

It remains to consider linear phi. At its single root the same boundary
argument gives g=A*w+B, h=C*w^2+D*w+E, with C=delta*phi, deg A<=2,
D(0)=1, B(0)=E(0)=0. Again A=0 is excluded by (phi B)'*hw=1.
If C!=0, the w^2 equation is

```
-phi'*A*C+phi*(2A'*C-A*C')=0,
```

hence phi*C=c*A^2. Since C=delta*phi, A=alpha*phi with alpha!=0.
The w equation is then

```
phi*(2delta*phi*B-alpha*phi*D)'=0.
```

The primitive is constant and vanishes at a root of phi, so
D=2delta*B/alpha. This contradicts D(0)=1, B(0)=0.
Thus C=0. The earlier relation phi^2*D=k*A and deg A<=2 force D=1
and A=phi^2/k. The constant equation in (5) now integrates to

```
phi*B-phi^3*E/k=u,
```

using B(0)=E(0)=0. Again phi cannot divide u. This excludes every
nonconstant univariate first quotient. The replay checks the universal
coefficient identities; squarefreeness, integration, and degree arguments are
written mathematics.

## Constant first quotient and polynomial inverse

If f=1, F1=a and the tail is a Keller pair over K(a) with both ordinary
degrees in (b,c) at most five. We use Nakai--Baba, *A generalization of
Magnus' theorem*, Osaka J. Math. 14 (1977), 403--409, Theorem 3(1)--(2),
printed page 408: a plane Keller pair is invertible if one degree is prime
or one degree is four. Degree one is elementary; this covers all degrees
one through five. This exact theorem was visually checked in the primary PDF
for KCUB15 and is sealed in `weighted-113/plane-degree5-import.json` together
with its original import and PDF hashes. It is not reproved here.

The tail fixes (b,c)=(0,0); over K[a] its linear matrix is
[[1,c7*a^2],[0,1]], which is unimodular. Therefore its unique formal inverse
has coefficients in K[a]. Over K(a) this formal inverse is the polynomial
inverse just obtained, so only finitely many coefficients are nonzero and
they are already in K[a]. The inverse is consequently polynomial globally,
including a=0. Uniqueness descends it from the algebraic closure. Thus every
remaining map is an automorphism, completing the theorem.

## Reproduction and preserved limits

Run from the repository root:

```
.venv/bin/python research/jacobian/cubic-reduction/verify_weighted_113.py
```

This uses independent sparse Fraction arithmetic, not SymPy or Singular.
The exact generators, Singular lifts, multipliers and source hashes are saved
under `weighted-113/`. Scripts `weighted_113.py`, `weighted_113_boundary.py`,
`weighted_113_lift.py`, `weighted_113_unit.py`, `weighted_113_quad.py`, and
`weighted_113_core.py` preserve the changed experiment and certificate exports.
The original broad coefficient basis timeout20 and real sparse/nonlinear-first
UNKNOWN30 in `weighted-sparse7/1_-1_-3/` remain unchanged. Changed full solves
and the hard quadratic field/clean queries also retain their timeout receipts.
They were superseded mathematically by the structural proof, not relabelled
as successful solver runs.

Together KCUB16--18 exclude all three additional initial weight classes
(1,-2,-3), (2,-1,-3), (1,-1,-3) through degree seven. KCUB15 classifies
(1,-1,-2) through that degree, where the established raw13 map remains.
Other weights, higher degrees, polynomial coordinates, global raw/normalized
sparsity, ordinary cubic dimension below12 and cubic-linear dimension below89
remain open. The next experiment must change one of these boundaries; it must
not repeat the unchanged time-limited solves in a now closed class.
