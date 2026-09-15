# Q80 genus-one k0 cannot have mixed abscissa valuations

On the literal direct11952 Q80 MW17 parent, a polynomial genus-one
common-quartic pair cannot have one integral and one negative abscissa
valuation at131. This includes arbitrary scalars, higher contacts, infinity
and unpointed genus-one bases. No denominator bound is assumed.

Together with the [quartic gate](Q80_GENUS_ONE_K0_QUARTIC_GATE_2026-09-14.md)
and [pole-degree theorem](Q80_GENUS_ONE_K0_POLE_BOUNDARY_2026-09-14.md), this
leaves only `k=0,j=2` with **both original abscissas negative**. That case
remains UNKNOWN. No infinite rank-at-least19 source is constructed.

## 1. The integral section has exactly one possible reduction

Retain the full chart

```
f(t,x_i)=x_i^3+A*x_i+B=eta*D_0*r_i^2,
deg x_i,deg r_i<=4,
D_0 squarefree binary quartic, gcd(D_0,Delta)=1,
gcd(D_0,x_1-x_2)=1.                                              (1)
```

The retained denominator proof puts a mixed pair at
`D_0,bar=q_0^2`, where `q_0=t^2+62*t+88`, after making D_0 monic.
The negative section has even abscissa valuation and eta is a unit by
the local contact proof in the pole-degree theorem. Denote the integral
section by P and the negative one by R.

The reduction of P is an ordinary section of the constant twist, with
ordinate proportional to q_0*r_bar. It cannot pass through the I1 node,
since differentiating a section through the node would contradict the
unit base derivative. Hence x_bar equals the simple root
`S=73+6*t mod q_0`. Divisibility of f(x_bar) by q_0^2 and the unit simple-root
derivative give a unique residue lift modulo q_0^2:

```
x_0=9+29*t+8*t^2+120*t^3,
x_bar=x_0+c*q_0^2,           c in F_131.                           (2)
```

All131 choices are checked, including every constant scalar squareclass.
Exactly one makes `f(x_bar)/q_0^2` a scalar times a square:

```
c=82,
X_T,bar=60+83*t+52*t^2+70*t^3+82*t^4,
Y_T,bar=59+50*t+77*t^2+41*t^3+66*t^4+33*t^5+3*t^6,
f(X_T,bar)/q_0^2=9*(95+98*t+83*t^2+80*t^3+t^4)^2.                 (3)
```

In particular eta has square residue and therefore is a square in Q_131.
This conclusion does not assume a rational point on the covering curve.
Absorb its square root into the r_i and now use `eta=1`. Choose the sign
of P so that `q_0*r_P,bar=Y_T,bar`.

The pair in (3) is point732 in the complete retained norm-four roster.
It is the reduction of the actual inherited section T with word

```
(0,1,1,-1,0,0,-1,1,0,-1,0,0,0,0,-1,0,0).                         (4)
```

The literal parent height Gram gives this word height4. Since all parent
fibres are irreducible, T misses O and has polynomial coordinates X_T,Y_T
of degrees at most4 and6. Their coefficients are integral at131: the
word reduces to the finite section (3) on the good coefficient Gauss
fibre. The retained word identities and a fresh exact Sage calculation
over F_131(t) identify its reduction, including the ordinate sign.

## 2. The nearby integral solutions form one formal parameter

The polynomial Y_T has a unique monic quadratic factor H_T over Z_131
reducing to q_0. Indeed q_0 occurs simply in Y_T,bar and is coprime to
the remaining degree-four factor. Thus the exact split solution

```
x=X_T,       D=H_T^2,       r=Y_T/H_T                             (5)
```

has the same reduction as P. This is only the centre of a formal
deformation calculation; it is not an admissible squarefree cover.

Fix the leading coefficient of D to be1. There are14 variable coefficients:
five for x, five for r and four for D. The13 coefficient equations are
`f(x)-D*r^2=0`. At the reduced solution (3), write

```
q=q_0,  D=q^2,
r_0=Y_T,bar/q=23+32*t+118*t^2+109*t^3+3*t^4,
F_x=3*X_T,bar^2+A_bar.
```

Their linearization is

```
F_x*delta x - 2*q^2*r_0*delta r - r_0^2*delta D.                  (6)
```

The13-by14 matrix has rank13. Removing the column for the coefficient
x_4 gives a square minor of determinant29 in F_131. The unique tangent
normalized by `delta x_4=3` is

```
delta x=r_0,
delta r=42+23*t+51*t^2+106*t^3+20*t^4,
delta D=21+78*t+110*t^2+29*t^3.                                  (7)
```

These are exact polynomial identities over F_131. At q_0,

```
r_0=38+81*t,     F_x=48+93*t,     delta D=11+2*t,                 (8)
```

all nonzero. The complete census and these ranks, tangents and units
have independent integer-polynomial replay.

The unit minor gives a unique integral formal solution through (5) after
choosing

```
kappa=(x_4-X_T,4)/3.
```

One can construct its coefficients successively by solving with that
invertible13-by13 matrix over Z_131. The resulting power series converge
for kappa in p*Z_131. Multivariate Hensel uniqueness shows that every
integral solution with the specified reduction lies on this formal
branch. Only this necessary local parametrization is used; neither a
rational parametrization nor a positive global family is inferred.

Let K/Q_131 be the unramified quadratic extension, and work in either
q_0 residue disc. Let tau be the root of H_T there. Hensel factorization
of D(kappa) into its two residue clusters gives a monic quadratic local
factor over O_K[[kappa]]. Its centre c(kappa) and discriminant satisfy

```
c(kappa)-tau in kappa*O_K[[kappa]],
disc_local(D(kappa)) = kappa*u(kappa),       u integral unit.      (9)
```

For the second assertion, at kappa=0 the local factor is `(t-tau)^2`.
Its discriminant derivative modulo p is
`-4*delta D(theta)/q_0'(theta)^2`, a unit by (8). Here theta is a root
of q_0 in the residue field; q_0' accounts for the other cluster's
unit factor. Thus (9) retains every higher-order deformation term.
For the actual distinct branch roots a,b in this disc, kappa is nonzero
and

```
v(a-b)=v(kappa)/2,       v((a+b)/2-tau)>=v(kappa).                 (10)
```

The roots a,b need not lie in K.

## 3. The negative section forces a particular deformation valuation

The local notation is the same as in the
[nodal orientation proof](Q80_GENUS_ONE_K2_NODAL_ORIENTATION_2026-09-14.md):

```
f=(x-rho)*(U^2-V),      U=x+rho/2,
V(alpha)=0,            V'(alpha) a unit,                          (11)
```

where alpha is the actual nodal base value in K and rho is the simple
root function. The negative section R chooses node roots at both local
branches by the pole-degree theorem. Preparation gives

```
x_R-rho=p^(-2*ell)*h(t)*(t-gamma)^2,       h integral unit,
gamma,beta in K,      gamma!=beta,       ell>=1,                  (12)
```

where gamma is the double simple-root contact and beta the other ordinate
zero. The valuation calculations in the earlier local proof give

```
v(beta-alpha)=2*ell,       v(gamma-alpha)=ell,
v(a-alpha)=v(b-alpha)=ell,
v(a-b)=3*ell/2,
v((a+b)/2-alpha-2*(gamma-alpha))>ell.                             (13)
```

Here only gamma,beta, not a,b, must lie in K. To make the use of (13)
with a ramified branch field explicit, differentiate (12) at beta.
Then `v(U'(beta))=-ell`, `v(U(beta))=ell`, and the term `2*(U')^2`
dominates H'' for `H=U^2-V`, so beta is exactly a double zero. Scaling
`t=gamma+p^ell*z` puts beta in one root cluster and a,b in the opposite
one. This proves their valuations, their separation greater than ell,
and the last congruence in (13).

Over a finite extension containing a,b, label the two near-node roots on
their smaller common disc by choosing a square root of V at a. The
polynomial abscissa difference has valuation `v(a-b)-ell`; a fixed node-root
label has difference valuation `v(a-b)-ell/2`. They cannot agree, so R
switches labels. The opposite-label difference has valuation ell/2,
forcing `v(a-b)=3*ell/2`. This repeats only the valuation part of the
old proof; its k2 agreement contradiction and K-rational branch assumption
are not imported.

Comparing (10) and (13) gives `v(kappa)=3*ell`. Since the branch centre
has distance of valuation ell from alpha, while it differs from tau
with valuation at least3*ell,

```
ell=v(tau-alpha)=v(Y_T(alpha)).                                  (14)
```

The last equality uses the simple H_T contact and the unit r_0 in (8).
This determines the denominator valuation from the literal inherited
section; it is not a chosen cutoff.

## 4. The exact contact gives an impossible leading scalar

The field is represented by

```
K residue field = F_131[theta]/(theta^2+62*theta+88).
```

We compute the actual section (4) at the lifted discriminant root alpha,
not at an arbitrary lift of theta. At the nodal fibre write e for its
double cubic root. Then `e=-3*B(alpha)/(2*A(alpha))`, and its smooth
group has coordinates

```
u=Y/(X-e),       X=u^2-2*e,       Y=u*(u^2-3*e).
```

The producer evaluates the word by multiplying
`(u+sqrt(3*e))/(u-sqrt(3*e))` with certified p-adic precision. The
independent checker works exactly in
`(Z/131^2)[theta]/(theta^2+62*theta+88)`: it Hensel-lifts alpha once,
checks the literal basis points there, and uses projective addition

```
(a:b)+(c:d)=(a*c+3*e*b*d : a*d+b*c).                            (15)
```

It uses neither square roots nor approximate p-adic arithmetic. All
inverted denominators are certified units. Both methods give

```
u_T(alpha)/p = 4+128*theta mod p,
Y_T(alpha)/p = 63+6*theta mod p,
Y_T'(alpha)  = 21+32*theta mod p.                                (16)
```

Thus `v(Y_T(alpha))=1`, so (14) gives ell=1 and the negative abscissa
valuation is exactly -2. The contact displacement is

```
g=(tau-alpha)/p mod p
 = -(63+6*theta)/(21+32*theta)
 = 91+113*theta.                                                (17)
```

For the primitive negative abscissa, the degree-four identity forces

```
p^2*x_R mod p = Lambda*q_0^2,       Lambda in F_131^*.             (18)
```

The centre congruence in (13) and (9), with ell=1, imply
`(gamma-alpha)/p=g/2 mod p`. In (12), the reduction of h(alpha) is
`Lambda*q_0'(theta)^2`. Evaluate at beta: x_R-rho reduces to3N,
where `N=29+128*theta` is the double root. Since beta-alpha has valuation2,
this gives

```
Lambda*q_0'(theta)^2*(g/2)^2 = 3*N,
Lambda = 12*N/(q_0'(theta)^2*g^2) = 69+40*theta.                  (19)
```

The last element does not belong to F_131. It contradicts (18). Changing
the sign of T changes both numerator and derivative in (17), so the
displacement and contradiction are unchanged. All mixed pairs are excluded.

## 5. A constraint on the remaining two-negative case

Return to the general unit scalar eta in (1); the square-scalar conclusion
above used an integral section. Apply (13) separately to two negative
original sections. The quantity ell is the valuation of the same branch
centre's displacement from alpha, so their valuations coincide at
`v(x_i)=-2*ell`. Their rescaled gamma displacements coincide by the last
congruence in (13). Evaluating (12) at each beta therefore forces the same
primitive abscissa scalar. Consequently

```
p^(2*ell)*x_1 mod p = p^(2*ell)*x_2 mod p = Lambda*q_0^2.
```

The primitive coefficient identities then give equal squares for the two
ordinate scalars. After changing one section's sign, which preserves
orthogonality, they too agree:

```
p^(3*ell)*r_1 mod p = p^(3*ell)*r_2 mod p = mu*q_0^2,
Lambda^3=eta_bar*mu^2.
```

Thus the first terms cancel in both differences. The common ell remains
unbounded here; the integral-contact calculation giving ell=1 cannot be
applied to this case.

## 6. Evidence and remaining scope

The [frozen input](../artifacts/generated-results/elkies-k3-q80-genus-one-k0-orthogonal-contact-v1/input.json)
binds the literal parent, retained norm-four word, nodal/pole premises,
producer, checker and preflights. The
[independent receipt](../artifacts/generated-results/elkies-k3-q80-genus-one-k0-orthogonal-contact-v1/independent-replay.json)
checks all131 coefficient values, the rank13 minor and tangent, and the
nodal word modulo131^2, including the nonrational scalar. The finite
producer used0.1291 CPU seconds and the independent replay0.0085 seconds,
under30 CPU seconds and2GiB each. Both frozen runs passed.

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 research/elkies-k3/scripts/verify_q80_k0_mixed_exclusion.py
```

The formal deformation, its exhaustiveness for nearby integral solutions,
the ramified branch valuation comparison and the unbounded exclusion are
written proofs. Finite replay does not independently or formally verify
those analytic arguments, and no external review is claimed. The existing
large quartic census was reused, not restarted.

Only the orthogonal `k=0,j=2` case with both original abscissas negative
remains for this Q80 polynomial chart. Both sum/difference pole quadratics
still reduce to q_0 and their branch values can be ramified over K. The
old genus-zero norm-eight singular-polynomial replay gap is neither used
nor upgraded. Other actual MW17 parents, higher-height rational abscissas
and the required infinite rank-at-least19 construction remain open.
