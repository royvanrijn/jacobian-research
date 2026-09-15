# Q80 genus-one k2 is empty by nodal branch orientation

On the literal direct11952 alternate-Q80 MW17 parent, all three polynomial
genus-one common-quartic allocations `k=2,j=0,1,2` are empty over Q_131,
and hence over Q. This includes every coefficient denominator, literal
scalar, higher contact and unpointed base allowed by the chart. The proof
closes the [nodal denominator boundary](Q80_GENUS_ONE_K2_NODAL_BOUNDARY_2026-09-14.md)
left by the earlier theorem; it does not strengthen that frozen certificate.

This component leaves `k=0` and `k=1` outside its exclusion. The later
[cubic reciprocity proof](Q80_GENUS_ONE_K1_CUBIC_RECIPROCITY_2026-09-14.md)
closes k1; five genus-one allocations `k=0,j=0,...,4` remain open.
Other actual MW17 parents and higher-height rational abscissas also remain
open. The required [MW17 parent with two independent gains on an infinite
quadratic base](CORRELATED_QUADRATIC_GAINS_2026-09-12.md) is unconstructed.
The old genus-zero norm-eight singular-pencil replay gap is not used here.

## 1. Exact remaining hypotheses

Use the complete [polynomial chart](COMMON_QUARTIC_BRANCH_STRATA_2026-09-14.md):

```
f(t,x_i)=x_i^3+A(t)*x_i+B(t)=eta*g(t)*d(t)*r_i(t)^2,
deg x_i<=4, deg r_i<=4,
g,d monic quadratic, g*d squarefree, gcd(g*d,Delta)=1,
g=gcd(g*d,x_1-x_2),   x_1!=x_2.                              (1)
```

The last equality fixes agreement at both g roots and disagreement at both
d roots. The characteristic-zero branch fibres are smooth. The prior
boundary theorem proves, in the original t coordinate, that

```
p=131, q_0=t^2+62*t+88,    g_bar=d_bar=q_0,
x_i=p^(-m_i)*X_i, r_i=p^(-n_i)*R_i,
X_i,bar=lambda_i*q_0^2, R_i,bar=mu_i*q_0^2,
m_i>0, lambda_i*mu_i!=0.                                  (2)
```

With `v_p(eta)=epsilon` in `{0,1}`, it also gives
`2*n_i=3*m_i+epsilon`. Its integral coefficient cases, denominator
specialization argument and finite norm tables are retained premises.

Let K be the unramified quadratic extension of Q_131 and normalize its
valuation by `v(p)=1`. The two roots of each of g and d lie in K, one in
each of the two residue discs cut out by q_0. Work in either one. Write
a for the g root and b for the d root in that disc. Both are K-rational,
distinct and at smooth fibres.

Both sections choose **node roots** at a and b after reduction. This
assertion does not require integral abscissa coefficients: at an integral
branch value, a selected root of the integral monic cubic is integral.
The [degree-two reciprocity theorem](Q80_DEGREE_TWO_RECIPROCITY_2026-09-14.md)
applies to the rational sections themselves. The unused root at the d
branches must have zero norm character. The nodal root multiset has codes
`(13412,13412,0)`, with zero only at its simple root. Thus both sections
choose the double root at d. Each section's character product over all
four branches then forces code13412 at g, again the double root. These
norm evaluations specialize on the smooth cubic root curve, even at the
ramification above the nodal fibre; no coefficient reduction of x_i is
being substituted into a character formula.

## 2. An analytic equation at the node

The reduced discriminant has q_0 as a simple factor and its cubic has a
double nonzero root N and a distinct simple root S:

```
N=29+128*t,       S=73+6*t=-2*N modulo q_0.                  (3)
```

Hensel lifting gives a unique actual nodal base value alpha in the chosen
disc, and an integral analytic simple-root function rho(t), with
`rho(alpha)=-2*e`, where e is the double cubic root. Put `u=t-alpha`.
Then rho belongs to `O_K[[u]]`, rho and `3*rho^2+A` are units throughout
the open residue disc, and

```
f(t,x)=(x-rho(t))*(U(t)^2-V(t)),
U=x+rho/2,       V=-A-3*rho^2/4.                            (4)
```

The identity

```
-4*A^3-27*B^2=4*V*(3*rho^2+A)^2                            (5)
```

shows that V has a simple zero at alpha, with V'(alpha) a unit.
All Taylor coefficients are integral; V' is a unit throughout the disc.
Consequently `v(V(t))=v(t-alpha)` there and the difference quotient of V
at any two distinct points of this disc is a unit. The finite certificate
checks the actual residues of rho, its implicit derivative and V' from
the literal rational A,B. It does not replace rho'(alpha) by the derivative
of the degree-one polynomial S representing its residue value.

In `F_131[theta]/(q_0(theta))`, the exact residues are
`rho'=51+97*theta`, `V'=118+113*theta`,
`3*rho^2+A=48+93*theta` and `(4*A^3+27*B^2)'=78`.

We use the classical Weierstrass preparation theorem over the complete
DVR O_K: a series whose reduction has order n is an integral unit times
a distinguished monic polynomial of degree n. In particular it has exactly
n zeros in the open geometric residue disc, counted with multiplicity.
See [Berger, section1, classical statement and Corollary1.2](https://perso.ens-lyon.fr/laurent.berger/articles/article33.pdf).
Units and their derivatives have integral Taylor coefficients throughout
this disc. This is the only power-series preparation imported below.

## 3. One section must switch node roots

Fix either section and suppress its index. By (2), both
`p^m*(x-rho)` and `p^n*r` have Weierstrass degree two in this disc.
Every zero of x-rho has even multiplicity: it cannot be a or b, which
select node roots, and the other cubic factor in (4) is a unit there.
The identity (1) therefore makes its multiplicity twice that of r.
Thus x-rho has just one double zero gamma. Its distinguished quadratic
is a square, so gamma lies in K (p is odd), and

```
x-rho=p^(-m)*h(t)*(t-gamma)^2,    h in O_K[[t-alpha]]^*.     (6)
```

At gamma, r has multiplicity one. Its distinguished quadratic therefore
has one further simple K-rational zero beta, distinct from gamma. The
factor `H=U^2-V` vanishes to order at least two at beta. This initially
allows beta to be a branch value and H to have order three; that possibility
will be excluded by differentiation, without assuming ordinary contacts.

At any zero of H in the disc, x reduces to the node root. Hence x-rho is
a unit there. Applying (6) at beta gives

```
m=2*v(beta-gamma)=2*ell,       ell>=1.                       (7)
```

Here the valuation is an integer because beta,gamma lie in K. Differentiate
(6). At beta its term `2*p^(-2*ell)*h*(beta-gamma)` has valuation `-ell`;
the term involving h' has valuation at least zero, as do derivatives of
rho. Therefore `v(U'(beta))=-ell`. Since `H'(beta)=0` and V' is a unit,

```
v(U(beta))=ell,
v(beta-alpha)=v(V(beta))=2*ell,
v(gamma-alpha)=ell.                                        (8)
```

The second derivative settles all higher contacts:

```
H''=2*(U')^2+2*U*U''-V'',
v(2*(U')^2)=-2*ell,
v(2*U*U'')>=-ell,           v(V'')>=0.                      (9)
```

The first term has strictly smallest valuation. Thus H'' is nonzero and
beta is exactly a double zero of H. In particular beta is neither a nor b:
at either branch (1) would give odd order `1+2*ord(r)` for H.

Now scale `t=gamma+p^ell*z`. Every zero of H has x-rho a unit, so (6)
places it on the annulus `v(t-gamma)=ell`. The scaled H has integral
coefficients and reduction

```
H(gamma+p^ell*z)_bar=(h_bar(gamma)*z^2+3*rho_bar/2)^2.       (10)
```

There are two distinct nonzero z roots, since p is odd and rho is a unit.
The root supplied by beta lies in the residue field of K, so the other
does too. Preparation in each z residue disc gives exactly two zeros.
The double zero beta exhausts its disc. Equation (1), together with the
already accounted-for double zero gamma of x-rho, leaves a and b as the
other two zeros of H. Hence a and b occupy the same other z disc. It follows
that

```
v(a-beta)=v(b-beta)=ell,
v(a-alpha)=v(b-alpha)=ell,
N_ab=v(a-b)>ell.                                            (11)
```

Since `U(a)^2=V(a)` and a lies in K, ell is even. In particular (7) already
forces epsilon zero through `2*n=3*m+epsilon`, without a rational point
assumption on the covering curve.

Compare the two branch evaluations of x using (6). The difference quotient
of `h(t)*(t-gamma)^2` at a,b has valuation ell: its leading term is a unit
times `(a-gamma)+(b-gamma)`, whose valuation is ell because a,b are in the
same nonzero z residue class and p is odd. The term from the difference
of h has valuation at least `2*ell`. The integral rho term cannot cancel
the term of valuation `-ell` after multiplication by `p^(-2*ell)`. Thus

```
v(x(a)-x(b))=N_ab-ell.                                      (12)
```

On the smaller disc `v(t-a)>ell`, the two near-node roots can be labelled
analytically as

```
e_+(t)=-rho(t)/2+s(t),    e_-(t)=-rho(t)/2-s(t),   s(t)^2=V(t).
```

Indeed V(a) is a nonzero square in K, and V(t)/V(a) differs from1 by an
element of positive valuation there. Hensel's square root with value1
labels s uniquely after choosing s(a). For either fixed label, the unit
difference quotient of V and `v(s(a)+s(b))=ell/2` give

```
v(e_+(a)-e_+(b))=v(e_-(a)-e_-(b))=N_ab-ell/2.               (13)
```

The integral rho difference has larger valuation. Equations (12) and (13)
are incompatible because ell is positive. Therefore **the section chooses
opposite node-root labels at a and b**. Comparing opposite labels also
gives `N_ab=3*ell/2`; this last equality is a useful control, not a further
assumption in the contradiction.

## 4. Two sections cannot have k2

For each section, (11) identifies ell with the same intrinsic quantity
`v(a-alpha)`. The smaller disc and the pair of analytic node-root labels
are consequently common to both sections. They agree at a by definition
of g. Each must switch labels at b by section3, so they agree at b as well.
This contradicts the definition of d in (1). One of the two conjugate
nodal residue discs already suffices. All three j allocations are excluded.

This proof uses no coefficient-height, denominator or p-adic precision
bound. The scalar is arbitrary and the covering genus-one curve need not
be pointed. Branch at a characteristic-zero singular fibre, extra poles
or larger-degree abscissas would change the hypotheses and are not included.

## 5. Exact local controls and replay boundary

An explicit local I1 model tests the orientation and contact assertions:

```
f(t,x)=(x+2)*((x-1)^2-t),   rho=-2, V=t,
c=z^2-3, s=4*z*c, beta=c^2, gamma=-3*c*(z^2+1),
x=-2+(t-gamma)^2/s^2,
Q=(t-beta)^2+4*s*z*(t-beta)+6*s^2*(z^2-1),
r=(t-gamma)*(t-beta)/s^3,
f(t,x)=Q*r^2,              disc(Q)=-8*s^2*c.                (14)
```

These are exact rational identities, independently expanded in two
variables. At p=131, z=30954 has `v(c)=2`. Its two Q_131 branch roots have
valuation2 and separation valuation3; their x values differ with valuation1
and have opposite leading node offsets. At z=38, `v(c)=1` and the quadratic
discriminant has odd valuation3, so its branch points require ramification.
At z=13793 the contact valuation is even but the residual discriminant
unit is a nonsquare over F_131: the branches lie in its unramified quadratic
extension. Even parity alone therefore does not assert Q_131 splitting.
The toy is a local control, not a K3/MW17 or infinite-base construction.

The [frozen input](../artifacts/generated-results/elkies-k3-q80-genus-one-k2-orientation-v1/input.json)
binds the literal parent, earlier boundary evidence and both new programs.
The [Sage producer](scripts/certify_q80_genus_one_k2_orientation.sage) verifies
the symbolic identities, actual nodal first derivatives and three controls.
The [independent standard-library checker](scripts/verify_q80_genus_one_k2_orientation.py)
uses sparse integer polynomial arithmetic, a separate finite-field
implementation and exact rational Hensel balls. Its checks establish the
finite hypotheses and controls; the all-denominator analytic proof above
is written mathematics, not a machine proof or external review.
The producer used0.0171 CPU seconds; the independent replay used0.0038 CPU
seconds and0.0052 wall seconds, each under20 CPU seconds and2GiB.
All nine targeted failure/control tests pass. The earlier complete norm
and coefficient censuses are reused with their frozen hashes and receipts.

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 research/elkies-k3/scripts/verify_q80_genus_one_k2_orientation.py
python3 -m unittest discover -s research/tests -p 'test_q80_genus_one_k2_orientation.py' -q
```

All runs have explicit CPU and memory caps. A discovery run stopped after
its first symbolic identity when Sage's nested-ring discriminant invoked an
incompatible variable ordering and its fallback exceeded the available
mapping space. That source and failure receipt are preserved. The fixed
quadratic discriminant formula and constrained BLAS startup were used before
freezing. A separate preview-loader path error stopped before its input
could be frozen and was corrected in the harness. No failed frozen run is
relabelled as a pass.
