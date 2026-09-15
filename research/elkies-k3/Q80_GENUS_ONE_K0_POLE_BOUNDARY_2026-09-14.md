# Only the orthogonal Q80 genus-one k0 allocation survives

On the literal direct11952 Q80 MW17 parent, the polynomial genus-one
common-quartic allocations `k=0,j=0,1,3,4` are empty over Q_131, hence
over Q. The proof includes arbitrary coefficient denominators and scalar
classes, higher contacts, infinity and unpointed genus-one bases.

The sole remaining polynomial genus-one allocation on Q80 is `k=0,j=2`.
It requires branch reduction `q_0^2`, at least one negative abscissa
valuation, and two degree-two pole divisors, one for each sum/difference,
both reducing to `q_0=t^2+62*t+88`. This boundary is **UNKNOWN**, in both
the mixed integral/negative and the two-negative case. No infinite
rank-at-least19 source is constructed.

## 1. The retained boundary and its local model

The [quartic norm and contact theorem](Q80_GENUS_ONE_K0_QUARTIC_GATE_2026-09-14.md)
already confines every k0 pair in the full chart

```
x_i^3+A*x_i+B = eta*D_0*r_i^2,
deg x_i,deg r_i<=4,
D_0 squarefree binary quartic, gcd(D_0,Delta)=1,
gcd(D_0,x_1-x_2)=1                                                (1)
```

to `D_0,bar=q_0^2` after normalizing a unit, with at least one negative
Gauss abscissa valuation. The quartic norm census, complete norm-four
contact roster and first lifts are retained evidence and are not repeated
here. Normalize `v_131(eta)=epsilon` in `{0,1}`.

Let K be the unramified quadratic extension of Q_131. In either nodal
residue disc, the actual bad base value alpha lies in K. The four roots
of D_0 give two geometric branch values a,b in this disc and two in its
conjugate. The values a,b can be ramified over K; neither is assumed to
belong to K. They are distinct and avoid alpha, since the characteristic-zero
branch fibres are smooth.

The exact local data from the
[nodal orientation proof](Q80_GENUS_ONE_K2_NODAL_ORIENTATION_2026-09-14.md#2-an-analytic-equation-at-the-node)
give an analytic simple cubic root rho and the factorization

```
f(t,x)=(x-rho)*(U^2-V),       U=x+rho/2,
rho in O_K[[u]],             V=u*v_0(u),
u=t-alpha,                  rho, v_0, 3*rho^2+A units.             (2)
```

The two near-node roots reduce to `N=29+128*t`, whereas rho reduces to
`S=73+6*t=-2*N mod q_0`. These are nonzero and distinct. In particular V
has exactly one zero, of multiplicity one, in the open geometric residue
disc. The imported analytic fact below is complete-DVR Weierstrass
preparation, as stated in
[Berger, section1 and Corollary1.2](https://perso.ens-lyon.fr/laurent.berger/articles/article33.pdf).
All applications are to bounded series with specified primitive reductions.

## 2. A pole-free negative abscissa cannot select the simple root twice

We first prove a local lemma. Suppose x belongs to `p^(-m)*O_K[[u]]`,
with m positive, and the primitive series `p^m*x` has reduction of order
`n>=2`. Suppose also that

```
f(t,x)=eta*D_0*r^2
```

for meromorphic r on the disc, and that x(a)=rho(a), x(b)=rho(b) at its
two branch values. The scalar eta is arbitrary and nonzero. Then such x
does not exist.

Put `H=U^2-V`. At a,b, H is a unit because `3*rho^2+A` is a unit.
At every other zero of H, both D_0 and x-rho are nonzero; equation (2)
therefore gives an even zero multiplicity for H. Zeros of x-rho cannot
also be zeros of H, by the same unit condition. Since H has no poles in
the disc, all its zeros have even multiplicity.

The primitive series `p^(2*m)*H` reduces to the square of the reduction
of `p^m*x`, and has Weierstrass degree2n. Its distinguished polynomial is
a square of a monic degree-n polynomial: the half-zero divisor is unique
and Galois invariant, so the square root polynomial belongs to O_K[u].
Its remaining factor is an integral unit. After at most an unramified
quadratic extension, that unit has an integral power-series square root
(p is odd). Thus H has a bounded analytic square root Z on the whole
disc, with `p^m*Z` integral.

Choose the sign of Z so that its primitive reduction equals that of
`p^m*U`. Then `p^m*(U+Z)` reduces to twice that same series and has
Weierstrass degree n. In particular U+Z has at least two zeros in the
disc. Both U+Z and U-Z are analytic, with no poles, but

```
(U+Z)*(U-Z)=V.                                                    (3)
```

The right side has only one simple zero. This is impossible. The argument
uses geometric zeros and imposes no splitting condition on a,b.

The pole hypothesis matters. Formula (3) is a zero-count contradiction
only when both factors have no poles in the disc; a possible cancellation
against a pole is not discarded in the surviving case below.

## 3. Every original negative section selects node roots

For an original negative abscissa, write

```
x=p^(-m)*X, r=p^(-n)*R,
2*n=3*m+epsilon,
X_bar^3=eta_unit,bar*q_0^2*R_bar^2.                               (4)
```

X,R are primitive integral binary forms of degree4. Parity and the degree
budget force `X_bar=lambda*q_0^2` and `R_bar=mu*q_0^2`, with lambda,mu
nonzero. In each nodal disc, `p^m*(x-rho)` and `p^n*r` consequently have
Weierstrass degree two.

Every zero of x-rho outside the simple-root branch selections has even
multiplicity. Its total zero degree is two. The number of simple-root
branch selections is therefore either zero or two. Two are impossible
by section2 with n=2. Hence every original negative section chooses a
near-node root at both a and b, and likewise in the conjugate disc.

This also proves that eta is a unit. With no simple-root branch selection,
x-rho has one double zero gamma in K, and

```
x-rho=p^(-m)*h(t)*(t-gamma)^2,       h integral unit.               (5)
```

The distinguished quadratic is a square, so gamma lies in K. At gamma,
r has exactly a simple zero. Its degree-two distinguished polynomial
has one further simple zero beta, also in K, distinct from gamma. At beta,
H vanishes and x-rho is a unit, because a zero of H in this disc reduces
to the node. Evaluating (5) gives

```
m=2*v(beta-gamma).
```

Thus m is even. Equation (4) forces epsilon zero. This step does not
require a,b in K, an ordinary nodal contact, or a rational point on the
cover. It uses only gamma,beta in K. No branch-orientation contradiction
is asserted for k0.

## 4. Heights determine the available pole divisors

On `C:w^2=eta*D_0`, put `P_i=(x_i,w*r_i)` and `Q_-=P_1-P_2`,
`Q_+=P_1+P_2`. The pullback elliptic surface has chi=4 and only
irreducible fibres. The retained
[intersection calculation](COMMON_QUARTIC_BRANCH_STRATA_2026-09-14.md#2-the-factors-determine-the-height-matrix)
gives

```
height(P_i)=8,        <P_1,P_2>=4-2*j,
height(Q_-)=8+4*j,    height(Q_+)=24-4*j.                          (6)
```

At each branch point, the two P_i are different nonzero2-torsion points,
so both Q_- and Q_+ are the third nonzero2-torsion point. In particular,
neither Q meets O above a branch value.

For either Q, the height formula is `height(Q)=8+2*(Q.O)`. The intersection
divisor with O is invariant under the covering involution, since Q is
anti-invariant and negation preserves its intersection with O. It avoids
the branch points and therefore descends to a divisor Z on the original
parameter line, with

```
(Q.O)=2*deg Z,        height(Q)=8+4*deg Z.                         (7)
```

Locally the abscissa has twice the pole order recorded in Z. Thus, with
s=deg Z, homogeneous rational coordinates can be written

```
x_Q=N/H^2,   r_Q=M/H^3,
deg H=s,     deg N=4+2*s,     deg M=4+3*s,                         (8)
```

where H records the poles, including infinity, and the indicated degrees
are binary degrees. Normalize H to be primitive integral. Formula (8)
follows from the minimal Weierstrass coordinate bundles; it does not
assume that Q has polynomial coordinates or miss O.

Choose Q to have the smaller height in (6). Then

| j | Smaller height | Degree s of its pole divisor |
|---:|---:|---:|
|0|8|0|
|1|12|1|
|2|16|2|
|3|12|1|
|4|8|0|

For s<=1 the reduction of H cannot vanish at a root of the irreducible
quadratic q_0. Thus H is an integral analytic unit throughout both nodal
residue discs. A rational pole at infinity does not alter this conclusion.

## 5. The two possible valuation patterns exclude s<=1

Complete the coefficient Gauss valuation on Q_131(t) and extend it to the
covering function field. Since eta is now a unit and D_0 reduces to q_0^2,
this gives a good-reduction elliptic curve over that valuation field,
possibly after the constant unramified extension. Reduction is a group
homomorphism. A negative abscissa reduces to O, and every original
integral section reduces to a finite point.

**Both original abscissas negative.** Both P_i reduce to O, so the
nonzero Q also reduces to O and its abscissa has negative Gauss valuation.
At every branch, the two P_i choose the two different near-node roots
by section3. Hence Q chooses the simple root rho at all four branches.

Write `x_Q=p^(-m)*N_0/H^2` with N_0 primitive integral. With a similar
primitive normalization of M, the leading coefficient identity from
(8) gives

```
N_0,bar^3=eta_bar*q_0^2*M_0,bar^2.                               (9)
```

The multiplicity of q_0 in N_0,bar is positive and even, hence at least
two. For s<=1 the degree bound `deg N_0<=6` makes it exactly two.
As H is a unit on the nodal discs, the primitive local abscissa has
reduction of order two there. Section2 now contradicts its two simple-root
branch selections.

**One original abscissa integral.** Denote its section by P and the
negative section by R. The integral P reduces, on the constant twist,
to an ordinary section with ordinate `sqrt(eta_bar)*q_0*r_bar` after an
unramified constant extension. Such a section cannot meet an I1 node:
differentiating its equation there would give `0=f_t`, while f_t is a
unit. It therefore reduces to S and chooses rho at both local branches.
The negative R chooses a near-node root. Consequently either P+R or
P-R chooses the other near-node root at those branches.

Under the coefficient Gauss valuation, the abscissa of either sum/difference
reduces to the abscissa of P, because R reduces to O. Thus x_Q is integral
in that valuation. If s<=1, its primitive denominator H is a unit on the
nodal discs, so its bounded local reduction is the same as that of P.
Its branch values therefore reduce to S, contradicting its near-node
root selections. This also covers the case where the chosen difference
reduces to -P, since negation fixes the abscissa.

Both valuation patterns are impossible for s<=1. Equations (6)-(7) exclude
all four allocations `j=0,1,3,4`.

## 6. What the orthogonal case must still supply

For j=2, both Q_+ and Q_- have height16 and both pole divisors have
degree two. If either primitive H had reduction coprime to q_0, the
mixed-case argument would be unchanged. In the two-negative case, (9)
would give primitive local order at least two (now possibly four), and
the general lemma in section2 would still contradict the simple-root
branch choices. Thus both denominators necessarily satisfy

```
H_+,bar proportional to q_0,       H_-,bar proportional to q_0.    (10)
```

After making them monic, H_+ and H_- are integral quadratics with simple
irreducible reduction q_0. Their roots lie in K, one in each nodal residue
disc. Their divisors are not asserted to be disjoint. This is a condition
on the poles of the sum and difference, not a factorization of the
quartic branch divisor D_0: its local branch points may still be ramified
over K. The original two sections still have height matrix `diag(8,8)`
in this one surviving allocation.

No coefficient solution, independent pair, rational base point or positive
MW17 construction follows from (10). Other actual MW17 parents and
higher-height rational abscissas also remain possible. The earlier
genus-zero norm-eight singular-polynomial replay gap is not used or upgraded.

## 7. Assurance

The [source-bound receipt](../artifacts/generated-results/elkies-k3-q80-genus-one-k0-pole-boundary-v1/receipt.json)
and [frozen input](../artifacts/generated-results/elkies-k3-q80-genus-one-k0-pole-boundary-v1/input.json)
are checked by [the premise checker](scripts/verify_q80_k0_pole_boundary.py).
This is a written local-analytic and height proof. Its finite premises are
the previously replayed literal nodal model, the complete k0 integral gate
and the full branch-stratum height identity. The accompanying source-bound
receipt checks those identities and the finite pole-degree table; it does
not mechanically verify preparation, analytic square-root existence,
Gauss reduction or the divisor descent argument. No new large census or
bounded denominator search is used, and no formal or external verification
is claimed.
