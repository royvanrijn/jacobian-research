# The Q80 polynomial genus-one common-quartic chart is empty

On the literal direct11952 alternate-Q80 MW17 parent, there is no pair
of distinct polynomial abscissas satisfying

```
x_i^3+A*x_i+B = eta*D*r_i^2,
deg x_i,deg r_i<=4,
D a squarefree binary quartic, gcd(D,Delta)=1, eta!=0.            (1)
```

All coefficients are rational. The statement includes arbitrary
denominators and scalar classes, higher contacts, infinity and genus-one
bases without a rational point. It concerns this polynomial chart with
branching at smooth parent fibres. It does not exclude higher-height
rational abscissas, covers branching at singular fibres, or other parents.

The new part excludes the last `k=0,j=2` case with both original abscissas
negative at131. A meromorphic square-root factorization forces the sum
of the two sections to have branch separation incompatible with either
original section. No denominator search or new finite-field census is used.
The required infinite rank-at-least19 construction remains open.

## 1. The last case and the sum's leading term

The [quartic gate](Q80_GENUS_ONE_K0_QUARTIC_GATE_2026-09-14.md),
[pole-degree proof](Q80_GENUS_ONE_K0_POLE_BOUNDARY_2026-09-14.md) and
[mixed-valuation exclusion](Q80_GENUS_ONE_K0_MIXED_VALUATION_2026-09-14.md)
leave only an orthogonal pair of height-eight sections with

```
p=131,          q=t^2+62*t+88,           D_bar=q^2,
v(x_1)=v(x_2)=-2*ell,                    ell>=1,
p^(2*ell)*x_i mod p = Lambda*q^2,
p^(3*ell)*r_i mod p = mu*q^2,
Lambda^3=eta_bar*mu^2,                   v(eta)=0.                (2)
```

Here v is the coefficient Gauss valuation, normalized by v(p)=1;
Lambda and mu are nonzero elements of F_131. The sign of one section
has been chosen to align the ordinate reductions, as proved in section5
of the mixed-valuation note. Sign changes preserve orthogonality. No
bound on ell is imported from that note's separate integral-seed argument.

Write `C:w^2=eta*D`, `P_i=(x_i,w*r_i)`, and `Q=P_1+P_2`. Its height is16.
The descended pole divisor has degree two, with a monic integral equation
H satisfying `H_bar=q`. Thus

```
x_Q=N/H^2,             deg N<=8,       deg H=2.                  (3)
```

Both P_i select the two different near-node roots at each branch point;
Q selects the third, simple root. In particular Q has no pole at a branch.

The completed coefficient Gauss field has good elliptic reduction, after
an unramified constant extension if needed. The formal parameter `z=-x/y`
at O has integral formal group law and `x=z^(-2)+O(z^2)`. Equations (2)
give the same leading term for z(P_1) and z(P_2), of valuation ell.
Since2 is a unit, z(Q) has twice that leading term. Consequently

```
v(x_Q)=-2*ell,
p^(2*ell)*x_Q mod p = (Lambda/4)*q^2,
p^(2*ell)*N mod p = (Lambda/4)*q^4.                              (4)
```

These are reductions of rational functions, followed by clearing the
primitive denominator H. They do not assume that x_Q is polynomial.

## 2. The original sections force a tighter branch cluster

Let K/Q_131 be unramified quadratic. In either q residue disc, let alpha
be the actual nodal base value and write `u=t-alpha`. The retained literal
local equation is

```
f(t,x)=(x-rho)*(U^2-V),       U=x+rho/2,
rho in O_K[[u]] a unit,
V=u*v_0(u),                 v_0 and 3*rho^2+A units.              (5)
```

Let a,b be the two roots of D in this disc. They may initially lie in a
ramified extension of K. The negative polynomial contact calculation
in the [nodal proof](Q80_GENUS_ONE_K2_NODAL_ORIENTATION_2026-09-14.md)
and section3 of the mixed-valuation note gives

```
v(a-alpha)=v(b-alpha)=ell,            v(a-b)>ell.                 (6)
```

Only this strict inequality is needed here, not the stronger equality
`v(a-b)=3*ell/2`. Its origin is useful: preparation gives
`x_i-rho=p^(-2*ell)*h_i(t)*(t-gamma_i)^2` with h_i a unit.
The remaining ordinate contact beta_i is a double zero of `U^2-V`.
After scaling `t=gamma_i+p^ell*z`, beta_i consumes one quadratic cluster
and a,b occupy the opposite cluster, with the same nonzero residue.
Thus their separation is strictly greater than ell. This argument uses
the two contacts gamma_i,beta_i in K, and does not require a,b in K.

## 3. A single pole forces two simple contacts at a different scale

We prove a local lemma with an arbitrary positive integer L. Suppose a
meromorphic abscissa x on the disc satisfies the section equation in (1),
selects rho at both a,b, and has a unique double pole xi in the disc.
Suppose xi is K-rational and, writing its primitive pole denominator
locally as `H=(t-xi)*h_0` with h_0 an integral unit,

```
p^(2*L)*H^2*x is integral with reduction of order4.               (7)
```

Then the only two branch roots in the disc have separation valuation L.
All zeros and poles in this argument are counted geometrically.

Put `J=U^2-V`. At a,b, x=rho and J=`3*rho^2+A` is a unit. A zero of
x-rho cannot be a zero of J, by the same unit condition. At every zero
of J, the identity `f=eta*D*r^2` therefore gives an even multiplicity.
The only pole of J is xi, of order four.

The integral series `p^(4*L)*H^4*J` has primitive reduction of order eight,
which is the square of the reduction in (7). Its distinguished polynomial
has only even zeros, so it is the square of a monic degree-four polynomial.
The remaining integral unit has a square root after an unramified constant
extension. This is the same complete-DVR preparation used in the earlier
pole proof; see [Berger, section1 and Corollary1.2](https://perso.ens-lyon.fr/laurent.berger/articles/article33.pdf).
It supplies a meromorphic square root Z of J, with at most a double pole
at xi and no other poles.

Choose its sign so that `p^(2*L)*H^2*Z` has the same reduction as
`p^(2*L)*H^2*U`. Define

```
F=U+Z,             G=U-Z,              F*G=V.                    (8)
```

The primitive numerator `p^(2*L)*H^2*F` then has Weierstrass degree four.
The denominator H^2 has degree two on this disc. Therefore

```
degree(zeros(F))-degree(poles(F))=2.                             (9)
```

This counts cancellations as well; it does not treat F as pole-free in
advance. At xi one of F,G has a double pole, since their sum has one and
their product has no pole. If xi=alpha the other has a triple zero;
the two divisor degrees would be -2 and3, contradicting (9). Hence
xi differs from alpha.

For xi different from alpha, the factor with a pole has order -2 at xi
and the other has order2. Away from xi, the only possible zero of either
factor is the single zero alpha of V. Assigning that zero to one factor
gives possible divisor degrees -2,-1,2,3. Equation (9) forces F to have
its double zero at xi and no zero at alpha; G has the pole at xi and
the simple zero at alpha. In particular

```
F=p^(-2*L)*h(t)*(t-xi)^2,              h an integral unit.         (10)
```

To see the stated integral normalization, its previously prepared
degree-four numerator must be `(t-xi)^4` times an integral unit.
Cancel H^2=`(t-xi)^2*h_0^2`. Thus (10) follows with no uncontrolled
meromorphic or unbounded unit factor.

The equation `F^2-3*rho*F+V=0` has two integral analytic root functions
on the whole disc. Write them as W_0 and W_1, with

```
W_0 mod p = 3*rho mod p,       W_0 a unit,
W_1=V/W_0,                    W_1 vanishing simply at alpha.     (11)
```

The derivative at either residue root is a unit, so these functions
exist by Hensel lifting. Consider just `F=W_0`, and substitute
`t=xi+p^L*z` into (10). The equation reduces to

```
h(xi)_bar*z^2-W_0(xi)_bar=0.                                    (12)
```

It has two distinct nonzero roots over a finite residue extension,
and both lift in the corresponding unramified extension. Call the
lifted base values b_+,b_-. They are distinct from xi and satisfy

```
v(b_+-xi)=v(b_--xi)=v(b_+-b_-)=L.                               (13)
```

Both contacts are simple: the derivative of F has valuation -L there,
whereas the derivative of W_0 is integral. Also, by (8),

```
x-rho = (F-W_0)*(F-W_1)/(2*F).                                  (14)
```

At b_+,b_-, F and F-W_1 are units. Thus x-rho, and hence f(t,x),
has a simple zero. Neither value is a pole. The section identity
`f=eta*D*r^2` forces both values to be roots of D: an even-order
ordinate factor cannot account for a simple zero.

There are exactly two roots of D in the residue disc. They must therefore
be b_+,b_-, proving the lemma. In particular the two simple contacts
cannot be discarded in favour of a closer pair solving `F=W_1`.

## 4. Contradiction and complete allocation coverage

For Q in (3), H has simple reduction q, so it has exactly one simple
root xi in each nodal residue disc, lying in K. The actual double pole
does not cancel: H is the denominator of its pole divisor. Equation (4)
gives precisely (7) with L=ell. Q selects rho at both branch roots.
The lemma therefore gives

```
v(a-b)=ell,
```

contradicting (6). This excludes the two-negative orthogonal case over
Q_131. No disjointness condition on the sum and difference pole divisors
was used; the sum alone gives the contradiction. No point on C, local
branch splitting assumption, bound on ell or integral inherited seed
was used either.

The retained complete factor presentation has fifteen genus-one cases.
The following proof coverage now excludes all of them over Q:

| Agreement degree k | j values | Exclusion |
|---:|---|---|
|4|0|[Unramified descent and smooth genus-one injectivity](COMMON_QUARTIC_BRANCH_STRATA_2026-09-14.md#4-exclusion-of-the-entire-k4-genus-one-stratum-on-q80)|
|3|0,1|[Single-branch reciprocity](Q80_SINGLE_BRANCH_RECIPROCITY_2026-09-14.md)|
|2|0,1,2|[Nodal orientation](Q80_GENUS_ONE_K2_NODAL_ORIENTATION_2026-09-14.md)|
|1|0,1,2,3|[Cubic reciprocity](Q80_GENUS_ONE_K1_CUBIC_RECIPROCITY_2026-09-14.md)|
|0|0,1,2,3,4|[Integral quartic gate](Q80_GENUS_ONE_K0_QUARTIC_GATE_2026-09-14.md), [four pole-degree exclusions](Q80_GENUS_ONE_K0_POLE_BOUNDARY_2026-09-14.md), [mixed exclusion](Q80_GENUS_ONE_K0_MIXED_VALUATION_2026-09-14.md), and this proof|

The full-chart conclusion is stated over Q: the k4 descent uses the
retained global arithmetic theorem. The new last-case contradiction
itself is over Q_131. The binary degree bounds and retained coordinate
patch theorem include all degree drops and infinity cases.

## 5. Assurance and remaining objective

The [frozen premise input](../artifacts/generated-results/elkies-k3-q80-genus-one-polynomial-closure-v1/input.json)
binds this note, the earlier literal equation, checked local data, source
proofs and completed receipts. The [premise checker](scripts/verify_q80_genus_one_polynomial_closure.py)
rechecks the finite nodal premises, validates the old pole/mixed receipts
and enumerates the fifteen proof allocations. Its
[receipt](../artifacts/generated-results/elkies-k3-q80-genus-one-polynomial-closure-v1/receipt.json)
is a check of those premises and source bytes.

```
python3 research/elkies-k3/scripts/verify_q80_genus_one_polynomial_closure.py check
```

The meromorphic factorization, formal-group reduction, forced simple
contacts and unbounded contradiction are written proofs. They are not
independently or formally verified by the finite checker, and no external
review is claimed. Prior quartic/cubic censuses are reused, not rerun.
The frozen intermediate UNKNOWN boundaries are retained as records of
their component scopes; this note supplies their successor conclusion.

The earlier genus-zero norm-eight singular-polynomial independent-replay
gap is not used or upgraded. Other actual MW17 parents, higher-height
rational abscissas and the required single quadratic cover with infinitely
many rational points and two certified new independent sections remain
open. This exclusion is not a positive rank-at-least19 construction.
