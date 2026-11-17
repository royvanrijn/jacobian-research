# A projective pencil multiplicity bound and critical-point expansions

The [reduced-pencil theorem](REDUCED_PENCILS_AND_RATIONAL_MONODROMY.md)
closed one-exception parameters and all sparse two-positive-term outputs.
This note permits several nonreduced pencil members through one
degree-free ramification inequality. It also gives an exact coefficient
test for parameters affine in the second target coordinate, and closes
an all-degree family with mixed critical values and repeated leading
zeros. Claims `HC4MYIH0/1/2` are recorded in
[MATH_STATUS.json](../../MATH_STATUS.json).

Only the retained Petravchuk--Iena closed-rational equivalence and
projective Riemann--Hurwitz theorem are used beyond the existing algebraic
inputs. The [source audit](sources/pencil-multiplicity-context.json)
separates those proof imports from orbifold literature used as context.

## 1. A degree-free bound for projective pencil multiplicities (HC4MYIH0)

Let `B=C[x_1,...,x_n]`, `L=Frac(B)`, and let `N,M` be coprime and
algebraically independent. In particular `n>=2`. For finite `lambda`
put `F_lambda=N-lambda M`, and put `F_infinity=M`.
Let `Gamma` be a finite set containing every non-square-free projective
pencil member. For each `lambda in Gamma`, assume:

- `F_lambda` has an irreducible factor of multiplicity one;
- every factor multiplicity is at most an integer `M_lambda>=2`.

Define

```
sigma=sum_(lambda in Gamma) (1-1/M_lambda).    (1)
```

If `sigma<2`, then `C(N/M)` is relatively algebraically closed in `L`.
One may use the exact largest multiplicity for each nonreduced member;
larger supplied upper bounds are valid but give a weaker test.
All degrees and all numbers of source variables are unrestricted.

### 1.1 Proof

Suppose `N/M=Phi(h)` with `h=H_0/H_1` reduced and `deg Phi=d>1`.
As in [HC4MYIG0](REDUCED_PENCILS_AND_RATIONAL_MONODROMY.md), distinct
projective point forms `ell_c(H_0,H_1)` are coprime. None is constant:
a constant point form would express both `H_0,H_1` as affine polynomials
in one of them, and would make both homogeneous substitutions defining
`N,M` polynomials in one element. This contradicts their independence.

For a point `c` above `lambda`, with ramification index `e_c`, its
nonconstant point form occurs to power `e_c` in `F_lambda`. It follows
that `e_c<=M_lambda`. A multiplicity-one factor of `F_lambda` must come
from a point with `e_c=1`, and with reduced point form at that factor.
Thus at every `lambda in Gamma`, the fibre of `Phi` has an unramified
point. At a value outside `Gamma`, every point has index one.

Choose one unramified point in the fibre above `lambda`. The other
local degrees add to `d-1`. Since each is at most `M_lambda`,

```
sum_(c above lambda)(e_c-1)
    <= (d-1)(1-1/M_lambda).                  (2)
```

Indeed each remaining term is `e_c(1-1/e_c)`, at most
`e_c(1-1/M_lambda)`. Summing (2) and using the retained
Riemann--Hurwitz formula for the finite separable map `Phi:P1->P1`
gives

```
2d-2 <= (d-1) sigma.
```

Therefore `sigma>=2`, a contradiction. The explicit retained
Petravchuk--Iena Lemma 2(1) now gives relative algebraic closedness.
The proof uses maximum multiplicities, not their lcm, and does not
invoke an orbifold-cover classification.

### 1.2 The equality boundary is genuine

If `sigma=2` and a composition of degree `d>1` exists, equality holds
in every estimate above. At every exceptional value the ramification
profile of `Phi` must be

```
(1, M_lambda, ..., M_lambda),
d=1 mod M_lambda.                            (3)
```

In particular there is exactly one unramified point there. Every value
in `Gamma` must actually be a branch value of `Phi`. If the `M_lambda`
are exact integer maxima, the possible numerical lists at equality are

```
(2,2,2,2), (3,3,3), (2,4,4), (2,3,6).        (4)
```

This list follows arithmetically. Every summand in (1) is at least
one half and less than one, so there are three or four terms. Four
forces all two. For three terms, the reciprocals sum to one; sorting
the denominators gives the three displayed triples. No realization or
Galois-group classification is inferred from this numerical list.

The strict inequality cannot be replaced by `sigma<=2`. Put

```
N=s^3+t^3, M=t(s^2+t^2),
N/M=Phi(s/t), Phi(h)=(h^3+1)/(h^2+1).         (5)
```

These `N,M` are coprime and algebraically independent. The finite
members correspond to the cubic
`h^3-lambda h^2+(1-lambda)`, with discriminant

```
(1-lambda)(4lambda^3+27lambda-27).
```

It has four distinct zeros. At each, the profile is `(2,1)`; a triple
root would require the derivative root to have higher multiplicity,
whereas the displayed discriminant has simple zeros. At infinity the
three poles of `Phi` are distinct and simple. Thus the pencil has
four nonreduced members with maximum two and a simple factor each,
giving `sigma=2`, while its rational function is explicitly composite.
The examples in (5) are pencil controls, not Keller counterexamples.

### 1.3 Keller pullback preserves the certificate

The retained [HC4MYHW0](ABELIAN_CLASSGROUP_OUTPUT_OBSTRUCTION.md) proves
that distinct target irreducibles have coprime square-free nonconstant
pullbacks under a complex Keller map. Consequently every projective
member has exactly the same factor multiplicities after pullback,
with some distinct factors possibly splitting further. Its maximum
multiplicity and the presence of a multiplicity-one factor are unchanged.
Dominance preserves numerator-denominator independence. Hence every
target certificate (1) with `sigma<2` remains valid on the source,
in all Keller dimensions.

By [HC4MYIG1](REDUCED_PENCILS_AND_RATIONAL_MONODROMY.md), the full group
of every finite Galois cover of this parameter field is preserved.
The retained filtered class-monodromy criterion admits these parameters
as well, with every global filtration, affine height-one fibre and
torsion-free class-group hypothesis retained over the full target.
There is no general preservation assertion for arbitrary covers of the
whole multi-parameter target.

## 2. An exact coefficient test for affine-in-V parameters (HC4MYIH1)

Take target polynomials

```
q=(Q(T)V+C(T))/D(T),
Q!=0, D nonconstant and square-free,
gcd(Q,C,D)=1.                                (6)
```

The numerator and denominator are coprime and algebraically independent.
The infinity member is `D`, so it is reduced. For finite `lambda`, put

```
G_lambda=gcd(Q,C-lambda D), chosen monic.
```

Then

```
Q V+C-lambda D=G_lambda * L_lambda(T,V),      (7)
```

where `L_lambda` is primitive and linear in `V`, hence irreducible, and
coprime to `G_lambda`. Thus a finite member is nonreduced exactly when
`G_lambda` is non-square-free. Every finite member has the irreducible
factor `L_lambda` with multiplicity one. Its maximum multiplicity is
the largest exponent in `G_lambda`, or one if that gcd is square-free.
This counts the nonreduced part only, not distinct factors introduced
by Keller pullback.

Here is a finite symbolic description of all exceptional values. Define

```
E={ C(c)/D(c) : ord_c(Q)>=2, D(c)!=0,
                 (C'D-CD')(c)=0 }.          (8)
```

This is exactly the finite nonreduced-value set. The vanishing conditions
say that `C-lambda D` has order at least two at a multiple zero of `Q`.
If `D(c)=0` and `Q(c)=0`, coprimality forces `C(c)!=0`, so that point
does not give a finite exceptional value. In the square-free `D` setting,
the candidate zeros can also be found from
`gcd(Q,Q',C'D-CD')`; no such candidate lies at `D=0`.

Compute the exact maxima `M_lambda` from (7). If

```
sum_(lambda in E)(1-1/M_lambda)<2,           (9)
```

then (6) remains relatively algebraically closed after every Keller
pullback, with every finite parameter-cover group preserved. This
supplies an all-degree coefficient certificate, not a search over
possible composition degrees. In particular it holds with at most two
exceptional values of any multiplicities, or at most three values each
having maximum two. The equality lists (4) are not asserted closed.

The earlier one-exception tests remain useful for nonsquare-free
denominators and other cases outside (6). Square-freeness and
nonconstancy of `D` are hypotheses of this coefficient corollary;
neither may be silently discarded.

## 3. Outputs obtained by expansion at critical points (HC4MYIH2)

Let `R in C[Y]` have degree `m>=2`. Choose nonconstant square-free
`H,D in C[T]` with `gcd(H,D)=1`, and polynomials `beta,B`, such that

```
H divides R'(beta(T)).
```

Let

```
S={ R(beta(c)) : H(c)=0 }.
```

If `|S|<=3`, the polynomial output

```
P(T,z)=D(T)[R(beta(T)+H(T)z)-R(beta(T))]/H(T)^2+B(T)     (10)
```

stays irreducible in `P(f_1,Z)+f_2` under every complex Keller base
change in any dimension `n>=2`, and its full target splitting group is
preserved. For a plane Keller pair every carried closed stable embedding
through `(T,xy-P)`, `r>=1`, therefore forces an inverse. All polynomial,
source and embedding degrees and all numbers of leading levels are
unrestricted, subject to the displayed certificate. Its leading
coefficient is `r_m D H^(m-2)`; for `m>=4` the `H` levels are repeated.

First (10) is polynomial: its linear Taylor term is `H R'(beta)z`,
divisible by `H^2`, and every higher Taylor term contains `H^2`.
Set

```
q=H(T)^2(V+B(T))/D(T)-R(beta(T)).             (11)
```

Its reduced numerator and denominator satisfy (6), with `Q=H^2` and
`C=H^2 B-D R(beta)`. At a zero `c` of `H`, the parameter value is
`-R(beta(c))`. Its derivative in `T` vanishes there, since
`R'(beta(c))=0`. Thus (8) gives exactly `E=-S`. For each such value,
the content in (7) has order exactly two at its corresponding `H` zeros:
`Q` has order two and `C-lambda D` has at least that order. Therefore

```
M_lambda=2 for every lambda in E,
sigma=|S|/2<=3/2<2.                          (12)
```

The parameter (11) is regular after Keller pullback by Section 2. The
root equation multiplied by `H(f_1)^2/D(f_1)` becomes

```
R(beta(f_1)+H(f_1)Z)+q(K)=0.
```

The equation `R(Y)+q` is irreducible over `C(q)` and all its finite
splitting monodromy is preserved, by the retained IG1 field argument.
The affine change of root variable is invertible over the source field.
Since the full target field is purely transcendental over `C(q)`, this
also preserves the full target splitting group. The plane embedding
consequence is exactly the retained
[HC4MYHZ2 filtered output proof](LEADING_FIBRE_UNIT_GATE.md), with all
affine height-one fibres still included.

### 3.1 Unbounded families and a repeated-leading quartic control

For every `m>=3`, take

```
R(Y)=Y^(m-1)(Y-1), alpha=(m-1)/m,
beta(T)=T, H=T(T-alpha), D=T-1.
```

Here `R'(T)=m T^(m-2)(T-alpha)` is divisible by `H`, `H,D` are
square-free and coprime, and `S={0,-alpha^(m-1)/m}` has two values.
Thus (10) closes this family in every `m`, with arbitrary polynomial
`B(T)`, and leading coefficient `(T-1)[T(T-alpha)]^(m-2)` of unbounded
multiplicity. This family has all its Taylor coefficients; it is not
limited to two positive-degree terms.

For a small exact control take

```
R(Y)=Y^4+Y, beta=T, H=4T^3+1, D=T, B=0.
P=T[H^2 z^4+4T H z^3+6T^2 z^2+z].           (13)
```

The three zeros `c` of `H` have distinct critical values `3c/4`, so
`sigma=3/2`. The leading coefficient `a=T H^2` shares its multiple
levels with `b_3=4T^2 H`, failing the older leading triple gcd. The
new parameter certificate proves full irreducibility anyway, without
bounding any other degrees or classifying its source fibres.

## 4. Exact boundary and replay

The bound (1) requires an actual multiplicity-one factor in each
exceptional projective member and no omitted parameter point. Algebraic
independence of `N,M` supplies the latter. Maximum multiplicities cannot
be replaced by gcds, and a denominator fibre cannot be omitted from
the sum. At equality the four numerical lists (4) and congruence (3)
are necessary patterns, not a proof of existence or exclusion. Example
(5) shows that composition can occur there.

Coefficient families whose certificate has sum at least two, projective
fibres without a simple factor, dependent numerator-denominator pencils
not covered by IG0, general coefficient variation, uncovered output
geometries, singleton-zero and unit-ideal plane endpoints and unrestricted
HC4/JC2 remain open. No remaining coefficient pattern is asserted a
Keller realization, and no nonlinear source Hessian invariance is used.

Run [the exact replay](verify_pencil_multiplicities.py):

```sh
timeout 60s .venv/bin/python research/jacobian-followup/hc4-triangular-graphs/verify_pencil_multiplicities.py
```

The [receipt](receipt-pencil-multiplicities.json) checks finite identities
and countercontrols accompanying the all-degree arguments above. The
[validation](validation-pencil-multiplicities.json) separates scoped
validation from unrelated whole-repository failures. Formal verification,
independent replay, external review and novelty priority are not claimed.
