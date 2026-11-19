# Primitive norm sumsets: local criteria and exact higher sums

Let `P_Q = {Q(x,y)>0 : x,y in Z, gcd(x,y)=1}`. A target `N` is
**pairwise locally admissible** for `Q_1+Q_2` if, at every prime `p`, it
has a solution over `Z_p` with neither coordinate pair in `p Z_p^2`.
Primitivity is imposed separately in the two blocks. This note proves
finite local criteria, an effective universal four-summand bound, and
exact stabilization of every higher sumset for principal norms in all
imaginary quadratic orders. The higher-order exception lists have proved
finite cutoffs; the binary conjectures remain separate.

## 1. A finite local decision theorem, including conductor primes

Write `Q_i=a_i x^2+b_i xy+c_i y^2`, with nonzero discriminant `D_i`, and
put

```
B_i = [[2a_i,b_i],[b_i,2c_i]],    det(B_i)=-D_i,
t_p = max(v_p(D_1),v_p(D_2)),    h_p=2t_p+1.
```

**Theorem PL1.** Pairwise primitive `Z_p` solubility is equivalent to a
pairwise primitive solution modulo `p^h_p`. For odd `p` not dividing
`D_1 D_2`, every target is soluble. Consequently pairwise local
admissibility is determined by a finite residue test at the primes
dividing `2D_1D_2`, even for nonmaximal orders. This is an algorithmic
complete classification; the modulus is sufficient and need not be minimal.

**Proof.** For a primitive block `v`, the identity
`adj(B_i) B_i v = det(B_i) v` shows that some entry of its gradient has
valuation at most `v_p(D_i)`. Given a primitive root modulo `p^h_p`, choose
integer lifts and a coordinate whose derivative has valuation `s<=t_p`.
Keep the other three coordinates fixed. The resulting univariate
polynomial `f` satisfies `v_p(f(z)-N)>=2t_p+1>2s`. The elementary strong
form of Hensel's lemma gives a root congruent to `z` modulo `p^(s+1)`.
Both blocks remain primitive. Necessity is reduction modulo `p^h_p`.
The assertion at the remaining odd primes follows from the positive
root counts in PL2 below. Only the real condition remains outside the
finite primes; positive definite forms represent every positive real
target with both blocks nonzero. QED.

Here the Hensel step can be obtained directly by Newton iteration:
`z_new=z-(f(z)-N)/f'(z)`. The correction has valuation at least `s+1`,
so the derivative keeps valuation `s`, and the residual valuation strictly
improves. This also explains why primitivity survives.

Define the normalized coordinate density

```
beta_p(N)=lim_{h->infinity} p^(-3h)
  #{v mod p^h : Q_1(v_1)+Q_2(v_2)=N, v_1,v_2 nonzero mod p}.
```

**Theorem PL2.** The density stabilizes at the explicit finite level
`h_p=2t_p+1`: its value equals the root count modulo `p^h_p` divided
by `p^(3h_p)`, and normalized counts agree at every higher level.
Consequently the density depends only on `N mod p^h_p`. On locally
admissible targets, `beta_p(N)>=p^(-2t_p-3)`, uniformly in `N`.
At an odd prime not dividing
`D_1D_2`, put `epsilon_i=(D_i/p)`. Then

```
p does not divide N:
 beta_p(N)=1-(2+epsilon_1 epsilon_2)/p^2
                +(epsilon_1+epsilon_2)/p^3;
p divides N:
 beta_p(N)=1+epsilon_1 epsilon_2/p
              -(2+epsilon_1+epsilon_2+epsilon_1 epsilon_2)/p^2
              +(1+epsilon_1+epsilon_2)/p^3.
```

In particular the good-prime product is positive and is
`>=c_epsilon N^(-epsilon)` for every fixed positive `epsilon`.

**Proof.** Put `t=t_p` and `r=t+1`. Partition the pair-primitive domain
into balls `v_0+p^r Z_p^4`. On each ball the minimum gradient valuation
`s=min_j v_p(partial_j F(v_0))<=t` is fixed. Taylor expansion gives

```
F(v_0+p^r u)=F(v_0)+p^(r+s) [L(u)+p^(r-s) F(u)],
```

where `L` is an integral linear form with at least one unit coefficient.
Because `r-s>=1`, fixing the other three variables makes the bracket a
bijection of `Z_p` in the variable with that unit coefficient. Its
reduction modulo `p` is affine with unit slope, and each root lifts
uniquely. Hence normalized Haar measure on the ball pushes forward to
uniform measure on `F(v_0)+p^(r+s) Z_p`. The ball has measure `p^(-4r)`
and the image coset has measure `p^(-r-s)`, so its density contribution
is `p^(-3r+s)` on that coset and zero elsewhere. Since
`r+s<=2t+1`, this proves exact stabilization at `h_p`, as well as
existence and local constancy of the density.

For the lower bound start at an exact root and choose its minimum
gradient valuation `s<=t`. The same argument on the smaller-radius
ball with `r=s+1` gives the positive contribution
`p^(-3(s+1)+s)=p^(-2s-3)>=p^(-2t-3)`.

At a good odd prime the number of primitive binary roots of `Q_i=t`
modulo `p` is `(p-1)(1+epsilon_i)` for `t=0`, and `p-epsilon_i` for
`t!=0`. To obtain these formulas, diagonalize the binary form over `F_p`
and use `sum_x (x^2-a / p)=-1` for `a!=0`; the latter identity counts
`u^2-v^2=a` by its `p-1` ordered nonzero factor pairs. Convolving the
two binary counts gives the displayed formulas. Every primitive root
has a nonzero derivative, so it has exactly `p^(3(h-1))` lifts modulo
`p^h`. The formulas are positive for all four sign choices and `p>=3`.
Away from target primes the factors are `1+O(p^-2)`; at target primes
they are at least `1-C/p`. Discard finitely many small primes and use
`1-C/p>=p^(-epsilon)` for sufficiently large `p` to obtain the stated
product bound. QED.

These are coordinate densities. They are not leading constants for
unweighted counts of distinct norm values.

## 2. The maximal-order valuation rules

For an imaginary quadratic maximal order `O_K`, an element is primitive
at `p` exactly when it is outside `p O_K`. Unique ideal factorization
then gives the following necessary valuation rules for its positive norm:

| Rational prime | Allowed valuation in a primitive norm |
|---|---|
| split | any nonnegative integer, using only one prime above `p` |
| inert | zero |
| ramified | zero or one |

Locally the split rule permits every `p`-adic norm value. At an inert
prime the values are exactly `Z_p^*`: the residue-field norm is
surjective, and a nonzero residue root lifts by Hensel. At an odd
ramified prime write a local integral basis as `1,pi`, where
`pi^2=p u` for a unit `u`. A primitive norm has either valuation zero,
with square unit residue, or valuation one, with unit residue in
`-u` times the squares. These conditions are sufficient: fix one
coordinate and lift the other from its nonzero square root. Thus the
ramified valuation rule alone omits a unit squareclass condition.

The maximal-order bound must not be imposed at conductor primes.
For example `x^2+4y^2` is the norm in the conductor-two Gaussian order,
and `(x,y)=(0,1)` and `(2,1)` are primitive with norms `4` and `8`.
The ramified valuation can therefore be two or three. PL1 includes
these cases without assuming that the localized order is maximal.

## 3. Explicit self-sum local classifications

**Corollary PL3.** The exact all-prime pairwise local conditions are:

| Form | Local condition on a positive target `N` |
|---|---|
| `x^2+2y^2`, discriminant `-8` | `16` does not divide `N` |
| `x^2+xy+2y^2`, discriminant `-7` | `49` does not divide `N` |
| `x^2+xy+((p+1)/4)y^2`, `p=11,19,43,67,163` | `N` even and `p^2` does not divide `N` |
| `x^2+5y^2`, discriminant `-20` | `N mod 8` lies in `{2,3,4,6,7}` |

**Proof at odd ramified primes.** For `p=3 mod 4`, put `a=2x+y`,
`b=y`. This is an invertible local change at `p`, and
`4Q=a^2+p b^2`. If `p` divides a self-sum, anisotropy of
`a_1^2+a_2^2` modulo `p` forces both `a_i=0 mod p`.
Primitivity forces both `b_i` to be units. If `p^2` also divided the
target, the same anisotropy would force both `b_i=0 mod p`, a
contradiction. Conversely, a unit target has a representation by the
first two squares with nonzero gradient; choose primitive second
coordinates where needed and lift. For a target of valuation one,
set both `a_i=0` and represent `4N/p` as `b_1^2+b_2^2` with both
`b_i` units. For each nonzero residue `t`, the number of such ordered
pairs is `p-1-2(t/p)>0` when `p>=7`. Hensel gives the required exact
solution. This explains the lower bound on `p`: at `p=3`, a target
divisible by three must be `6 mod 9`, not merely have valuation one.

For `x^2+5y^2` at five, use that `-1` is a square in `Z_5`.
The two `x`-squares can represent every `5`-adic target; if an `x`
coordinate vanishes modulo five, choose its `y` coordinate a unit.
The other nonzero `x` derivative lifts; at a target divisible by five
choose both `x` coordinates nonzero using isotropy. This establishes
primitive local solubility for every target at five.

**Proof at two.** For an odd fundamental discriminant, two is split
when `D=1 mod 8` and inert when `D=5 mod 8`. For the displayed forms,
this follows directly from `Q=x^2+xy+c y^2`: the primitive norm values
are all `Z_2` when `c` is even, and exactly the odd units when `c` is
odd. Each corresponding root modulo two has an odd partial derivative
and hence lifts. Two odd norms sum to every even target.

For `x^2+2y^2`, primitive values modulo sixteen are
`{1,2,3,6,9,11}`. Their two-fold sumset is precisely the fifteen
nonzero residues modulo sixteen. Every such root lifts to every
higher power: if an `x` coordinate is odd, changing it by
`2^(k-1)` toggles the `2^k` bit for `k>=4`; otherwise a `y`
coordinate is odd, and changing it by `2^(k-2)` changes its norm by
`2^k y+2^(2k-3)`, toggling exactly that bit. Primitivity is preserved.
For `x^2+5y^2`, primitive values modulo eight are `{1,5,6}` and
their sumset is `{2,3,4,6,7}`. Every root lifts by changing an odd
coordinate by `2^(k-1)` for `k>=3`, since both diagonal coefficients
are odd. All remaining odd primes are covered by PL2. QED.

### The complete maximal-order self-sum criterion

**Corollary PL3All.** Let \(D<0\) be a fundamental discriminant
and \(Q_D\) the principal integral-basis norm of its maximal order.
A positive target \(N\) is a sum of two primitive \(Q_D\)-values
over every \(\mathbb Z_p\) if and only if all applicable conditions
below hold:

* For every odd \(p\mid D\) with \(p\ge7\) and
  \(p\equiv3\pmod4\), require \(p^2\nmid N\).
* If \(3\mid D\), put \(\kappa=-D/3\bmod3\), a nonzero residue.
  Require \(3\nmid N\) or \(N\equiv6\kappa\pmod9\).
* At two, impose the condition in this table.

| Fundamental discriminant | Exact condition at two |
|---|---|
| \(D\) odd, \(D\equiv1\pmod8\) | none |
| \(D\) odd, \(D\equiv5\pmod8\) | \(N\) even |
| \(D=-4m\), \(m\) odd and squarefree, \(m\equiv1\pmod4\) | \(N\bmod8\in\{2,3,4,6,7\}\) |
| \(D=-8k\), \(k\) odd and squarefree, \(k\equiv1\pmod4\) | \(16\nmid N\) |
| \(D=-8k\), \(k\) odd and squarefree, \(k\equiv3\pmod4\) | none |

Odd ramified primes \(p\equiv1\pmod4\) impose no sum condition.
There are no other local conditions. In particular the table covers
every allowed higher target valuation, without a valuation cutoff.

A sufficient exact admissibility modulus is the product of the
dyadic modulus \(1,2,8\) or \(16\) in the table, the factor \(9\)
if \(3\mid D\), and \(p^2\) for each ramified
\(p\ge7\), \(p\equiv3\pmod4\). This modulus need not be minimal.

**Proof at odd primes.** Write \(Q=x^2+bxy+cy^2\), with
\(b=0\) or \(1\). At an odd ramified \(p\), put
\(\kappa_p=-D/p\in\mathbb Z_p^\times\). Then
\[
                 4Q=(2x+by)^2+p\kappa_p y^2.
\]
The exact primitive norm image consists of unit values with square
residue, and valuation-one values whose unit quotient is in
\(\kappa_p(\mathbb F_p^\times)^2\). Necessity follows from the
displayed equation and \(v_p(D)=1\). Sufficiency follows from
\[
             Q(t,0)=t^2,\qquad Q(-bt,2t)=p\kappa_p t^2,
\]
and square-root lifting for units. These are full \(p\)-adic value
sets, rather than necessary residue restrictions.

Put \(R=(\mathbb F_p^\times)^2\), and write \(\chi_p\) for the
Legendre symbol. For \(t\ne0\),
\[
 C_p(t):=\#\{(a,c)\in R^2:a+c=t\}
            =\frac{p-2-\chi_p(-1)-2\chi_p(t)}4.             \tag{PL3A.1}
\]
Indeed \(X^2+Y^2=t\) has \(p-\chi_p(-1)\) roots: the identity
\((X-Y)(X+Y)=t\) gives
\(\sum_X\chi_p(X^2-t)=-1\). Remove its four axis roots when
\(t\) is a square, then divide by four.

Every unit target is a sum of two primitive norms. If its residue
is square, use the ramified value \(p\kappa_p\) and the unit
square \(N-p\kappa_p\). If it is nonsquare, (PL3A.1) is positive;
choose two unit-square residues, fix the first unit norm, and
lift the second square root. For \(p\equiv1\pmod4\) and
\(p\mid N\), the unit values \(1\) and \(N-1\) suffice,
since \(-1\) is square modulo \(p\). This includes every
target valuation.

For \(p\equiv3\pmod4\), two unit norms cannot sum to a
multiple of \(p\). Such a sum therefore has two valuation-one
summands. Their quotients lie in \(\kappa_p R\) and cannot
cancel modulo \(p\), excluding \(p^2\mid N\). If \(N=pM\)
has valuation one, solubility is equivalent to
\(C_p(M/\kappa_p)>0\). For \(p\ge7\) this holds for every
nonzero residue, since its minimum is \((p-3)/4>0\).
For \(p=3\), the only element of \(R\) is \(1\), so the
condition is \(M/\kappa_p\equiv2\pmod3\), precisely
\(N\equiv6\kappa_p\pmod9\). The exact norm images already
proved supply every higher digit in these cases.

At unramified odd primes, a split local algebra permits every
primitive norm value, using the two factors \(1,n\). An inert
local algebra permits exactly the units; surjectivity of the
finite-field norm and nonsingular lifting prove sufficiency.
Two units sum to any target when \(p\ge3\), by choosing a
nonzero residue different from \(N\bmod p\). Thus these
primes impose no sum condition.

**Proof at two.** If \(D\) is odd, write
\(Q=x^2+xy+c y^2\). For even \(c\), primitive vectors
attain both residues modulo two; for odd \(c\), all three
nonzero vectors have norm one. At a primitive vector the
gradient has an odd component. Consequently every selected
root lifts to every value with its parity. These are,
respectively, all of \(\mathbb Z_2\) when
\(D\equiv1\pmod8\), and exactly its odd units when
\(D\equiv5\pmod8\).

Recall that an odd \(2\)-adic unit is square if and only if it
is \(1\bmod8\). Sufficiency follows directly by lifting an
odd root modulo \(2^h\): changing it by \(2^{h-1}\) toggles
the next square bit, for every \(h\ge3\).

For \(D=-4m\), with \(m\) odd, one necessarily has
\(m\equiv1\pmod4\) and \(Q=x^2+m y^2\). Opposite-parity
coordinates give exactly the full unit classes \(1\bmod4\):
choose \(y=0\) or \(2\) and solve for \(x\) using the
square-unit criterion. Two odd coordinates give exactly the
full class \(m+1\bmod8\): fix \(y=1\) and solve
\(x^2=n-m\). Hence the complete primitive value image is
\[
           \{n:n\equiv1\pmod4\}\ \cup\
           \{n:n\equiv m+1\pmod8\}.
\]
Its residues modulo eight are \(\{1,5,2\}\) or
\(\{1,5,6\}\). In either case their self-sum is exactly
\(\{2,3,4,6,7\}\bmod8\).

For \(D=-8k\), with \(k\) odd, the form is
\(Q=x^2+2k y^2\). Odd norms have \(x\) odd, giving exactly
the full classes \(1,1+2k\bmod8\), by fixing \(y=0\)
or \(1\). For an even primitive norm, write \(x=2z\);
then \(y\) is odd and
\[
                         n/2=2z^2+k y^2.
\]
These quotients are \(k\) or \(k+2\bmod8\). Conversely,
fixing \(z=0\) or \(1\) gives a square-unit equation for
\(y\) in each class. Thus the exact even norm classes are
\(2k,2k+4\bmod16\). The complete residue algebra is:

| \(k\bmod8\) | Primitive norm residues modulo sixteen | Self-sum |
|---|---|---|
| 1 | \(\{1,2,3,6,9,11\}\) | all except 0 |
| 3 | \(\{1,6,7,9,10,15\}\) | all |
| 5 | \(\{1,3,9,10,11,14\}\) | all except 0 |
| 7 | \(\{1,2,7,9,14,15\}\) | all |

Every value class used here contains all its higher \(2\)-adic
lifts. Adding two full residue classes gives a full target
residue class. Thus the finite residue algebra proves the
all-depth assertion; no extrapolation from a finite search is
needed. This completes the local classification. \(\square\)

By [AdditiveClassDisappearance](PRIMITIVE_NORM_ANALYTIC_BRIDGE.md#additive-disappearance-of-class-obstructions),
every sufficiently large target in this explicit local set is
a sum of two positive primitive norms. The exceptional locally
admissible targets form a finite set depending on \(D\); this
corollary does not enumerate it. For nonmaximal orders, PL1--PL2
supply the conductor-sensitive local test instead of PL3All.

**Global counterexamples.** Local admissibility alone is not an exact
all-target theorem. For discriminant `-7`, the target `7` is locally
admissible, but the primitive positive values at most six are exactly
`{1,2,4}` and have no sum seven. For discriminant `-19`, the target
`4` is locally admissible, while the only primitive positive value
below five is one. For `x^2+5y^2`, both three and four are locally
admissible but cannot be sums: its only primitive positive value below
five is one. These facts follow by completing the square and checking
the resulting bounded `y` range. They do not assert infinitely many
exceptions or disprove eventual local-global statements.

### Explicit ramified factors and the discriminants -8 and -11

**Corollary PL3D (explicit local densities).** For the self-sum of
`Q=x^2+xy+((p+1)/4)y^2`, where `p=3 mod 4` is an odd prime,
write `chi_p` for the Legendre symbol. Its ramified coordinate density is

```
v_p(N)=0: beta_p(N)=1+1/p-2(1+chi_p(N))/p^2;
v_p(N)=1: beta_p(N)=(p-1-2chi_p(N/p))/p^2;
v_p(N)>=2: beta_p(N)=0.
```

At two the density is `9/8` for even `N` and zero for odd `N` if
`p=3 mod 8`. If `p=7 mod 8`, it is `5/8` for even `N` and `1/2`
for odd `N`. For `x^2+2y^2` the factor at two is exactly

```
beta_2(N)=0       if N=0 mod16;
          5/4    if N=4 mod8;
          1/2    otherwise.
```

These formulas, together with PL2, give every local factor of the
coordinate-count main term, without an unevaluated limiting density.
For example the factor at eleven is `128/121` or `132/121` on units,
according as the target is a square or nonsquare modulo eleven. At
valuation one it is respectively `8/121` or `12/121`, using the square
class of `N/11`; at deeper valuations it vanishes.

**Proof.** At the odd ramified prime, make the unit-determinant
change `a=2x+y,b=y` and multiply the equation by four. The form becomes
`a_1^2+a_2^2+p(b_1^2+b_2^2)`. For a unit target, the equation
`a_1^2+a_2^2=4N` has `p+1` roots modulo `p`. Each normally allows
`p^2` choices of the two `b` coordinates. When an `a_i` vanishes,
primitivity removes the choice `b_i=0`. There are two such roots for
each vanishing coordinate exactly when `chi_p(N)=1`, and they cannot
both vanish. Thus the root count is
`p^2(p+1)-4p*1_(chi_p(N)=1)`. Every root is nonsingular, proving the
first formula after division by `p^3`.

For valuation one, anisotropy of `a_1^2+a_2^2` forces both `a_i`
into `p Z_p`. Substitute `a_i=p c_i` and divide the equation by `p`.
The two `b_i` must be units and their squared sum is `4N/p` modulo
`p`. There are `p-1-2chi_p(N/p)` such ordered pairs, each nonsingular.
The substitution contributes `p^-2` and dividing the equation
contributes `p`; the resulting density is the displayed quotient by
`p^2`. The same argument excludes higher target valuations.
At two for odd discriminant all primitive roots have a unit gradient.
The primitive binary residue counts are `(0,3)` for `p=3 mod8`, and
`(2,1)` for `p=7 mod8`. Convolving and dividing by `2^3` gives the
stated factors.

For discriminant `-8`, PL2 reduces the exact density to the finite
count modulo `2^7`. The checker enumerates all primitive binary
residues at that modulus, convolves their multiplicities, and verifies
the displayed factor for every one of the 128 target classes. This is
a finite certificate for the unrestricted formula, using the proved
stabilization theorem; the next-level check independently challenges
its normalization. QED.

For the self-sum of any fixed binary norm form of discriminant `D`,
write `chi(q)=(D/q)` at good odd primes. The ratio of the target factor
to the nontarget factor simplifies to

```
(q+3)/(q+2)   if chi(q)=1;
(q-1)/(q-2)   if chi(q)=-1.
```

Indeed the two factors for `chi=1` are
`(1-1/q)^2(1+3/q)` and `(1-1/q)^2(1+2/q)`; for `chi=-1` they are
`(1+1/q)^2(1-1/q)` and `(1+1/q)^2(1-2/q)`. The target ratios follow solely from the local density formulas.

## 4. Exact ideal-class criterion for individual primitive norms

Let `O_K` be maximal. A positive integer `n` is a primitive principal
norm if and only if:

1. No inert prime divides `n`, and every ramified prime has exponent
   at most one.
2. For each split prime power `p^e || n`, one can choose one of the two
   conjugate prime ideals `pfrak,pfrak_bar` so that the product of the
   chosen `pfrak^e` and the required ramified prime ideals has class one.

Indeed a primitive principal ideal contains at most one prime above
any split rational prime, contains no inert prime, and has ramified
exponents at most one. Conversely such a product ideal of class one
is generated by an integral element of norm `n`; its ideal
factorization prevents every rational prime from dividing its two
coordinates. This is an exact ideal factorization criterion, not a
splitting-semigroup criterion. Orientation is chosen once for each
whole prime power; choosing orientations independently for each copy
would incorrectly permit nonprimitive elements.

For a class group of exponent two, orientations have the same class,
so the condition is simply a parity condition in that finite group.
For `Q(sqrt(-5))` the class group has order two (also [Milne, version 3.08, Example 4.6](https://www.jmilne.org/math/CourseNotes/ANT.pdf#page=73)): Minkowski's bound is
`(2/pi)sqrt(20)<3`, so classes are generated by primes over two;
that ramified prime is not principal (there is no element of norm two)
and has square `(2)`. Thus the class condition here is of genus type;
one must not describe `x^2+5y^2` as an example of a hidden odd-order
class obstruction.

**Corollary PL4 (classical class-two specialization).** An individual
positive integer is a primitive `x^2+5y^2` norm if and only if it is
primitively represented over every `Z_p`. More explicitly, this is
equivalent to: no inert prime divides `n`, `v_2(n),v_5(n)<=1`, and
`n=1 mod 4` if odd, or `n=6 mod 8` if even.

**Proof.** Write `R` for the nonprincipal prime above two. A split odd
prime ideal `P` of norm `p` has class one only when `p=1 mod 4`, since
an odd value of `x^2+5y^2` is `1 mod 4`. If `P` has class `R`, then
`R P` is principal of norm `2p`. Its generator has both coordinates
odd, so `2p=6 mod 8` and `p=3 mod 4`. Since there are only two
classes, these necessary alternatives classify all split odd prime
ideals. The prime above five is principal. Hence the oriented ideal
product for `n=2^e m` with `m` odd has class one exactly when
`m=(-1)^e mod 4`. This is precisely the displayed dyadic condition.
Local primitive representation gives the inert and ramified valuation
rules and that dyadic condition. The ideal-class criterion then supplies
a global primitive generator. The converse is immediate by localization.
This proves sufficiency without assuming a general integral Hasse
principle for binary forms. QED.

The class-character test here is quadratic and of genus type. Its
sufficiency does not extend to arbitrary class groups.

A direct example of a genuine individual-norm local-global failure is
`Q=x^2+xy+6y^2` (discriminant `-23`), `n=2`. At two, two splits;
at twenty-three, two is a nonzero square; at every other prime it is
a unit norm or a split norm. Hence it is primitively represented
locally everywhere. Globally
`Q=(x+y/2)^2+23y^2/4`, so `|y|>=1` gives `Q>=23/4>2`, while `y=0`
only gives squares. No global representation exists. An individual
norm counterexample does not by itself produce arbitrarily large
exceptions in a binary additive sumset.

## 5. Effective additive stability in every quadratic order

The following higher-summand results are elementary corollaries of the
classical three- and four-square theorems. They do not require the
quaternary asymptotic, class number one, or maximality of the order.
They apply to **principal** norm forms; an arbitrary ideal-class form
with leading coefficient greater than one is outside this statement.

Use the integral basis of a quadratic order to write its principal norm
as `Q=x^2+bxy+c y^2`, with `b in {0,1}` and `c>=1`. Put `P=P_Q`.
The classical inputs are Lagrange's four-square theorem and the
Legendre--Gauss three-square theorem, as stated on PDF pages 1 and 3 of
[Goldmakher--Pollack, *Refinements of Lagrange's Four-Square Theorem*](https://www.pollack-math.net/4squares-arxiv.pdf).
Gauss's three-triangular-number theorem follows from the latter: represent
`8m+3` by three squares, all of which must be odd, then divide by eight.

**Theorem PL5 (effective universal bound and exact stabilization).**

* If `b=0`, then `4P` contains every `N>=4c`.
* If `b=1`, then `3P` contains every `N>=3c` with `N=c mod 2`.
  If `c` is odd, every primitive norm is odd. In that case `4P`
  contains every even `N>=4c`. If `c` is even, `4P` contains every
  `N>=4c`, of either parity.
* Consequently every `N>=4c` is a sum of at most four elements of `P`.
* If `b=1` and `c` is odd, then for **every** `k>=3` there is the exact
  set identity `kP=(k-3)+3P`, where the right side denotes translation.
* In all other cases, for every `k>=4`, there is the exact identity
  `kP=(k-4)+4P`.

Thus every sumset in the last two assertions is determined by a finite
exhaustive calculation: enumerate `3P` below `3c` in the odd-norm case,
or `4P` below `4c` otherwise. This gives a finite classification algorithm
for all higher sumsets in every imaginary quadratic order. The bound four
for at-most-four sums is sharp over these orders.

**Proof.** If `b=0`, write `N-4c=x_1^2+...+x_4^2` by Lagrange and use
`N=sum_i Q(x_i,1)`. Each pair has a coordinate equal to one.

If `b=1`, `Q(x,1)=2T_x+c`, with `T_x=x(x+1)/2`. Three triangular
numbers represent `(N-3c)/2`, giving the asserted `3P` coverage. Adding
`Q(0,1)=c` gives `4P` at parity `N=0 mod 2` when `c` is odd, and at
even targets when `c` is even. In the latter case, adding
`Q(1,0)=1` instead gives every odd `N>=3c+1`; this includes every
odd `N>=4c`. When `c` is odd, reducing `Q` modulo two shows that a
primitive pair always has odd norm.

Every value of `P` other than one is at least `c`. For `b=0` this is
immediate unless `y=0`, in which case primitivity forces `|x|=1`.
For `b=1`, complete the square:
`Q=(x+y/2)^2+(c-1/4)y^2`. At `|y|=1` the minimum integer value is
`c`; at `|y|>=2` the lower bound is `4c-1>=c`.

In the odd-norm case, let `k>=3`. Adding `k-3` ones proves
`(k-3)+3P` is contained in `kP`. Conversely, a represented target has
parity `k`. If `N>=3c+k-3`, the translated target is covered by the
three-triangular construction. If `N<3c+k-3=k+3(c-1)`, a
`k`-summand representation has at most two summands different from one,
since three such summands would already give `N>=3c+k-3`.
Removing `k-3` ones leaves an element of `3P`. For `c=1`, the lower
range has `N<k` and contains no representation, so the same conclusion
holds. This proves the first translation identity. The proof for
`k>=4` in the remaining cases is identical: above `4c+k-4` use the
four-summand coverage; below it there are at most three nonunit
summands, so enough ones can be removed.

For sharpness take `Q=x^2+8y^2`, a principal norm form of discriminant
`-32`. Its primitive norm residues modulo eight are `{0,1,4}`.
At most three such residues cannot sum to seven modulo eight.
Thus infinitely many targets need four summands, while the construction
above covers every `N>=32` with exactly four. QED.

## 6. Exact higher sumsets in the nine class-number-one fields

**Corollary PL6.** For the Gaussian norm, `3P_-4` is exactly
`{N>=3 : N!=2 mod 8}`, and `kP_-4={N>=k}` for all `k>=4`.
For discriminants `-7` and `-8`, `kP_D={N>=k}` for every `k>=3`.

For `D=1-4c` with `c=1,3,5,11,17,41`, let `E_D` be the set in the
following table. Then for every `k>=3`, exactly,

```
kP_D = {N>=k : N=k mod 2, N not in E_D+(k-3)}.
```

| `D` | `c` | Complete exceptional set `E_D` for three summands |
|---|---:|---|
| `-3` | 1 | empty |
| `-11` | 3 | empty |
| `-19` | 5 | `{5}` |
| `-43` | 11 | `{5,7,9,11,17,21}` |
| `-67` | 17 | `{5,7,9,11,13,15,17,23,27,29,33,45}` |
| `-163` | 41 | `{5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,47,51,53,57,59,61,65,67,69,71,75,77,79,81,93,111,117,121}` |

In particular the four-summand exception set in each odd-norm row is
exactly `E_D+1`, among even targets `N>=4`. The binary classifications
for -3, -7, -8 and -11 are also [complete](PRIMITIVE_NORM_EXACT_CLASSIFICATIONS.md).
The binary cases -19, -43, -67 and -163 remain open.

**Proof for `D=-7`.** The three-triangular construction with `c=2`
covers even `N>=6`. Every nonnegative `m` is `T_a+T_b+2z^2`:
represent `8m+2` by three squares. There are two odd square roots and
the remaining root is divisible by four; writing them as
`2a+1,2b+1,4z` gives the identity. Therefore every odd `N>=11` is

```
N=11+2m=Q(a,1)+Q(b,1)+Q(2z-1,2).
```

All three pairs are primitive. The remaining targets from three through
ten are sums of three members of `{1,2,4,7}`: respectively
`1+1+1, 1+1+2, 1+2+2, 2+2+2, 1+2+4, 2+2+4, 1+1+7, 2+4+4`.
Adding ones supplies all higher `k`.

**Proof for `D=-8`.** For every positive odd `m`, the form
`A=x^2+u^2+2v^2` represents `m`: represent `2m` by three squares,
choose its even coordinate as `2v`, and use half the sum and half the
difference of the other two coordinates. Also
`B=x^2+2u^2+2v^2` represents every positive odd `m!=7 mod 8`:
represent `m` by three squares, choose an odd coordinate as `x`, and
take half the sum and difference of the remaining coordinates, which
have the same parity.

For even `N>=6`, apply `A` to `N-5` and write
`N=Q(x,1)+Q(u,1)+Q(1,v)`. For odd `N>=5` with `N!=3 mod 8`,
apply `B` to `N-4` and write
`N=Q(x,1)+Q(1,u)+Q(1,v)`. If `N=3 mod 8` and `N>=11`,
then `N-6=5 mod 8` is a sum of three squares, so use
`N=sum_i Q(x_i,1)`. Finally three and four are `1+1+1` and
`1+1+2`. All pairs in these formulas contain a coordinate one.

**Proof for `D=-4`.** Reduction modulo eight excludes `N=2 mod 8`
from `3P`, since the primitive norm residues are `{1,2,5}`.
If `N>=3` lies in neither `2 mod 8` nor `3 mod 4`, then `N-3`
passes the three-square criterion, so `N=sum_i Q(x_i,1)`.
Every nonnegative `m` also has the form `a^2+b^2+z(z+1)`:
apply the three-square theorem to `4m+1`, whose representations have
one odd coordinate and two even coordinates. For `N=3 mod 4`,
`N>=7`, set `m=(N-7)/4` and use

```
N=Q(2a,1)+Q(2b,1)+Q(2z+1,2).
```

These pairs are primitive. The target three is `1+1+1`.
Lagrange with `c=1` gives `4P={N>=4}`; adding ones proves the
remaining assertions.

**Proof of the finite table.** PL5 covers every odd `N>=3c` and
excludes every even target. The checker exhaustively enumerates all
primitive positive norm values below `3c`, convolves their set three
times, and records exactly the omitted odd targets in `[3,3c)`.
Completing the square bounds `|y|<=sqrt(4(3c-1)/|D|)`; for each such
`y`, `|x|<=sqrt(3c-1)+|y|/2` suffices. Thus the finite enumeration
has a proved completeness bound. The largest cutoff is only 123.
The exact translation identity PL5 supplies every `k>=3`, without
any additional scan or analytic input.

## 7. Minimal asymptotic additive orders

Distinguish the least `h` for which `hP` contains every sufficiently
large integer from the least `h` for which `h(P union {0})` does so.
The latter permits at most `h` positive norm summands.

**Corollary PL7.** The minimal orders in the class-number-one fields are:

| Discriminant | Exactly `h` positive summands | At most `h` positive summands |
|---|---:|---:|
| `-4` | `4` | `3` |
| `-7,-8` | `3` | `3` |
| `-3,-11,-19,-43,-67,-163` | no such `h` | `4` |

All these assertions follow from elementary results above and local
obstructions, except the Gaussian at-most-three upper bound. That entry
also uses the eventual two-summand theorem from
[PRIMITIVE_NORM_ANALYTIC_BRIDGE.md](PRIMITIVE_NORM_ANALYTIC_BRIDGE.md)
on targets `N=2 mod 8`, which `3P_-4` excludes. No exact two-summand
classification is being asserted here.

**Proof.** PL6 gives the exact-three and exact-four upper bounds and
parity exclusions. For Gaussian at-most-three coverage, exact `3P`
covers every residue except two modulo eight, and the eventual
pair-primitive theorem supplies the remaining class. For minimality,
multiples of eight belong to neither Gaussian `P` nor `2P`.
Multiples of sixteen, or of forty-nine, belong to neither set for
`D=-8,-7`. These rule out at most two summands. Gaussian `3P`
misses the class two modulo eight, proving its exact-four order is
minimal. In the odd-norm fields, infinitely many even targets are
excluded from `2P` by the ramified local condition: use multiples of
`p^2` for `D=-p`, and multiples of nine for `D=-3`. They lie in
neither `P` nor `3P`, by parity. Thus at most three summands cannot
suffice, while PL5 gives the at-most-four bound. QED.

## 8. Completely specified coordinate main terms

**Corollary PL8.** For the self-sum of a fixed primitive positive binary
norm form of discriminant `D`, combine PL2 and PL3D with the
[quaternary theorem](PRIMITIVE_NORM_ANALYTIC_BRIDGE.md). Put
`chi(q)=(D/q)` for odd primes not dividing `D`. Then

```
R_D^*(N)=C_D N B_D(N) product_(q|N, q good) T_q
                                      + O_D(N^(1-1/1000)),
C_D=(4 pi^2/|D|) product_(q odd, q not dividing D)
                         (1-3/q^2+2chi(q)/q^3),
B_D(N)=product_(q dividing 2D) beta_q(N),
T_q=(q+3)/(q+2) if chi(q)=1, and (q-1)/(q-2) if chi(q)=-1.
```

The constant product converges absolutely and positively. For `D=-8`
and `D=-11`, every factor of `B_D` is given explicitly in PL3D. 

**Proof.** At a good prime the target and nontarget factors for `chi=1`
are `(1-1/q)^2(1+3/q)` and `(1-1/q)^2(1+2/q)`; for `chi=-1` they are
`(1+1/q)^2(1-1/q)` and `(1+1/q)^2(1-2/q)`. Their ratios are precisely
`T_q`. The real factor for a self-sum is `4 pi^2/|D|`, since the binary
Gram determinant is `|D|`. Substitution in the general coordinate
asymptotic proves the formula. QED.

This counts integral coordinate representations. It is not an
unweighted value-pair asymptotic or a numerical effective cutoff.

## Verification and scope

[verify_primitive_norm_local_structure.py](verify_primitive_norm_local_structure.py)
checks the good-prime density formulas by exhaustive finite-field
counting, the displayed dyadic and additive-order residue tables,
representative ramified residue conditions, exact density stabilization
at the PL2 level for four forms (including a conductor prime), the
explicit factors of PL3D, and the
finite global counterexamples. Its JSON
receipt is
[primitive_norm_local_structure_checks.json](primitive_norm_local_structure_checks.json).
The unrestricted local claims depend on the lifting proofs above;
finite tests do not replace them. No independent review or historical priority claim is asserted. Global
analytic assertions here are explicitly cited corollaries; finite tests
do not establish their analytic inputs.

The ideal-form correspondence and background are also developed in
[Milne, Algebraic Number Theory, version 3.08, Theorem 4.29, printed page 82](https://www.jmilne.org/math/CourseNotes/ANT.pdf#page=84).
The explicit arguments above specify the primitive restrictions needed
here, rather than importing an unrestricted representation statement.
