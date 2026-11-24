# Exact row-merger formulas and a positive-measure obstruction

This note distinguishes the actual factorial row multiplier from positive
averaging, and gives an exact failure of regularity preservation in a larger
positive-measure class. It does not give an ordinary factorial counterexample.

## 1. Beta, Jacobi, and confluent divided-difference identities

For two nonnegative integer support rows `u,v`, let
`U=sum_j u_j c_j d/dc_j`, `V=sum_j v_j c_j d/dc_j` be the commuting Euler
operators. On `c^k`, write `P=u.k`, `Q=v.k`. The exact row-merger operator is

\[
              \mathcal M(c^k)=\binom{P+Q}{P}c^k.
 \tag{1}
\]

It is invertible on every polynomial coefficient space. Its inverse, rather
than the merger itself, has the following Beta representation:

\[
 \mathcal M^{-1}H(c)=\int_0^1
        (U+V+1)H(c_j t^{u_j}(1-t)^{v_j})\,dt.
 \tag{2}
\]

Indeed the integral on a monomial is
`(P+Q+1) B(P+1,Q+1)=P!Q!/(P+Q)!`. Integration here can be interpreted
coefficientwise over the rationals; no analytic convergence is involved.
If `u+v=d(1,...,1)`, the Euler factor is the scalar `dm+1` on a degree-`m`
form. In general it cannot be dropped. Nor does vanishing of `H` at one
coefficient vector imply vanishing on the entire coefficient orbit occurring
inside the integral.

The faithful correlated-Gaussian interpolation retained in
[UNRESTRICTED_SFC_MECHANISMS.tex](UNRESTRICTED_SFC_MECHANISMS.tex) has multiplier

\[
 R_{P,Q}(\rho)=\sum_{j=0}^{\min(P,Q)}\binom Pj\binom Qj\rho^j.
 \tag{3}
\]

Thus its operator is `sum_j rho^j U^(falling j)V^(falling j)/(j!)^2`.
For `P<=Q`, the exact Jacobi representation is

\[
 R_{P,Q}(\rho)=(1-\rho)^P
       P_P^{(0,Q-P)}\!\left(\frac{1+\rho}{1-\rho}\right).
 \tag{4}
\]

This follows directly from the finite Jacobi expansion

\[
 P_P^{(0,Q-P)}(x)=\sum_{j=0}^{P}\binom Pj\binom Q{P-j}
        \left(\frac{x-1}{2}\right)^{P-j}
        \left(\frac{x+1}{2}\right)^j.
\]

The normalization agrees with [DLMF 18.5.8](https://dlmf.nist.gov/18.5#E8)
and its Rodrigues formula. The Jacobi weight is positive on `(-1,1)`
([DLMF 18.3](https://dlmf.nist.gov/18.3)); therefore when `P>0`, the
zeros of (3) are simple and negative. This entrywise fact does not establish
regularity of a system of transformed coefficient forms.

At the endpoints, `R_(P,Q)(0)=1` and
`R_(P,Q)(1)=binom(P+Q,P)` by Vandermonde. Equivalently,

\[
 (z^{P+Q})[\underbrace{1,\ldots,1}_{Q+1\text{ nodes}}]
       =\binom{P+Q}{P},
 \tag{5}
\]

the confluent divided difference of order `Q`. Its order depends on `k`;
it is not one fixed finite-order divided-difference operation on all degrees.

There is no probability measure on a positive quadrant whose mixed moments
are `binom(P+Q,P)` for all `P,Q`. Such a measure would have
`E X=E X^2=E Y=E Y^2=1`, forcing `X=Y=1` almost surely, while (1) requires
`E XY=2`. Equivalently its proposed moment value of `(X-Y)^2` is `-2`.
Thus the forward merger is not multiplication by independent positive
auxiliary variables. The Gamma/Beta coupling and the correlated-Gaussian
coupling must not be replaced by that nonexistent positive Mellin kernel.

A valid scalar special case is immediate: if `u=a(1,...,1)` and
`v=b(1,...,1)`, the merger multiplies the degree-`m` moment form by
`binom((a+b)m,am)`. It then preserves its entire zero set. Nonuniform
symmetric batch mergers are instead supplied by
[SFC_CYCLIC_POLYNOMIAL_TILTS.md](SFC_CYCLIC_POLYNOMIAL_TILTS.md); neither statement is an arbitrary merger
preservation theorem.

## 2. Rational atoms with an exact second-and-third-moment zero

Let `(A,B)` have the following three atoms:

| A | B | Probability |
|---:|---:|---:|
| 1 | 2 | `121/296` |
| 1 | -2 | `121/296` |
| `11/3` | 0 | `27/148` |

All odd-in-`B` moments vanish, and

\[
 EA=\frac{55}{37},\qquad EA^2=EB^2=\frac{121}{37},
 \qquad EA^3=3EAB^2=\frac{363}{37}.
 \tag{6}
\]

Consequently `E(A+iB)^2=E(A+iB)^3=0`.

For a strict polynomial Gram form, add an independent vector `(xi,eta)`
uniform on the disk of radius `delta=1/100`. Set `A'=A+xi`, `B'=B+eta`.
Rotational invariance gives `E(xi+i eta)^m=0` for every `m>=1`, so this
addition preserves all holomorphic moments of `A+iB`. If
`v=delta^2/4=1/40000`, write

\[
 H_2=EA'^2=EB'^2=\frac{121}{37}+v,\qquad
 H_3=EA'^3=3EA'B'^2=\frac{363}{37}+\frac{165}{37}v.
 \tag{7}
\]

Put

\[
                U=(A'+B'/3)/2,\qquad V=(A'-B'/3)/2.
 \tag{8}
\]

The atomwise bound `A-|B|/3>=1/3` gives `U,V>=4/25>0` after the noise.
Their distribution is exchangeable and has a positive density on open
ellipses. Thus `E|H(U,V)|^2>0` for every nonzero complex polynomial `H`.
All their mixed moments are rational and finite. The phase

\[
                  (1+3i)U+(1-3i)V=A'+iB'
\]

has its moments of orders two and three equal to zero.

## 3. A collision-free decorated merger destroys regularity

Let `X,Y,Z` be independent unit exponential variables, independent of `(U,V)`.
The source features are `XU,YV`; the merged features are `ZU,ZV`.
With coefficient variables `x,y`, the source moment forms at orders two and
three are exactly

\[
 P_2=\frac{H_2}{9}S,
 \qquad S=5x^2+4xy+5y^2,
\]
\[
 P_3=\frac{H_3}{18}(x+y)Q,
 \qquad Q=15x^2-2xy+15y^2=3S-14xy.
 \tag{9}
\]

They have no common nonzero projective zero. Indeed `S` is nonzero on each
coordinate axis, `S(x,-x)=6x^2`, and `S=Q=0` would force `xy=0`.
The normalized homogeneous resultant is exactly

\[
                       \operatorname{Res}(S,(x+y)Q)=29400.
 \tag{10}
\]

The merger multiplies the coefficient of `x^k y^(m-k)` by `binom(m,k)`.
The resulting endpoint forms are

\[
 \widetilde P_2=\frac{H_2}{9}D,\qquad
 \widetilde P_3=\frac{H_3}{18}\,3(x+y)D,
 \qquad D=5x^2+8xy+5y^2.
 \tag{11}
\]

They have the common nonzero zero `(x,y)=(1+3i,1-3i)`.

This is a genuine integer decoration merger with distinct endpoint columns:
the columns of `XU,YV` in variables `X,Y,U,V` are `(1,0,1,0)` and
`(0,1,0,1)`. Merging the first two rows gives `(1,1,0)` and `(1,0,1)`,
which remain distinct. No feature is constant. Both source and endpoint
feature distributions have strictly positive polynomial Gram forms: condition
on the positive exponential variables and use the open support of `(U,V)`.

The entire physical Gaussian correlation path is explicit. Its normalized
forms are

\[
 S_\rho=5(x^2+y^2)+4(1+\rho)xy,
\]
\[
 (x+y)Q_\rho=(x+y)\bigl(15(x^2+y^2)+(-2+26\rho)xy\bigr).
\]

Since `Q_rho=3S_rho+14(rho-1)xy`, their exact resultant is

\[
               9800(3-2\rho)(\rho-1)^2.
 \tag{12}
\]

It is positive throughout `0<=rho<1` and vanishes at the physical merger
endpoint `rho=1`. Strict feature-Gram positivity persists there.

**Scope.** The joint law of the extra variables `U,V` is not the independent
unit-exponential law. Thus this is not a counterexample to ordinary factorial
row merging or SFC. It refutes preservation based only on positive-measure
moment coupling, strict Gram positivity, integer decorations, distinct columns,
and the genuine correlated-Gaussian merger path. A successful ordinary proof
must use the factorial structure of the other rows as well.

## 4. The obstruction survives uniform private amplification

The same example gives a useful all-order source statement. For every
integer `D>=12`, the features `X^D U,Y^D V` satisfy the two-term nonzero
window at every positive starting order. Merging the private variables to
one `Z` produces `Z^D U,Z^D V`, whose second and third moments still vanish
at the coefficient pair above.

For the source statement, exchangeability gives common positive pure moments
`H_m=EU^m=EV^m`, and Hölder gives

\[
                     EU^kV^{m-k}\le H_m.
\]

The proof of balanced private amplification in
[SFC_BALANCED_PRIVATE_AMPLIFICATION.md](SFC_BALANCED_PRIVATE_AMPLIFICATION.md), Sections 1–3,
uses precisely this inequality, so its suppression and Turán argument apply
to these features too. Explicitly, after normalization by `H_m(Dm)!`, the
pure terms are `x^m+y^m`; every nonpure coefficient is bounded by
`binom(m,k)^(1-D)`. On `max(|x|,|y|)=1`, the error is at most
`B(m)=(m+1)m^(-11)` for `m>=2`, and zero for `m=1`.

The cited Turán bound supplies a power sum in each pair at least
`1/[8e^2(n+1)^2]>1/[72(n+1)^2]` in modulus. For `n>=2`, the relative
error is at most

\[
         72(n+1)^3/n^{11}\le243/256<1;
\]

the left side decreases in `n`. For `n=1`, its bound is `27/64<1`.
Thus the source pair is regular for every starting order. The endpoint
moments merely acquire the nonzero factor `(Dm)!`, so the exact two zeros
from Section 2 persist for every `D`.

This disproves a broad positive-measure de-amplification rule even after
the explicit `6s` amplification threshold. It leaves the ordinary integer
factorial de-amplification obligation open.

## 5. Exact controls

[verify_sfc_row_merging_obstructions.py](verify_sfc_row_merging_obstructions.py) checks the Beta and independent Rodrigues/Jacobi
identities, the nonexistent forward moment kernel, all rational atom/noise
moments needed for (9)–(11), the exact resultant (12), and finite controls of
the private-amplification inequalities. The open-density and all-order
arguments are the written proofs, not consequences of finite sampling.


Replay from the repository root:

```sh
python research/factorial/image-rigidity/verify_sfc_row_merging_obstructions.py
```

The receipt is
[sfc_row_merging_obstructions_verification.json](sfc_row_merging_obstructions_verification.json).
The checker freezes the SHA-256 of this note and the three linked retained
proof inputs, rejects disabled assertions, and exercises both the optimized-mode
and altered-note rejection controls. Its finite cases do not establish the
all-order Turán consequence or the positive-density assertion; those use the
written arguments and the precise retained amplification input.

Internal mathematical review: the root agent separately read and reconstructed
the Beta/Jacobi formulas, rational atoms and noise, source regularity, exact
correlation resultant, and the all-order `D>=12` amplification extension via
Hölder. The source and endpoint support columns and the ordinary-factorial
scope limitation were checked explicitly. This is internal mathematical review,
not external peer review or formal verification; no novelty claim is made.
