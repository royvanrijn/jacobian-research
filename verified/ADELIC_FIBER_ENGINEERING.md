# Adelic engineering of complete Keller fibers

This note combines the
[exact real-sheet spectrum](REAL_FIBER_SPECTRUM.md) with the
[effective finite-field Chebotarev theorem](../extended-geometry/FINITE_FIELD_CHEBOTAREV.md).
It gives simultaneous archimedean and nonarchimedean control of one complete
fiber of the explicit all-degree Keller map.

## 1. Simultaneous realization theorem

For every `N>=3`, let

\[
 F_N:\mathbb A^3_{\mathbb Q}\longrightarrow\mathbb A^3_{\mathbb Q}
\]

be the explicit map attached to the rational seed `H_N` in the
[all-degree rational-fiber theorem](ALL_DEGREE_RATIONAL_FIBERS.md).  Its
inverse pencil on `C=1` is

\[
 E_{s,t}(W)=H_N(W)-sW+t,\qquad (A,B,C)=(t,s,1).       \tag{1}
\]

Fix the following data:

1. an allowed real-root count
   \[
   j\in\{N,N-2,\ldots,N\bmod 2\};
   \]
2. a finite set `S` of rational primes outside the explicit bad integer of
   `H_N`; and
3. for every `p in S`, a partition `lambda_p` of `N` for which there is a
   squarefree pencil witness over `F_p`.

Then there is `(s,t) in Q^2` such that:

* `E_(s,t)` is squarefree and has exactly `j` real roots;
* for every `p in S`, its reduction is squarefree of factor-degree partition
  `lambda_p`; and
* the target `y=(t,s,1)` has a complete regular Keller fiber of degree `N`
  whose real points and residue degrees have exactly these prescribed data.

The finite-field Chebotarev theorem guarantees the local witness in item 3
for every partition simultaneously whenever `p` is good and its residue-field
size satisfies the effective inequalities displayed there.  Thus arbitrary
finite collections of partitions can be imposed at all sufficiently large
good primes.  Smaller good primes are also allowed whenever direct search
finds the required witnesses.

## 2. Constructive weak approximation

Let `Omega_j` be a nonempty real pencil chamber with `j` real roots, supplied
by the real-sheet theorem.  For every `p in S`, choose a finite-field witness

\[
 (\bar s_p,\bar t_p)\in U_{H_N}(\mathbb F_p)
\]

of type `lambda_p`.  Its residue box

\[
 \Omega_p=\{(s,t)\in\mathbb Z_p^2:
 (s,t)\equiv(\bar s_p,\bar t_p)\pmod p\}             \tag{2}
\]

is nonempty and `p`-adically open.  Weak approximation gives a rational point
in

\[
 \Omega_j\times\prod_{p\in S}\Omega_p.              \tag{3}
\]

In this situation weak approximation has an elementary constructive proof.
Choose a bounded rational rectangle

\[
 I_s\times I_t\subset\Omega_j
\]

and put `M=prod_(p in S) p`.  The Chinese remainder theorem gives residues
`r_s,r_t mod M` representing all selected local pairs.  Choose

\[
 D=1+kM
\]

with `k` sufficiently large.  The grids

\[
 \left\{\frac{r_s+mM}{D}:m\in\mathbb Z\right\},\qquad
 \left\{\frac{r_t+mM}{D}:m\in\mathbb Z\right\}       \tag{4}
\]

have mesh `M/D`, so they meet `I_s` and `I_t` once `D` is large enough.
Because `D=1 mod p`, the selected rational coordinates reduce to
`(bar s_p,bar t_p)` at every `p`.  Formula (4) is a deterministic witness
generator once the real rectangle and the finite-field witnesses are known.
It is implemented as
`jcsearch.chebotarev.constructive_weak_approximation_lift`.

## 3. The complete fiber algebra and local splitting

Membership in `Omega_j` makes `E_(s,t)` squarefree over `Q`.  In the quotient

\[
 K_{s,t}=\mathbb Q[W]/(E_{s,t}),                     \tag{5}
\]

the derivative `E'_(s,t)` is therefore a unit.  The reconstruction formulas
on `C=1`, interpreted modulo `E_(s,t)`, identify the complete fiber scheme
with

\[
 \boxed{F_N^{-1}(t,s,1)\simeq\operatorname{Spec}K_{s,t}.}       \tag{6}
\]

At a selected good prime, the coefficients are integral, the leading
coefficient is a unit, and the reduction is squarefree.  After dividing by
the leading coefficient,

\[
 \mathbb Z_p[W]/(E_{s,t})
\]

is finite etale over `Z_p`.  If

\[
 \lambda_p=(d_{p,1},\ldots,d_{p,r_p}),
\]

then

\[
 K_{s,t}\otimes\mathbb Q_p
 \simeq\prod_{i=1}^{r_p}K_{p,i},\qquad
 [K_{p,i}:\mathbb Q_p]=d_{p,i},                     \tag{7}
\]

and every factor `K_(p,i)/Q_p` is unramified.  Thus the prescribed partitions
are actual local decomposition types of the fiber algebra, not only modular
fingerprints.

## 4. Number fields of every signature

Adjoin one more sufficiently large good prime `p_0` and prescribe the
one-part partition

\[
 \lambda_{p_0}=(N).                                 \tag{8}
\]

The reduction of `E_(s,t)` is irreducible at `p_0`; hence `E_(s,t)` is
irreducible over `Q`.  Therefore (5) is a degree-`N` number field.  Its real
embeddings are the real roots of `E_(s,t)`, so its signature is

\[
 \boxed{(r_1,r_2)=\left(j,\frac{N-j}{2}\right).}     \tag{9}
\]

Consequently every degree `N>=3` and every degree-`N` signature occurs as a
complete fiber field of the single explicit map `F_N` for that degree.
Arbitrarily many further unramified splitting conditions may be imposed at
the same time.  No Hilbert-irreducibility input is needed: the single local
`N`-cycle in (8) supplies global irreducibility.

## 5. An explicit quartic witness

For `N=4`, the all-degree seed is

\[
 H_4(W)=\frac14W^4-\frac32W^3+\frac54W^2.
\]

Take

\[
 s=-\frac{308}{103},\qquad t=\frac{617}{309}.        \tag{10}
\]

These numbers lie in a zero-real-root chamber near `(s,t)=(-2.999,2)`.  The
inverse polynomial is

\[
 E(W)=\frac14W^4-\frac32W^3+\frac54W^2
      +\frac{308}{103}W+\frac{617}{309}.             \tag{11}
\]

Equivalently, the fiber field is generated by a root of the primitive
integral polynomial

\[
 P_4(W)=309W^4-1854W^3+1545W^2+3696W+2468.          \tag{11a}
\]

Its reductions have the following monic factorizations, up to a nonzero
scalar:

\[
 \begin{aligned}
 E(W)\bmod7&=W^4+W^3-2W^2-3,
     &&\text{irreducible},\\
 E(W)\bmod11&=(W+2)(W-4)(W^2-4W+5).
 \end{aligned}                                      \tag{12}
\]

Hence (11) defines a quartic field of signature `(0,2)`, the prime `7` is
inert, and `11` has unramified splitting type `(2,1,1)`.  The corresponding
complete Keller target is

\[
 (A,B,C)=\left(\frac{617}{309},-\frac{308}{103},1\right).       \tag{13}
\]

This small-prime example uses direct local witnesses; it is not an assertion
that `7` and `11` satisfy the deliberately coarse uniform Chebotarev bound.

## 6. Boundaries of the statement

The theorem prescribes squarefree, hence unramified, local types.  Ramified
local algebras require separate congruence conditions on discriminant strata.
It constructs number fields occurring in the two-parameter pencil; it does
not assert that every number field of a given degree and signature occurs.
Finally, all fiber statements use `C=1`; the separate `C=0` boundary
decomposition is irrelevant here.

## Verification

Run

```bash
.venv/bin/python scripts/verify_adelic_fiber_engineering.py
```

The checker verifies the constructive approximation grid, the primitive
integral polynomial, the exact quartic real-root count, both squarefree
modular factorization types, global irreducibility, and the displayed target
conversion.
