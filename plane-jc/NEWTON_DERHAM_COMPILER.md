# Newton-chain to de Rham compiler

## Status

The repository now has an executable middle end for weighted-Wronskian
blocks:

```text
derived Newton bands
    -> WeightedWronskianIR
    -> superelliptic character data
    -> supported primitive solve
    -> certified de Rham obstruction coordinates
    -> Gauss--Manin connection
```

The implementation is
[`cas/newton_derham_compiler.py`](cas/newton_derham_compiler.py), built on
[`cas/superelliptic_derham.py`](cas/superelliptic_derham.py).  It does not yet
derive Laurent bands from a corner chain.  An incomplete chain is rejected
with a list of the missing front-end data.

The preceding arrow from an already published Laurent polygon to exact band
supports and bracket layers is now implemented in
[`LAURENT_BAND_FRONTEND.md`](LAURENT_BAND_FRONTEND.md) and
[`cas/laurent_band_frontend.py`](cas/laurent_band_frontend.py).  The remaining
gap is specifically the theorem-heavy arrow from a corner chain to an
exhaustive normalized polygon list.

## 1. The `(72,108)` obstruction certificate

For

\[
 A=t+a_2t^2+\cdots+a_7t^7+t^8,
 \qquad
 L_A(D)=2AD'-3A'D,
\]

the supported triangular solve constructs

\[
 D_*=d_2t^2+\cdots+d_{12}t^{12}
\]

and proves the polynomial identity

\[
 L_A(D_*)-t^2
 =C_0t^{13}+C_1t^{14}+\cdots+C_5t^{18}.       \tag{1}
\]

The six `C_i` are exactly the six archived first-block compatibility
equations.  Dividing (1) by `2y^5`, where `y^2=A`, gives

\[
 \left[\frac{t^2dt}{2y^5}\right]
 =-\frac12\sum_{i=0}^5 C_i
   \left[\frac{t^{13+i}dt}{y^5}\right].        \tag{2}
\]

Thus it remains only to prove that the six displayed tail classes form a
basis of compact `H^1_dR`.

They are second-kind differentials.  At a finite branch point, `t` is a
power series in `y^2`, so `t^k dt/y^5` has no `dy/y` term.  At infinity,
`y` has order `-4`; the exponents `13,...,18` make the differentials regular.

If a linear combination of them were exact, its anti-invariant primitive
could be written `D/y^3`.  The finite pole orders force `D` to be a
polynomial, and regularity at infinity forces `deg D<=12`.  Projecting
`L_A(D)` to degrees `0,...,12` gives a triangular matrix with diagonal

\[
 2k-3,\qquad 0\le k\le12.
\]

Its determinant is

\[
 \prod_{k=0}^{12}(2k-3)=41247931725\ne0.       \tag{3}
\]

Hence `D=0`, so the six tail classes are independent.  The curve has genus
three, so they are a basis.  Equations (1)--(3) are the missing explicit
certificate that the six supported equations are six de Rham obstruction
coordinates.  It avoids a large symbolic modular inverse and introduces no
parameter denominator.  The cohomological interpretation is, as usual,
restricted to the squarefree locus of `A`.

## 2. Compiler IR

`NewtonChainIR` records corners, multiplicities, source, status, missing
front-end data, and source-level reconciliation notes.  A complete
`WeightedWronskianIR` additionally records:

- the normalized `A`, `R`, and coordinate `t`;
- the covering exponent `a` and primitive weight `b`;
- the allowed support of `D` and its full valuation bounds;
- endpoint normalizations and residual scaling weights.

Compilation returns the genus, character `b mod a`, affine and compact
dimensions, residue functional, compact basis, supported primitive,
compatibility vector, exact residual identity, tail coordinates, and two
fingerprints.  The local-system fingerprint deliberately omits `R` and
primitive support; the supported-problem fingerprint includes them.

Fingerprints are filters, not proofs of reuse.  The executable
`certify_local_system_reuse` routine instead requires an explicit base map and
nonzero coordinate scalings and checks the polynomial identity

\[
 A_{\rm candidate}(\lambda t)=\mu^a A_{\rm reference}(t).
\]

It also checks the covering exponent and character and forbids a purported
base map from depending on the fiber coordinate.  A successful certificate
identifies the base-changed cyclic curve families, hence the corresponding
pulled-back character Gauss--Manin local systems on the locus where the two
coordinate scales are nonzero.  It deliberately does not compare `R`,
primitive support, or residual scaling: those belong to the stronger
supported obstruction problem.

## 3. Gauss--Manin layer

For parameters `u_j` and a compact basis representative

\[
 e_i=\left[p_i(t,u)\frac{dt}{y^r}\right],
\]

the implementation uses

\[
 \nabla_{\partial_{u_j}}e_i
 =\operatorname{Red}_r\left[
 \frac{A\,\partial_{u_j}p_i-(r/a)p_i\,\partial_{u_j}A}
      {y^{a+r}}dt
 \right].                                           \tag{4}
\]

The ordinary Hermite reducer supplies the columns of the connection matrix.
The regression family

\[
 y^2=t^3+ut+v
\]

produces two exact `2x2` matrices, verifies zero curvature, and verifies that
their poles are supported on `4u^3+27v^2=0`, the discriminant divisor.
A degree-eight slice `y^2=t^8+ut^7+t` also produces a dense `6x6`
genus-three matrix whose denominators are scalar multiples of
`46656u^7+823543`, its discriminant up to sign.

The obstruction section need not be horizontal.  Reusing a local system and
reusing a plane-JC obstruction problem are therefore different claims.

## 4. Repeated-tail audit for `(96,144)`

The 2017 Section 7 table contains

\[
 (8,40)\longrightarrow(8,28)\longrightarrow(11/4,7),
 \qquad(m,n)=(3,2).
\]

That table supplies no Laurent coordinate change, band supports, leading
Wronskian, fixed predecessor coefficients, or residual scaling.  The compiler
records all five items as missing and refuses to assign the old genus-three
block to the new chain.

There is also a source-level issue that must be reconciled first.  The 2016
paper *The two-dimensional Jacobian conjecture and the lower side of the
Newton polygon*, in the remark following Proposition 3.29, says that
`B_0=(8,28), B_1=(8,40)` leads to the impossible last lower corner `(8,4)`
and can be discarded.  The impossibility of points of the relevant form is
proved in Proposition 3.29, but the application to this `B_0,B_1` pair is
stated as a straightforward remark.  The later table nevertheless retains
the length-two row.

Consequently `(96,144)` is currently useful as a compiler-reuse thought
experiment, but it should not be called the first live frontier case until
the two source statements are reconciled.  No `(96,144)` Laurent system can
be faithfully derived from the table alone.

If the row survives that reconciliation, the required reuse certificate is
a base map `phi`, an equivariant curve-family isomorphism, and a gauge/exact
comparison

\[
 C_{96}\simeq\phi^*C_{72},
 \qquad
 [\omega_{96}]=G\phi^*[\omega_{72}].
\]

Primitive support and residual scaling must then be compared separately.
The new exact curve-identity checker implements the first of these gates for
monomial coordinate scalings.  Since `(96,144)` still has no derived `A`, it
cannot pass even that gate; matching fingerprints may not be substituted for
the missing identity.

## 5. The `(75,125)` architecture row

The first numerical frontier is family `F2` at `j=1`:

\[
 (5,20)\longrightarrow(7/5,2),\qquad(m,n)=(3,5).
\]

The 2017 chain theorem additionally forces the Puiseux translation, two
consecutive edges, and the transformed bracket `X^4`; the terminal type-I
edge has an exact rank-zero normalization.  This partial source certificate
is encoded in [`cas/f2_75_125_frontend.py`](cas/f2_75_125_frontend.py) and
derived in [`F2_75_125_DERIVATION.md`](F2_75_125_DERIVATION.md).  It does not
classify the lower Laurent boundary.  The 2014 polynomial-system paper treats
the `j=0` member with degrees `(50,75)` and ratio `2:3`; its two modified coefficient systems
do not derive the `3:5` bands needed at `(75,125)`.  The compiler records three
missing inputs and rejects this row too.  In particular, neither the covering
exponent, character, genus, nor obstruction rank can yet be inferred
faithfully from the family table.

This makes `(75,125)` the clean next front-end derivation problem once the
published normal-form machinery is formalized.  It is not yet a test vector
for the de Rham middle end.

## 6. Reproduction

Run:

```bash
python3 plane-jc/cas/test_superelliptic_derham.py
python3 plane-jc/cas/test_newton_derham_compiler.py
```

The second command certifies the `(72,108)` tail basis and confirms that the
local-system checker accepts exact coordinate rescalings, rejects a distinct
curve with the same coarse fingerprint, and confirms that the incomplete
`(96,144)` and `(75,125)` records cannot be compiled.
