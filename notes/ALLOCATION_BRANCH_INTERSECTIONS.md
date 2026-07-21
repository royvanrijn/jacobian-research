# Hensel products and allocation-branch intersections

Fix a collision polynomial

\[
M_0=\prod_{\rho=1}^{\ell}(W-r_\rho)^{m_\rho}
\]

with distinct roots, and two allocations in the same component,

\[
A=((i_\rho,j_\rho))_\rho,
\qquad
A'=((i'_\rho,j'_\rho))_\rho,
\]

where

\[
2i_\rho+3j_\rho=2i'_\rho+3j'_\rho=m_\rho,qquad
\sum i_\rho=sum i'_\rho,qquad
\sum j_\rho=sum j'_\rho.
\]

There are unique integers `k_rho` such that

\[
i_\rho-i'_\rho=3k_\rho,qquad
j_\rho-j'_\rho=-2k_\rho.                         \tag{1}
\]

The global degree constraints give

\[
\sum_\rho k_\rho=0.                              \tag{2}
\]

Thus every transfer in one direction is compensated at other collision
roots.

## Local transfer blocks

For `k>=1`, let `Z_k` be the completion along `U=S^3,V=S^2` of

\[
U^2=V^3,qquad \deg U=3k,quad\deg V=2k,quad\deg S=k.
\]

The cases currently known are

\[
\mathfrak Z_1=k[[s,\epsilon]]/(\epsilon^2)
\]

and

\[
\mathfrak Z_2=
k[[p,q,X,Y]]/(X^3,2XY-pX^2,Y^2-qX^2).
\]

They are finite flat of ranks `2` and `4` over their reduced factor spaces.

## Hensel-product theorem for strong equality

Let `D^=` denote the completed correspondence of the two allocations with the
strong equation

\[
M_A=M_{A'}.
\]

Then

\[
\boxed{
D^=\simeq R_{\mathrm{com}}
\widehat\otimes
\bigwidehat\otimes_{\rho:k_\rho\ne0}
\mathfrak Z_{|k_\rho|}.}                         \tag{3}
\]

Here `R_com` is a regular complete local ring containing the root-position
parameters and the deformations of factor blocks on which the two allocations
agree.

### Proof

Choose pairwise coprime discs around the distinct roots `r_rho`.  Formal
Hensel factorization uniquely decomposes every monic deformation of `M_0` as
a product of monic cluster polynomials of degrees `m_rho`.  Equality of the
two global monic products is therefore equivalent to equality of the
corresponding cluster polynomials root by root.

If `k_rho=0`, the two allocations of that cluster agree and its correspondence
is the diagonal, hence smooth.  If `k_rho!=0`, canceling the common double and
triple factors leaves

\[
U_\rho^2=V_\rho^3
\]

with degrees `3|k_rho|,2|k_rho|`; this is `Z_|k_rho|`.  Independence of the
Hensel clusters gives the completed tensor product (3).

The incidence equation `Phi=0` has nonzero derivative in a root-position
direction on a dense open subset of every exact collision stratum.  Formal
implicit-function elimination removes one smooth parameter from `R_com` and
does not alter the Artin multiplication of the transfer blocks.

## The global affine equalizer theorem

The normalized coefficient map forgets the constant and linear coefficients.
The actual branch correspondence `D^aff` is consequently defined by

\[
M_A-M_{A'}=\lambda W+\mu.                         \tag{4}
\]

There is a closed immersion

\[
D^=\hookrightarrow D^{\mathrm{aff}}.             \tag{5}
\]

Equation (4) is global: unlike strong equality, it does not formally split
root by root. It is nevertheless rigid.

### Product formulation

Let \(U_\alpha,V_\alpha\) be monic local transfer factors of degrees
\(3h_\alpha,2h_\alpha\), and let \(e_\alpha\) be signed multiplicities
satisfying

\[
                 \sum_\alpha e_\alpha h_\alpha=0.                \tag{6}
\]

Put

\[
\begin{array}{ll}
 Q_+=\prod_\alpha U_\alpha^{e_\alpha^+},&
 Q_-=\prod_\alpha U_\alpha^{e_\alpha^-},\\[2mm]
 R_+=\prod_\alpha V_\alpha^{e_\alpha^+},&
 R_-=\prod_\alpha V_\alpha^{e_\alpha^-}.
\end{array}                                                       \tag{7}
\]

The two allocation polynomials after common-factor cancellation are

\[
                   M_+=Q_+^2R_-^3,\qquad M_-=Q_-^2R_+^3.          \tag{8}
\]

Thus the global affine equalizer is defined by

\[
                   M_+-M_-\in\langle1,W\rangle.                  \tag{9}
\]

Separate affine conditions on \(Q_+-Q_-\) and \(R_+-R_-\) are not the
allocation correspondence: the square and cube factors occur crosswise in
(8). If an exponent \(|e_\alpha|>1\) is concentrated at one Hensel cluster,
its effective local block is \(Z_{|e_\alpha|h_\alpha}\); it becomes a tensor
power of \(Z_{h_\alpha}\) only when the copies lie at distinct clusters.

**Theorem (arbitrary global transfer equalizer).** For every finite transfer
system (6) over a characteristic-zero coefficient algebra,

\[
                         \boxed{D^{\mathrm{aff}}=D^=.}            \tag{10}
\]

For distinct Hensel clusters with effective transfers \(k_\rho\), the
completed correspondence is canonically

\[
 \boxed{
 D^{\mathrm{aff}}\simeq D^=
 \simeq R_{\mathrm{com}}\widehat\otimes
       \bigwidehat\otimes_{\rho:k_\rho\ne0}Z_{|k_\rho|}.}         \tag{11}
\]

In particular it has transverse rank

\[
 \boxed{L(A,A')=\prod_{\rho:k_\rho\ne0}2^{|k_\rho|}
                   =2^{\sum_\rho|k_\rho|}.}                      \tag{12}
\]

### Universal Wronskian proof

Write the two global branch polynomials as

\[
 M=Q^2R^3,\qquad N=S^2T^3,
\]

where the factors are monic and

\[
 2\deg Q+3\deg R=2\deg S+3\deg T=n.
\]

Suppose `M-N=lambda W+mu` in the universal completed coordinate ring of
`D^aff`. Set

\[
 \mathcal W=NM'-MN'=N(M-N)'-(M-N)N'.             \tag{13}
\]

The derivative factorizations make `mathcal W` divisible by

\[
 F=QR^2ST^2,
 \qquad \deg F=n+{\deg R+\deg T\over2}.           \tag{14}
\]

The two triple-factor degrees have the same parity. Unless both vanish, their
sum is at least two, so `deg F>=n+1`; the right side of (8) has degree at most
`n`. Since `F` is monic, this forces `mathcal W=0` over any `Q`-algebra,
including rings with nilpotents. The coefficients of degrees `n` and `n-1`
in

\[
 N\lambda-N'(\lambda W+\mu)=0
\]

then give `(1-n)lambda=0` and `-n mu=0`. Thus `lambda=mu=0`. If both triple
factors are absent, `(Q-S)(Q+S)` is affine; the second factor has degree
`n/2>=2` and unit leading coefficient, so the product is zero as well.

Every step takes place in the universal coordinate algebra and remains valid
with zero divisors and nilpotents: only monicity and the invertibility of
nonzero integers are used. Therefore the affine-equalizer ideal itself
contains \(\lambda,\mu\). The reverse closed immersion is tautological, so
\(D^{\mathrm{aff}}=D^=\) as schemes, not merely on geometric points.

For the product formulation, set

\[
 t=\sum_\alpha e_\alpha^+h_\alpha
  =\sum_\alpha e_\alpha^-h_\alpha.
\]

Then

\[
 \deg Q_+=\deg Q_-=3t,\qquad
 \deg R_+=\deg R_-=2t,\qquad n=12t,
\]

and the Wronskian divisor has degree \(14t>12t\). Thus the argument applies
uniformly to arbitrary exponents and any number of clusters.

This proof takes place before Hensel decomposition. It therefore kills the
two coefficients shared by all local blocks at once, rather than trying to
cancel separate rootwise affine errors. Combining (10) with (3) and C22 gives
(11)--(12), including every higher compatibility and every nilpotent in the
Boolean local relations.

The earlier length-four and length-sixteen calculations for `(1,-1)`,
`(2,-2)`, and `(2,-1,-1)` remain exact coordinate audits of this theorem; they
are no longer separate hypotheses or exceptional cases.

### Diagnostic higher-transfer configurations

A transfer vector must sum to zero. Thus “one \(Z_3\) plus one compensating
\(Z_1\)” is not literally possible: \(3-1\ne0\). The elementary
compensation is three distinct \(Z_1\) blocks. If one instead cubes a single
primitive \(Z_1\) contribution at one cluster, that cluster is the effective
block \(Z_3\).

The first corrected configurations are:

| requested configuration | effective transfer vector | strong/affine transverse factor | rank |
|---|---:|---|---:|
| one \(Z_3\) with elementary compensation | \((3,-1,-1,-1)\) | \(Z_3\widehat\otimes Z_1^{\widehat\otimes3}\) | \(64\) |
| two interacting \(Z_3\) blocks | \((3,-3)\) | \(Z_3\widehat\otimes Z_3\) | \(64\) |
| one \(Z_4\) with elementary compensation | \((4,-1,-1,-1,-1)\) | \(Z_4\widehat\otimes Z_1^{\widehat\otimes4}\) | \(256\) |
| mixed \(Z_3\) and \(Z_2\) | \((3,-2,-1)\) | \(Z_3\widehat\otimes Z_2\widehat\otimes Z_1\) | \(64\) |
| three nontrivial blocks | \((4,-3,-1)\) | \(Z_4\widehat\otimes Z_3\widehat\otimes Z_1\) | \(256\) |

All five have affine and strong equations equal by the same Wronskian degree
gap. Their special-fiber Hilbert series are respectively
\((1+x)^6,(1+x)^6,(1+x)^8,(1+x)^6,(1+x)^8\). The executable regression checks
the cross-coupled product identity, divisor degree, rank, and Hilbert
coefficients in every row. These calculations diagnose the theorem; none is
used as an inductive endpoint.

## Conductor

Let `B_A,B_A'` be the
two regular completed normalization branches and let

\[
B_A\twoheadrightarrow D,
\qquad
B_{A'}\twoheadrightarrow D
\]

be their quotient maps to the product algebra (3), after incidence
elimination. Then the completed component ring is

\[
\widehat{\mathcal O}_{C,\bar\xi}
=B_A\mathop\times_D B_{A'},                      \tag{9}
\]

and its conductor in `B_A\oplus B_A'` is the direct sum of the two
quotient kernels. Thus C22 and the global Wronskian determine both the
completed intersection and conductor for every allocation pair.

## The global `(2,-2)` theorem

Take transfer vector

\[
(k_1,k_2)=(2,-2),
\]

over a collision with two multiplicity-twelve roots.  Formula (3) gives the
strong transverse algebra

\[
\mathfrak Z_2\widehat\otimes\mathfrak Z_2
\]

of length sixteen.

Use the twenty nonleading coefficients of the four monic factors as ambient
coordinates.  Equality of normalized coefficients supplies the coefficients
of degrees `2,...,23`.  Their linearization has rank twelve, hence tangent
dimension eight.  Four tangent directions are the reduced `S_0,S_1`
factor-space directions; slice them away and use transverse coordinates
`e_0,...,e_3`.

Twelve equations and twelve normal coefficient variables have an invertible
linear block.  Formal implicit elimination through degree three leaves the
homogeneous local standard basis

\[
\begin{aligned}
e_0^2-9e_1^2,&\qquad e_0e_1+3e_1^2,\qquad e_1^3,\\
e_2^2,&\qquad e_2e_3,\qquad e_3^3.               \tag{10}
\end{aligned}
\]

The first line becomes `(X^3,XY,Y^2)` after taking
`X=e_0`, `Y=e_0+3e_1`; the second does so with `X=e_3`, `Y=e_2`.
Consequently (10) is the tensor product of two coincident `Z_2` fibers.  Its
standard basis has sixteen elements and Hilbert function `(1,4,6,4,1)`.

The initial ideal of the full affine correspondence contains (10), so its
transverse length is at most sixteen.  The strong Hensel product is already a
closed length-sixteen subscheme.  The equalizer criterion gives

\[
\boxed{
D^{\mathrm{aff}}_{(2,-2)}=D^=_{(2,-2)}
=\mathfrak Z_2\widehat\otimes\mathfrak Z_2.}     \tag{11}
\]

This is the first genuinely global higher-transfer branch-intersection
theorem.

## The mixed `(2,-1,-1)` theorem

Place the `k=2` block at `W=0` and the two compensating elementary blocks at
`W=1,2`.  The two sheets have base factors

\[
\begin{array}{c|cc}
&Q&R\\ \hline
A&W^6&((W-1)(W-2))^2\\
A'&((W-1)(W-2))^3&W^4.
\end{array}
\]

The same twenty-variable implicit reduction again has tangent dimension eight:
four reduced directions and four transverse directions.  Through cubic order,
the transverse local standard basis is

\[
\begin{aligned}
e_0^2-20e_1^2,&\qquad 2e_0e_1+9e_1^2,\qquad e_1^3,\\
e_2^2,&\qquad e_2e_3,\qquad e_3^3.               \tag{12}
\end{aligned}
\]

It has sixteen standard monomials.  Therefore the affine correspondence has
transverse length at most sixteen.  The strong Hensel product

\[
\mathfrak Z_2\widehat\otimes
\mathfrak Z_1\widehat\otimes\mathfrak Z_1
\]

is a closed length-sixteen subscheme, and hence

\[
\boxed{
D^{\mathrm{aff}}_{(2,-1,-1)}=D^=_{(2,-1,-1)}
=\mathfrak Z_2\widehat\otimes
 \mathfrak Z_1\widehat\otimes\mathfrak Z_1.}    \tag{13}
\]

The particular linear transverse slice used in (12) is not the canonical
product slice: its two-variable factors are flat degenerations of the local
blocks.  The length argument identifies the unsliced completed
correspondence, where the Hensel factors are canonical.

The two length-sixteen calculations predate the universal proof. The
higher-transfer diagnostics above replace the former “next untested case”
boundary: no global transfer vector remains to be classified separately.
