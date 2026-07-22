# A weighted-seed to Gaussian-moment bridge

This note turns every nonconstant weighted inverse seed into an explicit
four-real-Gaussian counterexample to the Gaussian Moments Conjecture.  The
construction is uniform in the seed and retains an exact branch of the pencil

\[
 H(W)-sW+t=0.
\]

It uses the Lagrange--Good viewpoint highlighted by Christopher D. Long, but
the correction below is a repository-derived construction.  It has not been
externally reviewed and is not claimed to be Long's formula or a consequence
stated in his paper.

## 1. Statement

Let `H in C[z]` be nonconstant with `H(0)=0`, and fix `lambda != 0`.  Put

\[
 h(z)=1+\lambda H(z),\qquad
 D(z)=h(z)-zh'(z),\qquad
 R(z)=\frac{D(z)-1}{z}.                              \tag{1.1}
\]

Because `h(0)=D(0)=1`, the quotient `R` is a polynomial.  In variables
`z,y`, define

\[
 \begin{aligned}
 L(z,y)&=D(z)y+zR(z),\\
 B(z,y)&=h(z)R(z)(1+y)-h(z)h'(z)(1+y)^2,\\
 \Phi_1(z,y)&=h(z),\\
 \Phi_2(z,y)&=-h(z)R(z)(1+y)+L(z,y)B(z,y).
 \end{aligned}                                      \tag{1.2}
\]

All four expressions are polynomials.

Let `X_1,Y_1,X_2,Y_2` be independent standard real Gaussians and set

\[
 Z_j=\frac{X_j+iY_j}{\sqrt2},\qquad
 W_j=\frac{X_j-iY_j}{\sqrt2}\quad(j=1,2).           \tag{1.3}
\]

Define the explicit polynomials

\[
 P_{H,\lambda}
 =W_1\Phi_1(Z_1,Z_2)+W_2\Phi_2(Z_1,Z_2),
 \qquad Q=Z_1.                                      \tag{1.4}
\]

### Theorem 1.1

For every `m>=1`,

\[
 \mathbb E(P_{H,\lambda}^m)=0,                      \tag{1.5}
\]

and

\[
 \mathbb E(QP_{H,\lambda}^m)
 =(m-1)![z^{m-1}]h(z)^m.                            \tag{1.6}
\]

The right side of (1.6) is nonzero for infinitely many `m`.  Thus
`(P_{H,lambda},Q)` is a counterexample to `GMC(4)` for every nonconstant
normalized seed `H` and every nonzero `lambda`.

This is a family statement, not merely a bounded computation.  The exact
checker supplies independent finite-range Wick regressions after the written
all-order proof.

## 2. The formal Gaussian--Lagrange input

Let `Phi=(Phi_1,...,Phi_r)` be a polynomial map and let

\[
 \boldsymbol g(u)=u\Phi(\boldsymbol g(u))
 \in u\mathbb C[[u]]^r                                \tag{2.1}
\]

be its unique formal fixed point.  For independent circular complex Gaussian
pairs `(Z_j,W_j)` and

\[
 P=W\mathbin\cdot\Phi(Z),
\]

circular Wick contraction gives

\[
 \mathbb E(W^\alpha A(Z))=\partial^\alpha A(0).
\]

Consequently the formal exponential generating function is

\[
 \mathbb E\!\left(A(Z)e^{uP}\right)
 =\sum_{\alpha\in\mathbb N^r}
   \frac{u^{|\alpha|}}{\alpha!}
   \partial^\alpha\!\left(A\Phi^\alpha\right)(0).
                                                               \tag{2.2}
\]

The [constant-term Gaussian--Lagrange lemma](FORMAL_GAUSSIAN_LAGRANGE_LEMMA.md)
proves, for arbitrary polynomial `Phi`, that

\[
 \boxed{
 \mathbb E\!\left(A(Z)e^{uP}\right)
 =\frac{A(\boldsymbol g(u))}
        {\det(I-uJ\Phi(\boldsymbol g(u)))}.}
                                                               \tag{2.3}
\]

Its standalone proof works in `C[[u,z_1,...,z_r]]`, proves coefficientwise
finiteness, justifies the formal residue change when `Phi(0) != 0`, and fixes
both the determinant orientation and the complex Gaussian derivative
convention.  That constant-term case is essential here: by (1.2),
`Phi_1(0,0)=1`.

Good's multivariable inversion theorem and Long's use of this architecture
retain their external provenance.  Formula (2.3), including the moving fixed
point caused by the constant term, is proved locally rather than used as a
black box.

All series in this note are formal.  No analytic integrability assertion for
`exp(uP)` is needed.

## 3. The determinant-cancelling auxiliary coordinate

Let `g(u)` be the unique solution of

\[
 g=u h(g),                                            \tag{3.1}
\]

and put

\[
 k(u)=\frac1{D(g(u))}-1.                              \tag{3.2}
\]

Since `D=1+zR`, one has

\[
 L(g,k)=D(g)\left(\frac1{D(g)}-1\right)+gR(g)=0.     \tag{3.3}
\]

It follows from (1.2) and (3.1) that

\[
 \begin{aligned}
 u\Phi_1(g,k)&=g,\\
 u\Phi_2(g,k)
 &=-\frac{u h(g)R(g)}{D(g)}
   =-\frac{gR(g)}{D(g)}
   =\frac1{D(g)}-1=k.
 \end{aligned}                                      \tag{3.4}
\]

Thus `(g,k)` is exactly the fixed point in (2.1).

The purpose of the apparently redundant term `LB` is to prescribe the normal
derivative without changing the value on the fixed graph `L=0`.  Since
`partial_y L=D`, restriction to that graph gives

\[
 \begin{aligned}
 \partial_y\Phi_2
 &=-hR+D\left(\frac{hR}{D}-\frac{hh'}{D^2}\right)\\
 &=-\frac{hh'}D.                                     \tag{3.5}
 \end{aligned}
\]

The first component is independent of `y`, so the relevant Jacobian matrix is
triangular.  Along `(g,k)`,

\[
 \begin{aligned}
 \det(I-uJ\Phi)
 &=(1-uh'(g))\left(1+u\frac{h(g)h'(g)}{D(g)}\right)\\
 &=\frac{D(g)}{h(g)}
   \frac{D(g)+g h'(g)}{D(g)}=1.                     \tag{3.6}
 \end{aligned}
\]

The final equality uses `D+zh'=h`.  This is the exact determinant
cancellation promised in the construction.

## 4. Moments and nonvanishing

Apply (2.3) and (3.6).  For every polynomial `A(z,y)`,

\[
 \mathbb E\!\left(A(Z_1,Z_2)e^{uP_{H,\lambda}}\right)
 =A(g(u),k(u)).                                      \tag{4.1}
\]

Taking `A=1` gives a generating function equal to one, which proves (1.5)
coefficient by coefficient.  Taking `A=z` gives

\[
 \sum_{m\ge0}\frac{u^m}{m!}
 \mathbb E(QP_{H,\lambda}^m)=g(u).                  \tag{4.2}
\]

Univariate Lagrange inversion applied to (3.1) yields

\[
 [u^m]g(u)=\frac1m[z^{m-1}]h(z)^m,
\]

which is (1.6).

It remains to check the Mathieu counterexample condition, which asks for
mixed moments that do not vanish eventually.  If `g` were a polynomial of
degree `e>=1` and `h` had degree `d>=1`, equation (3.1) would equate degrees
`e` and `1+de`, an impossibility.  Hence `g` is not a polynomial and has
infinitely many nonzero positive-degree coefficients.  Equation (4.2) proves
that the mixed moments are nonzero infinitely often.

## 5. Exact relation with the weighted inverse pencil

Equation (3.1) and `h=1+lambda H` give

\[
 H(g)-\frac1{\lambda u}g+\frac1\lambda=0.            \tag{5.1}
\]

Thus the Gaussian generating series follows the small formal root of the
original weighted pencil along the exact one-parameter slice

\[
 s=\frac1{\lambda u},\qquad t=\frac1\lambda.         \tag{5.2}
\]

As `u` tends formally to zero, `s` tends to infinity and `g=u+O(u^2)`.  The
bridge therefore marks one inverse branch near target infinity.  It does not
identify all roots simultaneously and does not by itself recover the full
discriminant, conductor, or node-pairing invariant.

Every admissible weighted seed in this repository satisfies `H(0)=0`, so the
construction applies to the canonical, split, and positive-dimensional
stable-moduli families.  It produces parameter families of Gaussian witnesses
whose coefficients vary algebraically with the seed.  No claim is made here
that distinct stable seed classes give inequivalent Gaussian witnesses; that
requires a separate equivalence theory for Gaussian pairs.

## 6. Exact recovery from the mixed moments

Put

\[
 M_m(H,\lambda)=
 \mathbb E\!\left(QP_{H,\lambda}^m\right).
\]

Equation (4.2) says

\[
 g(u)=\sum_{m\ge0}\frac{M_m(H,\lambda)}{m!}u^m,
 \qquad g=u h(g).                                    \tag{6.1}
\]

Here `M_0=E(Q)=0` and `M_1=1`, so `g(u)=u+O(u^2)` has a unique
compositional inverse
\(\eta(z)=g^{\langle-1\rangle}(z)\).  Substituting `u=eta(z)` into the
fixed-point equation gives

\[
 z=\eta(z)h(z),
 \qquad
 \boxed{\eta(z)=\frac{z}{h(z)},\quad
        h(z)=\frac{z}{\eta(z)}.}                     \tag{6.2}
\]

The reciprocal orientation in (6.2) is forced already by the linear example
`h(z)=1+az`, for which `g(u)=u/(1-au)` and
`eta(z)=z/(1+az)`.

### Theorem 6.1 — injective mixed-moment fingerprint

On the set of polynomials `h` with `h(0)=1`, the map

\[
 h\longmapsto
 \left(\mathbb E(QP_h^m)\right)_{m\ge0}              \tag{6.3}
\]

defined by the determinant-cancelling bridge is injective.  Explicitly, the
complete sequence determines `g` by (6.1), formal series reversion determines
`eta`, and (6.2) recovers `h` exactly.
The same coefficientwise reconstruction remains valid for formal `h`, though
then the bridge generally produces formal rather than polynomial `Phi`.

If `h` has degree at most `d`, only the finite list

\[
 M_1,M_2,\ldots,M_{d+1}                              \tag{6.4}
\]

is needed: these moments determine `g mod u^(d+2)`, hence
`eta mod z^(d+2)` and `z/eta mod z^(d+1)`, which contains every coefficient
of `h`.  Thus distinct polynomials `1+lambda H` have distinct mixed-moment
sequences.  For fixed nonzero `lambda`, this recovers `H`; if `lambda` is also
allowed to vary, the moments recover the product `lambda H`, not its two
factors separately without a normalization on `H`.  On the repository's
normalized admissible seed locus, however, `H'(1)=-1`, and therefore

\[
 \boxed{\lambda=-h'(1),\qquad
 H=\frac{h-1}{\lambda}.}                            \tag{6.5}
\]

Thus the map `(H,lambda) -> (M_m)_(m>=0)` is injective on normalized seeds
with `lambda!=0`.

The bound (6.4) treats all coefficients of an arbitrary polynomial `h` as
independent.  It is not minimal on the normalized seed slice.  The endpoint
conditions remove exactly two more coefficients, and the moment map itself
is triangular.

### Theorem 6.2 — optimal moment coordinates on the normalized seed space

Fix `N>=4` and `lambda!=0`, and write

\[
 H(z)=\sum_{j=2}^N c_jz^j,\qquad
 h(z)=1+\lambda H(z)=1+\sum_{j=1}^N a_jz^j.          \tag{6.6}
\]

On the normalized degree-`N` slice

\[
 H(0)=H'(0)=H(1)=0,\qquad H'(1)=-1,                 \tag{6.7}
\]

the `N-3` mixed moments

\[
 \boxed{M_3,M_4,\ldots,M_{N-1}}                    \tag{6.8}
\]

recover `H` uniformly.  More precisely, after factorial normalization
`mu_m=M_m/m!`,

\[
 \mu_{j+1}
 =\frac1{j+1}[z^j]h(z)^{j+1}
 =a_j+F_j(a_1,\ldots,a_{j-1}),                     \tag{6.9}
\]

where `F_j` is a universal polynomial.  Since `a_1=0`, the coordinates
`mu_3,...,mu_(N-1)` recover successively
`a_2,...,a_(N-2)`, hence `c_2,...,c_(N-2)`.  Put

\[
 A=-\sum_{j=2}^{N-2}c_j,
 \qquad B=-1-\sum_{j=2}^{N-2}j c_j.                \tag{6.10}
\]

The two endpoint equations then give

\[
 c_{N-1}=NA-B,\qquad c_N=B-(N-1)A.                 \tag{6.11}
\]

Thus (6.8) is a global polynomial coordinate system on the ambient
normalized affine slice.  Its Jacobian with respect to
`c_2,...,c_(N-2)` is triangular with determinant

\[
 \lambda^{N-3}                                     \tag{6.12}
\]

for the factorial-normalized moments, or a nonzero constant multiple of
this for the raw moments.  It restricts to a biregular coordinate chart on
the open degree-`N` admissible locus, including every ordinary
boundary-clean open used in the decorated-normalization theorem.

This number of scalar algebraic observations is optimal both generically
and uniformly.  Indeed, the normalized seed space has dimension `N-3`.  If
any `r` mixed moments—or any `r` algebraic scalar functions at all—gave a
generically injective rational map to affine `r`-space, the image would have
dimension `N-3`, forcing `r>=N-3`.  Uniform injectivity implies generic
injectivity, while (6.8) attains equality uniformly.

If `lambda` is allowed to vary, the parameter space has dimension `N-2` and
the optimal uniform list is

\[
 M_3,M_4,\ldots,M_N.                               \tag{6.13}
\]

These moments first recover `a_2,...,a_(N-1)`; the equation `h(1)=1`
recovers `a_N`, after which
`lambda=-h'(1)` and `H=(h-1)/lambda`.  The same dimension argument proves
minimality on the locus `lambda!=0`.

### Moment-image equations, sparse slices, and Padé form

The triangular inverse in (6.9) also gives explicit equations for the
nonminimal fingerprint.  Let

\[
 a_j=A_j(\mu_2,\ldots,\mu_{j+1})                   \tag{6.14}
\]

be the recursively recovered coefficient.  Inside the full jet space with
coordinates `mu_2,...,mu_(N+1)`, the fixed-`lambda` normalized fingerprints
are cut out by

\[
 A_1=0,\qquad \sum_{j=1}^N A_j=0,
 \qquad \sum_{j=1}^N jA_j=-\lambda,                \tag{6.15}
\]

together with the open conditions expressing exact degree and whatever
admissibility open is desired.  After projection to the minimal coordinates
(6.8), there are no residual closed equations: the ambient normalized slice
is affine `(N-3)`-space and admissibility only removes closed subsets.

There is also a determinantal presentation of the infinite degree bound.
Write

\[
 \eta=g^{\langle-1\rangle},\qquad
 e(z)=\frac{\eta(z)}z=\sum_{n\ge0}e_nz^n=\frac1{h(z)}.            \tag{6.16}
\]

For `n>=1`, set

\[
 T_n(e)=
 \begin{pmatrix}
 e_1&1&0&\cdots&0\\
 e_2&e_1&1&\ddots&\vdots\\
 \vdots&\vdots&\ddots&\ddots&0\\
 e_{n-1}&e_{n-2}&\cdots&e_1&1\\
 e_n&e_{n-1}&\cdots&e_2&e_1
 \end{pmatrix}.                                    \tag{6.17}
\]

The reciprocal-series identity gives

\[
 [z^n]\frac1{e(z)}=(-1)^n\det T_n(e)=a_n.          \tag{6.18}
\]

Consequently the moment sequence comes from a polynomial `h` of degree at
most `N` exactly when these Toeplitz--Hessenberg determinants vanish for all
`n>N`.  Equivalently, `e=1/h` satisfies the order-`N` constant-coefficient
recurrence determined by `h`, so its infinite Hankel matrix has rank at most
`N`.  This is the natural Padé/Prony avatar of the bridge.  It is useful when
the degree is unknown or moments are supplied as a long sequence; for known
normalized degree, the triangular coordinates (6.8) are strictly smaller.

Finally, if the coefficient support `S subset {2,...,N}` is known in advance,
has size `s`, and contains at least two indices, the corresponding normalized
slice has dimension `s-2`.  Order the support and use the moments
`M_(j+1)` for all but its two largest indices.  Triangular recovery finds the
earlier supported coefficients (the gaps are known zeros), and the two
endpoint equations recover the last two.  Hence `s-2` moments are again
uniformly sufficient and dimensionally minimal on every fixed-support
slice.  Recovering an *unknown* support uniformly is a separate sparse
identification problem; the fixed-support statement does not silently assume
that combinatorial step.

This is an injective realization at the level of exact moment sequences.  It
does not assert that the resulting polynomial pairs are inequivalent under
arbitrary transformations of the Gaussian variables or under any broader
notion of Gaussian-witness equivalence.

### Corollary 6.3 — minimal Gaussian fingerprints in positive-dimensional families

Fix `N>=4` and `lambda!=0`.  On the normalized admissible degree-`N` seed
space, the bridge is an algebraic family of explicit `GMC(4)` witnesses, and
the finite moment vector

\[
 \mathfrak f_N(H)=
 \bigl(M_3(H,\lambda),\ldots,M_{N-1}(H,\lambda)\bigr)             \tag{6.19}
\]

is injective and has the smallest possible number `N-3` of scalar algebraic
coordinates.  This is Theorem 6.2.  Each entry of (6.19) is polynomial in the
coefficients of `H`, by (1.6).

The decorated-normalization theorem supplies a nonempty ordinary
boundary-clean seed open whose stable-moduli image has dimension `N-3`.
Restricting the bridge to that open therefore gives

\[
 \boxed{\text{an explicit algebraic }(N-3)\text{-dimensional family of
 four-real-Gaussian GMC witnesses separated optimally by }N-3
 \text{ mixed moments}.}                                      \tag{6.20}
\]

The separation in (6.20) is separation of the seed parameters by a finite
observable vector.  It is deliberately not phrased as inequivalence of the
Gaussian polynomial pairs under an unspecified transformation group.
Christopher D. Long's direct witness uses three real Gaussian variables and
is therefore dimensionally smaller.  The point of (6.20) is instead the
positive-dimensional, finitely fingerprinted family; the optimality claim
concerns the number of scalar algebraic moment coordinates, not the number of
underlying Gaussian variables.

## 7. Three-real-variable obstruction through quadratic real chaos

For the more economical ansatz

\[
 P=Wh(Z)+v(Z)T^2
\]

with one circular complex Gaussian and one real Gaussian, the half-pair
formula forces

\[
 v(z)=
 -\frac{h(z)h'(z)(2h(z)-zh'(z))}{2D(z)^2}           \tag{7.1}
\]

if the pure-moment generating function is to equal one.  For Long's affine
choice `h=1+z`, one has `D=1`, and (7.1) is exactly

\[
 v=-\tfrac12(1+z)(2+z).
\]

The same denominator is forced for the **most general separated correction
of real-Gaussian degree at most two**,

\[
 P=Wh(Z)+c(Z)+a(Z)T+b(Z)T^2.                        \tag{7.2}
\]

Indeed, the one-variable Gaussian--Lagrange identity reduces pure-moment
cancellation to

\[
 \mathbb E_T\!\left[
  \exp\!\left(\frac z{h(z)}V(z,T)\right)
 \right]=\frac{D(z)}{h(z)},                         \tag{7.3}
\]

where `V=c+aT+bT^2`.  Exact Gaussian integration and the fact that the
exponential of a nonconstant rational function cannot be algebraic force

\[
 b=-\frac{hh'(2h-zh')}{2D^2},
 \qquad c=-\frac{za^2D^2}{2h^3}.                   \tag{7.4}
\]

For every nonlinear polynomial `h` with `h(0)=1`, the first expression in
(7.4) is not a polynomial.  To see this uniformly, `D` has degree `deg(h)`.
If every root of `D` were a root of `h`, it would be a multiple root; a root
of multiplicity `m` in `h` has multiplicity exactly `m-1` in `D`.  The total
of these multiplicities is strictly below `deg(h)`, a contradiction.  Hence
`D` has a root `beta` with `h(beta) != 0`.  At that root
`h'(beta)=h(beta)/beta != 0` and `2h(beta)-beta h'(beta)=h(beta) != 0`, so even
`D` does not divide the numerator in (7.4).

Thus no nonlinear weighted seed admits a polynomial correction of the form
(7.2).  This strictly strengthens the earlier seed-by-seed half-pair check,
but it is still not an obstruction to every polynomial in `(Z,W,T)`.
Higher powers of `T`, corrections involving `W`, and the obstruction to
importing a second Lagrange--Good fixed-point coordinate into an odd Gaussian
dimension are analyzed in the
[three-real-variable search note](THREE_REAL_WEIGHTED_GAUSSIAN_SEARCH.md).
The unrestricted three-real polynomial problem remains open.

## 8. Reproduction and status

Run

```bash
.venv/bin/python scripts/verify_formal_gaussian_lagrange.py
.venv/bin/python scripts/verify_weighted_gaussian_bridge.py
.venv/bin/python scripts/verify_gaussian_moment_fingerprint.py
.venv/bin/python scripts/search_three_real_weighted_ansatz.py \
  --z-degree 1 --w-degree 1 --moment-order 5
python3 scripts/audit_weighted_gaussian_bridge_independent.py
```

The first checker tests the standalone identity through `u^6` for a nonlinear
two-variable `Phi` with two nonzero constant terms.  The bridge checker then
verifies:

- the parameterized fixed point and pencil identity;
- the determinant cancellation symbolically;
- exact Wick moments through `m=12` for canonical cubic, quartic, and quintic
  seeds and a split quartic seed;
- the Lagrange coefficient formula; and
- the original half-pair divisibility regression for the audited nonlinear
  seeds (now subsumed by the all-nonlinear proof in Section 7).

The fingerprint checker directly reconstructs a symbolic quartic `h` by
exact compositional reversion and audits a concrete degree-five instance.  It
then verifies the optimal normalized coordinates `M_3,...,M_(N-1)` and their
Jacobian in degrees four through eight, the variable-scale `N-2` bound, and
the Toeplitz--Hessenberg determinant equations.

The three-real search checker constructs exact pure- and mixed-moment ideals
for finite polynomial supports involving `(Z,W,T)`.  Its unit-ideal results
are exact exclusions of the specified supports; nonunit truncations are
reported only as survivors, not witnesses.  The companion search note records
the all-orders separated-quadratic theorem, the odd-dimensional polarization
obstruction for the existing determinant architecture, and the remaining
higher-chaos frontier.

The second checker reconstructs the correction using only sparse `Fraction`
arithmetic, independently solves the fixed branch through order 14, and
replays the determinant and Wick moments without SymPy or project imports.

The theorem is **proved in the repository** and **exactly computationally
checked in bounded ranges**.  It is **not externally specialist-reviewed**.
The standalone formal Gaussian--Lagrange lemma is the priority target for such
review.
Christopher D. Long is credited for the Gaussian Lagrange--Good search
architecture; his direct three- and four-variable witnesses remain
independently authored external results and are not reclassified as weighted
seed constructions.
