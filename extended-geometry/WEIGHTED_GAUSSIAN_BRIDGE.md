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

This is an injective realization at the level of exact moment sequences.  It
does not assert that the resulting polynomial pairs are inequivalent under
arbitrary transformations of the Gaussian variables or under any broader
notion of Gaussian-witness equivalence.

## 7. Why Long's three-real-variable correction does not transfer directly

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

For the audited nonlinear canonical and split weighted seeds, the numerator
in (7.1) is not divisible by `D^2`; the required `v` is rational rather than
polynomial.  This is an exact obstruction to the direct half-pair ansatz, not
an obstruction to every possible three-real-variable construction.

## 8. Reproduction and status

Run

```bash
.venv/bin/python scripts/verify_formal_gaussian_lagrange.py
.venv/bin/python scripts/verify_weighted_gaussian_bridge.py
.venv/bin/python scripts/verify_gaussian_moment_fingerprint.py
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
- the half-pair divisibility obstruction for those nonlinear seeds.

The fingerprint checker directly reconstructs a symbolic quartic `h`
from `M_1,...,M_5` by exact compositional reversion and also audits a concrete
degree-five instance.

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
