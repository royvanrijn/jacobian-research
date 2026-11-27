# Simultaneous mirrored merger limits and permutation matchings

When finitely many mirrored factorial blocks approach complete
correlation at comparable rates, the scaled angular polynomial
is a finite free additive convolution of the individual limits.
It is also the expected weighted matching polynomial of a union
of independent colored permutation graphs. The limiting roots are simple
and real in every degree. The first four Euclidean pivots of every
consecutive pair are positive, uniformly in the number of blocks and
the starting degree. This proves consecutive interlacing through upper
degree six. The comparison in unrestricted degrees for arbitrary blocks
remains open.

## The simultaneous limit

For each block \(i\), let \(a_i\ge0\) be its smaller exponent
and \(d_i\ge1\) its positive gap, both integers. Set
\(b_i=2a_i+d_i\). Start with independent standard complex
Gaussians \(Z_i,W_i\), with variance one in each complex
coordinate, and let
\[
 X_i=|Z_i|^2,\quad
 Y_i=\left|\sqrt{1-\kappa_i\delta}\,Z_i+
             \sqrt{\kappa_i\delta}\,W_i\right|^2,
 \qquad \kappa_i\ge0 .
\]
At least one \(\kappa_i\) is positive. Take \(\delta>0\)
small enough for all square roots to be real. The two monomials
are
\[
 U_\delta=\prod_i (X_iY_i)^{a_i}X_i^{d_i},\qquad
 V_\delta=\prod_i (X_iY_i)^{a_i}Y_i^{d_i}.
\]
For \(m\ge1\) put \(P_m(z;\delta)=\mathbb E(U_\delta+zV_\delta)^m\).
Use \(y^{\underline{k}}\) for a falling factorial and set
\(h_i=\kappa_i d_i^2\).

**Theorem.** Locally uniformly in complex \(u\),
\[
 \frac{P_m(-e^{i\sqrt{\delta}u};\delta)}
 {(-i)^m\delta^{m/2}\prod_i\Gamma(mb_i+1)}
 \longrightarrow
 p_m(u)=\sum_{k=0}^{\lfloor m/2\rfloor}
  \frac{(-1)^k m^{\underline{2k}}}{k!}
        \mathbb E\left(\sum_i\frac{h_i}{G_i}\right)^k u^{m-2k},
                                                               \tag{1}
\]
where the \(G_i\) are independent, unit-scale Gamma variables
of shapes \(mb_i+1\). The polynomial \(p_m\) is monic, with
\(m\) distinct real roots.

For the limit, \((Y_i-X_i)/\sqrt\delta\) converges in every
finite \(L^p\) to \(\sqrt{\kappa_i X_i}\,H_i\), where the
\(H_i\) are independent real Gaussians of variance two and
independent of the \(X_i\). This follows by expanding the
squared modulus and using rotational symmetry of \(W_i\).
Let \(M=\prod_i X_i^{b_i}\). Factoring each difference
\(Y_i^{d_i}-X_i^{d_i}\) and telescoping the products gives
\[
 \frac{U_\delta-e^{i\sqrt\delta u}V_\delta}{-i\sqrt\delta}
   \longrightarrow
 M\left(u-i\sum_i\frac{\sqrt{\kappa_i}\,d_i H_i}{\sqrt{X_i}}\right)
\]
in every finite \(L^p\), uniformly for \(u\) in compact
sets. The apparent inverse powers on the right are
multiplied by \(M\), and \(b_i\ge1\). The finite polynomial
products and the exponential Taylor remainder justify
passage to the \(m\)-th moment by Hölder's inequality.
Conditioned on the \(X_i\), the real sum in this expression
has variance \(2\sum_i h_i/X_i\). Integrating its even
moments and tilting the independent \(X_i\) by \(M^m\)
proves (1).

To see simplicity, write \(\partial=d/du\). For each block
the constant-coefficient operator
\[
  \sum_{k=0}^{\lfloor m/2\rfloor}
          \frac{(-h_i)^k}{k!\,(mb_i)^{\underline{k}}}\partial^{2k}
\]
takes \(u^m\) to that block's polynomial. The operator product
takes \(u^m\) to \(p_m\) by the multinomial expansion in (1).
The product is precisely the degree-\(m\) symmetric additive
convolution. Every active single-block polynomial is a scaled
Gegenbauer polynomial of shape \(1+(b_i-1)m\), with simple
real roots. An inactive block is the identity \(u^m\).
[Single-block merger limit](SFC_PHYSICAL_MIRROR_MERGER_EDGE.md#the-angular-scaling-limit-for-an-ordinary-block)
[Symmetric additive convolution and Walsh's real-rootedness theorem](https://doi.org/10.1007/s00440-021-01105-w)

For completeness, the linear convolution with a fixed real-rooted
degree-\(m\) factor preserves real-rootedness also in degrees
below \(m\), by adjoining real roots that escape to infinity
and taking a coefficient limit. For every positive shift
smaller than the minimum gap of a simple real-rooted input,
that input and its shifted copy interlace. All their real
linear combinations are real-rooted. Convolution is linear,
commutes with translation, and preserves real-rootedness, so
the two output translates weakly interlace by the same
interlacing characterization. If the output had a repeated
root, a shift avoiding its finitely many root differences
would violate interlacing between those two equal roots.
Starting with an active block therefore gives simplicity.

For each fixed \(m\), the simple scaled roots account for all
\(m\) zeros at sufficiently small positive \(\delta\), by
Rouché's theorem on small disjoint disks in \(u\).
The scaled angular roots are real by the exact circle-location
theorem. Strict interlacing of \(p_n,p_{n+1}\), when available,
gives strict angular interlacing in a neighborhood of the
merger for that fixed pair of orders.
[Physical circle location](SFC_PHYSICAL_MIRRORED_ROOT_CIRCLE.md)

## A matching-polynomial representation

Give block \(i\) an independent uniform permutation on \(mb_i\)
slots, of which \(m\) are distinguished and identified with
common labeled vertices. Whenever the permutation maps
one distinguished vertex to a different one, put an undirected
edge copy between them, of color \(i\) and weight \(h_i\).
Opposite arrows give parallel copies. Discard loops and edges
involving nondistinguished slots. Then \(p_m\) is the expected
matching polynomial of the resulting weighted multigraph,
with the convention that a matching of size \(k\) contributes
\((-1)^k u^{m-2k}\) times the product of its edge weights.

In a matching with \(k_i\) edges of color \(i\), the prescribed
arrows have distinct starts and ends. The probability that
all occur is \(\prod_i[(mb_i)^{\underline{k_i}}]^{-1}\).
For fixed counts \(k_i\) with sum \(k\), the number of labeled
oriented colored matchings is
\(m^{\underline{2k}}/\prod_i k_i!\). This proves the identity
using \(\mathbb E G_i^{-k_i}=1/(mb_i)^{\underline{k_i}}\).

For the consecutive comparison, uniform permutations can be
coupled by inserting \(b_i\) slots of color \(i\), one
distinguished. Deleting those slots and stitching their
cycles recovers the lower uniform permutation. Stitching
changes edges between surviving vertices, so this coupling
does not identify the lower graph with a vertex deletion
of the upper graph. The usual matching-polynomial interlacing
under vertex deletion does not directly compare these two
averages.

Simplicity and adjacent interlacing of separate inputs also
do not suffice for a convolution argument. For example, take
\(a=99/100\), \(c=1/100\), and the strictly interlacing pair
\(f_3=u(u^2-a^2)\), \(f_4=(u^2-1)(u^2-c^2)\). Adding two copies
at each respective degree gives \(u(u^2-2a^2)\) and
\[
 u^4-2(1+c^2)u^2+2c^2+(1+c^2)^2/6 .
\]
The quartic is positive at \(u=\sqrt{2}\,a\), by direct
rational substitution. This precludes strict interlacing.
The example is an abstract convolution pair. Formula (1)
imposes additional inverse-Gamma constraints in the ordinary
merger limit.

## Exact adjacent slot shifts

There is an exact correction to derivative compatibility. For slot
parameters \(T_i\ge j\), define the monic polynomial
\[
 F_j(u;T,h)=\sum_{2|k|\le j}(-1)^{|k|}j^{\underline{2|k|}}
       u^{j-2|k|}\prod_i\frac{h_i^{k_i}}{k_i!T_i^{\underline{k_i}}}.
 \tag{2}
\]
The physical merger limit is \(p_j=F_j(u;jb,h)\). Put \(N=n+1\),
and \(E_i=h_i\partial_{h_i}\). Coefficient comparison proves
\[
 \frac1N p_N'(u)=F_n(u;Nb,h).
 \tag{3}
\]
For integer \(b_i\ge1\), its relation to the actual lower polynomial is
\[
 \boxed{\quad \frac1N p_N'
 =\prod_i\prod_{r=0}^{b_i-1}
       \left(1-\frac{E_i}{Nb_i-r}\right)p_n .\quad}
 \tag{4}
\]
Indeed, on a monomial of degree \(k_i\) in \(h_i\), the multiplier is
\[
 \frac{(Nb_i-k_i)^{\underline{b_i}}}{(Nb_i)^{\underline{b_i}}}
       =\frac{(nb_i)^{\underline{k_i}}}{(Nb_i)^{\underline{k_i}}}.
\]
Every eigenvalue in (4) is positive on the required finite degree space.
This statement alone does not order the roots.

For real \(b_i\ge1\), let the independent variables
\(Z_i\sim\operatorname{Beta}(nb_i+1,b_i)\). The complementary identity is
\[
 p_n(u)=\mathbb E\,F_n(u;Nb,(h_i/Z_i)_i).
 \tag{5}
\]
Its proof is the beta integral
\(\mathbb E Z_i^{-k}=(Nb_i)^{\underline k}/(nb_i)^{\underline k}\).
All these moments exist since \(k\le n/2<nb_i+1\).

The closest classical input is the degree-reduction identity for finite
free additive convolution: [Marcus, Spielman and Srivastava,
arXiv:1504.00350v2, Lemma 1.15 and its proof on printed page 9](https://arxiv.org/pdf/1504.00350v2).
It differentiates a fixed convolution factor and keeps its slot parameters
fixed. Equations (4) and (5) account for the additional change from
\(Nb_i\) to \(nb_i\). Thus that classical identity does not by itself
prove adjacent physical interlacing. This is a comparison of hypotheses,
not a claim that the general multiblock question is absent from the
literature.

## Four positive Euclidean pivots at every start

The following algebraic assertion allows any finite number of blocks,
real \(b_i\ge1\), and \(h_i\ge0\) with at least one \(h_i>0\).
It includes all the physical integer blocks above.

**Theorem.** Start Euclidean division with the monic polynomials
\(p_N,p_n\), replacing each negative remainder by its monic version.
For \(1\le k\le\min(4,n)\), the unnormalized negative remainder
has degree \(n-k\) and strictly positive leading coefficient \(\beta_k\).
In particular \(\beta_1=H\), \(\beta_2\ge H\), where
\(H=\sum_i h_i/b_i>0\). There is no bound on \(n\) or the block count.

Here is a complete finite polynomial certificate for this infinite
parameter statement. Put \(w_i=h_i/b_i\), and define
\[
 \begin{split}
 J_m&=\sum_i\frac{w_i^2}{mb_i-1},\qquad
 K_m=\sum_i\frac{w_i^3}{(mb_i-1)(mb_i-2)},\\
 L_m&=\sum_i\frac{w_i^4(5mb_i-6)}
                    {(mb_i-1)^2(mb_i-2)(mb_i-3)},\\
 c_m&=\frac{(m-1)(m-2)(m-3)}{2m},\quad
 e_m=\frac{\prod_{r=1}^5(m-r)}{6m^2},\quad
 f_m=\frac{\prod_{r=1}^7(m-r)}{24m^3},\\
 B_m&=c_m(H^2+J_m),\\
 C_m&=e_m(H^3+3HJ_m+4K_m),\\
 D_m&=f_m(H^4+6H^2J_m+3J_m^2+16HK_m+6L_m).
 \end{split}                                                     \tag{6}
\]
Use \(B_m\) only for \(m\ge2\), \(C_m\) for \(m\ge3\),
and \(D_m\) for \(m\ge4\), so all displayed denominators are positive.
The finite sums in (2) show that the first descending coefficients are
\[
 p_m=u^m-(m-1)Hu^{m-2}+B_mu^{m-4}-C_mu^{m-6}+D_mu^{m-8}+\cdots.
 \tag{7}
\]
A coefficient below degree zero is zero: the products in (6) supply
exactly that zero in their stated ranges. One can derive (6) directly
from (2) by grouping equal indices. Equivalently, the first four
cumulants of \(\sum_i h_i/G_i\), multiplied by \(m^r\) in order
\(r\), are \(H,J_m,4K_m,6L_m\), respectively.

Write the Laurent expansion at infinity as
\[
 \frac{p_n}{p_N}=u^{-1}+\mu_2u^{-3}+\mu_4u^{-5}
                        +\mu_6u^{-7}+\mu_8u^{-9}+\cdots.
\]
Then coefficient division gives
\[
 \begin{split}
 \mu_2&=H,\qquad \mu_4=nH^2-(B_N-B_n),\\
 \mu_6&=C_N-C_n+nH\mu_4-B_NH,\\
 \mu_8&=-(D_N-D_n)+nH\mu_6-B_N\mu_4+C_NH.
 \end{split}                                                     \tag{8}
\]
Set
\[
 A_2=\mu_4-H^2,\quad A_3=H\mu_6-\mu_4^2,\quad
 A_4=(\mu_4-H^2)\mu_8-\mu_6^2+2H\mu_4\mu_6-\mu_4^3.
 \tag{9}
\]
Expanding four levels of the formal continued fraction, or performing
four monic divisions using (7), gives
\[
 \beta_1=H,\qquad \beta_2=A_2/H,\qquad
 \beta_3=A_3/(HA_2),\qquad \beta_4=HA_4/(A_2A_3).
 \tag{10}
\]
These identities use no prior interlacing assumption. Each is used only
after the preceding positive denominator has been established.

For the second pivot put \(\alpha=c_N-c_n\) and
\(\delta_b=c_N/(Nb-1)-c_n/(nb-1)\). For \(n\ge2,b\ge1\),
\[
 \begin{split}
 \alpha&=n-\tfrac52+\frac3{nN},\qquad
 0\le\delta_b\le\delta_1=\tfrac12-\frac3{nN},\\
 \delta_1-\delta_b
   &=\frac{(b-1)(n-2)(bn+3b-2)}{2(bn-1)(bN-1)}.
 \end{split}
\]
The nonnegative lower bound follows directly from
\[
 \delta_b=\frac{(n-2)(n-1)(bn^2+5bn+3b-2n-3)}
                    {2nN(bn-1)(bN-1)}.
\]
Since \(\sum_i w_i^2\le H^2\), these formulas give
\(A_2=(n-1-\alpha)H^2-\sum_i\delta_{b_i}w_i^2\ge H^2\).

The remaining two signs have degrees four and six in the nonnegative
weights \(w_i\). To check the sign for arbitrarily many blocks, it
suffices to check every exponent partition of these fixed degrees.
The parameter at an index with weight exponent one drops out: the
nonlinear sums \(J,K,L\) use that index at least twice. Thus at most
two or three independent \(b_i\)'s occur in any coefficient, respectively.
For the coefficient of a particular monomial, factor its rational
function, and substitute \(n=3+x\) for \(A_3\), \(n=4+x\) for
\(A_4\), and \(b_i=1+y_i\) for its repeated indices.
The following table records the full numerator expansions. Every
nonzero coefficient of both numerator and denominator is positive,
and every numerator has the positive constant shown.

| Polynomial | Weight exponents | Numerator terms | Least coefficient | Constant |
|---|---|---:|---:|---:|
| \(A_3\) | 1111 | 5 | 62 | 3240 |
| \(A_3\) | 211 | 21 | 27 | 9936 |
| \(A_3\) | 22 | 81 | 22 | 59616 |
| \(A_3\) | 31 | 35 | 21 | 1728 |
| \(A_3\) | 4 | 49 | 12 | 864 |
| \(A_4\) | 6 | 144 | 144 | 4478976 |
| \(A_4\) | 51 | 99 | 312 | 5557248 |
| \(A_4\) | 42 | 405 | 2160 | 6652772352 |
| \(A_4\) | 411 | 117 | 2544 | 566009856 |
| \(A_4\) | 33 | 350 | 2556 | 6648293376 |
| \(A_4\) | 321 | 195 | 3228 | 5591669760 |
| \(A_4\) | 3111 | 55 | 3648 | 468115200 |
| \(A_4\) | 222 | 405 | 3544 | 221891788800 |
| \(A_4\) | 2211 | 117 | 4056 | 18543974400 |
| \(A_4\) | 21111 | 33 | 4558 | 1541376000 |
| \(A_4\) | 111111 | 9 | 5095 | 127656000 |

For an explicit short reconstruction of the five degree-four
coefficients, put
\[
 \begin{split}
 a&=(n-1)\alpha-\alpha^2-c_n+(e_N-e_n),\\
 q_b&=(n-1-2\alpha)\delta_b-c_n/(nb-1)
                         +3(e_N/(Nb-1)-e_n/(nb-1)),\\
 t_b&=e_N/((Nb-1)(Nb-2))-e_n/((nb-1)(nb-2)).
 \end{split}
\]
In the table's order they are
\(24a,12a+2q_b,6a+q_b+q_c-2\delta_b\delta_c,
4a+2q_b+4t_b,a+q_b-\delta_b^2+4t_b\).
For degree six, equations (6), (8), and (9) specify each coefficient
without any search bounds. Expansion and coefficient extraction produce
the eleven listed polynomials. Since some \(w_i>0\), the pure-power
coefficient already makes \(A_3>0\) and \(A_4>0\). Equations (10)
complete the proof. This is an exact computational polynomial-positivity
proof, rather than a claim inferred from sampled parameter values.

## Consequences and the remaining comparison

**Corollary.** The merger limit polynomials \(p_n,p_{n+1}\) strictly
interlace for every \(1\le n\le5\), for arbitrary finitely many blocks
in the parameter domain above.

For \(n\le4\), the available Euclidean steps reach a nonzero constant.
For \(n=5\), four steps bound the gcd degree by one. The two parents
have opposite parity, and the even parent's constant term in (1) is
nonzero. Any common root would therefore be nonzero and come with its
negative, contradicting that gcd bound.

At each fixed degree the existing convolution proof of simple real
roots also holds for real \(b_i\ge1\): the single-block Gegenbauer
shape is \(1+(b_i-1)m>0\). Move the nonnegative weights continuously
to a single active block, then move that block's \(b\) to one, always
retaining a positive weight. No common roots occur along this path at
the stated orders. At the endpoint the polynomials are the scaled
Chebyshev polynomials \(h^{m/2}U_m(u/(2\sqrt h))\), which strictly
interlace by their cosine zeros. Simplicity and the absence of a common
root preserve that order along the path.

Consequently, for any fixed finite physical block data and any one of
these five starts, the original correlated factorial moments strictly
interlace for all sufficiently small positive \(\delta\). This follows
from the locally uniform scaling limit and its simple zeros. The
permitted \(\delta\) depends on the block data and the start. This does
not prove interlacing at independence, at arbitrary starts, or a uniform
factorial tail theorem near an exceptional Gaussian contact.

At general \(n\ge4\), the four pivots give only
\(\deg\gcd(p_n,p_{n+1})\le n-4\). Continuing the positive Jacobi
fraction, or proving all of its Hankel determinants positive uniformly,
is still an obligation. Positivity of the sixteen polynomials above
does not imply positivity at later pivot orders. Even full merger-limit
interlacing would require a further argument to reach unrestricted
ordinary linked factorial arrays, full dimension-two SFC, or support
three.

## Exact replay of the pivot certificates

Run from the repository root, using Python 3 and SymPy 1.14.0:

```sh
.venv/bin/python research/factorial/image-rigidity/verify_multiblock_merger_euclidean_pivots.py \
  --output /tmp/multiblock-merger-pivots.json
```

The checker reconstructs all sixteen rational coefficient functions
symbolically, checks every coefficient of their shifted numerators and
denominators, and records canonical coefficient-list hashes in the
[receipt](multiblock_merger_euclidean_pivots_verification.json).
The fixed formulas are its deterministic exact inputs. It also checks
the second-pivot bound, compares fixed-slot differentiation and all four
pivot formulas with independent rational polynomial division at 64 pairs,
and rejects a deliberately misordered real-rooted pair. Those 64 controls
test the implementation; they are separate from the unrestricted symbolic
certificates. The current replay is internal. There is no formal
verification, independent replay, or external mathematical review.
