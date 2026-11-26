# Uniform zero-rank exclusion with a prime at the second order

Let \(f=ax^2+bxy+cy^2+k\) with \(abck\ne0\), and let
\(M_j=\Lambda_2(f^j)\), where \(\Lambda_2(x^i y^j)=i!j!\).
For every prime \(p\ge7\), the moments at \(p-1,p,p+1,p+2\)
cannot all vanish in characteristic zero.

The [negative](BINARY_EVEN_QUADRATIC_PRIME_FRONT_NEGATIVE_RANK_LIFT.md)
and [diagonal](BINARY_EVEN_QUADRATIC_PRIME_FRONT_DIAGONAL_LIFT.md)
lifts leave only the zero residue of \(t=b/2\) in the
normalization \(a=1\) from the
[prime-front reduction](BINARY_EVEN_QUADRATIC_PRIME_FRONT_REDUCTION.md).
Here the first digit in that residue forces both the
scaled rank correction and one of its factorial-polynomial
values to fixed rational residues. The next digit gives
an explicit nonzero rational obstruction.

## Two uniform layers beyond the initial reduction

Work over the algebraic \(p\)-adic closure with \(v(p)=1\).
Use the [zero-rank gate](BINARY_EVEN_QUADRATIC_PRIME_FRONT_ZERO_RANK_GATE.md)
at \(\ell=1\), and put
\[
 h=(p-1)/2,\quad C=(-1)^h,\quad d=c-t^2,\quad
 \kappa=k^p/p^2,\quad \alpha=p/k^h,\quad x=d/k,\quad z=t/\alpha .
\]
Then \(p=\kappa^h\alpha^p,\ k=\kappa\alpha^2\), \(x,z\)
are units, \(v(\alpha)=1/p\), and \(v(\kappa+2)\ge5/p\).
Write \(g_m(x)=\sum_i\binom mi(2i)!x^i\) and
\[
 H=g_h(x),\quad G=g_{h+1}(x),\quad G_1=g_{h+2}(x),\quad G_2=g_{h+3}(x).
\]
All expressions in polynomial values and derivatives below
are over the residue field, with integral lifts used when
forming expansions. This convention is valid at the retained
precision: the omitted coefficient errors have a full factor
\(p\), and the largest retained relative power of \(\alpha\)
is five, less than \(p\).

Let \(Q\) be any unit lift of \((-2)^h\bmod p\), and let
\(\theta=\overline{(\kappa+2)/\alpha^5}\).
The normalized moments have expansions
\[
 \begin{split}
 M_{p-1}/(-C Q\alpha^{p-1})&=H+\alpha F_-+o_v(1/p),\\
 M_p/(2p^2\alpha)&=E_0+\alpha F_0+\alpha^2 K_0+o_v(2/p),\\
 M_{p+1}/(4p^2\alpha^3)&=E_1+\alpha F_1+\alpha^2 K_1+o_v(2/p).
 \end{split}                                                 \tag{1}
\]
Here \(o_v(q)\) denotes terms of valuation strictly above \(q\).
The coefficients are
\[
 \begin{split}
 F_-&=C+z\,g_{h-1}(x),\\
 E_0&=z+2CG,\\
 F_0&=z^2+2Cz(G+2xG')+1+2x,\\
 K_0&=z^3+z(1+6x)+2Cz^2(G+5xG'+2x^2G'')+(2/3)C G_1,\\
 E_1&=z^3-z(1+6x)-(4/3)C G_1,\\
 F_1&=z^4-z^2(1+12x)-(4/3)Cz(G_1+2xG_1')
                                           -(1+4x+24x^2)/2,\\
 K_1&=z^5-z^3(1+20x)-z(1+12x+120x^2)/2\\
    &\qquad -(4/3)Cz^2(G_1+5xG_1'+2x^2G_1'')
                                  -(4/15)C G_2+\theta/2.
 \end{split}                                                 \tag{2}
\]

To verify (1)–(2), use the
[rank perturbation formula](BINARY_EVEN_QUADRATIC_PRIME_FRONT_NEGATIVE_RANK_LIFT.md#perturbation-formulas-at-the-negative-residue)
for \(f=(x_{\rm var}+t y)^2+d y^2+k\), where \(x_{\rm var}\)
denotes the original first variable. In order \(p+j\),
\(j=0,1\), pair each group \(0\le L\le j\) with the
group \(p+L\), and put \(R=2(j-L)\).
At the retained precision their sum, normalized by \(p^2\)
and after replacing the corresponding \(\kappa\) by \(-2\), is
\[
  2\binom jL R! k^L\sum_{i=0}^L\binom Li x^i
                         \sum_{m\ge R+1}(m+1)_{2i}t^m .
                                                               \tag{3}
\]
Indeed, the two order binomials reduce to \(\binom jL\),
the low inner binomials to \(\binom Li\), and
\((2p+R)!/p^2\equiv 2R!\).
Periodicity of the rising-factorial coefficients shows that
\(S_{2p+R,i}(t)\equiv(2i)!(1-t)^{-1-2i}\)
through \(t\)-degree \(2j+3<p\) in characteristic \(p\).
The later group has coefficient \(\kappa k^L R!\)
and exactly the initial terms \(S_{R,i}(t)\), which cancel
when \(\kappa=-2\).
The replacement error contributes only to \(K_1\):
at relative order \(\alpha^5\) it is \(2\theta\) before
division by \(4\), from \(L=0,i=0\) in the later group.

The critical group \(L=h+j+1\) contributes
\[
             \alpha^{2j+1}\beta_j
                  (1-t)^{-1}g_{h+j+1}(x/(1-t)^2),
 \qquad   \beta_0=4C,\quad \beta_1=-16C/3,                  \tag{4}
\]
to \(M_{p+j}/p^2\) at this precision.
This uses \(S_{p-1,i}\equiv(2i)!(t-1)^{p-1-2i}\),
\((p-1)!\equiv-1\), and
\[
      \binom{p+j}{h+j+1}/p
          \equiv \frac{Cj!}{\prod_{m=1}^{j+1}(h+m)} .
\]
The term from expanding \((1-t)^{p-1}\) through the
required \(t\)-degree equals the expansion of \((1-t)^{-1}\).
Next \(L=j+1\) has \(R=2p-2\), one factor \(p\) in
the order binomial and one in \(R!\).
It contributes \(\alpha^{2j+2}\) times
\[
 -\frac{(-2)^{j+1}}{j+1}
       \sum_{i=0}^{j+1}\binom{j+1}i x^i
                       ((2i)!+(2i+1)!\,t).
\]
The next critical group \(L=h+j+2\) gives
\(\alpha^{2j+3}\beta_j g_{h+j+2}(x)/(2j+3)\):
the consecutive order-binomial ratio is
\(h/(h+j+2)\), and the factorial ratio modulo \(p\)
is \(1/2\), together with the additional factor
\(-2\) from \(k/\alpha^2\).

These are all groups through relative degree \(2j+3\).
For \(j+1<L\le h+j\) the order binomial and
\(R!\) each have a factor \(p\), so the normalized
contribution starts at degree \(2j+4\).
For \(h+j+2<L<p\) the order binomial alone has
a factor \(p\) and the first possible degree is
\(2j+5\). If \(L\ge p\), inner indices outside the
low range \(i\le L-p\) have either an inner binomial
factor \(p\) or a rising-factorial factor \(p^2\).
Thus they lie above the displayed precision.
Finally, in \(M_{p-1}\) the critical order \(L=h\)
gives the leading \(H\) and the first expansion of
\((1-t)^{-1}H(x/(1-t)^2)\); its linear coefficient
is \(H+2xH'=g_{h-1}\) in characteristic \(p\).
The earlier \(L=0\) contributes the constant \(C\)
to \(F_-\). Positive earlier \(L<h\) have a factor
\(p\) in \(R!\) and at least two powers of \(\alpha\);
later \(L>h\) have at least \(p+1\) powers.
This proves (1)–(2).

## The first digit fixes two rational residues

Work first in characteristic \(p\).
The factorial recurrence and derivative identities give
\[
 xg_j'=j(g_j-g_{j-1}),\quad G=1+xg_{h-1},\quad
                         G_1=1+6xG-3xH.
\]
Put \(V=1+2x-4G^2\) and \(T=3GV-2+6xH\).
At a possible residue contact, \(H=0,\ z=-2CG,\ T=0\)
by the leading equations in (1). One has the polynomial identity
\[
          2xT'=3T+3H(12G^2-1-6x),                         \tag{5}
\]
as follows from \(2xG'=G-H\) and
\(2xH'=-H+(G-1)/x\).

Eliminate \(z\) through the first digit with
\(z=-2CG-\alpha F_0+o_v(1/p)\) from \(M_p\).
The constant layer of \(E_1\) at \(z=-2CG\) is
\((2C/3)T\); the coefficient from the change in \(z\)
is \(-\alpha D F_0\), where \(D=12G^2-1-6x\).
Here and below, the same expressions evaluated in the
valuation ring mean arbitrary integral lifts of the indicated
polynomials. Replacing a polynomial congruence by such a
lift incurs error of valuation at least one.
Choose the representative of the residue of \(x\) in the
finite unramified extension with the required residue field.
Vanishing reductions at this representative then give values
of valuation at least one.
If \(\bar G\ne1\), \(H'\) is a unit at the residue,
so \(M_{p-1}=0\) forces the change of \(x\)
from a representative to have valuation at least \(1/p\).
Then (5) gives \(v(T)>1/p\).
If \(\bar G=1\), \(H'=0\) but
\(H''=1/(4x^3)\) at the residue. This follows by
differentiating \(2xH'=-H+(G-1)/x\), and uses that
\(x\) is a unit. Taylor expansion of \(H=O(\alpha)\)
now gives a change in \(x\) of valuation at least \(1/(2p)\):
a smaller change would make the quadratic term uniquely
lowest. Equation (5) gives \(T'=T''=0\) at this
residue, so the cubic remainder again gives \(v(T)>1/p\).
Thus in both cases the \(M_{p+1}\) equation requires
\[
                              F_1-D F_0=0
\]
over the residue field after \(z=-2CG,\ H=0\).

Use \(3GV=2\), so \(G\ne0\) and
\(x=2G^2-1/2+1/(3G)\). Then the formulas (2) give
\[
     F_0=V,\quad
     F_1=16G^4+(16x-12)G^2+(32/3)G-1/2-2x-12x^2 .
\]
Substitution simplifies \(F_1-D F_0\) to \(-5/2+2/G\).
It follows that
\[
                         G=4/5,\qquad x=359/300             \tag{6}
\]
in the residue field. In particular the singular possibility
\(G=1\) is excluded. If \(p=359\), (6) contradicts that
\(x\) is a unit. In the other cases the displayed rational
representatives of \(x,G\) have denominators prime to \(p\).

## The second digit

We first determine \(\theta\) in (2) from the upper
weighted boundary \(p+2\). In the notation
\(a_2=(2p+4)!/(2p)!\), \(C_2(Z)=B_{p+2}(Z)/(a_2 Z^2)\),
the [two-band computation](BINARY_EVEN_QUADRATIC_PRIME_FRONT_ZERO_RANK_GATE.md#comparing-the-two-boundary-orders)
at \(Z=1\) gives
\[
         C_2(1)/p^2=\kappa+2+\mu_2 k^{h+3}/p+o_v(5/p),
                                \qquad \mu_2=-2C/45 .
\]
Here the only middle term through \(5/p\) is
index \(h\), and its coefficient is (2) of the
[diagonal lift](BINARY_EVEN_QUADRATIC_PRIME_FRONT_DIAGONAL_LIFT.md#boundary-and-perturbation-formulas).
The boundary equality is
\(C_2(1)=(c-t)B_{p+2}(c)/[(t-1)a_2]\).
Its right side divided by \(p^2\alpha^5\) reduces
to \(-z(1+4W+24W^2)/3\), where
\(W=\overline{c/k}=x-2G^2\).
This uses only the low exterior band
\(B_{p+2}(c)/k^{p+2}\), whose first three coefficients
give \(1+4W+24W^2\).
Thus at (6), where \(W=-1/12\) and \(z=-2CG\),
\[
                     \theta=-z(1+4W+24W^2)/3-16C/45
                            =4C/45.                        \tag{7}
\]

Choose rational representatives \((x_0,z_0)=(359/300,-8C/5)\).
At (6) the derivative \(H'=(G-1)/(2x^2)\)
is a unit. The equations for \(M_{p-1}\) and \(M_p\)
therefore determine the first corrections
\(\xi=\overline{(x-x_0)/\alpha}\),
\(\zeta=\overline{(z-z_0)/\alpha}\).
They satisfy
\[
 H'\xi+F_-=0,\qquad
           2CG'\xi+\zeta+F_0=0 .
\]
From (2) their values are
\[
        \xi=C\bigl(4xG-2x^2/(G-1)\bigr)=32669C/1800,
                         \qquad \zeta=-389/30.              \tag{8}
\]
Let \(v_1=(\xi,\zeta)\) and \(D=2-2/G\), the value
of \(12G^2-1-6x\) at (6). The leading gradients satisfy
\(\nabla E_1=D\nabla E_0\) there. Indeed the \(z\)
derivative gives \(D\), and the remaining \(x\)
derivative after eliminating \(z\) is \(T'\) times a
unit, which vanishes by (5).
Hence subtracting \(D\) times the normalized \(M_p\)
equation from the normalized \(M_{p+1}\) equation
removes the linear term in the unknown second correction.
At precision \(2/p\) its residue is
\[
 \left[K_1+\nabla F_1\,v_1+\tfrac12 v_1^{\mathsf T}
                     \nabla^2 E_1\,v_1\right]
  -D\left[K_0+\nabla F_0\,v_1+\tfrac12 v_1^{\mathsf T}
                     \nabla^2 E_0\,v_1\right].              \tag{9}
\]
The argument requires only the residues in (8), since
the difference of the leading gradients is zero and
any higher-valuation change in the first correction
enters (9) above the retained precision.

Evaluation of (9) gives \(-553C/72\).
One can use solely rational arithmetic in (2), (6)–(8),
with
\[
 \begin{gathered}
 G'=(G-H)/(2x),\qquad H'=-\bigl(H-(G-1)/x\bigr)/(2x),\\
 G_1=1+6xG-3xH,\quad G_1'=3(G_1-G)/(2x),\quad
 G_2=1+5x(4G_1-3G).
 \end{gathered}
\]
The second derivatives in (2) and (9) follow by
differentiating \(G',G_1'\).
The [exact replay](verify_binary_even_quadratic_prime_second_zero_rank.py)
performs this rational evaluation for both signs \(C\).
Thus (9) can vanish only at \(p=7\) or \(p=79\).
For these primes, the factorial coefficient sum at
the forced \(x=359/300\) gives respectively
\(g_3(5)=5\bmod7\) and \(g_{39}(22)=45\bmod79\),
contrary to the necessary \(H=0\).
The zero residue is therefore excluded at every \(p\ge7\).
