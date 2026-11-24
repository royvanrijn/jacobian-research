# Native-prime Hasse coprimality at arbitrary offsets

Let \(p\) be a prime, \(n,h\ge1\) integers and \(\gamma>0\) rational
and integral at \(p\). Write
\[
 T_{m,\beta,p}(t)=\sum_{k=0}^{m}\binom mk(\beta)_{pk}t^k,\qquad
 C=\gamma+p(n+1),\qquad \kappa=\frac{p(n+1)}{(C)_h}.
\]
Define \(\sigma_h(x)=[Y^2](x+Y)_h\). The actual pair
\(T_{n,\gamma,p},T_{n+1,\gamma+h,p}\) is coprime over
\(\mathbb Q[t]\) if either
\[
             v_p(\kappa)\ge-1\qquad\text{or}\qquad
             \sigma_h(\bar\gamma)\ne0\text{ in }\mathbb F_p .
\]
This extends the [native-prime criterion](PRIME_PLUS_ONE_EXPONENT_GAP_SFC.tex),
Theorem \(\mathrm{thm:criterion}\), on rational integral shapes by
allowing arbitrary \(h\). The reduction below uses its exact
upper identity and reduced rising-jet kernel. The conclusion
extends to every characteristic-zero field.

## 1. Reducing higher rising jets

Use the notation of the native-prime source,
\[
 \epsilon=(-1)^{p+1},\quad
 R=z^{pn}T_{n,\gamma,p}(\epsilon z^{-p}),\quad
 M=\gamma+pn,\quad E=M-zD_z,\quad K_j=(E)_jR.
\]
Its lower Euler equation and commutator identities are
\[
 (E-M)K_p=-\epsilon z^p(E-\gamma)R,\qquad
 Ez^p=z^p(E-p).
\]
For \(r\ge1\), divide \((E+p)_r\) at \(E=M\), writing
\[
 H_r(E)=\frac{(E+p)_r-(M+p)_r}{E-M}.
\]
This is an integral polynomial of degree \(r-1\).
The identities yield
\[
 K_{p+r}
   =(M+p)_r K_p-\epsilon z^p H_r(E-p)(E-\gamma)R.
 \tag{1}
\]
The last operator has degree \(r\). Expanding it in the monic rising
basis expresses (1) in jets of lower index, with integral
coefficients and a factor \(z^p\) on the correction.
Induction therefore expresses each \(K_j\), \(j>p\), at a root
\(R=0\) as a \(\mathbb Z_p[u]\)-linear combination of
\(K_1,\ldots,K_p\) with \(u=z^p\).
The specialization of that row at \(u=0\) is
\((M+p)_{j-p}e_p^{\mathsf T}\).

The lower monic reciprocal is integral. In fact \(v_p(u)\ge1\)
at a lower root, since the coefficient of index \(k\) in \(T_n\)
contains a rising product with \(k\) complete blocks modulo \(p\).
The native-prime reduced kernel theorem at this root gives a vector
\[
 k_j=\left.\frac{d}{dx}(x)_j\right|_{x=\bar\gamma},
 \qquad 1\le j\le p,\qquad k_p=-1
\]
spanning the kernel of the reduced cyclic concomitant \(Z(u)\).
The rank is \(p-1\). The full rising jet is nonzero, since the
lower equation has a regular leading coefficient at the chosen
\(z\ne0\). Normalize its components \(K_1,\ldots,K_p\) by an
actual component of minimum valuation, giving a primitive vector
\(w\) with \(\bar w\) a nonzero multiple of \(k\).

For \(j=p+r\), the reduced row in (1) evaluates on \(k\) as
\(-(\bar\gamma)_r\). This continues the derivative rule to all \(j\):
in characteristic \(p\), \((x)_{p+r}=(x^p-x)(x)_r\) and
\(\bar\gamma\in\mathbb F_p\).
Thus, with \(V_h(E)=((E)_h-(C)_h)/(E-C)\), let
\(v_h(u)\in\mathbb Z_p[u]^{1\times p}\) be the row expressing
\(V_h(E)R\) in the first \(p\) jets, omitting \(K_0=R\).
We obtain
\[
              \overline{v_h(0)}\,k
              =V_h'(\bar\gamma)=\sigma_h(\bar\gamma).
 \tag{2}
\]
The equality uses \(\bar C=\bar\gamma\) and the expansion of
\((x)_h\) around \(\bar C\).

## 2. Upper exclusion and descent

The exact upper identity in the native-prime source,
Lemma \(\mathrm{lem:upper}\), holds at all positive offsets and gives
\[
                    w_p=\epsilon\kappa\,u\,v_h(u)w
 \tag{3}
\]
at any common root.
Choose \(p-1\) independent rows \(I\) of the reduction of \(Z(0)\).
If \(v_p(\kappa)\ge-1\), set \(q=p\). Then \(\eta=u/q\) is integral
by the lower-root bound, and \(\lambda=\epsilon\kappa q\) is integral.
Otherwise, if \(\sigma_h(\bar\gamma)\ne0\), use
\(q=(\epsilon\kappa)^{-1}\) and \(\lambda=1\).
In this case (2) makes \(v_h(u)w\) a unit and (3) makes \(\eta=u/q\)
a unit. Both choices give \(v_p(q)>0\).
Form the polynomial
\[
 F(X)=\det\begin{bmatrix}
 Z_I(qX)\\ e_p^{\mathsf T}-\lambda Xv_h(qX)
 \end{bmatrix}\in\mathbb Z_p[X].
\]
At every common root the vector \(w\) lies in its matrix kernel at
\(X=\eta\). The rows \(I\) and the scaling \(q\) are fixed independently
of the root. The reduction is affine with nonzero constant term because
\(v_p(q)>0\) and \(k_p=-1\). If the linear coefficient vanishes, there
is no integral root. If the linear coefficient is nonzero, the
divided difference at any two distinct integral inputs is a unit.
Thus there is at most one integral root over a complete algebraic
closure. The map \(t\mapsto\eta=\epsilon/(qt)\) is
injective. All distinct roots of any rational common irreducible
factor occur in that closure. In characteristic zero at most one
such root forces the factor to be linear. But a common rational
root would be real, and the even-order member of the pair is a
strictly positive Gamma integral at every real \(t\). Rational
Bézout gives the stated conclusion.

## 3. Ordinary consequences

For the ordinary phase \(X^a(A+BX^p)\) use
\(\gamma=an+1\) and \(h=a\). Multiplication by
\(A^m(am)!\) gives the ordinary moment when \(A\ne0\) and
\(t=B/A\). A coefficient-zero phase is a monomial.

* For every prime \(p\), lower exponent \(a=2p\) is covered at
  all orders. In \(\mathbb F_p[Y]\) the product
  \((\bar\gamma+Y)_{2p}=(Y^p-Y)^2\) has second coefficient \(1\).
  Together with the
  [prime-power Hasse result](PRIME_POWER_HASSE_COPRIMALITY.md),
  this covers all \(a=2p,\ d=p^s\) with \(s\ge1\).
* For every prime \(p\ge5\), lower exponent \(a=p-1\) is covered
  at gap \(p\). For a residue \(x=\bar\gamma\), the product at
  offset \(p-1\) is
  \((Y^p-Y)/(Y+x-1)\). When \(x\ne1\), its second coefficient is
  \(1/(x-1)^2\). When \(x=1\), every factor in \((C)_{p-1}\)
  is a unit, so \(v_p(\kappa)\ge0\).
* If \(p\ge5\) and \(3\) is a nonsquare modulo \(p\), both offsets
  \(a=p-3\) and \(a=p+3\) are covered at gap \(p\).
  At \(p-3\), the missing product is
  \(D_x(Y)=\prod_{j=1}^{3}(Y+x-j)\).
  If \(D_x(0)=0\), every factor in \((C)_{p-3}\) is a unit.
  Otherwise the second coefficient of
  \((Y^p-Y)/D_x(Y)\) is \(D_x'(0)/D_x(0)^2\), where
  \(D_x'(0)=3(x-2)^2-1\) never vanishes.
  At \(p+3\), the product is \((Y^p-Y)(x+Y)_3\).
  Its second coefficient is \(-\big(3(x+1)^2-1\big)\),
  which likewise never vanishes.
