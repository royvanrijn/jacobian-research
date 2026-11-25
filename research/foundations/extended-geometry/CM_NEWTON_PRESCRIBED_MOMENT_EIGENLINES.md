# Prescribed sphere moments and vanishing Frobenius returns

At every fixed odd prime \(p\), each positive moment order prime to
\(p\) occurs in a family of supporting-face certificates with
unbounded connected cyclic-cover conductor. At the
certifying prime the primitive factor is nonordinary and the return
map on the character-one coherent line is zero. The certificate
instead detects a nonzero transition from a different character.
There is a smaller branch-exponent construction for some moment
orders at \(p\ge5\).

The coefficient calculation uses the usual cyclic-cover Hasse–Witt
matrix, as in [Bouw, §§4–5](https://doi.org/10.1023/A:1017513122376)
and [Lin–Mantovan–Singhal, §2.4](https://arxiv.org/html/2406.09632v1#S2.SS4).
The added ingredients are the sphere moment congruence, the
prescribed angular coefficients, and the identification of a
zero transition in the return orbit.

## A small branch-exponent angular polynomial

Fix a prime \(p\ge5\). Use a monic \(p\)-integral polynomial with
reduction
\[
 P_p(z)=\frac{(z^{p-1}-1)^2}{z-1}
       =(z-1)\prod_{\beta\in\mathbb F_p^\times,\ \beta\ne1}(z-\beta)^2
       \quad\text{in }\mathbb F_p[z],
\]
and with this same root multiplicity pattern before reduction.
For example lift the displayed distinct nonzero roots to their
integer representatives. Set \(b=2p-4\), so \(\deg P_p=b+1\).
The angular coefficient needed at moment order \(m\) is
\[
 c_m=[z^{bm}]P_p(z)^m
     =[T^m](1-T^{p-1})^{2m}(1-T)^{-m}.
\]
For \(1\le m<p-1\) it equals \(\binom{2m-1}{m}\) modulo \(p\).
In particular every prescribed \(m\ge1\) has \(c_m\ne0\) by
choosing \(p\ge5\) with \(p>2m-1\).

There are also unbounded \(m\) with \(c_m\ne0\) at a fixed \(p\).
Let \(G(z)=z^{-b}P_p(z)\). For \(m=1+p^t\) with \(p^t>b\),
over the residue field
\[
 \operatorname{CT}G(z)^{1+p^t}
 =\operatorname{CT}\bigl(G(z)G(z^{p^t})\bigr)
 =\bigl(\operatorname{CT}G(z)\bigr)^2=1.
\]
Indeed the first angular support lies in \([-b,1]\), so an
angle from it cannot cancel a nonzero multiple of \(p^t\).
More generally take
\[
 m=\sum_{\ell=0}^{v} m_\ell p^{t_\ell},\qquad
 t_0=0,\qquad 1\le m_\ell\le(p-1)/2,
\]
with gaps \(p^{t_{j+1}}>b\sum_{\ell\le j}m_\ell p^{t_\ell}\).
In the characteristic-\(p\) factorization by these digits,
the last factor with nonzero chosen angle cannot be canceled
by the earlier factors. The coefficient \(c_m\) is the product
of the nonzero \(c_{m_\ell}\). In all these cases \(p\nmid m\).

## Lift to phases and choose the source character

Fix any \(m\) with \(p\nmid m\) and \(c_m\ne0\). Put \(u=(p+1)/2\)
and choose an integer \(\ell\ge0\) and an odd integer \(s>1\)
with \(s\equiv1\bmod mp\) such that
\[
 \alpha=\frac{m p^\ell}{s}
 \quad\text{lies in}\quad
 I_p=\left(\frac{u}{p},\frac{u}{p}+\frac1{bp}\right).
\]
These choices exist. The interval of eligible \(s\) has endpoints
\(m p^\ell/(u/p+1/(bp))\) and \(m p^\ell/(u/p)\), and its
length tends to infinity with \(\ell\). It therefore contains a
member of the progression \(1\bmod 2mp\) for large \(\ell\).

Take arbitrarily large \(r\) with \(p^r\equiv1\bmod m\) and put
\[
                      n_r=\frac{s p^r-1}{m}.
\]
On the sphere \(t^2+2xy=1\) use the literal Laurent weight
\(W(x^it^j)=j+2s i\). Prescribe highest weight \(n_r\) and its
angular polynomial \(z^{-b}P_p(z)\), taking the stated
characteristic-zero lift of \(P_p\). Allow arbitrary
strictly smaller weights congruent to \(n_r\bmod s\), with
coefficients integral at the chosen prime of residue characteristic
\(p\). The weights and coefficients refer to the literal Laurent
expansion after substituting \(y=(1-t^2)/(2x)\).

These equality terms admit regular lifts. For a negative angle
\(-v\), the term of prescribed leading height is supplied by
\[
 (-1)^v x^{-v}t^{n_r}(1-t^{2s})^v
 =
 (-2)^v y^v t^{n_r}(1+t^2+\cdots+t^{2s-2})^v .
\]
All its smaller literal weights satisfy the congruence. For large
\(n_r\), the terms of positive angle also have nonnegative heights.
The equality lift has least sphere degree \(n_r+(2s-1)b\).
The full phase has that degree if the additional normal terms
have degree at most this value.

The [transition isolation lemma](CM_NEWTON_SELECTED_EIGENLINES.md)
at order \(m\) gives
\[
                s p^r\mathcal A(f^m)\equiv c_m\pmod p .
\]
Indeed an angular-constant height \(j\) satisfies
\(j+1\equiv0\bmod s\) and \(j+1\le s p^r\). The only such
denominator divisible by \(p^r\) is \(s p^r\).
Its coefficient comes from the equality terms and is \(c_m\).
Thus the prescribed moment is nonzero.

There is also a residue-level version of the support restriction.
Put \(d=\lfloor\log_p s\rfloor\). One may add arbitrary
strictly lower weights with coefficients in
\(p^{d+1}\mathcal O\), where \(\mathcal O\) is the local
coefficient ring. No congruence of their weights is required.
Every product using such a factor has a coefficient in
\(p^{d+1}\mathcal O\). Its interval denominator is at most
\(s p^r\), hence has \(p\)-valuation at most \(r+d\).
Its contribution after multiplication by \(s p^r\) is zero in
the residue field. The assertion is termwise and allows completely
arbitrary Laurent supports for these smaller terms, as long as
the resulting phase is regular on the sphere.

## A zero return in the selected orbit

For large \(r\), the auxiliary cover
\[
                        w^{n_r}=z^{n_r-b}P_p(z)
\]
is connected and tamely ramified at the certifying prime.
Its branch exponents are \(n_r-b,1,2,\ldots,2,n_r-1\), with
simple branch at \(z=1\) and doubled branches at the other nonzero
residue representatives. Their greatest common divisor is one,
and \(p\nmid n_r\). The primitive factor carries
\(\mathbf Q(\zeta_{n_r})\) of relative rank \(p-1\).
For a primitive character \(i\), the coherent dimension at
\(x=i/n_r\) is
\[
                  D(x)=1+\lfloor b x\rfloor-(p-2)\lfloor2x\rfloor .
\]
This follows by reducing \((-b,1,2,\ldots,2,-1)\) at that character.

For \(x\in I_p\), both \(D(x)\) and \(D(px-u)\) are one.
In fact \(1/2<x<1/2+1/b\) and \(0<px-u<1/b\).
The normalized one-step Frobenius polynomial at this \(x\) has
finite exponents
\[
 \lfloor p\{-bx\}\rfloor=1,\qquad
 \lfloor p\{x\}\rfloor=u,\qquad
 \lfloor p\{2x\}\rfloor=1 .
\]
Over \(\mathbb F_p\) it is consequently
\[
 z(z-1)^u\prod_{\beta\ne0,1}(z-\beta)
       =z(z-1)^{u-1}(z^{p-1}-1).
\]
The coefficient at degree \(p-1\) is zero. The factor
\(z(z-1)^{u-1}\) has degree \(u<p-1\) and no constant term.
Thus Frobenius on this directed pair of coherent lines is zero.

The primitive character represented by \(p^{r+\ell}\) has ratio
\[
               \frac{p^{r+\ell}}{n_r}
                =\frac{m p^\ell}{s-p^{-r}}\longrightarrow\alpha .
\]
This ratio is in \(I_p\), and in particular less than one, for all
sufficiently large \(r\). Hence that character is an actual
member of the character-one Frobenius orbit and supplies a zero
map in its return. This shows nonordinarity of the primitive factor.
The \(r\)-fold Frobenius from character \(s\) to character one
still sends the overlap class of \(w^s/z^s\) to the nonzero
multiple \(c_m[w/z]\). Thus neither bijectivity of full
coherent Frobenius nor a nonzero return on the detected line is
needed for this moment obstruction.

This result is local at the specified prime. It does not use
a density theorem or assert that the same return fails at
different primes of the characteristic-zero cover.

## Every prime-to-\(p\) moment at a fixed odd prime

Fix an odd \(p\) and a positive \(m\) with \(p\nmid m\). Choose a
power \(v=p^A>m\) and choose \(w\) with \(1\le w<v\) and
\(wm\equiv C\bmod v\). Use \(C=1\) for \(p\ge5\). For \(p=3\),
use \(C=1\) if \(m\equiv2\bmod3\) and \(C=2\) if
\(m\equiv1\bmod3\). These choices make \(w\equiv-1\bmod3\)
when \(p=3\). Replace the small-exponent polynomial by
\[
 P(z)=\frac{(z^{p-1}-1)^v}{(z-1)^w},\qquad
 b=v(p-1)-w-1
 \quad\text{over }\mathbb F_p .
\]
Use any fixed characteristic-zero lift of the distinct roots with
the indicated multiplicities \(v-w\) at 1 and \(v\) at every
other nonzero residue. The reciprocal coefficient is now
\[
 [T^m](1-T^{p-1})^{vm}(1-T)^{-wm}
   =[T^m](1-T)^{-C}\ne0\pmod p.
\]
Indeed both factors with exponents divisible by \(v\) have only
degrees divisible by \(v\) in characteristic \(p\), and \(v>m\).
The last coefficient is 1 for \(C=1\), and \(m+1\equiv2\bmod3\)
in the remaining case.

There is still an open interval giving a zero transition in a
character-one return orbit. Choose \(X\in\{1,\ldots,p-2\}\) such
that \(K=(-wX\bmod p)\) satisfies \(1\le K\le X\).
For \(p\ge5\), \(X=p-2\) works unless
\(w\equiv(p-1)/2\bmod p\). In the exception \(X=p-3\) gives
\(K=(p-3)/2\). For \(p=3\) take \(X=K=1\).
Put
\[
 x_0=\frac{X}{p}+\frac{1}{pv},\qquad c=\frac{v-w}{v}.
\]
At \(x_0\), the fractions \(\{v x_0\}\) and
\(\{(v-w)x_0\}\) are \(1/p\) and \((K+c)/p\).
Write
\[
 S(x)=(p-2)\{vx\}+\{(v-w)x\}-x .
\]
The coherent dimension is \(1+\lfloor S(x)\rfloor\), since
\(b=(p-2)v+(v-w)-1\). At \(x_0\),
\(pS(x_0)=p-2+K-X+c-1/v\), strictly between 0 and \(p\).
For \(x\) just to its right the source dimension is one. The
target ratio \(y=\{px\}\) approaches \(1/v\) from the right,
and \(S(y)\) approaches \(c-1/v\in[0,1)\). In the boundary
case \(c-1/v=0\), it starts increasing to the right, with
slope \(b>0\). Thus the target dimension is also one.

At such \(x\) the normalized finite exponents are
\[
 \lfloor p\{-bx\}\rfloor=1+X-K,\qquad
 \lfloor p\{(v-w)x\}\rfloor=K,\qquad
 \lfloor p\{v x\}\rfloor=1 .
\]
For the first identity, note that
\(p\{-b x_0\}=2+X-K-c+1/v\) and decreases to the right.
This also gives the asserted floor when \(c=1/v\).
The normalized Frobenius polynomial is
\[
 z^{1+X-K}(z-1)^{K-1}(z^{p-1}-1).
\]
Its degree-\(p-1\) coefficient is zero because the front factor
has degree \(X<p-1\) and vanishes at zero. This proves the desired
zero map on an open interval immediately right of \(x_0\).

Choose \(\ell,s\) using this interval in place of \(I_p\), and
apply the same lift with \(n_r=(s p^r-1)/m\). For large \(r\)
the order exceeds all the fixed interior exponents and the
highest Laurent angles have regular lifts. The connected cover
has the same angular coefficient \(c_m\ne0\), and \(p^{r+\ell}\)
lies in the zero-map interval for all sufficiently large \(r\).
This proves the prescribed-moment result for any \(p\)-unit order
with \(p\) fixed. In characteristic two the interval denominators
are all odd, so this form of prime-power isolation does not apply.

## Exact replay

The [standard-library checker](verify_cm_newton_prescribed_moment_eigenlines.py)
constructs allowed parameters and verifies the zero map and
the isolated coefficient at \((p,m)=(5,2),(5,26),(7,3),(11,5)\).
The connected orders are \(202,25024,590,2422\), with respective
selected entries \(3,1,3,5\). It checks the general multiplicity
construction at \((p,m)=(3,4),(3,5),(5,7),(5,11)\), with connected
orders \(218,1960,80158482,16031534\) and selected entries
\(2,1,1,1\). It also expands the second moment
of the first lift with a strictly lower off-congruence term.
The term \(5^3t^{172}\) preserves the predicted residue \(3\).
Replacing it by \(5^2t^{172}\) changes the normalized residue to
\(2\), since the cross term with height 202 has denominator
\(375=3\cdot5^3\). This exhibits the role of the precision level.
The replay concerns these explicit inputs. The parameter
arguments establish the unbounded families.

    python3 -I verify_cm_newton_prescribed_moment_eigenlines.py
