# Excluding the negative rank-one residue in prime-front quadratic windows

Let \(f=ax^2+bxy+cy^2+k\) with \(abck\ne0\), and let
\(M_j=\Lambda_2(f^j)\), where \(\Lambda_2(x^i y^j)=i!j!\).
Fix a prime \(p\ge7\). Suppose the moments at
\(p-1,p,p+1,p+2\), or those at \(p,p+1,p+2,p+3\),
all vanish. At each place over \(p\), use the normalization
\(a=1\) of the [prime-front reduction](BINARY_EVEN_QUADRATIC_PRIME_FRONT_REDUCTION.md).
Then \(b/2\) cannot reduce to \(-1\).
In particular, that reduction leaves only the residues \(0,1\)
for either prime placement. The argument is over the algebraic
\(p\)-adic closure, and applies to a characteristic-zero
solution by the algebraic reduction in the cited proof.

## The strength forced by the two weighted boundary equations

Assume \(t=b/2\) reduces to \(-1\) and put
\[
 h=(p-1)/2,\quad d=c-t^2,\quad s=v(k)=2/p,\quad
 \kappa=k^p/p^2,\quad r=v(c-1).
\]
Thus \(v(d)>0,\ r>0,\ \bar\kappa=-2\), where bars indicate
reduction in the residue field and \(v(p)=1\).
We may assume \(c\ne1\) by the
[equal-axis four-window theorem](BINARY_EVEN_QUADRATIC_SINGLE_FACE_STRONG_FACTORIAL.md#3-the-plane-ac).
Write
\[
 B_j(Z)=\sum_{i=0}^j\binom ji(2i)!\,Z^i k^{j-i},
 \qquad a_j=\frac{(2p+2j)!}{(2p)!},\qquad
 C_j(Z)=\frac{B_{p+j}(Z)}{a_j Z^j}.
\]
All \(a_j\) with \(0\le j\le3\) are \(p\)-units.
The recurrence
\[
 B_{m+1}(Z)=k^{m+1}+2(m+1)Z
           \bigl((2m+1)B_m(Z)-2m kB_{m-1}(Z)\bigr)
\]
follows by comparing coefficients. Define
\[
 J=\frac{4p(p+1)}{a_1},\quad
 b_j=\frac{4(p+j)(p+j+1)a_{j-1}}{a_{j+1}}\quad(j=1,2),
 \qquad K_j=\frac{k^{p+j+1}}{a_{j+1}} .
\]
Here \(v(J)=1\), the two \(b_j\) are units and
\(v(K_j)=2+(j+1)s\). The recurrence gives
\[
 C_0-C_1=J kB_{p-1}-K_0Z^{-1},
 \qquad
 C_{j+1}-C_j=-b_j k Z^{-1}C_{j-1}+K_j Z^{-j-1}.              \tag{1}
\]

Put \(g=(c-t)/(t-1)\), a unit with reduction \(-1\), and
\(\mathcal D_j(F)=F(1)-g c^j F(c)\). In the first window,
take \(\ell=1\), and in the second take \(\ell=2\).
The two terminal [weighted boundary cancellations](BINARY_EVEN_QUADRATIC_SINGLE_FACE_STRONG_FACTORIAL.md#4-the-exact-reduction-for-general-b-and-its-limitation)
are \(\mathcal D_\ell(C_\ell)=\mathcal D_{\ell+1}(C_{\ell+1})=0\).
Subtracting these zeros and using (1) gives
\[
 g c^\ell(c-1)C_\ell(c)
       =-b_\ell k\mathcal D_\ell(C_{\ell-1})+K_\ell(1-g).
                                                                  \tag{2}
\]
The unique lowest term of \(B_{p-1}\) on unit arguments
is at index \(h\), of valuation \(hs=1-1/p\).
The weighted difference \(\mathcal D_1(B_{p-1})\)
has the same valuation because its normalized reduction
is a unit times \(1-(-1)=2\). Consequently
\[
              v\bigl(\mathcal D_1(C_0-C_1)\bigr)=2+1/p .      \tag{3}
\]

For \(\ell=1\), (3) gives \(v(\mathcal D_1(C_0))=2+1/p\)
directly. For \(\ell=2\), write \(F=C_0-C_1\).
Equation (1) and \(\mathcal D_2(C_2)=0\) imply
\[
 \mathcal D_2(C_1)=b_1 k\mathcal D_1(C_0)-K_1(1-g),\quad
 \mathcal D_1(C_1)=\mathcal D_2(C_1)+g c(c-1)C_1(c).
\]
Initially \(v(\mathcal D_1(C_0))\ge2\) by its coefficient
table. Substituting the first equality in (2) at \(\ell=2\)
gives \(v((c-1)C_2(c))\ge2+4/p\).
Also \((1-b_1 k/c)C_1(c)=C_2(c)+b_1 kF(c)/c-K_1/c^2\).
The left coefficient is a unit and \(v(F(c))=2+1/p\).
Hence \(v((c-1)C_1(c))>2+1/p\).
The preceding difference equations now give
\(v(\mathcal D_1(C_1))>2+1/p\).
Combining with (3) shows \(v(\mathcal D_1(C_0))=2+1/p\)
also when \(\ell=2\). Returning to (2) gives in the two cases
\[
           v((c-1)C_\ell(c))=2+(2\ell+1)/p.                 \tag{4}
\]

For \(\ell=1,2\) the following expansion holds
coefficientwise for the Laurent polynomials:
\[
 \begin{split}
 C_\ell(Z)&=p^2 P_\ell(k/Z)(\kappa+2Z^p)+E_\ell(Z),\\
 P_\ell(Y)&=\sum_{j=0}^{\ell}
       \binom{\ell}{j}\frac{(2\ell-2j)!}{(2\ell)!}Y^j ,
 \qquad v([Z^i]E_\ell)\ge 2+(2\ell+1)/p.
 \end{split}                                                \tag{5}
\]
Indeed, the two exterior coefficient bands of \(B_{p+\ell}\)
are indexed by \(i=\ell-j\) and \(i=p+\ell-j\).
The binomial coefficients in these bands reduce to
\(\binom{\ell}{j}\) by expanding
\((1+X)^{p+\ell}\equiv(1+X^p)(1+X)^\ell\).
The residual factorials are \((2\ell-2j)!\) in the low
band and \(2(2\ell-2j)!\) after removing \(p^2\)
in the high band, by Wilson's theorem on the two full
blocks of nonmultiples of \(p\).
Any error in the two bands gains a factor \(p\).
In the middle band the binomial coefficient contributes
one to the valuation. For \(i\le h\) the smallest possible
valuation is \(1+(p+\ell-h)s=2+(2\ell+1)/p\).
For \(i>h\) the factorial contributes one more, and the
smallest valuation is \(2+(\ell+1)s\). This proves (5).

Both \(P_\ell(k)\) and \(P_\ell(k/c)\) have residue one.
The error term in (5) is strictly higher than \(C_\ell(c)\)
by (4). The same holds at \(1\), since the weighted
boundary ratio there is \(g c^\ell\) with residue \(-1\).
Therefore \(\kappa+2\) and \(\kappa+2c^p\) have opposite
nonzero leading residues and valuation \((2\ell+1)/p-r\).
Their difference is \(2(c^p-1)\). If \(r\ge1/(p-1)\) the
binomial expansion gives \(v(c^p-1)\ge 1+r\), impossible
here. Otherwise \(v(c^p-1)=pr\). Thus
\[
           r=\frac{2\ell+1}{p(p+1)},\qquad
           v(\kappa+2)=pr=\frac{2\ell+1}{p+1}.              \tag{6}
\]

## Perturbation formulas at the negative residue

Use the
[prime-terminal perturbation expansion](BINARY_EVEN_QUADRATIC_PRIME_TERMINAL_FOUR_WINDOW.md#the-off-diagonal-case)
for \(f=(x+ty)^2+d y^2+k\). The contribution of total
perturbation count \(L\), with \(i\) occurrences of \(d\),
to \(M_j\) is
\[
 \binom jL\binom L i k^{L-i}d^i R!\,S_{R,i}(t),\qquad
 R=2(j-L),\quad S_{R,i}(t)=\sum_{m=0}^R(m+1)_{2i}t^m .
                                                                  \tag{7}
\]
For \(i\le h\), \(S_{p-1,i}(t)\equiv(2i)!(t-1)^{p-1-2i}\)
modulo \(p\). For \(i>h\) every summand of \(S_{R,i}\)
is divisible by \(p\).
Write \(G_R(t)=S_{R,0}(t)=\sum_{m=0}^R t^m\).

For each \(i>0\), pairing the indices \(m\) and \(m+p\)
gives \(S_{2p-2,i}(-1)=0\bmod p\): the coefficients are
periodic modulo \(p\), the signs are opposite, and the
unpaired index \(m=p-1\) has a factor \(p\).
For \(i=1\) the coefficientwise congruences are
\[
 \begin{split}
 S_{2p-2,1}(t)&\equiv(1+t^p)\,2(t-1)^{p-3},\\
 S_{2p,1}(t)&\equiv(1+t^p)\,2(t-1)^{p-3}+2t^{2p}
 \pmod p .
 \end{split}                                                 \tag{8}
\]
These follow from the same periodicity without substituting
\(t=-1\). Whenever \(v(d)>r\), the relation \(c=t^2+d\)
gives \(v(t+1)=r\). In both cases of (6), \(r<1/(p-1)\)
and \(pr<1\). Then \(v(1+t^p)=v(t^{2p}-1)=pr\).
In particular, \(G_{2p}(t)-1,\ G_{2p+2}(t)-G_2(t)\),
and \(S_{2p,1}(t)-2\) have valuation at least \(pr\),
and \(v(S_{2p-2,1}(t))=pr\). The claims about the \(G\)'s
follow by summing the corresponding geometric series,
and those about \(S\) follow from (8).

## Prime at the second index

Here \(M_{p-1}=0\) first forces \(v(d)=s\).
Otherwise the unique pure term at \(L=h\) with the
smaller of the valuations \(v(d),s\) has valuation
\(h\min(v(d),s)<1\). Earlier factorials are divisible by
\(p\), and later terms strictly gain in powers of \(k,d\).
Thus \(v(d)>r\), and at \(L=h\) the leading group is a
unit times \(k^h g_h(X)\) modulo errors above one, where
\[
       X=d/[k(t-1)^2],\qquad
       g_i(X)=\sum_{m=0}^i\binom im(2m)!X^m .
\]
All other terms have valuation at least one, giving
\(v(g_h(X))\ge1/p\).

At \(M_p\), the terms \(L=0\) and \((L,i)=(p,0)\)
combine to a term of valuation at least \(2+pr\)
by (6) and (8). Here \(pr=3/(p+1)>s\).
For \(1\le L\le h\) the factorial and the order binomial
each contribute a factor \(p\). The lowest possible term
is at \(L=1\), of valuation \(2+s\), and its \(i=1\)
term is strictly higher by (8). For \(L=h+1\) the group
has leading form a unit times \(p k^{h+1}g_{h+1}(X)\),
of initial valuation \(2+1/p\). The \(i=h+1\) term can
be included in this expression modulo \(p\) after
extracting the displayed powers and a unit: its rising
factorial and \((2h+2)!\) are divisible by \(p\).
The error gains a factor \(p\).
For \(h+2\le L<p\), the binomial contributes one
and all terms have valuation at least \(2+3/p\).
At \(L=p,i>0\) the inner binomial or factorial gains at
least one, giving more than \(2+s\).

It follows that \(g_h(\bar X)=g_{h+1}(\bar X)=0\).
The former root is simple and the latter derivative
vanishes there, by
\[
 Xg_i'=i(g_i-g_{i-1}),\qquad
        g_{i+1}=1+2(i+1)X[(2i+1)g_i-2i g_{i-1}] .
\]
Lift the simple root in a finite unramified extension
with the residue of \(X\). The first valuation bound
places the lift at distance of valuation at least \(1/p\)
from \(X\). The next polynomial and its derivative have
valuation at least one at the lift. Taylor expansion
then gives \(v(g_{h+1}(X))\ge2/p\).
Consequently the \(L=h+1\) group of \(M_p\) lies above
\(2+s\). Its \(L=1,i=0\) term is uniquely lowest, contrary
to \(M_p=0\).

## Prime at the first index

Now \(pr=5/(p+1)>2s\). Let \(e=v(d)\).
The equation \(M_p=0\) first implies \(e\ge(1-s)/h\).
Otherwise \(e<s\), and \(L=h+1,i=h\) has valuation
\(1+s+he<2\), uniquely lowest. Earlier positive \(L\)
have an order binomial and a factorial factor \(p\).
Later \(L<p\) with \(i\le h\) have an extra power of \(k\)
relative to that valuation bound, and \(i>h\) adds a
factor \(p\). At \(L=p\) the constant has valuation two.
Its \(i>0\) terms with \(i<p\) either have a power of
\(k\) of degree at least \(h+1\), together with the
binomial factor \(p\), or inherit an additional
rising-factorial factor of \(p\). The pure \(i=p\)
term has factorial valuation two before counting \(d\).

This gives \(e>r\) and \(2e>s\), since
\(e\ge(1-s)/h\), \(p\ge7\), and (6) holds.
The \(L=0\) and \(L=p,i=0\) sum in \(M_p\) is now above
\(2+s\), as in the preceding case. At \(1\le L\le h\)
every \(L\ge2\) contributes above \(2+s\) by
\(2\min(e,s)>s\). At \(L=1,i=1\) equation (8) gains
\(pr>s\), and at \(i=0\) the valuation is \(2+s\).
For \(L=p,i>0\) the mixed terms have valuation at least
\(1+p\min(e,s)>2+s\) and the pure term at least
\(2+pe>2+s\). If \(e<s\), the critical \(L=h+1,i=h\)
has valuation \(1+s+he<2+1/p\), uniquely below both the
earlier groups and the later bounds. If \(e>s\), the
critical \(L=h+1,i=0\) has valuation \(2+1/p\)
and is uniquely lowest. Thus \(e=s\).
With \(X=d/[k(t-1)^2]\), the critical group of \(M_p\)
now gives \(v(g_{h+1}(X))\ge1/p\) by the same estimates.

At \(M_{p+1}\), first remove the following two
cancellations. The \(L=0\) and \(L=p,i=0\) sum has
valuation at least \(2+pr\) by (6) and
\(G_{2p+2}(t)-G_2(t)\) above. At \(L=1\) and \(L=p+1\),
retain the \(i=0,1\) terms. Upon division by \(p^2\)
their sum is
\((2+\kappa)(k+2d)\) modulo terms of valuation at least
\(s+pr\). This follows from (8), \(G_{2p}(t)-1\)
and the coefficient units
\((p+1)(2p)!/p^2\equiv2\) and \(p+1\equiv1\bmod p\).
Other terms at \(L=p+1\) gain an inner binomial or
factorial factor \(p\), and other terms at \(L=p\)
gain one as well. All these terms and the two
cancellations are strictly above \(2+2s\).

For \(2\le L\le h+1\), the factorial and the order
binomial supply a factor \(p^2\). Only \(L=2\) can
contribute at valuation \(2+2s\), and there only
\(i=0\) does: \(S_{2p-2,i}(-1)=0\bmod p\) for \(i>0\),
whereas \(G_{2p-2}(-1)=1\).
The critical \(L=h+2\) contribution has initial
valuation \(2+3/p\) proportional to \(g_{h+2}(X)\),
with a unit coefficient after extracting
\(p k^{h+2}\). All later \(L<p\) have valuation at
least \(2+5/p\). Thus \(g_{h+2}(\bar X)=0\).
The same unramified root lift, now using the simple
root of \(g_{h+1}\), gives \(v(g_{h+2}(X))\ge2/p\).
It suppresses the critical group above \(2+2s\).
The surviving term \(L=2,i=0\) contradicts \(M_{p+1}=0\).
