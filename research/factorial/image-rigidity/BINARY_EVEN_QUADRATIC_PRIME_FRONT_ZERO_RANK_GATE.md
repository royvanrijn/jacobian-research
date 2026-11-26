# Finite-field conditions at the zero cross coefficient in prime-front quadratic windows

Let \(f=ax^2+bxy+cy^2+k\) with \(abck\ne0\), and let
\(M_j=\Lambda_2(f^j)\) with \(\Lambda_2(x^i y^j)=i!j!\).
Fix a prime \(p\ge7\). Suppose \(M_j\) vanishes on one of
the two four-windows in the
[prime-front reduction](BINARY_EVEN_QUADRATIC_PRIME_FRONT_REDUCTION.md).
Take the window \(p-1,p,p+1,p+2\) when \(\ell=1\), and
\(p,p+1,p+2,p+3\) when \(\ell=2\).
Normalize \(a=1\) as there. If \(t=b/2\) reduces to zero, put
\(d=c-t^2,\ h=(p-1)/2\). Then
\[
                v(t)=1/p,\qquad v(d)=v(k)=2/p.             \tag{1}
\]
In particular \(X=\overline{d/k}\ne0\).
For the factorial polynomials
\[
             g_j(X)=\sum_{i=0}^j\binom ji(2i)!X^i,\qquad
             G=g_{h+1}(X),\quad H=g_h(X),
\]
the following equations are necessary in the residue field:
\[
\begin{array}{c|l}
 \ell=1 & H=0,\quad 3G(1+2X-4G^2)=2\\[2pt]
 \ell=2 & 15G(1+4W+24W^2)=8,\quad
          -12G^3+3G+6XG-2+6XH=0,\quad W=X-2G^2.
\end{array}                                                  \tag{2}
\]
The equations are over \(\mathbb F_p\). A gcd equal to one
excludes this residue for that prime. Nontrivial gcds are
possible, and (2) asserts only a necessary condition.
As in the prime-front reduction, the statements apply to
characteristic-zero solutions via algebraic specialization.

## Comparing the two boundary orders

Write \(s=v(k)=2/p\), \(e=v(d)>0\), \(r=v(t)>0\),
\(\rho=v(c)>0\), and \(\kappa=k^p/p^2\).
The initial reduction gives \(\bar\kappa=-2\).
Use the notation of the
[negative-residue boundary calculation](BINARY_EVEN_QUADRATIC_PRIME_FRONT_NEGATIVE_RANK_LIFT.md#the-strength-forced-by-the-two-weighted-boundary-equations):
\[
 B_j(Z)=\sum_{i=0}^j\binom ji(2i)! Z^i k^{j-i},\quad
 a_j=(2p+2j)!/(2p)!,\quad C_j(Z)=B_{p+j}(Z)/(a_j Z^j),
 \quad g=(c-t)/(t-1).
\]
Here \(\sigma=v(g)=v(c-t)>0\). The last two weighted
identities read \(C_j(1)=g B_{p+j}(c)/a_j\)
for \(j=\ell,\ell+1\). At unit arguments the
[Laurent coefficient expansion](BINARY_EVEN_QUADRATIC_PRIME_FRONT_NEGATIVE_RANK_LIFT.md#the-strength-forced-by-the-two-weighted-boundary-equations)
\[
 C_j(Z)=p^2P_j(k/Z)(\kappa+2Z^p)+E_j(Z),\qquad
 P_j(Y)=\sum_{b=0}^j\binom jb
                     \frac{(2j-2b)!}{(2j)!}Y^b
\]
applies for \(j\le3\) with every error coefficient at
valuation at least \(2+(2j+1)/p\). The same two-band
proof applies at \(j=3\), since \(2j<p\) and
\(2+(2j+1)/p\le3\).
For \(j=1,2\), the error evaluated at \(1\) has exactly
that valuation, from the unique middle coefficient
at index \(h\) in \(B_{p+j}\). Both exterior bands
gain a full \(p\) after their matched leading terms,
so their errors are strictly higher for \(j=1,2\).

It follows that
\[
 v\left(C_{\ell+1}(1)
       -\frac{P_{\ell+1}(k)}{P_\ell(k)}C_\ell(1)\right)
                                  =2+(2\ell+1)/p .           \tag{3}
\]
The factors \(P_j(k)\) are units. For nonunit \(c\),
if \(\rho<s\), \(B_{p+j}(c)\) has the unique lowest
index \(j\) with valuation \(2+j\rho\).
Indeed, the low exterior band has values
\(2+(j-i)s+i\rho\), uniquely smallest at \(j\).
The high exterior band has valuation at least
\(2+(p+j-i)s+i\rho>2+j\rho\).
For \(j<i\le h\) the difference from \(2+j\rho\) is
\(1-(i-j)s+(i-j)\rho>0\). For \(h<i<p\) the
factorial adds one to that lower bound.
If \(\rho\ge s\), all terms instead have valuation at
least \(2+js\) by the same bands, or by their value
at \(\rho=s\) and monotonicity in \(\rho\).
Consequently (3) gives
\[
 \rho<s\ \Longrightarrow\ \sigma+\ell\rho=(2\ell+1)/p,
 \qquad
 \rho\ge s\ \Longrightarrow\ \sigma\le1/p.                 \tag{4}
\]
In the first case \(C_{\ell+1}(1)\) is strictly above
\(2+(2\ell+1)/p\), since its valuation is that of
\(C_\ell(1)\) plus \(\rho>0\). In the second case it
is also strictly above that value, by the bound
\(2+\sigma+(\ell+1)s\). In both cases the expansion
at \(Z=1\) shows \(v(\kappa+2)>(2\ell+1)/p\).

## The valuations of the cross coefficient and rank correction

We use the
[perturbation formula](BINARY_EVEN_QUADRATIC_PRIME_FRONT_NEGATIVE_RANK_LIFT.md#perturbation-formulas-at-the-negative-residue)
for \(f=(x+ty)^2+d y^2+k\).
For an order \(j\), the term with \(L\) perturbations,
\(i\) of them from \(d\), has factor
\(\binom jL\binom L i k^{L-i}d^i R!S_{R,i}(t)\)
with \(R=2(j-L)\). At \(R=p-1\), the factor
\(S_{p-1,i}(t)\) for \(i\le h\) is a unit because
\(t-1\) is a unit. For \(i>h\) it is divisible by
\(p\) for every \(R\).

If \(\ell=1\), the order \(p-1\) forces \(e=s\):
any strict inequality makes the appropriate pure term
at \(L=h\) uniquely lowest at \(h\min(e,s)<1\).
Earlier terms have a factorial factor \(p\), and
later terms add a positive perturbation valuation.
If \(\rho<s\), the condition \(e=s\) implies
\(\rho=2r<s\) and \(\sigma=r\), so (4) would give
\(r=1/p\), a contradiction. Thus \(\rho\ge s\);
then \(r\ge1/p\), and (4) forces \(r=1/p\).
This proves (1) when \(\ell=1\).

If \(\ell=2\), first \(M_p=0\) implies
\[
                     e\ge e_0:= (1-s)/h .
\]
If \(e<e_0<s\), the term \(L=h+1,i=h\) has valuation
\(1+s+he<2\), uniquely lowest. Earlier positive \(L\)
have factorial and binomial factors contributing \(p^2\).
Later \(L<p\) with \(i\le h\) raise the bound by
an extra power of \(k\), while those with \(i>h\)
gain another factor \(p\). For \(L=p\), the constant
has valuation two. A mixed term with \(i\le h\)
has an inner binomial factor \(p\) and a power of
\(k\) of degree at least \(h+1\). One with \(i>h\)
has a rising-factorial factor \(p\) in addition to
the binomial, and the pure \(i=p\) term has factorial
valuation two.

Suppose \(\rho<s\). If \(2r<e\), then \(\rho=2r\),
\(\sigma=r\), and (4) forces \(r=1/p\), impossible.
Otherwise \(e\le2r\), so \(\rho\ge e\) and \(e<s\).
If \(r\ge e\), (4) gives \(3e\le5/p\).
For \(p>7\), \(3e_0>5/p\). For \(p=7\) equality
forces \(e=e_0\), and the unique critical term at
\(L=h+1,i=h\) has valuation exactly two, while the
\(L=0\) and \(L=p,i=0\) pair cancels at two, since
\(t\) reduces to zero and \(\bar\kappa=-2\).
All other terms are higher, by the same bounds.
Thus \(r<e\) is necessary, and \(\sigma=r\).
Now the \(L=0\) and \(L=p,i=0\) pair has valuation
exactly \(2+r\): its first remaining contribution
is the linear \(t\) from \(G_{2p}(t)\), since
\(v(\kappa+2)>5/p>r\).
The critical \(L=h+1,i=h\) term has valuation
\(1+s+he\). All other preceding positive \(L\)
are at least \(2+e>2+r\); subsequent \(L<p\)
gain over the critical bound as above. In the
last group, mixed terms have valuation at least
\(1+pe_0>2+s\) and the pure \(d\) power is higher.
Thus \(M_p=0\) forces \(r=s+he-1\).
If \(\rho=e\), substitute \(r+2e=5/p\) to obtain
\(e=s\). If \(\rho>e\), cancellation in \(c=t^2+d\)
requires \(2r=e\); substituting also gives \(e=s\).
Both contradict \(\rho<s\).

Therefore \(\rho\ge s\). If \(e<s\), this requires
\(2r=e,\ \sigma=r<1/p\). The same \(M_p\) comparison
applies, since \(r<e\) and \(v(\kappa+2)>5/p\);
it again gives \(e=s\). We have proved \(e\ge s\)
and \(\rho\ge s\), which force \(r=1/p\) by (4).
For both prime positions we now have
\(\sigma=r=1/p\), and the upper boundary gives the
stronger estimate
\[
                         v(\kappa+2)\ge (2\ell+3)/p .       \tag{5}
\]
This follows from its right side \(g B_{p+\ell+1}(c)/a_{\ell+1}\)
of valuation at least \(2+(2\ell+3)/p\).

## The two residue conditions

At this stage \(X=\overline{d/k}\) is defined, and may
provisionally be zero for \(\ell=2\).
Set \(C=(-1)^h,\ \alpha=p/k^h\), and
\(\zeta=\overline{t/\alpha}\). We have \(v(\alpha)=1/p\),
\(k/\alpha^2=\kappa\), and
\(\overline{c/k}=W=X-\zeta^2/2\).
At \(M_p/p^2\) use valuation \(1/p\).
The \(L=0\) and \(L=p,i=0\) pair contributes
\(2t\) at that valuation. At \(L=h+1\) the
critical group contributes
\(-2C(k^{h+1}/p)g_{h+1}(d/[k(t-1)^2])\).
The unit \(-2C\) is the product of
\(\binom p{h+1}/p\equiv2C\) and \((p-1)!\equiv-1\).
All other terms are higher: earlier positive \(L\)
give at least \(2+s\) before normalization, later
\(L<p\) at least \(2+3/p\), and \(L=p,i>0\)
gain an inner binomial or factorial factor \(p\).
Dividing the retained terms by \(\alpha\) and reducing
gives \(\zeta=-2C G\). In particular \(G\ne0\), since
\(t\) and \(\alpha\) have the same valuation.
Hence \(W=X-2G^2\).

For the lower boundary order \(p+\ell\), use valuation
\((2\ell+1)/p\) after dividing by \(p^2\).
At \(C_\ell(1)/p^2\), only the middle index \(h\)
contributes, by (5). Its coefficient relative to
\(k^{h+\ell+1}/p\) reduces to \(-2C/3\) for \(\ell=1\)
and \(-2C/45\) for \(\ell=2\).
Indeed,
\[
 \binom{p+\ell}{h}/p=
   \frac{(p+1)\cdots(p+\ell)}{h(h-1)\cdots(h-\ell)}
                  \binom{p-1}{h-\ell-1},
\]
and the remaining normalized factorial is
\((p-1)!/a_\ell\equiv-1/(2\ell)!\).
At the nonunit \(c\), \(B_{p+\ell}(c)/k^{p+\ell}\)
reduces to \(1+2W\) for \(\ell=1\) and
\(1+4W+24W^2\) for \(\ell=2\): exactly the low
exterior band can contribute. Moreover \(g/t\) reduces
to one. Converting the weighted equality
\(C_\ell(1)=g B_{p+\ell}(c)/a_\ell\) to the scale
\(p^2\alpha^{2\ell+1}\) now gives
\[
 \begin{array}{ll}
  \ell=1:& 3\zeta(1+2W)=-4C,\\
  \ell=2:& 15\zeta(1+4W+24W^2)=-16C.
 \end{array}                                                \tag{6}
\]
With \(\zeta=-2CG\), these are the respective boundary
equations in (2). For \(\ell=1\), \(M_{p-1}\) gives
\(H=0\) from its first critical group, so (2) follows.

For \(\ell=2\) use \(M_{p+1}\) at valuation \(3/p\)
after dividing by \(p^2\). Its \(L=0\) and \(L=p,i=0\)
sum contributes \(4t^3\) at this valuation. Indeed,
\(G_{2p+2}(t)-G_2(t)=t^3+\cdots\), and (5) makes
the remaining constant scale difference higher.
The \(L=1\) and \(L=p+1\) contributions with \(i=0,1\)
cancel at valuation \(s\). At valuation \(3/p\) their
sum is \(2tk+12td\), using
\(G_{2p}(t)=1+t+\cdots\) and
\(S_{2p,1}(t)=2+6t+\cdots\).
The omitted powers of \(t\) add at least \(2/p\),
the omitted binomials and factorials in the extreme
groups add at least one, and \(2\le L\le h+1\)
is at least \(2+2s\) before normalization.
The critical \(L=h+2\) group contributes
\(-4 C(k^{h+2}/p)g_{h+2}(d/[k(t-1)^2])/3\).
This uses \(\binom{p+1}{h+2}/p\equiv4C/3\).
Later \(L<p\) give at least \(2+5/p\) before
normalization. Dividing by \(\alpha^3\) and reducing
gives
\[
                  \zeta^3-\zeta-6\zeta X
                          -(4/3)C g_{h+2}(X)=0.
\]
The binomial recurrence for the \(g_j\) in the
[negative-residue lift](BINARY_EVEN_QUADRATIC_PRIME_FRONT_NEGATIVE_RANK_LIFT.md#prime-at-the-second-index)
gives \(g_{h+2}=1+6XG-3XH\) in characteristic \(p\).
Substitution proves the last equation of (2).
If \(e>s\), then \(X=0\) and \(G=H=1\).
The two equations at \(\ell=2\) would demand
\(p\mid1327\) and \(p\mid11\), respectively.
The inconsistency proves \(e=s\), completing (1).

## Exact finite-field checks and scope

The [exact replay](verify_binary_even_quadratic_prime_front_zero_rank_gate.py)
uses Euclidean gcd over \(\mathbb F_p[X]\) and gives the
following surviving monic factors in the stated finite
test ranges. Compute the \(g_j\) by their displayed
coefficients, substitute \(W=X-2G^2\), and take the
gcd of the two polynomials obtained by moving each
right side of (2) to the left. The two ranges contain
98 and 43 primes respectively.

| Prime position | Primes tested | Nonconstant gcd |
|---|---|---|
| Second | \(7\le p<550\) | \(p=19:\ X+2\) |
| First | \(7\le p<200\) | \(p=11:\ X+6;\quad p=37:\ X+9\) |

The listed nonzero roots satisfy the necessary equations.
Their lifting to characteristic-zero solutions is unresolved.
The derivation of (1)–(2) holds for every prime \(p\ge7\).
The residue \(1\) remains outside this argument.
