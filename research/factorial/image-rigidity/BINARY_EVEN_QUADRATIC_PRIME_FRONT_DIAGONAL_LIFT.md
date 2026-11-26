# Excluding the diagonal residue in prime-front quadratic windows

Let \(f=ax^2+bxy+cy^2+k\) with \(abck\ne0\), and let
\(M_j=\Lambda_2(f^j)\), where \(\Lambda_2(x^i y^j)=i!j!\).
Fix a prime \(p\ge7\). If the four moments on \(p-1,p,p+1,p+2\)
or \(p,p+1,p+2,p+3\) vanish, then in the normalization \(a=1\)
of the [prime-front reduction](BINARY_EVEN_QUADRATIC_PRIME_FRONT_REDUCTION.md),
\(t=b/2\) cannot reduce to \(1\).
We work over the algebraic \(p\)-adic closure, which suffices
by the algebraic specialization in the cited reduction. Use
\(\ell=1,2\), respectively, for the two windows.

Suppose otherwise. Put
\[
 h=(p-1)/2,\quad C=(-1)^h,\quad u=t-1,\quad d=c-t^2,\quad
 D=u^2+d=c+1-2t,\quad s=v(k)=2/p,\quad \kappa=k^p/p^2.
\]
Then \(\bar\kappa=-2,\ r=v(u)>0,\ e=v(d)>0\).
Bars denote residue reduction, and \(v(p)=1\).
The [single-face](BINARY_EVEN_QUADRATIC_SINGLE_FACE_STRONG_FACTORIAL.md)
and [critical boundary](BINARY_EVEN_QUADRATIC_CRITICAL_BOUNDARY_WINDOWS.md)
results allow us to assume \(u,d,D,c-1,c-t\ne0\).
Let \(G_m(t)=\sum_{i=0}^m t^i\), and let \(o_v(z)\)
denote a term of valuation strictly above \(z\).

## Boundary and perturbation formulas

Write \(a_j=(2p+2j)!/(2p)!\),
\(B_n(Z)=\sum_i\binom ni(2i)! k^{n-i}Z^i\), and
\(b_j=(2j+1)/p\). For \(j=\ell,\ell+1\), the
[weighted boundary equations](BINARY_EVEN_QUADRATIC_SINGLE_FACE_STRONG_FACTORIAL.md#4-the-exact-reduction-for-general-b-and-its-limitation)
are \((c-t)B_{p+j}(c)-uB_{p+j}(1)=0\).
The [two exterior bands](BINARY_EVEN_QUADRATIC_PRIME_FRONT_NEGATIVE_RANK_LIFT.md#the-strength-forced-by-the-two-weighted-boundary-equations)
give
\[
 \frac{B_{p+j}(Z)}{a_j Z^j}
     =p^2P_j(k/Z)(\kappa+2Z^p)+E_j(Z),\qquad
 P_j(Y)=\sum_{i=0}^j\binom ji\frac{(2j-2i)!}{(2j)!}Y^i .
                                                               \tag{1}
\]
All coefficients of \(E_j\) have valuation at least \(2+b_j\)
for \(j\le3\); the band errors have valuation at least three.
This uses \(2j<p\), also at \(j=3\). For \(j=1,2\) the
first error is uniquely the middle index \(h\), with Laurent
exponent \(h-j\) and value at \(1\) having leading coefficient
\(\mu_j\) relative to \(p k^{h+j+1}\), where
\[
                        \mu_1=-2C/3,\qquad \mu_2=-2C/45.     \tag{2}
\]
These are the [middle coefficients](BINARY_EVEN_QUADRATIC_PRIME_FRONT_ZERO_RANK_GATE.md#the-two-residue-conditions).

In \(M_n\) the term of \(L\) perturbations (with \(i\) occurrences
of \(d\)) in \(f=(x+ty)^2+d y^2+k\) is
\[
 \binom nL\binom L i k^{L-i}d^i R!S_{R,i}(t),\qquad
 R=2(n-L),\qquad S_{R,i}(t)=\sum_{m=0}^R(m+1)_{2i}t^m .
                                                               \tag{3}
\]
At \(R=p-1\) and \(i\le h\),
\(S_{p-1,i}(t)\equiv(2i)!u^{2(h-i)}\pmod p\);
for \(i>h\) every summand has a factor \(p\).
For \(n=h+L_0\) and \(L>L_0\), the same congruence
holds if \(i\ge L-L_0\), because the deleted terms contain
a factor \(p\). If \(i<L-L_0\) there are at least \(L_0+1\)
powers of \(k\). These are the
[prime-terminal group comparisons](BINARY_EVEN_QUADRATIC_PRIME_TERMINAL_FOUR_WINDOW.md#the-diagonal-transition-s2p).
For \(1\le L\le h\),
\[
 S_{2p-2L,L}(t)\equiv(1+t^p)S_{p-1,L}(t)\pmod p,              \tag{4}
\]
by pairing the two blocks and discarding the missing suffix of
terms with a factor \(p\). The same pairing at \(R=2p-2\) works
for every \(i>0\). Also
\(S_{2p,1}(t)\equiv(1+t^p)2u^{p-3}+2t^{2p}\pmod p\).

## The rank correction has higher valuation than \(u^2\)

If \(e\le2r\), the equation \(M_p=0\) first gives
\(e>(1-s)/h\). Otherwise \(e<s\), the term \(L=h+1,i=h\)
has valuation \(1+s+he\le2\), uniquely lowest in its group.
Every preceding positive \(L\) has an order binomial and a
factorial contributing \(p^2\). Subsequent \(L<p\) with
\(i\ge L-h-1\) and \(i\le h\) have lower bound
\(1+2hr+(L-i)s+i(e-2r)\), which gains at least \(s\)
over the critical minimum. The congruence errors and \(i>h\)
gain a further factor \(p\); \(i<L-h-1\) requires at least
\(h+2\) powers of \(k\) besides the order binomial.
In \(L=p\), the mixed terms with \(i\le h\) have
\(h+1\) or more powers of \(k\) and an inner binomial
factor \(p\); those with \(h<i<p\) have both inner
binomial and rising-factorial factors \(p\), and the pure
\(i=p\) term has factorial valuation two. The \(L=0\)
and \(L=p,i=0\) pair is above two because
\(G_{2p}(1)\equiv1\) and \(\bar\kappa=-2\). Thus the
critical term would survive.

Let \(b=b_\ell\), \(g=(c-t)/u\), and put
\[
 U_j(Z)=B_{p+j}(Z)/(a_jZ^j P_j(k/Z)),\qquad
 g'_j=g c^j P_j(k/c)/P_j(k).
\]
Then \(U_j(1)=g'_j U_j(c)\). In (1), division of \(E_j\)
by \(P_j(k/Z)\) preserves its coefficient bound in a
convergent Laurent series on units. Its difference between
\(1,c\) gains at least \(v(c-1)\), by factoring each
Laurent monomial difference.

If \(e<2r\), put \(q=e-r=v(g-1)\) and \(\rho=v(c-1)\).
Then \(v(g'_j-1)=q<\rho\). Indeed the extra unit factor
differs from one by valuation at least \(\rho\);
for \(q\le0\) this remains a strict gain after
multiplication by \(g\).
We also have \(v(c^p-1)>b+q\). For \(e\le r\)
use \(\rho\ge e>(1-s)/h>1/(p-1)\) and \(q\le0\).
For \(r<e<2r\), \(\rho=r\) and the binomial bound is
\(\min(pr,1+r)>b+r\), since
\((p-1)r>(p-2)/p\ge b\).
At least one \(U_j(c)\) has valuation at most \(2+b\),
because the lower error has exactly that value and the
upper error is strictly higher. But its difference from
\(U_j(1)\) is above \(2+b+q\), contrary to the weighted ratio.

For \(e=2r\), \(\rho=r\) and \(v(c^p-1)>b+r\).
Let \(\lambda=\overline{D/u^2}\); the coefficient of \(u\)
in \(g'_j-1\) reduces to \(\lambda+2j\).
If \(v(\kappa+2)<b\), the leading difference in each row is
strictly higher than \(u U_j(c)\), forcing both
\(\lambda=-2\ell\) and \(\lambda=-2(\ell+1)\).
If \(v(\kappa+2)>b\), the first lower error has Laurent
exponent \(h-\ell\); its leading difference gives
\(\lambda+2\ell=-2(h-\ell)\), so \(\lambda=1\),
contrary to \(v(d/u^2)=0\).
The remaining \(v(\kappa+2)=b\) forces
\(\lambda=-2(\ell+1)\) in the upper row.
The lower row then gives, at precision \(b\),
\[
         \kappa+2=(h-\ell-1)\mu_\ell k^{h+\ell+1}/p+o_v(b).
                                                               \tag{5}
\]
In fact its leading difference is
\(-2(h-\ell)u E_\ell(1)\), and the right side of the
difference equation is \(-2u[p^2(\kappa+2)+E_\ell(1)]\).

At \(M_p\), the critical \(L=h+1,i=h\) term has
valuation \(1+s+he\); the \(L=0\) and \(L=p,i=0\)
pair is above \(2+s\), by (5) and the geometric identity
\(G_{2p}(t)-1=t(t^{2p}-1)/u\). The only other term
at or below \(2+s\) is \(L=1,i=0\): before the critical
group this follows from (4) and \(he>(p-2)/p>s\);
later groups use the previous bounds, now up to \(2+s\).
Thus \(he=1\), \(r=1/(p-1)\), and \(e>s\).
Write \(A=d^h/p\). The two groups at \(s\) in \(M_p/p^2\)
give \(k(1+CA)\), so \(\bar A=-C\).

In \(M_{p+1}/p^2\) the top pair contributes \(6(\kappa+2)\)
through valuation \(5/p\): the other term uses
\(G_{2p+2}(t)-G_2(t)=t^3(t^{2p}-1)/u\), of valuation at
least one. The \(L=1,p+1\) constant and \(d\)-linear terms
combine to \((2+\kappa)(k+2d)\) plus errors above \(5/p\),
using \(G_{2p}-1\) and \(S_{2p,1}-2\), the latter of
valuation at least \((p-3)r>3/p\).
Other extreme terms gain an inner binomial or factorial
factor \(p\). Of the middle groups, the \(L=2,i=0\)
and \(L=h+2,i=h\) terms combine to \(k^2(1+CA)/2\)
up to errors above \(5/p\). The \(i=1,2\) terms at
\(L=2\) use the periodic identity after (4); other
early \(L\) start at \(6/p\), and decreasing \(i\)
at the critical group costs \(s\). The only subsequent
term at \(5/p\) is \(L=h+3,i=0\) from the exception
\(i<L-h-2\).

For \(\ell=1\) only \(6(\kappa+2)\) contributes at
\(3/p\), which is impossible. For \(\ell=2\), the
extra precision \(3/p\) in \(M_p/p^2\) has only the
exception \(L=h+2,i=0\), with coefficient
\(-2C/3\) times \(k^{h+2}/p\). Thus
\(1+CA=(2C/3)k^{h+1}/p+o_v(1/p)\).
The propagated pair in \(M_{p+1}/p^2\) contributes
\((C/3)k^{h+3}/p\), and its own exception contributes
\(-4C k^{h+3}/(15p)\), since
\(\binom{p+1}{h+3}/p\equiv-4C/15\) and
\((p-3)!G_{p-3}(1)\equiv1\).
Equation (5) contributes \(6(\kappa+2)=
(14C/15)k^{h+3}/p+o_v(5/p)\).
The sum \(C k^{h+3}/p\) has valuation \(5/p\), another
contradiction. We have proved \(e>2r\).

## The new strength constraints

Put
\[
 \Delta=e-2r,\quad \delta=d/u^2,\quad U=(p-1)r,\quad
 \psi=(c^p-1)(c-t)/D,\qquad
 R_i=\frac{(c-t)c^i-u}{D}.
\]
Then \(v(\delta)=\Delta>0\), and
\[
 R_i=1+\frac{(t+u\delta)(t+1+u\delta)}{1+\delta}
                               \frac{c^i-1}{c-1}.           \tag{6}
\]
The \(R_i\) are integral with residue \(2i+1\).
At \(\delta=0\), \(R_i=G_{2i}(t)\), and
\(R_{p+i}-R_i=c^i\psi\) is exact.
If \(U\ne1\), \(v(\psi)=\min(U,1)\), and if \(U=1\)
it is at least one. Also \(v(R_h)\ge\min(U,1,\Delta)\),
with first term \(\delta\) at precision \(\Delta\)
when \(\Delta<\min(U,1)\).

Apply the weighted equations to the paired exterior
bands after division by \(D a_j p^2\). They give
\[
         L_j(\kappa+2)+2c^j P_j(k/c)\psi,\qquad
         L_j=\sum_{i=0}^j [Y^{j-i}]P_j(Y)\,k^{j-i}R_i,       \tag{7}
\]
up to coefficient errors of valuation at least one. The
matrix of the two displayed linear forms reduces to rows
\((2j+1,2)\), with determinant \(-4\).
The middle terms have bounds
\[
    b_j+v(R_h)\ (i=h),\qquad b_j+s\ (i<h),\qquad
    b_j+1/p\ (i>h).                                         \tag{8}
\]
Index \(p-1\) alone attains the last bound.
Since \(v(R_h)>0\), inversion first gives
\(v(\kappa+2),v(\psi)>b\), hence \(U>b\).
Using (8) again gives
\[
                    v(\kappa+2),v(\psi)\ge b+\min(\Delta,1/p).
                                                               \tag{9}
\]
If \(\Delta\ne1/p\), the lower residual at this valuation
is unique at \(h\) for \(\Delta<1/p\) and at \(p-1\) otherwise.
The upper residual is strictly higher. If \(2\ell+3\)
is a unit, both values in (9) are then equalities.

At \(M_p\) the top pair \(L=0\) and \(L=p,i=0\)
is above \(2+s\) by \(U>b\) and (9).
At positive \(L\le h\) the only term at or below \(2+s\)
is \(L=1,i=0\). Mixed terms with \(i<L\) have a power
of \(k\) and a positive additional exponent unless
\(L=1,i=0\). For pure \(d\) terms (4) gives the bound
\(2+\min(U+L\Delta,1+Le)>2+s\).
The critical \(L=h+1\) group has initial bound
\(1+U+s+h\min(\Delta,s)\). When it is at most \(2+s\)
it is realized uniquely at \(i=h\), since \(\Delta<s\).
Later \(L<p,\ i\ge L-h-1,\ i\le h\) raise the
critical lower bound by at least \((L-h-1)s\).
If \(i<L-h-1\), the \(k\)-power gives at least
\(1+(h+2)s>2+s\); congruence errors and \(i>h\)
gain an additional \(p\) over a perturbation group.
For \(L=p\) the mixed terms with \(i\le h\) have bound
\(1+(h+1)s+h e\) if \(e<s\), or \(1+ps\) otherwise.
The first is above \(2+s\) by \(he\ge U>b\).
Those with \(i>h\) are at least \(2+p\min(e,s)>2+s\),
including the pure \(i=p\) term.
Thus \(M_p=0\) gives
\[
                 U+h\Delta=1,\qquad e=1/h>s,\qquad U<1.     \tag{10}
\]
In particular \(v(\psi)=U\).
At valuation \(s\) in \(M_p/p^2\) the surviving pair
again gives \(\bar A=-C\) for \(A=d^h/p\).

## Excluding the two strict branches

For \(J_1=M_p-pk M_{p-1}\) and
\(J_2=M_{p+1}-(p+1)k M_p/2\), compare terms with the
same \(R,i\) in (3). If \(L\) is the higher perturbation
count, the binomial coefficient ratios, including the
extra \(k\), are \(pk/(L-i)\) and \((p+1)k/(L-i)\)
respectively. Thus the earlier pure term at \(L=\ell,i=0\)
and the critical pure term at \(L=h+\ell,i=h\) cancel
exactly in \(J_\ell\).

For \(\Delta<1/p\) in \(J_\ell/p^2\), the next critical
index \(h-1\) is above \(b\), since its gain over the
canceled group is \(s-\Delta>1/p\).
The first pure exception from \(L=h+\ell+1,i=0\),
combined with its shifted counterpart in the lower
order, gives the nonzero term
\[
      \frac{C}{m_\ell}\frac{k^{\ell+1}}{\alpha},
       \qquad \alpha=p/k^h,\quad m_1=3,\quad m_2=15
                                                               \tag{11}
\]
at valuation \(b\). Its coefficient is found using
\((p-3)!G_{p-3}(1)\equiv1\) and
\(\binom{p-1}{h+1}\equiv-C,\ \binom p{h+2}/p\equiv-2C/3,\
\binom{p+1}{h+3}/p\equiv-4C/15\).
The top pairs are above this precision by (9)–(10).
In the earlier groups, \(e>s\) puts additional \(L\) at
valuation at least \((\ell+1)s>b\) after division by
\(p^2\), except for the mixed terms at the same \(L=\ell\)
and the groups considered next. The mixed terms at
\(L=\ell,\ i>0\) have bound \(U+(\ell-i)s+i\Delta>b\)
using the periodic identity at \(R=2p-2\), or the bound
from the factor \(p\) in its congruence error.
For \(J_2\), the \(M_{p+1}\) groups \(L=1,p+1\) combine
their constant and linear \(d\) terms to
\((2+\kappa)(k+2d)\) plus errors of valuation at least
\(\min(s+U,U+\Delta,1+s)\), by the two identities after (4) and
\(G_{2p}-1=t(t^{2p}-1)/u\).
The omitted extreme mixed terms gain a binomial or factorial
factor \(p\). After the critical group the bounds from (3)
raise the minimum by \(s\) at each step for
\(L-L_0\le i\le h\), where \(L_0\) is the critical count;
the cutoff \(i<L-L_0\) leaves only the pure exception
at this precision. Thus (11) survives, a contradiction.

If \(\Delta>1/p\), (9)–(10) imply \(p>4\ell+3\),
so \(v(\kappa+2)=U=(2\ell+2)/p\) by (8).
The upper row gives
\((2\ell+3)(\kappa+2)+2\psi=0\) at precision \(U\).
The top pair of \(M_{p+1}/p^2\) at this precision is
\(6(\kappa+2)+4\psi=-4\ell(\kappa+2)\).
Indeed \(G_{2p+2}-G_2=t^3(t^{2p}-1)/u\) and \(\psi\)
have the same leading value because \(U<1,\ e>2r\).
For \(\ell=1\), the other two terms at precision \(4/p\)
combine to \(k^2(1+CA)/2\), which cancels by \(M_p\);
the group comparisons above apply since the mixed
critical index \(h-1\) gains \(s-\Delta>0\).
Thus the top pair survives.
For \(\ell=2\), use \(J_2/p^2\): the critical index
\(h-1\) is uniquely lowest at \(4/p+s-\Delta<5/p\).
Its coefficient ratio \((p+1)k/3\) differs from the
subtraction multiplier \((p+1)k/2\). The top pair
has valuation \(6/p\), and the pure exception \(5/p\).
Earlier mixed and later terms are higher by the previous
estimates, using \(U=6/p\) and \(\Delta>1/p\).
This excludes both strict branches.

## The second correction at the remaining strength

We have
\[
 \Delta=1/p,\qquad U=(p+1)/(2p),\qquad e=1/h,\qquad
             v(\alpha)=v(\delta)=1/p .                      \tag{12}
\]
By (9), \(p\ge4\ell+3\).
At precision \(b\) in \(J_\ell/p^2\), the pure exceptional
contribution (11) and the critical index \(h-1\) give
\[
                    \overline{\delta/\alpha}=S,\qquad
             S=3C/8\ (\ell=1),\quad S=5C/16\ (\ell=2).      \tag{13}
\]
The two leading terms and their stricter expansions are
listed below. At \(p=4\ell+3\), the two weighted middle
terms in (8) at precision \(U=b+1/p\) cancel by (13).
Their coefficients relative to \(k^{\ell+1}\) are
\(\mu_\ell\delta/\alpha\) and \(1/[(\ell+1)(2\ell)!]\):
\(R_h\) has leading value \(\delta\) at precision \(1/p\),
and \(R_{p-1}\equiv-1\). The upper residual is higher.
Both rows (7) would then force \(v(\psi)>U\), impossible.
We may take \(p>4\ell+3\); the upper boundary gives
\[
                  v(\kappa+2)\ge\min(U,b+3/p,1)>b+1/p .
                                                               \tag{14}
\]

The normalized equation \(M_p/p^2=0\) through valuation
\(3/p\) is
\[
 -kG_{2p-2}(t)+CkA
   -\frac C8 k^2 A/\delta-\frac{2C}{3}k^2/\alpha
                      +o_v(3/p)=0 .                         \tag{15}
\]
The first two terms are \(L=1,i=0\) and \(L=h+1,i=h\).
The next two are \(L=h+1,i=h-1\) and \(L=h+2,i=0\).
The mixed coefficient uses
\(\binom p{h+1}/p\equiv2C,\ \binom{h+1}{h-1}\equiv-1/8\)
and \((p-1)!(p-3)!\equiv1/2\);
the last is the pure exception from (11).
The top pair is higher by (14).
The earlier mixed \(d\) term at \(L=1\) gains \(U+\Delta\)
by (4); other earlier groups start at \(2s\)
after normalization. Later groups with \(i\ge L-L_0\)
gain another \(s\), and the only pure exception at
\(3/p\) is the displayed \(L=h+2,i=0\).
Thus, with \(A=d^h/p\),
\[
        A=C G_{2p-2}(t)+\bigl(C/(4S)-4/3\bigr)\alpha
                                         +o_v(1/p).         \tag{16}
\]

For \(J_\ell/(p^2 k^{\ell+1})\), the terms at valuations
at most zero are as follows. The index refers to the
group in the higher moment. The Laurent factor in \(t\)
is retained where its linear term can contribute.
\[
\begin{array}{c|rr}
 & \ell=1 & \ell=2\\ \hline
 L=h+\ell,\ i=h-1 & CA/(8\delta)& CA/(48\delta)\\
 L=h+\ell+1,\ i=0 & -C G_{p-3}(t)/(6\alpha)&
                         -C G_{p-3}(t)/(30\alpha)\\
 L=h+\ell,\ i=h-2 &-C kA/(96\delta^2)&-C kA/(768\delta^2)\\
 L=h+\ell+1,\ i=h &CA/8&CA/48\\
 L=\ell+1,\ i=0 &1/4&1/24
\end{array}                                                  \tag{17}
\]
Use the binomial ratios in the preceding section to combine
the shifted lower group with its higher counterpart.
At \(h-1\) the lower coefficients, apart from \(u^2\)
and the \(k,d\) powers, are \(-C/4,-pC/8\);
the ratio differences are \(-pk/2,-(p+1)k/6\).
At \(h-2\) the lower coefficients apart from \(u^4\)
and the powers are \(C/64,pC/192\), and the ratio
differences are \(-2pk/3,-(p+1)k/4\).
At the later \(i=h\) group the lower coefficients
are \(-C/4,-pC/8\) and the differences are
\(-pk/2,-(p+1)k/6\). At the earlier pure group
the lower coefficients are \(-p/2,-p^2/4\) after
\(G_{2p-4}(1)\equiv-3\), and the respective
differences are again \(-pk/2,-(p+1)k/6\).
For the pure exceptions use the order-binomial
residues from (11) with \((p-3)!\equiv-1/2\),
retaining \(G_{p-3}(t)\).

The list is exhaustive through \(2+(\ell+1)s\) before
normalization. At the critical group indices below
\(h-2\) gain more than \(2/p\). At the next \(L\)
only \(i=h\) or the pure exception can contribute
by the bounds from (3); further groups are higher.
Congruence errors behind \(u^2,u^4\) have relative
gains at least \(1-4r>1/p\).
Before the critical group, \(e>s\) raises every mixed
term at \(L=\ell+1\), and the mixed terms at \(L=\ell\)
have bound \(U+(\ell-i)s+i\Delta>(\ell+1)s\)
from the periodic \(R=2p-2\) identity.
For \(J_2\) the \(L=1,p+1\) cancellation leaves
\((2+\kappa)(k+2d)\) and errors at least
\(\min(s+U,U+\Delta,1+s)\), all above \((\ell+1)s\)
by (14). The top pairs are higher
by (12) and (14), and the extreme mixed terms
gain an inner binomial or factorial factor \(p\).
For (13) at the lower precision \(2+b\),
only the first two rows are needed, and the same
bounds hold at \(p=4\ell+3\), using (10) for the top
pair instead of (14).

Since \(2r>1/p\), the geometric sums satisfy
\(G_{2p-2}(t)=-1+u+o_v(1/p)\) and
\(G_{p-3}(t)=-2+3u+o_v(1/p)\).
At valuation zero the last three rows of (17) sum to
\(-5/216\) for \(\ell=1\), or \(-7/1200\) for \(\ell=2\),
using \(k/\alpha^2\equiv-2,\ \bar A=-C\), and (13).
Substitute (16) into (17) and multiply by \(\delta\).
It follows that
\[
 \frac{\delta}{\alpha}
      =S(1+u/2)+
       \begin{cases}
          (53/192)\alpha,&\ell=1,\\
          (149/768)\alpha,&\ell=2
        \end{cases}
          +o_v(1/p).                                       \tag{18}
\]
The two correction coefficients are \(1/4+5/192\)
and \(1/6+7/256\): the first summand comes from (16)
and the second from the last three rows of (17).

Use the lower weighted equation after division by
\(D a_\ell p^2 k^{\ell+1}\), at precision \(1/p\).
Its middle terms are
\[
 \mu_\ell \frac{R_h}{\alpha}
       +\omega_\ell R_{p-1}
       +\frac{\mu_\ell}{2\ell+3}\frac{k}{\alpha}
       +o_v(1/p),\qquad
              \omega_\ell=-\frac1{(\ell+1)(2\ell)!}.          \tag{19}
\]
Here the third term is at \(h-1\): the factorial
coefficient ratio apart from the extra \(k\)
is \(-1/[2(2\ell+3)]\), and \(\bar R_{h-1}=-2\).
All other middle indices are higher by (8).
By (6),
\[
 R_h=(1-3u/2)\delta-\delta^2+o_v(2/p),\qquad
 R_{p-1}=-1+u+2\delta+o_v(1/p).
\]
Indeed the rank term \(G_{p-1}(t)\) of \(R_h\)
has valuation \(U>2/p\); its first \(\delta\)
coefficient is \(1-3u/2+o_v(1/p)\), and the
quadratic coefficient reduces to \(-1\).
For \(R_{p-1}\) the rank term is \(G_{2p-2}(t)\)
and the linear \(\delta\)-coefficient reduces to \(2\).

When \(p>4\ell+5\), both exterior main terms of (7)
are higher than \(b+2/p\) by \(U>b+2/p\) and the
upper row. The coefficient errors from the exterior
bands are also higher, as their valuation before
division by \(k^{\ell+1}\) is at least one.
Then (19), with \(\omega_\ell/\mu_\ell=S\), requires
\[
       \delta/\alpha=S(1+u/2)
                    +\bigl(2/(2\ell+3)-S^2\bigr)\alpha
                    +o_v(1/p).
\]
Its \(\alpha\)-coefficient differs from (18)
by \(1/60\) for \(\ell=1\) or \(1/168\) for \(\ell=2\).
Both differences are units at the primes in scope.

The remaining prime in this branch is \(p=13,\ell=2\),
where \(U=b+2/p=7/13\).
The upper weighted row gives \(7(\kappa+2)+2\psi=0\)
at precision \(U\). The normalized lower equation (19)
therefore has the extra top contribution
\(-2(\kappa+2)/k^3\).
Since \(\mu_2=-2/45\) and \(k/\alpha^2\equiv-2\),
the discrepancy \(+\alpha/168\) in (18) would require
\[
         \overline{\psi/\alpha^7}
                         =\mu_2/12=-1/270=9
                         \quad\text{in }\mathbb F_{13}.
\]
But \(\psi\) has leading value \(2u^{12}\),
and \(u^{12}=d^6/\delta^6,\ p=\kappa^6\alpha^{13}\).
Here \(\overline{d^6/p}=-1,\ \bar\kappa^6=-1\), and
\(\overline{\delta/\alpha}=5/16=6\) with \(6^6=-1\).
Thus \(\overline{\psi/\alpha^7}=-2=11\), a contradiction.
