# Four factorial moments of an even binary quadratic with a prime penultimate order

Let \(\Lambda_2(x^i y^j)=i!j!\), and let
\(f=ax^2+bxy+cy^2+k\) with \(abck\ne0\) over a
characteristic-zero field. For every prime \(p\ge5\), the four
moments at \(p-2,p-1,p,p+1\) cannot all vanish.
In particular, \(p=13\) excludes the complete window at
orders \(11,12,13,14\), irrespective of the residue at \(11\).

Use the
[weighted boundary identities](BINARY_EVEN_QUADRATIC_SINGLE_FACE_STRONG_FACTORIAL.md#4-the-exact-reduction-for-general-b-and-its-limitation)
and the reduction of the
[prime-terminal theorem](BINARY_EVEN_QUADRATIC_PRIME_TERMINAL_FOUR_WINDOW.md#reduction-at-the-terminal-prime).
In particular, a putative zero block gives an algebraic solution,
and we work over \(\overline{\mathbb Q}_p\) with \(v(p)=1\).
The normalization using \(M_p=M_{p-1}=0\) gives
\[
 a=1,\quad p=2h+1,\quad t=b/2,\quad u=t-1,\quad
 d=c-t^2,\quad D=d+u^2,\quad s=v(k)>0,\quad v(d)>0 .
\]
Here \(M_j=\Lambda_2(f^j)\).
We may assume \(c\ne1,\ c\ne t,\ u\ne0,\ D\ne0,\ d\ne0\) by the
[single-face](BINARY_EVEN_QUADRATIC_SINGLE_FACE_STRONG_FACTORIAL.md)
and [critical boundary](BINARY_EVEN_QUADRATIC_CRITICAL_BOUNDARY_WINDOWS.md)
theorems. Put \(B_j(Z)=\sum_i\binom ji(2i)!Z^i k^{j-i}\).
The necessary weighted boundary values
\[
       T_j=(c-t)B_j(c)+(1-t)B_j(1)
\]
are zero for \(j=p,p+1\). Write \(\overline q\) for residue
reduction of an integral element. In the expansions below,
\(o_v(A)\) denotes a term of valuation strictly greater than \(A\).
A rational coefficient replaced by its residue incurs valuation
at least one before multiplication by the displayed quantities.

## The common boundary strength

At unit arguments the lowest terms of \(B_p\) change from index
\(0\) to index \(p\) as \(s\) increases through \(2/p\).
For \(B_{p+1}\) the lowest terms change from index \(1\)
to index \(p+1\) at the same transition. Indeed, the terms
\(2\le i\le h\) of the latter have valuations \(1+(p+1-i)s\).
They exceed \(ps\) when \(s\le2/p\), and exceed two for
\(s\ge2/p\). For \(h<i\le p\), the valuation is at least
\(2+(p+1-i)s\). At the transition only the two named terms
are lowest in each polynomial.

If \(\bar t=0\), the weight \(1-t\) is a unit and \(c\) is
nonunit. Below the transition the constant term of \(B_p\)
forces the unit weight sum \(D\), and above it \(B_p(1)\)
strictly dominates \(B_p(c)\). Both contradict \(T_p=0\).
If \(\bar t\ne0,1\), a unique lowest monomial of index \(i\)
in either boundary identity forces \(\bar t^{2i+1}=1\).
Below the transition use index zero, and above it use indices
\(p,p+1\), whose odd exponents are coprime. If \(\bar t=1\),
the normalized quotient \(B_{p+1}/B_p\) reduces to a nonzero
constant times \(Z\) on the disk of \(1\), strictly injecting
that disk by the
[residue-disk divided-difference argument](BINARY_EVEN_QUADRATIC_PRIME_TERMINAL_FOUR_WINDOW.md#the-two-boundary-transitions).
This contradicts equality at \(1,c\).

It follows that \(s=2/p\). With \(\kappa=k^p/p^2\), the two
normalized reductions are
\[
        B_p(Z)/p^2\equiv \bar\kappa+2Z^p,\qquad
        B_{p+1}(Z)/p^2\equiv 2Z(\bar\kappa+2Z^p)\pmod{\mathfrak m}.
\]
If \(\bar t\notin\{0,1,-1\}\), the two reduced boundary
equations force the normalized lower polynomial to vanish
at both \(1,\bar c\), since its two coefficients in the
weighted equations have nonzero determinant. Frobenius
then gives \(\bar c=1\), contrary to \(\bar t^2=\bar c\).
For \(\bar t=0,-1\), the two reduced weights give
\(\bar\kappa=-2\) from \(T_p\). For \(\bar t=1\), the same
holds because otherwise the ratio reduces to \(2Z\)
on the disk of \(1\) and strictly injects it.

## The two off-diagonal residues

Use the
[homogeneous-perturbation formulas](BINARY_EVEN_QUADRATIC_PRIME_TERMINAL_FOUR_WINDOW.md#the-off-diagonal-case)
for \(f=(x+ty)^2+d y^2+k\). In a moment of order \(j\),
the term indexed by total perturbation \(L\) and \(d\)-power \(i\)
is
\[
 \binom jL\binom L i k^{L-i}d^i R!\,S_{R,i}(t),\qquad
 R=2(j-L),\quad S_{R,i}(t)=\sum_{\ell=0}^R(\ell+1)_{2i}t^\ell .
\]
For \(R=p-1\) and \(i\le h\),
\(S_{p-1,i}(t)\equiv(2i)!(t-1)^{p-1-2i}\pmod p\).
For \(i>h\) it vanishes modulo \(p\) for every \(R\).

When \(\bar t=0\) or \(-1\), the equation \(M_{p-2}=0\)
forces \(v(d)=s\): otherwise the unique pure term at
\(L=h-1\) using the smaller valuation has lowest valuation,
less than one. The first relevant groups of \(M_{p-2},M_{p-1}\)
then give \(g_{h-1}(\bar X)=g_h(\bar X)=0\), where
\[
 X=d/[ku^2],\qquad
               g_i(X)=\sum_{\ell=0}^i\binom i\ell(2\ell)!X^\ell .
\]
The first root is simple and the second has derivative zero
there, by the recurrence and derivative identities in the
[prime-terminal proof](BINARY_EVEN_QUADRATIC_PRIME_TERMINAL_FOUR_WINDOW.md#the-off-diagonal-case).
The other terms of \(M_{p-2}\) are of valuation at least
\(hs\), so \(v(g_{h-1}(X))\ge s\). Lift its simple root in
the finite unramified extension with the required residue
field. At the lift \(g_h\) and its derivative have valuation
at least one. Since \(2s<1\), Taylor expansion gives
\(v(g_h(X))\ge2s\). Thus the group \(L=h\) in \(M_{p-1}\)
is strictly above valuation one, as are all later groups
and the earlier groups with \(L>0\). The \(L=0\) term
has exact valuation one because \(G_{2p-2}(\bar t)=1\) for
\(\bar t=0,-1\), with \(G_m(t)=\sum_{\ell=0}^m t^\ell\).
This contradiction leaves \(\bar t=1\).

## Valuations in the diagonal case

Set \(r=v(u)>0\), \(e=v(d)>0\), and \(\rho=v(c-1)>0\).
At \(j=h+L_0\) and \(L=L_0\), the perturbation formula
gives an initial valuation
\[
             A_{L_0}=2hr+L_0\min(s,e-2r).                  \tag{1}
\]
Every error in this group has valuation at least one, strictly
greater than one for \(L_0>0\). For \(L<L_0\) the factorial
has valuation at least one, with a strict gain if \(L>0\).
For \(L>L_0\), put \(v_0=L-L_0\). If \(v_0\le i\le h\),
the lower bound for the term is
\[
                         2hr+(L-i)s+i(e-2r)                \tag{2}
\]
up to an error above valuation one. Every removed term in the
comparison of \(S_{R,i}\) with \(S_{p-1,i}\) contains a factor
\(p\), because its last rising-factorial index is at least \(p\).
If \(i>h\) the terms are above one. If \(i<v_0\) the lower
bound is \((L_0+1)s\); equality is possible only at
\(L=L_0+1,i=0\), where it holds for \(L_0=h-1\).
These are the original [diagonal group estimates](BINARY_EVEN_QUADRATIC_PRIME_TERMINAL_FOUR_WINDOW.md#the-diagonal-transition-s2p),
and apply before deciding the sign of \(e-2r\).

If \(e\le2r\), \(M_{p-1}=0\) forces \(he\ge1\). Otherwise
the pure \(d\)-term at \(L=i=h\) has valuation \(he\),
and all later terms with \(i\le h\) have bound at least
\(he+(L-h)s\) or \((h+1)s>1\). If \(e=2r\), every term
in (2) for \(M_{p-2}\), and its critical \(L_0=h-1\) group,
is at least \(he\ge1\). The isolated term \(L=h,i=0\)
at \(hs<1\) contradicts vanishing.

If \(e<2r\), \(v(D)=e\), and \(\rho\ge\min(r,e)\ge1/(p-1)\).
Write \(B_j[c,1]=(B_j(c)-B_j(1))/(c-1)\). Then
\[
 v(B_p[c,1])=2+1/p,\qquad
 v(B_p(c))=r+\rho+2+1/p-e .
\]
The first equality uses the unique lowest derivative contribution
at index \(h\). The top divided power has valuation at least
one because \((p-1)\rho\ge1\); all other nonconstant
coefficients have larger valuations. The second equality is
\(T_p=0\) in the form
\(D B_p(c)=-u(c-1)B_p[c,1]\).
The identity
\[
 B_{p+1}(Z)=k^{p+1}+A Z B_p(Z)-C k Z B_{p-1}(Z),\qquad
 A=2(p+1)(2p+1),\quad C=4p(p+1),
\]
follows from the binomial coefficients. After using \(T_p=0\),
the equation \(T_{p+1}=0\) becomes
\[
 0=D k^{p+1}+A(c-t)(c-1)B_p(c)-C k B_{p-1}(c)\Pi,\qquad
 \Pi=D+(c-1)\bigl(c-t+u B_{p-1}[c,1]/B_{p-1}(c)\bigr).      \tag{3}
\]
The last quotient has residue \(h\), since \(B_{p-1}\) has
unique lowest index \(h\). The correction to \(D\) in \(\Pi\)
has valuation at least \(2\min(e,r)>e\). The last term of
(3) has valuation \(2+1/p+e\). The first term is higher.
The middle term is also higher because
\(v(c-t)+2\rho+r\ge3\min(e,r)+r>2e\).
Thus \(e>2r\) is necessary.

Put \(q=\min(e,s+2r)\). In the distinct case the \(M_{p-1}\)
equation now forces \(hq=1\), by (1), (2) and the unique
\(L=0\) valuation-one term when \(hq>1\). Similarly
\(M_{p-2}\) forces \(2r+(h-1)q=hs\). Therefore
\[
                 e=1/h,\qquad e-2r=1/p,\qquad
                 r=(p+1)/(2p(p-1)).                         \tag{4}
\]
In the tied case \(e=s+2r\), put \(U=2hr\).
The \(M_{p-1}\) equation requires \(hs+U\le1\).
If strict, \(M_{p-2},M_{p-1}\) first force a common reduced
root of \(g_{h-1},g_h\). The error after the \(g_{h-1}\)
group is at least \(hs\), giving \(v(g_{h-1}(X))\ge s-U>1/p\).
The same unramified lift gives \(v(g_h(X))\ge2(s-U)\).
This pushes the \(g_h\) group above one, leaving the
\(L=0\) contribution uncancelled in \(M_{p-1}\).
Thus the tied case requires \(U=1/p\) and \(r=1/[p(p-1)]\).

For this tied possibility use
\[
 L=\binom{p-1}h(p-1)!k^h,\quad
 H=\frac{(c-t)u(c-1)^2}{D^2},\quad
 R_i=\frac{(c-t)c^{i+1}+(1-t)}{D},\quad
 \Delta_i=\frac{c^i-1}{c-1}.
\]
In (3), eliminate \(B_p(c)\) via \(T_p=0\) and divide by
\(D\). The contributions at valuation \(2+1/p\) from
the index \(h\) in \(B_p[c,1]\) and \(B_{p-1}\) cancel:
their coefficients relative to \(p kL\) are
\(-A H\Delta_h/(p-h)\) and \(-4(p+1)R_h\), with
residues \(8,-8\). Here \(\bar H=4,\bar R_h=2\).
The top term \(-A H(2p)!\Delta_p\) has exactly this valuation,
because \(\Delta_p=(c-1)^{p-1}\) plus terms of valuation
at least one. All other terms are strictly higher.
This excludes the tied case.

## Resolving the distinct case

Assume (4), and put
\[
       \alpha=p/k^h,\qquad \delta=d/u^2,\qquad
       U_0=d^h/p,\qquad V_0=u^2d^{h-1}/k^h,\qquad C_0=(-1)^h .
\]
Then \(v(\alpha)=v(\delta)=1/p\), \(k=\kappa\alpha^2\),
and \(\delta=\alpha U_0/V_0\).
First expand the two earlier moments to their first correction:
\[
 \begin{split}
 U_0&= C_0G_{2p-2}(t)-\alpha+o_v(1/p),\\
 V_0&=G_{p-3}(t)+o_v(1/p).
 \end{split}                                                 \tag{5}
\]
Here are all the contributions at the retained orders. In
\(M_{p-1}/p\) the leading terms are
\(-G_{2p-2}(t)+C_0 U_0\).
The \((L,i)=(h,h-1),(h+1,0)\) terms are respectively
\(-C_0 V_0 k/(4\alpha)\) and \(-C_0 k/\alpha\)
at valuation \(1/p\). Their sum at that order is \(C_0\alpha\)
using \(\overline{V_0}=-2,\bar\kappa=-2\).
In \(M_{p-2}/k^h\) the leading terms are
\(C_0(V_0-G_{p-3}(t))/4\).
The \((L,i)=(h-1,h-2),(h,h),(0,0)\) terms are respectively
\(-C_0 V_0 k/(32\delta), C_0 U_0\alpha/4,\alpha/2\)
at valuation \(1/p\). The leading equations give
\(\overline{U_0}=-C_0\) and \(\overline{V_0}=-2\), so
\(\overline{\delta/\alpha}=C_0/2\). The three terms give \(-\alpha/4,-\alpha/4,\alpha/2\)
and cancel. In each list earlier positive \(L\) below the
critical group contributes a factorial factor \(p\) and at
least one power of \(k\) or \(d\). Inside the critical group
each successive index below its maximum adds \(1/p\).
Above the critical group (2) becomes
\(2hr+(L-i)s+i/p\), minimized at the largest possible \(i\)
and strictly increasing with \(L\) once \(i=h\).
For \(i<v_0\), only the named pure exceptional term can
contribute, since another power of \(d\) adds \(e>1/p\)
or an additional power of \(k\) adds \(s>1/p\).
This verifies that no other indices enter (5).
Consequently, since \(2r>1/p\),
\[
 \delta=C_0\alpha(1/2+u/4)+\alpha^2/2+o_v(2/p),              \tag{6}
\]
using \(G_{2p-2}(t)=-1+u+o_v(1/p)\) and
\(G_{p-3}(t)=-2+3u+o_v(1/p)\).

Next apply \(T_p=0\) to improve \(\kappa=-2+o_v(1/p)\).
Modulo terms strictly above \(1/p\), \(B_p(c)/p^2\) is
\(\kappa+2+4C_0\alpha\), since the only coefficient
beyond the constant and top affected at this order is at
\(h\), and \(c^p-1\) is above \(1/p\).
At this order \(B_p[c,1]/p^2=-2C_0\alpha+o_v(1/p)\).
The top divided-difference term is already higher since
\((p-1)r>1/p\). Thus the exact relation
\(B_p(c)=-u(c-1)B_p[c,1]/D\) gives the claim.

Divide (3) by \(DpkL\) and again eliminate \(B_p(c)\).
Let \(\epsilon=(2p-2)!/L\), and let \(\nu\) be the ratio
of the \(h-1\) coefficient of \(B_{p-1}\) to \(L\).
Thus
\(\epsilon=C_0\alpha+o_v(2/p)\) and
\(\nu=-k/2+o_v(2/p)\).
The only indices of \(B_{p-1}\) in this precision are
\(h,p-1,h-1\). Indeed their valuations relative to \(L\)
are \(0,1/p,2/p\) and each next coefficient away from
those indices adds at least \(s=2/p\). For \(i<p\)
the coefficient of \(B_p\) is \(pk/(p-i)\) times
the coefficient of \(B_{p-1}\), with \(p-i\) a unit.
The reduced equation consequently has the shape
\[
 0=\frac{\kappa\alpha}{\binom{p-1}h(p-1)!}
      +F_h+\epsilon F_{p-1}+\nu F_{h-1}
      -A H\,\frac{(2p)!}{p kL}\Delta_p+o_v(2/p),              \tag{7}
\]
where \(F_i=-A H\Delta_i/(p-i)-4(p+1)R_i\).
All \(R_i\) are integral by
\(R_i=1+(c-t)(c-1)\Delta_{i+1}/D\).
The terms in \(\Delta_p\) other than \((c-1)^{p-1}\)
have valuation at least one.

For completeness the required expansions of these rational
functions are
\[
 \begin{split}
 F_h&=(-12+6u)\delta+20\delta^2+o_v(2/p),\\
 \epsilon F_{p-1}&=(4-16\delta)\epsilon+o_v(2/p),\\
 \nu F_{h-1}&=8\nu+o_v(2/p).
 \end{split}                                                 \tag{8}
\]
They follow by writing \(c=t^2+u^2\delta\) and
\[
 H=\frac{(t+u\delta)(t+1+u\delta)^2}{(1+\delta)^2},
 \qquad
 R_i=1+\frac{(t+u\delta)(t+1+u\delta)}{1+\delta}\Delta_{i+1}.
\]
At \(\delta=0\), \(F_h\) reduces modulo \(p\) to
\(-4(t^2+t+1)u^{p-1}\), using
\(t(t+1)G_{p-2}(t)+G_{p+1}(t)=(t^2+t+1)u^{p-1}\)
in characteristic \(p\). This has valuation greater than \(2/p\).
The coefficient of \(\delta\) through first order in \(u\)
is \(-12+6u\), and the coefficient of \(\delta^2\) at
\(u=0\) is \(20\). The constant-in-\(\delta\) value of
\(F_{p-1}\) is \(4+O(u^2)\) modulo \(p\), and its coefficient
of \(\delta\) at \(u=0\) is \(-16\). The constant term of
\(F_{h-1}\) is \(8\). Higher terms exceed \(2/p\) because
\(2r>1/p\) and \(v(\delta)=1/p\).

By (6)–(8), the terms except the explicit top contribution
in (7) sum to \(-\alpha^2+o_v(2/p)\). Indeed
\((-12+6u)C_0\alpha(1/2+u/4)+6C_0\alpha\) has valuation
greater than \(2/p\), and the coefficient of \(\alpha^2\)
is \(-6+5-8+8=-1\). The top contribution has valuation
\((p-1)r-1/p=(p-1)/(2p)\) in (7). This exceeds \(2/p\)
when \(p\ge7\), proving the contradiction there.

When \(p=5\), \(h=2,\ C_0=1\), and the top contribution
is also at \(2/p\). At the leading units,
\(p=\kappa^2\alpha^5\) and \(u^4=d^2/\delta^2\), so
\(u^4=-16\alpha^3+o_v(3/p)\) by (5)–(6).
Also \((2p)!/(pkL)=1/\alpha+o_v(-1/p)\),
\(AH\) has residue \(8\), and \(\Delta_p\) has the leading
term \((2u)^4\). Hence this term of (7) contributes
\(3\alpha^2+o_v(2/p)\) modulo five. The sum in (7)
has nonzero coefficient \(2\) at valuation \(2/p\).
This excludes the last case.
