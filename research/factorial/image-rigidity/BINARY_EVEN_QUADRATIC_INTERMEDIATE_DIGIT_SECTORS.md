# Intermediate two-adic sectors for even binary quadratics

Normalize the nonzero constant coefficient to one. Put
\[
 f=ax^2+bxy+cy^2+1,\qquad M_j=\Lambda_2(f^j),\qquad
 A=v_2(4a),\quad C=v_2(4c),\quad Q=v_2(4b).
\]
The valuation may be extended to any valued characteristic-zero coefficient field, with \(v_2(2)=1\).

**Theorem.** No four consecutive moments starting at a positive index vanish in either of the following sectors:
\[
        \frac14<\min(A,C)<\frac12,
             \qquad\text{or}\qquad
        A,C\ge\frac12,\quad Q>\frac14.                         \tag{1}
\]
The first sector imposes no restriction on \(b\) and does not bound the larger of \(A,C\). The second sector uses only the original interior moments. The first combines the exact coupled boundary equations with the original interior moments.

## Coefficient valuations and a useful integer ratio

Use \(s(d)\) for the binary digit sum. At order \(j\), let \(i,\ell,h\) count the occurrences of \(a,b,c\), and put \(k=i+\ell+h\). For any nonconstant term the valuation is at least
\[
             iA+\ell Q+hC-s(\ell)+s(j-k)-s(j).                \tag{2}
\]
In fact the difference between the exact valuation and (2) is the sum of the two binary carry counts in the additions \(i+\lfloor\ell/2\rfloor\) and \(h+\lfloor\ell/2\rfloor\). This follows from the multinomial formula and the factorial valuations \(v_2(d!)=d-s(d)\).

There is also a direct comparison at each fixed \(k\). After taking out \(\binom jk(k!)^2\), the remaining coefficient on \(a^i b^\ell c^h\) is
\[
                   \frac{(2i+\ell)!(2h+\ell)!}{k!\,i!\,\ell!\,h!},
                                                               \tag{3}
\]
which is an integer. To see this use the factorial floor formula at each prime power. After replacing \(i,h,\ell\) by their remainders at that power, the terms from the quotients cancel. The lookup is then the sum of the floors of \(2x+z,2y+z\), for \(x,y,z\in[0,1)\), compared with the floor of \(x+y+z\). Their average is \(x+y+z\); the larger floor alone is at least its floor. This proves integrality at every prime power. Thus when \(Q\le A,C\), the pure \(b^k\) term has valuation no larger than any other term at total degree \(k\), and is uniquely least at that degree if the inequalities are strict.

For digit estimates, if \(j\equiv u\bmod 8\) and \(k=8d+w\) with \(0\le w<8\), write
\[
 \delta_u(k)=
 \begin{cases}
 s(u-w)-s(u)-s(d),&w\le u,\\
 s(8+u-w)-s(u)-s(d+1),&w>u.
 \end{cases}
\]
Then \(s(j-k)-s(j)\ge\delta_u(k)\) whenever \(k\le j\). Split \(j\) and \(j-k\) into their last three binary digits, and apply subadditivity of the digit sum to the upper parts. The analogous identity with last two digits will be used as well.

## The open axis interval

Write \(g_j(a)=\Lambda_1((1+a x^2)^j)\). In \(g_j(Z/4)\), the valuation of its coefficient at \(Z^i\) is \(s(j-i)-s(j)\). If \(r=v_2(Z)\) lies in \((1/4,1/2)\), the unique least term is \(i=j\bmod 4\). For a direct proof write \(j=4T+u\), \(i=4d+w\). If \(w\le u\), the valuation of the term relative to index \(u\) is at least
\[
                 (4r-1)d+(w-u)r+s(u-w).
\]
It is strictly positive unless \(d=0,w=u\), since each nonzero integer \(u-w\le3\) has digit sum strictly greater than \((u-w)r\). If \(w>u\), the relative lower bound is \((4r-1)d+(w-u)r+s(4+u-w)-1>0\).

Suppose first that \(1/4<A<C\) with \(A<1/2\). Among the two final boundary orders \(n+2,n+3\), the odd one has a negative least valuation for \(g_j(a)\). All terms of positive degree for \(g_j(c)\) raise their valuation strictly relative to those for \(g_j(a)\), and the constant has valuation zero, so \(g_j(a)\) dominates \(g_j(c)\) strictly. The [weighted boundary equation](BINARY_EVEN_QUADRATIC_SINGLE_FACE_STRONG_FACTORIAL.md#4-the-exact-reduction-for-general-b-and-its-limitation) at that order forces
\[
    \varepsilon=b/(2a)-1=\frac{c-a}{a}\frac{g_j(c)}{g_j(a)+g_j(c)},
                        \qquad v_2(\varepsilon)>0 .           \tag{4}
\]
The case \(\varepsilon=0\) is included with valuation \(+\infty\).

Let \(N\) be the multiple of four in the proposed window. The perturbation from replacing \(b\) by \(2a\) at \(M_N\) has positive valuation. In (2), use \(v_2(8a)=A+1\) for a term in the comparison plane, and add \(v_2(\varepsilon)\) for the changed cross factor power. The remaining lower bound is at least \(k A+s(N-k)-s(N)\), because \(\ell-s(\ell)\ge0\) and \(C>A\). It is nonnegative, with strict inequality for \(k>0\). One can check this by the last two digits with \(N\equiv0\bmod4\): for \(k=4d+w\), the bound is at least \((4A-1)d\) if \(w=0\), and \((4A-1)d+wA+s(4-w)-1>0\) otherwise.

The moment on the comparison plane \(b=2a\) is a unit at \(N\). The [exact plane weight identity](BINARY_EVEN_QUADRATIC_COUPLED_TWO_ADIC_SECTORS.md#a-dominant-axis-and-approximation-by-its-single-face-plane) gives, for total degree \(d=i+h\), the term
\(\binom Nd(d!)^2(4a)^i c^h\binom{2h}{h}\).
In a lower bound one may replace \(v_2(c)\) by \(c_0=\min(v_2(c),A-1)\). Since \(A-c_0\ge1\), the function \(h\mapsto h(A-c_0)-s(h)\) is nondecreasing. Thus the smallest valuation bound at fixed \(d\) is attained with all \(h=d\), and is at least \(d(2+c_0)+s(N-d)-s(N)\). Here \(2+c_0=\min(C,1+A)>1/4\), so it is positive for \(d>0\) by the same last-two-digit calculation. The constant term is one. This and (4) exclude the unequal-axis case.

It remains to consider \(A=C=r\in(1/4,1/2)\). Put \(j=n+2\). The two boundary equations force the same value of \(g_{j+1}/g_j\) at \(a,c\). On the annulus of valuation \(r\) for \(Z=4a\), the quotient has leading monomial exponent \(+1\) for \(j\equiv0,1,2\bmod4\), and \(-3\) for \(j\equiv3\bmod4\). More precisely, choose a scale of valuation \(r\); divide each polynomial by its unique leading coefficient times its leading monomial. The resulting Laurent polynomials on the unit annulus are \(1\) plus terms of uniformly positive valuation. A quotient of two such polynomials has the same structure by geometric expansion. Differences of that correction at unit points \(X,Y\) have valuation strictly greater than \(v_2(X-Y)\), by factoring each difference of Laurent monomials. The same is true after multiplying the correction by \(X^{+1}\) or \(X^{-3}\).

The \(+1\) map is therefore injective. For exponent \(-3\) it is injective on each unit residue class, since \((X^{-3}-Y^{-3})/(X-Y)\) is a unit when the residues agree. If \(a\ne c\), equality is possible only when \(j\equiv3\bmod4\) and the residue of \(c/a\) is a nontrivial cube root of one. We can take \(a\ne c\), since \(a=c\) is a [known four-window plane](BINARY_EVEN_QUADRATIC_SINGLE_FACE_STRONG_FACTORIAL.md#3-the-plane-ac).

Now \(N=n+3\) is divisible by four. If \(N\equiv4\bmod8\), the unique least nonconstant term of \(g_N(Z/4)\) on the annulus is index four, of valuation \(4r-1<1\). For \(i=4d+w\), nonzero \(w\) gives lower bound \(d(4r-1)+w r+s(4-w)-1>4r-1\); for \(w=0,d\ge2\) the bound is at least \(d(4r-1)\). The residues of the fourth powers of the two scaled axes are distinct. Thus the sum \(g_N(a)+g_N(c)\) has valuation \(4r-1\), since its constant part is two. The numerator \(2(a g_N(a)+c g_N(c))\) has valuation \(r-1\), because \(g_N(a),g_N(c)\) reduce to one and the scaled axes have different residues. The weighted equation therefore gives \(v_2(b)=-3r\), or \(Q=2-3r\).

At \(n=N-3\equiv1\bmod8\) the pure linear \(n b\) term is uniquely least. To verify it, subtract its valuation \(-3r\) in (2). The resulting lower bound is
\[
               (k-4\ell+3)r+2\ell-s(\ell)+\delta_1(k).
\]
At \(r=1/4\) this is at least \((k+3)/4+\delta_1(k)\), nonnegative by writing \(k=8d+w\) and using \(s(d)\le d\), \(s(d+1)\le d+1\). At \(r=1/2\) it is \((k+3)/2-s(\ell)+\delta_1(k)\), which is strictly positive except at the selected term. For \(k\ge9\), use \(s(\ell)\le\lfloor\log_2(k+1)\rfloor\le k/3\) and \(\delta_1(k)\ge-\lceil(k-1)/8\rceil-1\), obtaining at least \(k/24-1/4>0\). For \(k=2,\ldots,8\) the lower bounds are respectively \(5/2,1,3/2,1,5/2,1,3/2\); at \(k=1,\ell=0\) it is one. Affine interpolation gives strict positivity at the actual \(r\). The constant term is above the pure \(n b\) term as well.

If \(8\mid N\), every nonconstant term of \(g_N(Z/4)\) has valuation greater than one. Write \(i=8d+w\): if \(w=0\) its valuation is at least \(d(8r-1)>1\), and otherwise it is at least \(d(8r-1)+w r+s(8-w)-1>1\). Thus the axis sum has valuation one, and the numerator has valuation \(r-1\) as before. Hence \(Q=r\). Use \(j=N-1\equiv7\bmod8\), where the pure \(b^7\) term has valuation \(7r-6\). From (3) and \(\delta_7(k)=-s(k)\), the margin of any term at total degree \(k\) over this term is at least
\[
                             (k-7)r+6-2s(k).
\]
For \(0\le k<7\) use the endpoint \(r=1/2\), which gives a nonnegative bound and hence a strictly positive margin in the open interval. For \(k=8d+w\ge8\), the endpoint \(r=1/4\) bounds the margin below by \(17/4+w/4-2s(w)\ge0\), again strict for the actual \(r\). At \(k=7\), a term with \(\ell<7\) improves (2) over the pure term by at least one since \(s(\ell)\le2\). This completes the first sector of (1).

The unequal-axis proof also allows \(A=1/4<C\) or \(A=1/2<C\) if an odd boundary value at \(a\) has its lowest possible valuation. In particular this holds if the residue of \(4a/t\), for \(t^4=2\) or \(t^2=2\), respectively, is different from one. At the endpoint only the corresponding binary-bit face can tie the leading term, and its two unit terms have residue quotient \((4a/t)^4\) or \((4a/t)^2\). The unequal-axis argument then uses only nonnegativity for the perturbation bound and strict positivity for the comparison plane, as before.

## Axis coefficients of valuation at least \(-3/2\)

Assume \(A,C\ge1/2\), and write \(q=Q>1/4\). Let \(J\) be the multiple of four in the window. Under the assumptions, we show that \(M_J\) is nonzero or exhibit a neighboring nonzero moment in that same window.

### When \(J\equiv4\bmod8\)

For \(1/4<q\le1/2\), compare at each fixed \(k\) with the pure \(b^k\) term by (3). In order \(J\) the pure \(b^3\) term has valuation \(3q-2\). The margin of the pure term at a different positive \(k\) is bounded below by
\[
                         (k-3)q+2-s(k)+\delta_4(k).            \tag{5}
\]
For \(k=1,2\) this is \(2-2q,1-q\), respectively. For \(k\ge3\) use \(q=1/4\). If \(k=8d+w,\ w\le4\), the resulting bound is at least \((w-3)/4+1-s(w)+s(4-w)\); for \(w>4\) it is at least \((w-3)/4-s(w)+s(12-w)\). These are nonnegative. Only \(w=3,7\) can give zero, and all \(k>3\) improve strictly for the actual \(q>1/4\). Mixed terms at \(k=3\) are strictly higher because \(Q<A,C\) if \(q<1/2\), and because \(s(\ell)\le1\) for \(\ell<3\) at \(q=1/2\). The constant is higher.


### When \(J\equiv4\bmod8\)

For \(1/4<q\le1/2\), compare at each fixed \(k\) with the pure \(b^k\) term by (3). In order \(J\) the pure \(b^3\) term has valuation \(3q-2\). The margin of the pure term at a different positive \(k\) is bounded below by
\[
                         (k-3)q+2-s(k)+\delta_4(k).            \tag{5}
\]
For \(k=1,2\) this is \(2-2q,1-q\), respectively. For \(k\ge3\) use \(q=1/4\). If \(k=8d+w,\ w\le4\), the resulting bound is at least \((w-3)/4+1-s(w)+s(4-w)\); for \(w>4\) it is at least \((w-3)/4-s(w)+s(12-w)\). These are nonnegative. Only \(w=3,7\) can give zero, and all \(k>3\) improve strictly for the actual \(q>1/4\). Mixed terms at \(k=3\) are strictly higher because \(Q<A,C\) if \(q<1/2\), and because \(s(\ell)\le1\) for \(\ell<3\) at \(q=1/2\). The constant is higher.

For \(1/2\le q<2/3\), compare instead with the bounds at \(q=1/2\) and \(2/3\). The margin over \(3q-2\) increases with \(q\) if \(\ell\ge3\). If \(\ell\le2\) and \(k>0\), the bound at \(q=2/3\) is already positive, since \(k/2+\delta_4(k)\ge1\) and \(\ell/6-s(\ell)\ge-5/6\). The first inequality follows by writing \(k=8d+w\) in the definition of \(\delta_4\). Thus \(b^3\) remains uniquely least. For \(q>2/3\) every nonconstant term is strictly positive in valuation, by the same comparison and monotonicity in \(\ell q\). The constant is then uniquely least.

At \(q=2/3\) use \(j=J+1\) if \(J=n\), \(j=J-3\) if \(J=n+3\), and \(j=J-1\) otherwise. These have residues \(5,1,3\bmod8\), respectively. The pure \(b^1,b^1,b^3\) terms, respectively, are uniquely least. Here is a finite digit comparison for this assertion. At \(j\equiv u\bmod8\), use the bound
\[
                           k/2+\ell/6-s(\ell)+\delta_u(k).
\]
For \(k\ge16\), writing \(k=8d+w\) gives a lower bound \(-1\): use \(s(\ell)\le\lfloor\log_2(k+1)\rfloor\le d+3\), \(s(d),s(d+1)\le d\), and \(d\ge2\). For smaller positive \(k\) the following are the minima of that bound minus the valuation of the selected term, omitting the selected pair \((k,\ell)\).
\[
\begin{array}{c|c|c|c|c}
 u & (k,\ell) & \text{selected valuation}& 1\le k\le7&8\le k\le15\\ \hline
 1&(1,1)&-4/3&5/6&2\\
 3&(3,3)&-2&2/3&8/3\\
 5&(1,1)&-4/3&1/3&2
\end{array}
\]
Each entry is obtained by substituting \(0\le\ell\le k\) and the last three digits in \(\delta_u(k)\). The constant is higher as well.

### When \(8\mid J\)

At \(q\ge1/2\), compare all terms with \(A=C=Q=1/2\), and then with the pure terms by (3). At \(k=8d+w\), the bound \(k/2-s(k)+\delta_0(k)\) is strictly positive for \(k>0\). For \(w=0\) it is at least \(2d\). For \(w>0\) it is at least \(2d+w/2-s(w)+s(8-w)-1>0\). Hence \(M_J\) is a unit.

For \(1/4<q<1/2\), the pure \(b^k\) term is always uniquely least at fixed degree \(k\) by (3). If \(J\ne n\), use \(j=J-1\equiv7\bmod8\); the pure \(b^7\) comparison of the preceding section applies with \(r=q\), so this moment is nonzero.

If \(J=n=8V\), use \(j=J+1=8V+1\). The valuation of the pure \(b^1\) term is \(q-2\), and that of the pure \(b^7\) term is \(7q-4+v_2(V)\). For \(1/4<q\le1/3\), compare every pure term first with \(7q-4\) using \(\delta_1(k)\). At \(k=8d+w\) the lower bound for this margin is
\[
 (8q-2)d+
 \begin{cases}
 (w-7)q+4-s(w)-w,&w\le1,\\
 (w-7)q+2-s(w)+s(9-w),&w>1 .
 \end{cases}
\]
This is positive except for \(k=7\) and possibly \(k=1\) at \(q=1/3\). For \(q>1/3\), the margin over \(q-2\) is increasing for each \(k>1\), so \(b^1\) wins. For \(q<1/3\) and \(V\) odd, \(b^7\) wins. If \(V\) is even and \(q\le1/3\), compare with \(q-2\). For \(k=2,\ldots,7\) the bounds gain at least one over \(\delta_1(k)\), since \(s(V-1)-s(V)=v_2(V)-1\ge0\). This makes the comparison with \(q-2\) strict. For \(k\ge8\), use \(s(d+1)\le d\) as \(d\ge1\); the lower bound becomes
\[
 (8q-2)d+
 \begin{cases}
 (w-1)q+2-s(w)-w,&w\le1,\\
 (w-1)q+1-s(w)+s(9-w),&w>1 ,
 \end{cases}
\]
which is strictly positive. In all these cases the constant term lies above the selected term.

The remaining equality is \(q=1/3,\ V\) odd. Use order \(J\) itself. The pure \(b^7\) term has valuation \(-2/3\), and every other pure term has nonnegative valuation. To verify this, the bound \(k/3-s(k)+\delta_0(k)\) at \(k=8d+w\) is at least \(2d/3\) if \(w=0\), and at least \(2d/3+w/3-s(w)+s(8-w)-1\) otherwise. Only \(w=7,d=0\) can be negative. This finishes the second sector of (1).

The two sectors neither cover every two-adic realization of the coefficients nor prove the full nonsymmetric four-window.

## A uniform extension to deeper dyadic intervals

Let \(d=2^e\ge4\), and suppose \(r\in(1/d,2/d)\). The same boundary calculation has a unique least term at \(i=j\bmod d\) in the polynomial with normalized argument valuation \(r\). In the proof for \(d=4\), replace the four-digit blocks by blocks of size \(d\). If \(i=d v+w\) and \(w\le u=j\bmod d\), the margin over \(i=u\) is at least
\[
                  (dr-1)v+(w-u)r+s(u-w),
\]
which is strictly positive away from \(i=u\): each bit of \(u-w\) is at most \(d/2\). For \(w>u\), the margin is at least \((dr-1)v+(w-u)r+s(d+u-w)-1>0\).

If the smaller normalized axis valuation is \(r\) and the other axis valuation is strictly larger, the same odd-boundary dominance and comparison-plane argument work whenever the four-window contains an index divisible by \(d\). At such an index the normalized axis value with valuation greater than \(1/d\) has uniquely least constant term, by the same block inequality. This is the only replacement needed in the unequal-axis proof.

If both normalized axis valuations belong to \((1/d,2/d)\), equality of the boundary ratios forces them to coincide and the Laurent argument applies. The leading ratio exponent is \(+1\) unless \(n+2\equiv-1\bmod d\), when it is \(1-d\). Thus a collision requires \(d\mid n+3\), and the scaled residue of \(c/a\) must be a nontrivial \((d-1)\)-st root of unity. If \(2d\mid N=n+3\), every nonconstant term of \(g_N\) has valuation strictly greater than one. Indeed at \(i=2d v+w\) use the same last-block bound with size \(2d\): the lower bound is \(v(2dr-1)\) when \(w=0\), and \(v(2dr-1)+wr+s(2d-w)-1\) otherwise. In the second expression \(wr+s(2d-w)>2\): this is immediate for \(w>d\), is \(dr+1>2\) at \(w=d\), and for \(0<w<d\) one has \(s(2d-w)\ge2\). The boundary then gives \(Q=r\). At \(J=N-1\equiv-1\bmod 2d\), the pure term \(b^{2d-1}\) is uniquely least. Indeed the digit loss of order \(J\) at total degree \(k\) is at least \(-s(k)\), with equality for \(k=2d-1\). For \(k<2d-1\), missing bits each have size at most \(d\), so the margin over the chosen term is strictly positive by \(dr<2\). For \(k\ge2d\), write \(k=2d v+w\) and use \(2dr>2\) to obtain at least \(v(2dr-2)>0\) in addition to the nonnegative missing-bit contribution. At the chosen degree any mixed term has \(s(\ell)<s(2d-1)\) in (2). Thus this entire half of every deeper dyadic interval is excluded. The case \(v_2(n+3)=e\) for \(d>4\) is not covered by this extension.

## Exact checks

The [finite checker](verify_binary_even_quadratic_intermediate_digit_sectors.py) verifies 3,375 factorial ratio integrality instances, 1,101 axis leading terms, 45 deeper dyadic least terms, 1,538 open-annulus margins, 104 high-axis order bounds, 6,429 endpoint neighbor comparisons, and 165 pure-cross switching cases. The inequalities in the proof apply uniformly beyond these fixtures.
