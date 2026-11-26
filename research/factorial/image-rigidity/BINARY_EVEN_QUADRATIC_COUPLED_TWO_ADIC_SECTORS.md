# Coupled two-adic exclusions for even binary quadratics

Let \(v=v_2\) be a valuation extending the ordinary two-adic valuation, normalized by \(v(2)=1\). Normalize the nonzero constant coefficient to one, and put
\[
 f=ax^2+bxy+cy^2+1,\qquad M_j=\Lambda_2(f^j),\qquad
 B_j(t)=\Lambda_1((1+t x^2)^j).
\]
Here \(\Lambda_2(x^i y^j)=i!j!\). Swapping \(a,c\) preserves every moment. The valuations below are those of the normalized coefficients.

**Theorem.** No four consecutive moments starting at \(n\ge1\) can vanish in either sector
\[
\begin{array}{ll}
 v(a)\ge-1,\quad v(c)\ge-1, &\text{or}\\[2mm]
 -10/3\le v(a)\le-2,\quad v(c)\ge-1, &\text{or its exchange}.
\end{array}                                                     \tag{1}
\]
There is no valuation restriction on \(b\) in either sector. In a deeper mixed sector with \(v(a)<-10/3,\ v(c)>-1\), a four-zero block forces
\[
          r=v(4a)=1-\tau,\qquad \tau=v(n+3)\ge3.
                                                               \tag{2}
\]
If \(\tau\ge4\), it further requires
\[
                     (n+3)/2^\tau-1\ \ge\ 2^{(\tau-3)2^\tau}.
                                                               \tag{3}
\]
These are assertions for the original moments, with the two forced boundary cancellations used in conjunction with interior vanishing.

## The coupled boundary equations

The polynomial Stein and Euler identities of the [single-face calculation](BINARY_EVEN_QUADRATIC_SINGLE_FACE_STRONG_FACTORIAL.md#4-the-exact-reduction-for-general-b-and-its-limitation) give, from a four-zero block,
\[
 (c-b/2)B_j(c)+(a-b/2)B_j(a)=0,\qquad j=n+2,n+3.                \tag{4}
\]
They involve no division by the quadratic discriminant or by \(a+c-b\). For reference, the weighted derivative
\[
 (c-b/2)\partial_x f+(a-b/2)\partial_y f=2(ac-b^2/4)(x+y)
\]
is the source of (4). At those indices the full moment and the Euler expression for the preceding \((x+y)\)-moment vanish. It is essential here that (4) accompanies the other interior zeros.

Write \(s_2(u)\) for the binary digit sum of a nonnegative integer \(u\). Legendre's formula gives
\[
 v\!\left(\binom ji(2i)!\right)=2i+s_2(j-i)-s_2(j),\qquad
 v((2i)!/2^i)=v(i!).                                          \tag{5}
\]

## The integral-axis sector

If \(v(t)\ge-1\), set \(T=2t\). Equation (5) gives \(B_j(t)\equiv1+jT\pmod{2}\). For an even \(j\) it also gives
\[
       B_j(t)\equiv1+jT+6\binom j2 T^2\pmod{4}.                 \tag{6}
\]
In fact the cubic contribution is divisible by four, since \(\binom j3\) is even, and every contribution with \(i\ge4\) has valuation at least three. Congruences here mean divisibility in the valuation ring and also apply to ramified fields.

Assume \(v(a),v(c)\ge-1\). Put \(A=2a,C=2c\). Equation (4) reads
\[
                   b(B_j(a)+B_j(c))=A B_j(a)+C B_j(c).
\]
The right side is integral. If \(v(b)<-1\), both sums on the left have valuation greater than one. The odd index therefore forces \(A+C\equiv0\pmod2\). At the even index, (6) then gives \(B_j(a)+B_j(c)\equiv2\pmod4\), using \(A^2+C^2=(A+C)^2-2AC\). This contradiction proves \(v(b)\ge-1\).

With all three valuations at least \(-1\), every \(M_j\) with \(4\mid j\) is a unit. To see this use the multinomial term with respective exponents \(i,\ell,s,u\) of \(a,b,c,1\), summing to \(j\). The first factorial floors for \((2i+\ell)!(2s+\ell)!\) offset the factor valuations up to a loss of at most one, possible only for odd \(\ell\). If any count is odd, the multinomial contributes at least \(v(j)\ge2\). If all counts are even but some count is \(2\bmod4\), it contributes at least \(v(j)-1\ge1\) and there is no factorial loss. If all counts are multiples of four, any nonconstant term has a factorial argument at least four, so a higher factorial floor gives strictly positive valuation. The constant term is one. Every four-term window contains a multiple of four, which proves the first sector.

## A dominant axis and approximation by its single-face plane

Assume \(v(a)\le-2\), \(v(c)\ge-1\), and write \(r=v(4a)\le0\). At every positive index \(j\), the unique lowest term of \(B_j(a)\), by (5), is its leading term. It follows that
\[
                     v B_j(a)=j r-s_2(j)<0.
\]
Every \(B_j(c)\) is integral, and the even ones are units by (6). Set \(\varepsilon=b/(2a)-1\). Solving (4) gives
\[
 \varepsilon=\frac{c-a}{a}\,\frac{B_j(c)}{B_j(a)+B_j(c)},
 \qquad
 v(\varepsilon)=-jr+s_2(j)+v B_j(c),\quad j=n+2,n+3.             \tag{7}
\]
The denominator is nonzero by strict dominance of \(B_j(a)\). The even equation makes \(\varepsilon\ne0\), hence \(B_j(c)\ne0\) also at the odd boundary index. In particular \(v(\varepsilon)>0\).

If \(n\) is even, the two boundary indices are even followed by odd. In (7) the digit sum grows by one, the even \(B_j(c)\) is a unit, and the odd one is integral. The two values for \(v(\varepsilon)\) cannot coincide. Thus \(n\) is odd. Write \(j_b=n+3\) for its even boundary index and \(\tau=v(j_b)\). Comparing the two equations in (7) now gives
\[
            -r=\tau-1+v B_{j_b-1}(c),\qquad
            v(\varepsilon)=-j_b r+s_2(j_b).                    \tag{8}
\]

Let \(M_j^0=\Lambda_2((a x^2+2axy+c y^2+1)^j)\) be the single-face comparison. The actual perturbation satisfies
\[
           v(M_j-M_j^0)\ \ge\ (j-j_b)r+s_2(j_b)-s_2(j).         \tag{9}
\]
For the proof, a term with \(\ell>0\) cross factors gains \(v(\varepsilon)\) from \((1+\varepsilon)^\ell-1\) after comparison with \(b=2a\). For counts \(i,\ell,s,u\) of the other terms as above, the unperturbed valuation is bounded below by
\[
   jr-s_2(j)+\ell-s_2(\ell)+s(v(c)-v(a))+s_2(u)-ur .
\]
This follows by expanding the multinomial valuation and using
\(s_2(2i+\ell)\le s_2(i)+s_2(\ell)\) and its \(s\) analogue in the two factorials. All additions to \(jr-s_2(j)\) are nonnegative. Apply (8) to obtain (9).

There is a useful exact formula for the plane terms. For counts \(h,s,u\) of \(a,c,1\), put \(d=h+s\), so \(j=d+u\). The term is
\[
         \binom jd (d!)^2 (4a)^h c^s \binom{2s}s .              \tag{10}
\]
Indeed
\[
        \Lambda_2\bigl(y^{2s}(x^2+2xy)^h\bigr)
                  =4^h h!(h+s)!(2s)!/s!.
\]
One may verify this over the reals by \(w=x+y,\ t=y/(x+y)\) in the product exponential integral. The radial integral contributes \((2h+2s+1)!\); the other contribution is \(\int_0^1 t^{2s}(1-t^2)^h\,dt\). Its beta value gives the stated integer identity, hence (10) is a polynomial identity over characteristic zero. In particular the term valuation is
\[
   v\binom jd +2v(d!)+h r+s v(c)+s_2(s)
       \ \ge\ (2+\lambda)d-s_2(d)-[s_2(j)-s_2(j-d)],
       \qquad \lambda=\min(r,-1).                              \tag{11}
\]

## Excluding the mixed sector

Suppose \(-4/3\le r\le0\). For \(j\equiv0\bmod4\) and \(d=4u+w\), \(0\le w<4\), the digit-sum loss in (11) is at most \(s_2(u)\) for \(w=0\), and at most \(s_2(u+1)-s_2(4-w)\) otherwise. This follows by writing \(j=4V\) and subtracting \(d\), then using subadditivity of the digit sum on the high part. At \(\lambda=-4/3\), inserting these bounds in (11) gives the following lower bounds:
\[
 \begin{array}{c|cccc}
 w&0&1&2&3\\ \hline
 \text{lower bound}&2u/3&(2u+2)/3&(2u+1)/3&2u/3 .
 \end{array}
\]
Here \(s_2(u)\le u,\ s_2(u+1)\le u+1\) were used. Thus \(M_j^0\) is a unit if \(r>-4/3\). At \(r=-4/3\) the only possible nonconstant term of valuation zero is at \(d=3\); it has strictly positive valuation if \(8\mid j\), since \(s_2(j)-s_2(j-3)=-v(j/4)\).

If \(n\equiv3\bmod4\), take \(j=n+1\). Since \(j_b=j+2\), we have \(s_2(j_b)=s_2(j)+1\). Equation (9) makes the perturbation strictly positive in valuation. Therefore a unit \(M_j^0\) already contradicts the block. The only remaining subcase in this direction is \(r=-4/3\) and \(j\equiv4\bmod8\), so \(j_b\equiv6\bmod8\). At \(j_b\) the pure \(h=d=2\) term of (10) has valuation \(2+2r=-2/3\), since \(\binom{j_b}2\) is odd. It is unique at that lowest value. Indeed \(d=0,1\) have larger valuations, and the other terms with \(d=2\) increase because \(v(c)>r\). For \(d\ge3\), write \(d=8u+w\). Subtracting \(2+2r\) in (11) gives a lower bound of \((10/3)u\) plus the following constants:
\[
       \begin{array}{c|cccccccc}
          w&0&1&2&3&4&5&6&7\\ \hline
           &2/3&1/3&0&2/3&4/3&1&2/3&7/3 .
       \end{array}
\]
For \(w=2\), \(d\ge3\) forces \(u\ge1\), so the bound is still positive. To obtain the table, the digit loss is at most \(s_2(u)+2-s_2(6-w)\) for \(w\le6\), and \(s_2(u+1)-1\) for \(w=7\). Since (9) at \(j_b\) gives a nonnegative perturbation valuation, \(M_{j_b}\ne0\).

If \(n\equiv1\bmod4\), then \(\tau\ge2\) and (8) with \(-r\le4/3\) forces \(\tau=2\) and \(-4/3\le r\le-1\). Hence \(n\equiv1\bmod8\). Use \(j=n\), for which the pure \(h=d=1\) term of (10) has valuation \(r<0\). The constant and \(s=d=1\) terms have larger valuation. For \(d\ge2\), subtract \(r\) in (11); the result is at least the bound obtained at \(r=-4/3\), since \(\lambda=r\) and the coefficient of \(r\) in the difference is \(d-1\). Write \(d=8u+w\). The lower bound is \((10/3)u\) plus the following constants:
\[
       \begin{array}{c|cccccccc}
          w&0&1&2&3&4&5&6&7\\ \hline
           &4/3&0&8/3&4/3&3&5/3&3&2 .
       \end{array}
\]
Here \(w=1\) forces \(u\ge1\). To check the table, the digit loss is at most \(s_2(u)+w\) for \(w\le1\), and \(s_2(u+1)+1-s_2(9-w)\) otherwise. Thus the term of valuation \(r\) is unique. The indices \(n=1\bmod8\) and \(j_b=n+3=4\bmod8\) have equal digit sum. Equation (9) therefore bounds the perturbation at \(n\) by \(-3r>0\). This contradicts \(M_n=0\) and proves (1).

## An additional restriction in the deep strict mixed sector

If \(v(c)>-1\), all \(B_j(c)\) are units, even at odd indices, by the reduction modulo two above. In the remaining deep mixed sector equation (8) therefore gives \(r=1-\tau\). Since \(r<-4/3\), this forces \(\tau\ge3\). Take \(j=j_b-1=2^\tau w-1\), \(w\) odd, and \(L=2^\tau\). Divide \(M_j^0\) by \((j!)^2(4a)^j\) and group (10) by omitted degree \(e=j-d\). The inner sum in the grouped factor is
\(\sum_{s=0}^{j-e}\binom{2s}s(c/(4a))^s\), a unit since \(v(c/(4a))>0\). Hence the exact grouped valuation is
\[
             V(e)=(\tau-3)e+s_2(e)+s_2(j)-s_2(j-e).
\]
In a block \(e=Lh+b\) with \(0<b<L\), the low binary digits of \(j-Lh\) are all ones. It follows that
\(V(Lh+b)-V(Lh)=(\tau-3)b+2s_2(b)>0\). At a block start,
\[
 V(Lh)=(\tau-3)Lh+s_2(h)+s_2(w-1)-s_2(w-1-h).
\]
If \(\tau\ge4\) and \(w-1<2^{(\tau-3)L}\), then every term with \(h>0\) has positive valuation because the digit sum of \(w-1-h\) is at most \((\tau-3)L\). Thus the grouped constant term is uniquely lowest. The plane moment has valuation
\((3-\tau)j-2s_2(j)<0\). Equation (9) at \(j=j_b-1\) gives a nonnegative perturbation valuation, a contradiction. This proves (2) and (3).

The inequalities leave other valuation sectors and the unrestricted nonsymmetric four-term theorem open.

## Exact checks

The [finite checker](verify_binary_even_quadratic_coupled_two_adic_sectors.py) compares 484 plane weights with direct factorial sums, 70 coupled recurrences with actual moments, 75 perturbation bounds with actual moments, and 164 boundary congruences. It also checks 6,724 carry comparisons in the three index classes and 4,464 deep digit-block equalities. These are finite checks of the identities and inequalities used in the proof above.
