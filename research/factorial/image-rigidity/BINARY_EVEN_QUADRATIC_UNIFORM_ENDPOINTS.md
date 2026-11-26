# Uniform endpoint exclusion for even binary quadratics

Normalize the constant to one in
\[
 f=ax^2+bxy+cy^2+1,\qquad M_j=\Lambda_2(f^j),
 \qquad g_j(y)=\Lambda_1((1+yx_1^2)^j).
\]
Let \(v=v_2\) with \(v(2)=1\). Suppose \(a\ne c\), \(d=2^e\ge4\), and \(v(4a)=v(4c)=1/d\). Extend the valued field to choose \(t^d=2\), and let \(U,V\) be the residues of \(4a/t,4c/t\).

**Theorem.** Suppose \(U,V\ne1\). If \(M_{N-3},\ldots,M_N\) vanish at positive indices, then \(v_2(N)=e+1\) and \(S=U+V\) is a nontrivial cube root of unity. The two axis residues are distinct. Also, whether or not both axis residues differ from one, a collision with \(d\mid N\) cannot have exactly one of \(U,V\) equal to one.

The condition \(U=1\), as well as the other residue conditions, is independent of the choice of \(t\). Any two such choices differ by a root of unity of \(2\)-power order, which reduces to one in characteristic two.

## Axis polynomials at the endpoint

The coefficient of \(Z^i\) in \(g_j(Z/4)\) has valuation \(s(j-i)-s(j)\), with \(s\) the binary digit sum. Write \(j=dT+u\) with \(0\le u<d\). In \(g_j(tX/4)\), the uniquely least index for even \(T\) is \(i=u\). For odd \(T\) there is one additional least index \(u+d\). All other indices lie strictly higher.

Indeed, for \(i=dq+w\) and \(w\le u\), the valuation difference from the \(u\)-th term at unit \(X\) is
\[
       q+s(T-q)-s(T)+(w-u)/d+s(u-w).
\]
The upper-block contribution is at least \(q-s(q)\), and the low-block contribution is strictly positive if \(w<u\) since \(s(u-w)\ge2(u-w)/d\). The upper-block term can vanish only for \(q=0,1\). For \(q=1\) the contribution is \(v_2(T)\). For \(w>u\), the lower bound is
\[
                      q-s(q+1)+(w-u)/d+s(d+u-w)>0 .
\]
This proves the asserted minima.

Divide \(g_j(tX/4)\) by \(t^u2^{-s(u)}\). Its reduction at unit \(X\) is then
\[
                  X^u \quad\text{or}\quad X^u(1+X^d),               \tag{1}
\]
according as \(T\) is even or odd. The normalized polynomial differs from the corresponding characteristic-zero polynomial in (1) by terms with coefficients of strictly positive valuation.

## The boundary quotient

The [weighted boundary equations](BINARY_EVEN_QUADRATIC_SINGLE_FACE_STRONG_FACTORIAL.md#4-the-exact-reduction-for-general-b-and-its-limitation) at \(N-1,N\) imply
\[
                 g_N(a)g_{N-1}(c)=g_N(c)g_{N-1}(a)                 \tag{2}
\]
since the two constant weights are not both zero when \(a\ne c\). When \(U,V\ne1\), both normalized axis values in (1) are units at both arguments. If \(d\nmid N\), the reduced quotient of the normalized polynomials is \(X\). If \(d\mid N\), it is
\[
 L_+(X)=X^{1-d}(1+X^d)\quad\text{for }N/d\text{ odd},\qquad
 L_-(X)=\frac{X^{1-d}}{1+X^d}\quad\text{for }N/d\text{ even}.        \tag{3}
\]

Each of the three reduced rational maps has nonzero derivative on a nonzero residue other than one. The derivatives in (3) are \(1+X^{-d}\) and \(X^{-d}/(1+X^d)\) in characteristic two. The exact normalized quotient is injective within each such residue class. To see this, the difference of two evaluations of any error polynomial above has valuation strictly greater than the argument difference, by factoring the difference of each monomial. The errors have positive valuation, while the displayed leading denominators are units on the class. Subtracting the corresponding characteristic-zero rational map from the exact quotient preserves this stronger difference bound by the product and inverse rules. The leading rational map itself has divided difference a unit at two points with this same residue, since it reduces to its derivative there. Its leading difference therefore cannot vanish after adding the error.

Thus \(U,V\ne1\) and (2) require distinct residues and \(d\mid N\). When \(d\mid N\), one of the two normalized boundary polynomials has a factor \(1+X^d\) in its reduction and the other has none. Consequently, if exactly one axis has scaled residue one, the two sides of the normalized determinant (2) have different zero status in the residue field. This excludes the mixed exceptional residues in the wrapped case.

## The case \(N/d\) odd

Let \(R=L_+(U)=L_+(V)\), a nonzero residue. The weighted boundary equation at \(N\) gives
\[
                \frac{4b}{2t}
                =\frac{(4a/t)g_N(a)+(4c/t)g_N(c)}{g_N(a)+g_N(c)} .
\]
The denominator reduces to \(U^d+V^d\ne0\). Because \(U(1+U^d)=R U^d\) and \(V(1+V^d)=R V^d\), the quotient reduces to \(R\), giving \(v(4b)=1+1/d\).

Use the original moment at \(j=N-2\). Write \(i,\ell,m\) for the counts of \(a,b,c\), with \(k=i+\ell+m\). By the [exact multinomial digit formula](BINARY_EVEN_QUADRATIC_INTERMEDIATE_DIGIT_SECTORS.md#coefficient-valuations-and-a-useful-integer-ratio), the term valuation is
\[
 k/d+s(j-k)-s(j)+\ell-s(\ell)
        +\operatorname{car}(i,\lfloor\ell/2\rfloor)
        +\operatorname{car}(m,\lfloor\ell/2\rfloor).
\]
At \(j=N-2\) the high block \(\lfloor j/d\rfloor\) is even. The first three terms uniquely minimize at \(k=d-2\) by the axis calculation. The extra cross count penalty is nonnegative, with zero only at \(\ell=0,1\), and both carry counts vanish for these two values. With a common factor removed, the least-term residue is
\[
 S_{d-2}+R S_{d-3},\qquad S_m=\frac{U^{m+1}+V^{m+1}}{U+V}.
\]
Each unit integer coefficient in the normalized group reduces to one, and the cross term contributes the additional factor \(4b/(2t)\). Using \(R U^{d-2}=U^{-1}+U^{d-1}\) and the same equality at \(V\), the residue is \(1/(UV)\ne0\). Hence \(M_{N-2}\) is nonzero.

## The case \(N/d\) even

Put \(D=2d\) and \(T=N/D\). In \(g_N(tX/4)\) the nonconstant terms of valuation at most one are the two indices \(d,D\) when \(T\) is odd. Both have valuation one. There are none when \(T\) is even. To check this, write \(i=Dq+w\), \(0\le w<D\). If \(w=0\), the lower bound \(2q-s(q)\) is greater than one for \(q\ge2\); the exact value at \(q=1\) is \(1+v_2(T)\). If \(w>0\), then \(w/d+s(D-w)\ge2\), with equality only at \(w=d\). Combining this with \(2q-s(q+1)\ge q-1\) excludes \(q\ge1\). For \(q=0\) the actual upper-block contribution is \(v_2(T)-1\), giving equality only at \(w=d\) and \(T\) odd.

It follows that \((g_N(a)+g_N(c))/2\) reduces to one for even \(T\), or to \(1+S^d+S^{2d}\) for odd \(T\). The weighted boundary numerator after normalization is \((4a/t)g_N(a)+(4c/t)g_N(c)\), reducing to \(S\ne0\). If the denominator after division by two is nonzero in the residue field, then \(v(4b)=1/d\), with the residue of \(4b/t\) equal to \(S\) for even \(T\). In that case \(S\ne1\). In fact \(L_-(X)=X/[X(1+X)]^d\), so if \(U+V=1\) the denominators there agree at \(U,V\), and the two quotient values would differ.

It remains to show that \(v(4b)=1/d\) with these residue conditions gives a nonzero \(M_{N-1}\). At each total degree \(k\) a lower bound for all term valuations is the valuation of the pure cross term, namely
\[
                       k/d-s(k)+s(N-1-k)-s(N-1).
\]
This follows from the [integer coefficient ratio](BINARY_EVEN_QUADRATIC_INTERMEDIATE_DIGIT_SECTORS.md#coefficient-valuations-and-a-useful-integer-ratio) at fixed total degree. Write \(N-1=DH+D-1\) with \(H=T-1\), and \(k=Dq+w\) with \(0\le w<D\). The pure term valuation becomes
\[
          2q-s(q)+s(H-q)-s(H)+w/d-2s(w).
\]
The low-block contribution is uniquely minimized at \(w=D-1\) since every low bit has size at most \(d\). The high-block contribution is at least \(2q-2s(q)\), strictly positive for \(q\ge2\). At \(q=1\) it equals \(v_2(H)\) when \(H\ge1\). Therefore the minimizing degree is \(D-1\), with one more equal degree \(2D-1\) if \(H\) is odd. Each of these degrees has all bits set, so every other monomial at that degree is strictly higher by the exact digit valuation formula, using \(s(\ell)<s(k)\) for \(\ell<k\). If \(H\) is even, the sole least term proves nonvanishing. If \(H\) is odd, then \(T\) is even. The two pure terms have normalized residue ratio \(S^D\ne1\), since their exponents differ by \(D\) and \(S\ne1\). Their sum cannot cancel.

Thus a collision with \(U,V\ne1\) requires \(T\) odd and \(1+S^d+S^{2d}=0\). The last expression is \((1+S+S^2)^d\) in characteristic two, which gives exactly the asserted cube-root condition.

## Exact checks

The [finite replay](verify_binary_even_quadratic_uniform_endpoints.py) checks 521 endpoint axis minima, 45 even-wrap axis checks, 45 pure-cross minima, 216 odd-wrap and 150 even-wrap residue pairs, and 595 denominator identities. The proof establishes the bounds uniformly in \(d\) and \(N\).
