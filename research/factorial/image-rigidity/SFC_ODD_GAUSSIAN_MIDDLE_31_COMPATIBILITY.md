# The normalized valuation-one obstruction at the middle prime 31

The entire lower-row unit contributes at valuation one. Its contribution
changes the leading-row compatibility from 17 to 20 modulo 31; it does not
cancel it. The proof below accounts for all its children by cancelling the
lower digit factors before taking the last digit quotient.

**Theorem (Middle31).** For every integer \(e\ge1\), the three quadratic
Gaussian forms
\[
 G_m(A,B,C)=\mathbb E(A+BX_m+CX_m^2)^m,
 \qquad X_m\sim\mathcal N(0,1/m),
 \qquad m=6\cdot31^e-1,\ 6\cdot31^e,\ 6\cdot31^e+1
\]
have no common nonzero projective zero over any characteristic-zero field.
For \(e\ge3\), their fully normalized final-slope compatibility through
relative valuation one is \(20\ne0\) modulo 31 under the projection
\((c_+,c_-)\mapsto3c_+-7c_-\).

The cases \(e=1,2\) belong to the already established finite prefix in
[the local-window theorem, Sections 4 and 6](SFC_ODD_GAUSSIAN_LOCAL_WINDOW_TESTS.md#6-exact-tests-at-prime-divisors-of-the-prescribed-orders).
The rest of this proof is uniform in \(e\ge3\). It does not infer an
all-exponent statement from the finite controls at exponents two and three.

## 1. Reduction and precision

Put \(p=31\), \(M=p^e\), \(J=M/p\), \(a=(M-1)/2\),
\(\ell=(M-1)/30\), and \(L=\ell/M\). Use \(v(p)=1\).
The middle-order reduction in
[Section 5 of the local-window proof](SFC_ODD_GAUSSIAN_LOCAL_WINDOW_TESTS.md#5-odd-starts-with-a-large-prime-power-divisor-in-the-middle-order)
forces the chart \(B=1\), \(x=A/B\), \(h=C/B\), with
\(v(h)=e/2-L/2\). Its cubic comparison excludes \(v(x)\ge v(h)\).
The remaining Newton slope is \(v(x)=L/2\). Set
\[
 Y=x^2,\qquad Z=(-1)^\ell 2Y^M/p^\ell.
\]
The two pure odd rows, after removing their common linear factors, have
coefficients
\[
 U_r^+=\frac{(2(6M+1))^r(3M)_{\underline r}}{(2r+1)!},
 \quad 0\le r\le3M,
 \qquad
 U_r^-=\frac{(2(6M-1))^r(3M-1)_{\underline r}}{(2r+1)!},
 \quad 0\le r<3M.
\]
Normalize the upper row at \(r=0\) and the lower row by
\(U_a^-Y^a\). All weights below are relative to these minima.

The pure-row initial polynomials are
\[
 A_0(Z)=1+3Z/2+Z^2/4+Z^3/120,
 \qquad B_0(Z)=1-Z/3+Z^2/60.
\]
Their reductions have exactly one common root, \(z=7\), with
\(A_0'(7)=7\), \(B_0'(7)=3\) modulo 31. This follows either by
division or from the remainder \(160(7Z-18)\) of the integer-scaled
polynomials in the local-window proof.

Terms containing \(h\) have relative weight at least
\(v(h)-v(x)=e/2-L>1\), by the coefficientwise integral-ratio argument
in that proof. Thus only the pure arrays matter through weight one.
Write \(O_{>1}\) for an error of strictly greater relative valuation.
A coefficient-unit congruence modulo \(p\) at a positive weight has
error greater than one. At weight zero, coefficients must instead be
known modulo \(p^2\). We retain those errors explicitly in Section 4.

All associated-graded identities below mean identities of the actual
Laurent monomials in \(p,Y\), with each coefficient reduced after
stripping its power of \(p\). Multiplication and division by a unit are
performed with the weight cutoff in place; only finitely many terms
occur. No convergence of a \(p\)-adic exponential is asserted.

## 2. Complete upper inventory and the two lower shapes

For \(r=tM+b\), \(0\le b<M\), Legendre's formula gives the upper weight
\[
 w_+(b)=\frac{s_p(b)-b/M}{30}
   +\sum_{d=1}^e\mathbf1_{\{0<b\bmod p^d<(p^d-1)/2\}}.       \tag{1}
\]
If the sum of indicators is positive, the weight exceeds one. Otherwise
the digits after the first nonzero digit are all at least 15. If that
first digit is above position zero, it is at least 16. For \(e\ge3\),
any such string with at least two nonzero digits has digit sum at least
31, and hence weight greater than one. A string starting at position
zero has digit sum at least \(15e>30\). Consequently the complete
upper inventory through one is
\[
 b=0,\qquad b=vJ\quad(16\le v\le30),
 \qquad w_+(vJ)=v/p.                                      \tag{2}
\]
There is no upper group of exact weight one.

Grouping indices \(b+tM\), the factorial top digits give the shapes
\[
 H(Z)=1+Z/3+Z^2/60=\tfrac23 A_0'(Z),
 \qquad C(Z)=1-Z+Z^2/12=B_0(Z)+2ZB_0'(Z).                 \tag{3}
\]
The upper groups in (2) have shape \(H\). Every lower group with
\(b\ge a\) has shape \(B_0\), and every lower group with \(b<a\)
has shape \(C\). Indeed \(3M-1-(b+tM)\) has top digit \(2-t\);
\(2(b+tM)+1\) has top digit \(2t+1\) or \(2t\), respectively.
Wilson's digit-factorial formula supplies the signs and the scale in
\(Z\). The top factorial arguments are below \(p\).

Let \(D\) be the whole lower unit formed by the groups \(b\ge a\),
and \(K\) the scalar sum of groups \(b<a\), each normalized by
\(U_a^-Y^a\). The group at \(b=a\) has scalar one. Let \(E\)
denote the scalar sum of the fifteen upper secondary groups. With
\(A_e,B_e\) the exact leading-group polynomials modulo \(p^2\),
the normalized rows are therefore
\[
 F_+=A_e+EH+O_{>1},\qquad
 F_-=DB_e+KC+O_{>1}.                                    \tag{4}
\]
Using \(B_e\) instead of \(B_0\) in the positive-weight groups changes
only \(O_{>1}\). This also explains why no positive-weight coefficient
error is missing in (4).

## 3. Cancel every lower digit before dividing the last digit

The exact lower coefficient recurrence, with \(r=a+c\), is
\[
 \frac{U^-_{a+c+1}}{U^-_{a+c}}
 =\frac{(6M-1)(5M-1-2c)}{(M+2c+1)(M+2c+2)}.             \tag{5}
\]
For \(0\le c<a\), stripping powers of \(p\) gives the same valuation
and residue as \(1/(2c+2)\). For \(-a\le c<-1\) the same comparison
holds, while the inverse ratio at \(c=-1\) is
\(M(M-1)/((6M-1)(5M+1))\), with valuation \(e\) and unit residue one.
The compared nonzero integers have absolute value below \(M\), so
their valuations are below \(e\); this justifies every comparison.
Thus, coefficientwise in the associated graded ring,
\[
 D=\sum_{c=0}^a\frac{Y^c}{2^c c!},\qquad
 K=\sum_{d=1}^a M(-2)^{d-1}(d-1)!Y^{-d},                 \tag{6}
\]
with only terms of the required weight retained.

Put \(\ell_i=(p^i-1)/30\) and
\[
 X_i=\frac{(-1)^{\ell_i}Y^{p^i}}{2^{p^i}p^{\ell_i}},
 \qquad w(X_i)=\frac{1-p^{i-e}}{30},\qquad X=X_{e-1}.
\]
In particular \(w(X)=1/p\). Wilson's formula
\(\operatorname{unit}(c!)=(-1)^{v(c!)}\prod c_i!\bmod p\)
turns each coefficient of the first sum in (6) into
\(\prod X_i^{c_i}/c_i!\).

Here is a complete cutoff argument for the negative powers. Write
\(d=sJ-c\), \(0\le c<J\). Since \(d\le a\), one has \(1\le s\le16\).
The weight of the corresponding term of \(K\) is exactly
\[
 1-s/p+g(c),\qquad g(c)=\frac{s_p(c)-c/M}{30}.             \tag{7}
\]
For \(s=16\), the degree restriction requires \(c\ge(J+1)/2\).
Then \(s_p(c)\ge16\) and
\(g(c)>(16-1/p)/30>16/p\); the weight exceeds one.
For \(1\le s\le15\), (7) is at least \(16/p\), and a retained
term requires \(g(c)\le s/p\le15/p\). In particular
\(s_p(c)\le15\). The complement digits of
\(d-1=(s-1)J+(J-1-c)\), together with
\((30-j)!\equiv(-1)^{j+1}/j!\), show that its scalar is the
scalar at \(d=sJ\) times the scalar \(Y^c/(2^c c!)\).
The valuation identity is
\[
 e+v((sJ-c-1)!)=1+s\ell_{e-1}-v(c!).
\]
This verifies both the units and the powers of \(p\) in the factorization.

Define the finite lower-digit unit
\[
 D_{\rm lo}=\prod_{i=0}^{e-2}\left(\sum_{j=0}^{15}X_i^j/j!\right),
 \quad D_{\rm hi}(X)=\sum_{j=0}^{15}X^j/j!,
 \quad K_{\rm hi}=\frac p2\sum_{s=1}^{15}(-1)^{s-1}(s-1)!X^{-s}.
\]
Only weights through the relevant cutoff are used in these products.
The preceding calculation proves \(K=D_{\rm lo}K_{\rm hi}+O_{>1}\).
Since \(w(K)\ge16/p\), only weights at most \(15/p\) of \(D\)
can affect \(K/D\) through one. At that cutoff every index in \(D\)
has digit sum at most 15 and automatically is at most \(a\). Hence
\(D=D_{\rm lo}D_{\rm hi}\) to precisely the needed precision, and
\[
                  K/D=K_{\rm hi}/D_{\rm hi}+O_{>1}.       \tag{8}
\]
This cancels all lower digit children uniformly in their number.

Let \(\mathcal H=X^p/p\). From \(\ell=1+p\ell_{e-1}\) one obtains
\(\overline{\mathcal H}=-\bar Z/4\). The upper factorial digits give
\[
 E=\mathcal E(X),\qquad
 \mathcal E(X)=3\sum_{v=16}^{30}
       \frac{4^v X^v}{(31-v)!(2v-31)!}\quad\text{in }\mathbb F_{31}[X].
\]
The complete remaining calculation is the fixed polynomial identity
\[
 \frac{\sum_{s=1}^{15}(-1)^{s-1}(s-1)!X^{31-s}}
      {\sum_{j=0}^{15}X^j/j!}
       =-\mathcal E(X)/6+12X^{31}\pmod{(31,X^{32})}.       \tag{9}
\]
All denominators are units. Multiplication by the denominator verifies
every coefficient of this finite identity; the accompanying checker
performs exactly that check, including zero coefficients.
Equations (8)--(9) now give
\[
                         K/D=E/(3Z)+6p+O_{>1}.           \tag{10}
\]
Replacing \(\mathcal H\) by \(-Z/4\) in a term of weight at least
\(16/p\) incurs error greater than one. For the last term the
cancellation \(X^{31}/\mathcal H=p\) is exact. The term \(6p\)
is the boundary contribution that would be lost by stopping below one.

## 4. Every weight-zero coefficient error modulo the prime square

Write \(W=(p-1)!\bmod p^2\), \(h_{15}=\sum_{j=1}^{15}j^{-1}\bmod p\),
\(\kappa=-2^{30}/W\), and \(\eta=\kappa(1+ph_{15})\bmod p^2\).
Direct arithmetic gives \(h_{15}=19\), \(\kappa=776\), \(\eta=404\).
For every \(e\ge2\), the leading rows satisfy
\[
 A_e(Z)=A_0(\kappa Z),\qquad B_e(Z)=B_0(\eta Z)\pmod{p^2},
\]
or, in increasing coefficient order,
\[
 A_e=[1,203,628,612],\qquad B_e=[1,506,542].              \tag{11}
\]

Here is the uniform factorial-unit calculation behind (11). For
\(F(n)=n!/p^{v(n!)}\), and \(0\le t<p\), grouping nonmultiples
of \(p\) into blocks gives
\[
 F(tp^e)=(-1)^{t(e-1)}W^t t!\pmod{p^2}.
\]
Indeed a full block of units modulo \(p^2\) has product \(-1\),
and a block of length \(p\) has product \(W\), since
\(\sum_{j=1}^{p-1}j^{-1}=0\bmod p\). For half-digit factorials,
\[
 \frac{F(a+2M)}{F(a+(2-j)M)}
 =(-1)^{j(e-1)}W^j(1+jph_{15})\,2_{\underline j}
          \pmod{p^2},\qquad j=0,1,2.
\]
At levels at least two the two remainders modulo \(p^2\) agree;
their full-block counts differ by \(j\). At level one, the tail is
15 and its ratio is \(1+jph_{15}\); the top factorial ratio is
\(2_{\underline j}\). These are all levels of the recursion.
Substitution in the displayed \(U^\pm\) formulas proves (11):
\(\ell\equiv e\bmod2\), \(6M\equiv0\bmod p^2\), and
\(2^{p^e}\equiv2^p\bmod p^2\) remove the exponent dependence.
No finite-exponent extrapolation is used.

In particular,
\[
 A_e(7)/p=29,\quad B_e(7)/p=10,\quad
 (B_e'A_e-A_e'B_e)(7)/p=17\pmod p.                     \tag{12}
\]

## 5. Compatibility and the nonzero transverse residue

A hypothetical common zero in (4) has \(v(Z-7)\ge16/p\): the
upper leading polynomial has a simple root modulo \(p\), and all its
other terms have weight at least \(16/p\). Thus \(B_0(Z)\) also
has valuation at least \(16/p\). Since \(w(E)\ge16/p\), products
\(E^2\) and \(EB_0(Z)\) have valuation greater than one.
Divide the lower row by \(D\), use (3) and (10), and set
\(Z_*=Z+2E/3\). The exact value equations through one become
\[
 0=A_e(Z_*)+O_{>1},\qquad
 0=B_e(Z_*)+6pC(7)+O_{>1}.                             \tag{13}
\]
The derivatives of \(A_e,B_e\) agree with those of \(A_0,B_0\)
modulo \(p\); their errors times \(E\) lie above one. This justifies
the common shift and its entire Taylor remainder.

The first equation forces \(v(Z_*-7)\ge1\). Put
\(Z_*=7+p\zeta\), with \(v(\zeta)\ge0\), and divide (13) by \(p\).
Since \(C(7)=11\bmod31\), the equations are
\[
 29+7\bar\zeta=0,\qquad 10+3\bar\zeta+4=0.
\]
Their transverse compatibility is
\[
 3\cdot29-7(10+4)=20\ne0\pmod{31}.                    \tag{14}
\]
The leading contribution is 17 and the complete row-unit contribution
is \(-7\cdot4=3\). This contradicts a common zero for every \(e\ge3\).
Together with the two finite-prefix cases it proves Middle31.

## Evidence and closest classical inputs

The proof uses elementary Legendre and Wilson block arguments, including
all coefficient errors at the prime square. The closest general arithmetic
input is [Granville's prime-power factorial method, proof of Theorem 1](https://www.cecm.sfu.ca/organics/papers/granville/paper/binomial/html/node3.html).
The shapes in (3) are terminating Kummer polynomials, consistent with
[DLMF 13.2.1--13.2.2](https://dlmf.nist.gov/13.2).
Those sources supply the classical methods, not the normalized compatibility
(9)--(14), which is derived here. No novelty or historical-priority conclusion
is inferred from this comparison.

Run the exact finite audit from the repository root:

```sh
python3 research/factorial/image-rigidity/verify_sfc_odd_gaussian_middle_31.py
```

It verifies (9), (11)--(14), the grouped shapes, and the digit-factor
coefficients against full independent coefficient recurrences at \(e=3\).
It also recomputes the weight-one boundary directly from the earlier
complete normalized inventory. The retained
[receipt](sfc_odd_gaussian_middle_31_verification.json) records the actual run.
The all-exponent conclusion follows from the written cutoff and block
arguments. This is internal exact evidence, without independent mathematical
review, external review or formal verification. The unrestricted odd-start
Gaussian window remains open.
