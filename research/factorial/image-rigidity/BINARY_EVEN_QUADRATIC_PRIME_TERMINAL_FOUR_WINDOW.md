# Four factorial moments of an even binary quadratic ending at a prime

Let \(\Lambda_2(x^i y^j)=i!j!\), and let
\(f=ax^2+bxy+cy^2+k\) with \(abck\ne0\) over a
characteristic-zero field.

**Theorem.** For every prime \(p\ge5\), the four factorial moments
\(\Lambda_2(f^{p-3}),\Lambda_2(f^{p-2}),\Lambda_2(f^{p-1}),
\Lambda_2(f^p)\) do not all vanish.

The proof uses the two exact weighted boundary cancellations and
the original interior moments. At the terminal prime there are only
two possible valuation transitions in the boundary ratios. The
interior moments exclude each transition.

## Reduction at the terminal prime

Write \(M_j=\Lambda_2(f^j)\) and
\[
 B_j(Z)=\sum_{i=0}^j\binom ji(2i)!\,Z^i k^{j-i}.
\]
The [coupled weighted boundary recurrence](BINARY_EVEN_QUADRATIC_SINGLE_FACE_STRONG_FACTORIAL.md#4-the-exact-reduction-for-general-b-and-its-limitation)
makes the two boundary expressions
\[
 (c-b/2)B_j(c)+(a-b/2)B_j(a)=0,\qquad j=p-1,p               \tag{1}
\]
necessary for a four-zero block. We can exclude the already proved
cases \(a=c,\ b=2a,\ b=2c,\ a+c=b,\ b^2=4ac\) by the
[single-face theorem](BINARY_EVEN_QUADRATIC_SINGLE_FACE_STRONG_FACTORIAL.md)
and the [critical boundary theorem](BINARY_EVEN_QUADRATIC_CRITICAL_BOUNDARY_WINDOWS.md).

If there were a characteristic-zero solution, the polynomial
equalities and the required nonvanishing conditions would have a
solution over \(\overline{\mathbb Q}\) by the weak Nullstellensatz.
Apply an embedding into \(\overline{\mathbb Q}_p\), choosing
\(v(p)=1\). Scale the phase so all four coefficients are integral
and at least one is a unit. From Frobenius,
\(M_p\equiv k^p\bmod\mathfrak m\), where \(\mathfrak m\) is the
maximal ideal of the valuation ring. Thus \(v(k)>0\).
In the reduction of \(M_{p-1}\), only the monomial
\(x^{p-1}y^{p-1}\) survives from the homogeneous phase to the
\((p-1)\)-st power. Its coefficient is
\((b^2-4ac)^{(p-1)/2}\) in characteristic \(p\). Indeed, on
factoring the quadratic as \(a(x-r y)(x-s y)\), that coefficient
is \(a^{p-1}\sum_{i=0}^{p-1}r^is^{p-1-i}
=a^{p-1}(r-s)^{p-1}\) by the binomial coefficients modulo \(p\)
and Frobenius. Polynomial identity extends the formula to \(a=0\).
Consequently \(\bar b^2=4\bar a\bar c\). At least one axis
coefficient is a unit, since otherwise \(\bar b=0\) as well.
Exchange the axes if needed and scale once more to make \(a=1\).
Use the notation
\[
 p=2h+1,\quad t=b/2,\quad d=c-t^2,\quad u=t-1,\quad
 D=1+c-2t=d+u^2,\quad s=v(k)>0 .
\]
We have \(v(d)>0\), and (1) becomes
\[
        D B_j(c)-u\bigl(B_j(1)-B_j(c)\bigr)=0.             \tag{2}
\]
All bars below are reductions in the residue field.

## The two boundary transitions

For a unit argument \(Z\), the two polynomials in (2) have the
following strictly lowest terms, except at the indicated ties:
\[
\begin{array}{c|cc}
 &\text{small }s&\text{large }s\\ \hline
 B_{p-1}(Z)
    & \binom{p-1}{h}(p-1)!\,k^h Z^h\quad(s<1/h)
    & (2p-2)! Z^{p-1}\quad(s>1/h)\\
 B_p(Z)
    & k^p\quad(s<2/p)
    & (2p)! Z^p\quad(s>2/p).
\end{array}                                               \tag{3}
\]
At each tie exactly the two displayed terms have the same
valuation. For the first row use \(v\binom{p-1}i=0\),
with factorial valuation zero for \(i\le h\) and one for
\(h<i\le p-1\). For the second use
\(v\binom pi=1\) at \(0<i<p\), while the final factorial
has valuation two. The remaining candidates have strictly
higher valuation even at \(s=2/p\).

If \(\bar t=0\), then \(v(c)>0\), and the weight of \(B_j(1)\)
in (1) is a unit while the other weight is not. For \(s<2/p\)
the constant term uniformly dominates both \(B_p\)'s, giving
a unit weight sum. For \(s>2/p\) the leading term of \(B_p(1)\)
has smaller valuation than every term of \(B_p(c)\).
At \(s=2/p\) use \(B_{p-1}(1)\) instead, whose unique lowest
term is at \(i=h\); every term of \(B_{p-1}(c)\) has strictly
higher valuation. This excludes \(\bar t=0\).

Hence \(c\) is a unit. By (1), the projective boundary ratios
\([B_p(1):B_{p-1}(1)]\) and \([B_p(c):B_{p-1}(c)]\)
are equal. Away from the two transitions, the quotient in
(3), after removing its constant overall multiplier, reduces
to \(Z^q\) for \(q=-h,h+1,1\), respectively. For \(\bar t=1\)
this implies strict injectivity on the residue disk of \(1\),
since \(q\) is a \(p\)-unit. More explicitly expand in
\(Z-1\). The error from the strictly higher coefficients in
(3) has all coefficients of strictly positive valuation,
whereas \(Z^q\) has unit linear coefficient. Differences
on the open residue disk therefore have the valuation of the
difference of the arguments. This disallows \(1\ne c\).

If \(\bar t\ne1\) and \(\bar t\ne0\), the weights in (1)
reduce to \(\bar t(\bar t-1)\) and \(1-\bar t\).
When a lowest term of degree \(i\) is unique, cancellation
requires \(\bar t^{2i+1}=1\).
For \(s<2/p\) use the constant term at \(p\).
For \(2/p\le s<1/h\) use index \(h\) at \(p-1\), giving
\(\bar t^p=1\). For \(s>1/h\) use both top indices,
giving \(\bar t^{2p+1}=\bar t^{2p-1}=1\). All give
\(\bar t=1\), which is excluded in this case. Thus the only
remaining off-diagonal case is
\[
                  s=1/h,\qquad \bar t^{2p+1}=1 .           \tag{4}
\]

## The off-diagonal case

Expand \(f=(x+ty)^2+d y^2+k\) in its last two terms.
The term with total perturbation count \(L\) and \(i\)
occurrences of \(d y^2\) contributes
\[
 \binom jL\binom L i k^{L-i}d^i R!\,
       S_{R,i}(t),\qquad
 R=2(j-L),\quad
 S_{R,i}(t)=\sum_{\ell=0}^{R}(\ell+1)_{2i}t^\ell.           \tag{5}
\]
Here \((\ell+1)_{2i}\) is a rising factorial. Formula (5)
follows directly by expanding \((x+ty)^R\).
For \(0\le i\le h\) one has the polynomial congruence
\[
             S_{p-1,i}(t)\equiv (2i)!(t-1)^{p-1-2i}
                                  \pmod p.               \tag{6}
\]
Indeed the sum is the \(2i\)-th derivative of
\(t^{2i}(t^p-1)/(t-1)\). In characteristic \(p\),
divide \(t^{2i}\) by \(t-1\); only the proper fraction
contributes after differentiation. If \(i>h\), the rising
factorial in every summand of (5) is divisible by \(p\).

At \(p=5\), \(j=p-3=h\) gives a unit term with \(L=i=0\),
by (6) and \(\bar t\ne1\); every term with \(L>0\)
has positive valuation. Thus we may take \(p\ge7\) below.
For \(j=p-3,p-2\), the first possible term below valuation
one is at \(L=j-h\). Terms with smaller \(L\) have \(R\ge p\).
By (6), if \(v(d)\ne s\), the pure term at this first
index \(L\) using the smaller of \(v(d),s\) is uniquely
lowest and nonzero. Hence \(v(d)=s\). Put
\(X=d/(k(t-1)^2)\) and
\(g_i(X)=\sum_{\ell=0}^i\binom i\ell(2\ell)!X^\ell\).
At \(L=j-h\) the term in (5) is a unit multiplier times
\(k^L g_L(X)\), up to an error of valuation at least
\(1+Ls\). All other terms have valuation at least
\(Ls+s\). Thus \(g_{h-2}(\bar X)=g_{h-1}(\bar X)=0\).

The first of these polynomials has a simple root there
and the second has derivative zero there. In fact
\(Xg'_i=i(g_i-g_{i-1})\), and the recurrence
\[
 g_{i+1}=1+2(i+1)X[(2i+1)g_i-2i g_{i-1}]
\]
prevents three consecutive zeros, including here in
characteristic \(p\). It also gives \(g_h(\bar X)=1\).
Lift the simple root of \(g_{h-2}\) by Hensel's lemma
in the finite unramified extension of \(\mathbb Q_p\)
with residue field containing \(\bar X\). From \(M_{p-3}=0\) the distance from \(X\) to that lift
has valuation at least \(s\). Taylor expansion of
\(g_{h-1}\) about the lift gives \(v(g_{h-1}(X))\ge2s\),
since the value and derivative there are zero modulo \(p\)
and \(2s\le1\). This suppresses the term at \(L=h-1\)
in \(M_{p-2}\) strictly beyond valuation one.

Write \(G_r(t)=\sum_{\ell=0}^r t^\ell\).
At valuation one, \(M_{p-1}\) has only \(L=0,h\).
At the same valuation \(M_{p-2}\) has only \(L=0,h\),
by the suppression just proved. In the latter formula
\(S_{p-3,i}\) equals \(S_{p-1,i}\) modulo \(p\) for \(i\ge1\),
since both removed tail terms contain a factor \(p\).
At \(i=0\) it differs by \(-t^{p-2}-t^{p-1}\).
Divide the two equations by \(p\) and reduce. Setting
\(\kappa=\overline{k^h/p}\), \(C=(-1)^h\), and omitting
bars within the residue equations, this gives
\[
 \begin{split}
  0&=-G_{2p-2}(t)-\kappa C\,G_{p-1}(t),\\
  0&=-\tfrac16 G_{2p-4}(t)
        -\tfrac14\kappa C\bigl(G_{p-1}(t)-t^{p-2}-t^{p-1}\bigr).
 \end{split}                                               \tag{7}
\]
The constants follow from Wilson's theorem and
\(\binom{p-1}{h}\equiv C,\ \binom{p-2}{h}\equiv C/2\).
All uses of \(g_h=1\) are in the residue field.

Put \(z=t^p\), so \(t=z^{-2}\) by (4).
The numbers \(G_{2p-2}\) and \(G_{p-1}\) are nonzero,
because \(z\ne1\) and \(z^{2p+1}=1\).
Using
\[
 \frac{G_{2p-4}}{G_{2p-2}}=z^4+1,\qquad
 \frac{t^{p-2}+t^{p-1}}{G_{p-1}}=-z(z+1)(z^2+1)
\]
in (7) gives \(F(z)=z^4+3z^3+3z^2+3z+1=0\).
Since \(z^p\) is also a root and its square is \(z^{-1}\),
the palindromic \(F\) has a common root with \(F(Z^2)\).
Their resultant is \(1936=16\cdot11^2\).
For a short certificate, one has
\[
 \begin{split}
 22={}&(-12Z^7-3Z^6-30Z^5-5Z^4-27Z^3+3Z^2-36Z-1)F(Z)\\
      &+(12Z^3+39Z^2+39Z+23)F(Z^2).
 \end{split}
\]
Thus only \(p=11\) is possible. There \(z\ne1\) would
have order \(23\) and lie in an extension of degree at most
four, from its quartic equation. None of
\(11,11^2,11^3,11^4\) is \(1\bmod23\).
This excludes the off-diagonal case.

## The diagonal transition \(s=1/h\)

Now \(v(u)>0,\ v(c-1)>0\). Under \(s=1/h\), \(B_p(Z)\)
has its strictly lowest term at \(p\) on the unit disk.
Equation (2) at \(p\) gives
\[
 D=-u(c-1)W,\qquad
 W=\frac{B_p(c)-B_p(1)}{(c-1)B_p(c)} .
\]
The valuation of \(W\) is positive. After scaling by the
leading coefficient of \(B_p\), its numerator divided
difference has reduction \((Z^p-1)/(Z-1)\) at \(Z=1\),
which is zero in characteristic \(p\); higher coefficients
are strictly positive in valuation. Since \(c-1=2u+D\),
we obtain \(v(D)>2\rho\) for \(\rho=v(c-1)=v(u)>0\).

The ratio \(B_{p-1}(Z)/B_p(Z)\), after a constant
normalization, reduces to \(J Z^{h-p}+E Z^{-1}\), with
\(J=\binom{p-1}{h}(p-1)!k^h/p\) and \(E=(2p-2)!/p\)
read in the residue field. Its derivative at \(1\) must
vanish for (1) to hold with \(c\ne1\), by the residue-disk
injectivity argument above. Thus \(E=hJ\) and
\[
                   \overline{k^h/p}=\frac{1}{h(-1)^h}.
                                                               \tag{8}
\]

For the interior, expand the homogeneous phase as
\[
       Q=ax^2+bxy+cy^2=(x+y)(x+c y)-Dxy,\qquad a=1 .
\]
Use \(z=y/(x+y)\) to obtain the exact homogeneous moment
\[
 \Lambda_2(Q^j)
 =(2j+1)!\sum_{i=0}^j\sum_{\ell=0}^{j-i}
 \binom ji\binom{j-i}{\ell}(-D)^i(c-1)^\ell
       \frac{(i+\ell)!\,i!}{(2i+\ell+1)!}.                \tag{9}
\]
The identity follows from the elementary beta integral,
or by monomial expansion. For \(h\le j\le p-1\),
a term of (9) can contribute at valuation below one
only when \(2i+\ell\ge p-1\), since the numerator prefactor
has valuation one and the denominator has at most one
factor of \(p\). Such a term has valuation at least
\((p-1)\rho+i(v(D)-2\rho)\).

If \((p-1)\rho<1\), the unique lowest term of \(M_{p-1}\)
is \(j=p-1,i=0,\ell=p-1\) in (9), with valuation
\((p-1)\rho\). In fact for \(i\ge1\) the bound is strictly
larger, for \(j<h\) the power of \(k\) is at least \(h+1\),
and for \(h\le j<p-1\) a power of \(k\) gives a positive
extra valuation. Its coefficient before \((c-1)^{p-1}\)
is \((2p-1)!/p\), a unit. Thus \(M_{p-1}=0\) requires
\((p-1)\rho\ge1\).

Consequently the terms with \(i>0\) in \(M_{p-2}\) can
be dropped to valuation strictly beyond one: those at
\(j=p-2\) with \(2i+\ell<p-1\) have at least one factor
of \(p\) and an additional positive \(v(D)\); the others
use the strict bound \(v(D)>2\rho\), or an extra power of
\(k\). For \(j<h\) the power of \(k\) is at least \(h\),
and a \(D\) factor again increases valuation.
The homogeneous moments of the remaining plane \(D=0\) are
\((2j+1)!\sum_{\ell=0}^j c^\ell/(j+1)\).
At valuation one in \(M_{p-2}\) only \(j=p-2,h-1\)
contribute. Divide by \(p\), use \(\bar c=1\), Wilson's
theorem and (8); the residue is
\[
             \frac12+\frac{(-1)^h}{2}\frac{1}{h(-1)^h}
                         =-\frac12\ne0 .
\]
This excludes \(s=1/h\).

## The diagonal transition \(s=2/p\)

Use (5) and (6) with \(v(u)>0\), \(v(d)>0\).
For \(j=h+L_0\), with \(L_0=h-2,h-1,h\), let
\[
          q=\min(v(d),s+2v(u)),\qquad
          A_{L_0}=2(h-L_0)v(u)+L_0q .
\]
At \(L=L_0\) the terms of (5) have least valuation at
least \(A_{L_0}\) up to an error of valuation at least one.
If the two entries defining \(q\) are distinct, the lowest
term giving \(A_{L_0}\), whenever \(A_{L_0}<1\), is unique
and has nonzero leading factor. At a tie the group equals
a unit multiplier times \(k^{L_0}u^{2h}
g_{L_0}(d/(ku^2))\), up to valuation at least one.
At \(L_0=h\) the errors have valuation strictly above one.
These statements follow from (6) before evaluating
\(t=1+u\); its error coefficients are multiples of \(p\).

Here are the bounds on the other terms. For \(L<L_0\)
one has \(R\ge p\), so the terms have valuation at least
one, strictly above one if \(L>0\). For \(L>L_0\), set
\(v_0=L-L_0>0\). If \(i\ge v_0\), every removed tail term
of \(S_{p-1,i}\) when replacing it by \(S_{p-1-2v_0,i}\)
has a factor of \(p\). Indeed its index \(\ell\) is at
least \(p-2v_0\), giving \(\ell+2i\ge p\).
The terms with \(i\le h\) then have the bound
\[
  (L-i)s+i v(d)+2(h-i)v(u)
\]
up to errors of valuation strictly above one; for \(i>h\)
the entire term has valuation above one.
If \(i<v_0\), the power of \(k\) is at least \(L_0+1\).
It has valuation strictly above \((L_0+1)s\) unless
\(L=L_0+1,\ i=0\). This exceptional term has valuation
exactly \((L_0+1)s\) for \(L_0=h-2,h-1\):
the truncated sum \(G_{p-3}(t)\) has residue \(p-2\ne0\)
and the factorial and binomial coefficients are units.
It lies above one for \(L_0=h\).

The higher groups with \(i\ge v_0\) gain strictly
positive valuation over \(A_{L_0}\). To see this first note
that at \(p-1\), \(B_{p-1}(Z)\) has its strictly lowest term
at \(h\), so the expression \(W\) in (2) has residue
\(h\ne0\). Thus \(v(D)=2v(u)\) and
\(D/u^2\equiv-2h\equiv1\), using \(D=-u(2u+D)W\) as before.
In particular \(v(d)>2v(u)\). Since \(i\le L\), the
displayed bound for the higher group is at least
\[
 A_{L_0}+(L-L_0)\min\bigl(s,v(d)-2v(u)\bigr).
\]
This also gives a gain of at least \(s\) at the tie
\(v(d)=s+2v(u)\).

For \(j=p-1\), the \(L=0\) term of (5) has valuation
one: \((2p-2)!\) has valuation one and
\(G_{2p-2}(t)\) reduces to \(2p-1\ne0\).
If the two candidates for \(q\) are distinct,
\(M_{p-1}=0\) now forces \(hq=1\).
The \(p-2\) and \(p-3\) equations each force
\(A_{L_0}=(L_0+1)s\), because otherwise the unique
lowest group or the unique truncated term dominates.
With \(q=1/h\), subtract the two equations to get
\(q-2v(u)=s\). Substituting into either equation
gives \(hq=(h+1)s\), contradicting
\(hq=1<(h+1)2/p\).

It remains to consider the tie. Write \(U=2h v(u)>0\).
Then \(A_{L_0}=L_0s+U\).
The terminal interior equation \(M_{p-1}=0\) requires
\(hs+U\le1\), so \(U\le1/p<s\).
The other two equations force
\(g_{h-2}(\bar X)=g_{h-1}(\bar X)=0\), with
\(X=d/(ku^2)\), because their truncated term is
higher than the first group's least valuation.
When \(h=2\), \(g_0=1\) already contradicts this.
Otherwise set \(\gamma=s-U>0\). From \(M_{p-3}=0\),
the first \(g\)-value has valuation at least
\(\gamma\); all other terms are at least
\((h-1)s\). As in the off-diagonal case, lift its
simple root in an unramified extension and use the
derivative of \(g_{h-1}\).
This gives \(v(g_{h-1}(X))\ge2\gamma\), since
\(2\gamma\le1\). In \(M_{p-2}\) the first group is
now strictly above \(hs\), as are the higher groups
with \(i\ge v_0\) and all the earlier groups.
The unique truncated term of valuation \(hs\) survives.
This last contradiction proves the theorem.

## A prime in the penultimate order

For a possible four-zero block at \(p-2,p-1,p,p+1\), the
same projective normalization forces \(a=1,\ \bar t=\bar c=1\)
after exchange of axes, together with
\[
                        v(k)=2/p,\qquad
                        \overline{k^p/p^2}=-2.             \tag{10}
\]
These are necessary conditions in every place over \(p\).
The [prime-penultimate theorem](BINARY_EVEN_QUADRATIC_PRIME_PENULTIMATE_FOUR_WINDOW.md)
excludes the remaining diagonal case by using the next
interior equation and the higher boundary terms.

The equations at \(p,p-1\) again give \(v(k)>0\), \(v(d)>0\).
The boundary cancellations are now at \(p,p+1\).
For unit \(Z\), the strictly lowest term of \(B_{p+1}(Z)\)
is \(2(p+1)k^p Z\) for \(s<2/p\), and
\((2p+2)!Z^{p+1}\) for \(s>2/p\), with just those two terms
tied at equality. In fact the binomial coefficient is a unit
exactly at indices \(0,1,p,p+1\). Terms at \(2\le i\le h\)
have valuation \(1+(p+1-i)s\), and all the other interior
competitors lie above the stated minima. Away from the tie
\(B_{p+1}/B_p\) reduces, after a constant normalization,
to \(Z\) times a nonzero constant. This excludes distinct
unit arguments, by direct residue comparison on different
disks and residue-disk injectivity on the same disk.
When \(c\) is nonunit, the single equation at \(p\)
already disallows both sides of the tie, as in the terminal
proof. Hence \(s=2/p\).

At equality, set \(\kappa=\overline{k^p/p^2}\). The reductions
after division by \(p^2\) are
\[
      B_p(Z)/p^2\equiv\kappa+2Z^p,\qquad
      B_{p+1}(Z)/p^2\equiv2 Z(\kappa+2 Z^p).
\]
Wilson's theorem gives \((2p)!/p^2\equiv2\), and this
also gives the leading coefficient at \(p+1\).
If \(\bar t\ne0,\pm1\), both weights in (1) are units and
\(\bar c=\bar t^2\ne1\). The two residue equations force
both \(\kappa+2\) and \(\kappa+2\bar c^p\) to vanish, which
is impossible. At \(\bar t=0,-1\), the equation at \(p\)
forces \(\kappa=-2\). At \(\bar t=1\), if
\(\kappa\ne-2\) the quotient is defined on the residue disk
and has nonzero reduced linear coefficient, again giving
injectivity. Thus \(\kappa=-2\) there as well.

To exclude \(\bar t=0,-1\), use the first two interior
equations. In (5) the first group below valuation one at
\(j=p-2,p-1\) is \(L=h-1,h\), respectively, since
\(hs<1\) and \(\bar t\ne1\). The same unique-term argument
using (6) forces \(v(d)=s\) and
\(g_{h-1}(\bar X)=g_h(\bar X)=0\), for \(X=d/[k(t-1)^2]\).
The first root is simple, by the consecutive recurrence.
The remaining terms of \(M_{p-2}\) have valuation at least
\(hs\), so \(v(g_{h-1}(X))\ge s\). An unramified lift
and \(g_h'(\bar X)=0\) give \(v(g_h(X))\ge2s\).
Thus the \(L=h\) group in \(M_{p-1}\) is strictly past
valuation one, as are all later groups. The \(L=0\) term
has valuation one because \(G_{2p-2}(t)\) reduces to one
at both \(0\) and \(-1\). It cannot cancel. This proves
the necessary diagonal condition (10).

## Exact checks

The [replay](verify_binary_even_quadratic_prime_terminal_four_window.py)
compares (5) and (9) with independent native multinomial
expansions and checks the boundary valuation tables, the
finite-field rising sums, the critical residue in (8), and
the quartic certificate. It also verifies at \(p=13\)
that \(g_{h-2},g_{h-1}\) can share an off-diagonal candidate
parameter modulo \(p\), illustrating the need for the
higher terms in (7). The proof applies to all primes
regardless of the finite range of this replay.
