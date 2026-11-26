# Unbalanced factorial binomials at prime-power orders

Suppose both orientations of an exponent gap occur, with \(R>S\ge1\)
increasing and decreasing rows, after reversing the variable if
needed. At upper orders \(Q=p^q\) there are two exact separation
criteria.

- For unit gaps with positive smaller exponents \(A_i\) on the
  increasing side and \(B_j\) on the decreasing side, assume
  \(p>\max A_i+\max B_j\) and
  \(p-1\ge(R-S)(\max A_i-1)\). The common factor has degree at
  most \(\min A_i-1\). In particular, if \(\min A_i\le2\), all
  powers of each such prime give coprimality. A residual test
  independent of \(q\) applies without this minimum restriction.
- For gap-two rows \((1,3)\) and \((3,1)\), assume \(p>3\).
  Coprimality holds when \(p\equiv5\bmod6\) and \(q\) is even,
  and when \(p\equiv1\bmod6\) and \(q\equiv0\bmod(p-1)\).
  At each fixed \(q,R,S\), an explicit nonzero integer excludes
  finitely many primes and certifies all remaining ones.

The Newton polygon product rule used here is
[Filaseta and Rocks, Theorem 2.1](https://people.math.sc.edu/filaseta/papers/ColloquiumMathPaper2016.pdf),
with the polygon read from the constant coefficient to the leading
coefficient. A common factor can use only common slopes, and at
each common slope its residual polynomial divides both residual
polynomials after any monomial factors have been removed.
At the odd prime powers considered here the lower order \(Q-1\)
is even. Its real moment is positive, since the integrand
\((U+zV)^{Q-1}\) for real \(z\) is nonnegative and
\(U+zV\) is nonzero away from a set of Lebesgue measure zero
under the independent Gamma representation. A common rational
factor therefore has even degree.

## Unit gaps

Write \(n=Q-1\), \(w_0=(R-S)/(p-1)>0\), and
\(k_0=\min A_i-1\). Let \(v(k!)=v_p(k!)\) and let
\(s_p(k)\) denote the digit sum. For the upper order,
after subtraction of scalar valuation terms, the interior
ordinates are
\[
 a_k=(R-S)v(k!)-(S-1)(q-v_p(k)),\qquad 1\le k<Q,
\]
with endpoint ordinates \(a_0=0\),
\(a_Q=(R-S)v(Q!)\).
Indeed, shifting a factorial argument by a multiple of \(Q\)
preserves the valuations of all its short increments as long as
the arguments are below \(pQ\). This reduces the valuation to
those of \(k!, (Q-k)!\) and the binomial coefficient.
The binomial valuation is \(q-v_p(k)\).
After subtracting \(w_0 k\) the unique minimum is at \(k=Q-1\):
among interior indices it alone has digit sum \(q(p-1)\), and it
is a \(p\)-unit. Its ordinate is \(-(R-1)q\), strictly below
the endpoints \(0,-w_0\).
Hence all earlier polygon slopes are strictly less than \(w_0\),
and the last slope at width one is \((R-1)q>w_0\).

For the lower order, first shift each row ending at \((Q-1)A_i\)
or \((Q-1)B_j\) to the corresponding multiple of \(Q\).
Subtract the scalar row terms and \(S v(n!)\).
Because \(s_p(k)+s_p(n-k)=s_p(n)\) at \(n=Q-1\), the
ordinate after subtracting \(w_0 k\) is
\[
                            -w_0s_p(k)-U(k).
\]
Here \(U(k)\) is the valuation lost by removing the \(A_i\) and
\(B_j\) factorial tails. Each tail has length less than \(p\).
A left loss occurs only for \(r=k\bmod p<A_i\); all active
left rows then lose \(t=v_p(k-r)\), interpreted as \(q\) when
\(k-r=0\). A right loss occurs only for \(r\ge p-B_j\).
The left and right bands are disjoint by the prime-size
hypothesis. The right loss is at most \(Sq\), and
\(s_p(k)\le q(p-1)\), giving a shifted ordinate at least
\(-Rq\) in the right band and outside both bands.
In the left band with \(t<q\),
\(s_p(k)\le(q-t)(p-1)+r\); so the shifted ordinate is at least
\(-Rq+S(q-t)-w_0r\ge-Rq\).
Here \(w_0r\le1\) by the second size hypothesis.
If \(t=q\), then \(k=r<p\).
For \(k\le k_0\) all left rows are active and the unshifted
ordinate is exactly \(-Rq\). For \(k>k_0\), at least one
left row is inactive and the shifted ordinate is at least
\(-Rq+q-w_0r\ge-Rq\). At \(k=n\) the shifted ordinate is
exactly \(-Rq\).

It follows that the lower polygon consists of a horizontal
segment of width \(k_0\) followed by the segment from \(k_0\)
to \(n\) of slope \(w_0n/(n-k_0)\).
Indeed the shifted line from \(k_0\) to \(n\) lies at or below
\(-Rq\), and the points after \(k_0\) lie at or above that
number. This positive slope is at least \(w_0\) and strictly
below \((R-1)q\), since
\[
 \frac{w_0 n}{n-k_0}\le\frac{R-S}{p-\min A_i}< (R-1)q .
\]
Only the horizontal slope can be shared. This proves the
degree bound, and the even-degree condition gives the
assertion for \(k_0\le1\).

At \(q=1\) the second size hypothesis can be omitted. Then
\(v(k!)=0\) up to \(k=n\), and the lower ordinates are the
negatives of the active row counts on the two disjoint end
bands. Their left minimum has width \(k_0\), and the right
endpoint ordinate is \(-S\). Any point below that right
endpoint has \(k<\max A_i\). Consequently the largest lower
slope is at most \((R-S)/(p-\max A_i)<R-1\).
The only positive upper slope is \(R-1\), while all lower
slopes are nonnegative. Thus the common slopes again consist
at most of the horizontal one.

For a stronger residue test put \(z=(-1)^{S+1}t\).
On the lower horizontal face the residual polynomial is
\[
 L_A(z)=\sum_{k=0}^{k_0}(-1)^{Rk}
       \frac{\prod_i(A_i-1)^{\underline{k}}}
            {\prod_j(B_j+1)_k}\,z^k,
\]
where an underline denotes a falling factorial and
\((B_j+1)_k\) is rising. This follows by taking the
successive factorial ratios from the constant coefficient,
using \(\binom{Q-1}k\equiv(-1)^k\).
The upper horizontal face is supported on \(1\le k<p\)
when \(S\ge2\) and on \(0\le k<p\) when \(S=1\). On
\(1\le k<p\), its units after a common normalization are
\((k!)^{R-S}k^{S-1}\) in \(z\), since the binomial unit is
\((-1)^{k-1}/k\), each increasing tail gives \(k!\), and each
decreasing tail gives
\((-1)^{k-1}/((B_j+1)(k-1)!)\) after restoring its one
factor \(Q\). Add a constant \(B_1+1\) when \(S=1\);
otherwise remove the monomial \(z\) in the upper residue.
A gcd of degree at most one between these two modular
polynomials certifies coprimality for every \(q\) at that
prime, by the even-degree condition and the polygon product rule.

## The gap-two face

Now take \(R\) rows \((1,3)\) and \(S\) rows \((3,1)\), and let
\(D_0=(p-1)/2\), \(s_0=(R-S)/D_0\) and
\(P_r=(p^r-1)/(p-1)\). The same upper valuation calculation
replaces \(k!, (Q-k)!\) by \((2k)!, (2(Q-k))!\).
The valuation of \(\binom{2Q}{2k}\) is again \(q-v_p(k)\),
by the factorial-floor formula because \(p>3\).
The upper ordinate after subtracting \(s_0 k\) and the
scalar baseline is therefore
\[
     -\frac{R-S}{p-1}s_p(2k)-(S-1)(q-v_p(k))
                           \qquad(0<k<Q).
\]
For an even integer \(0<v<2Q\) one has
\(s_p(v)\le q(p-1)\). If \(v<Q\), equality means \(v=Q-1\).
If \(v\ge Q\), its lower \(q\) digits cannot all be \(p-1\),
by parity, and equality means those digits are one below their
maximum at exactly one place. Thus the shifted upper minimum
is attained exactly at
\[
                     k_r=Q-(p^r+1)/2,\qquad 0\le r\le q .
\]
All these indices are \(p\)-units and the two endpoints of the
original polynomial lie higher than this minimum.

For the lower order, put \(l=2k,\ v=2(n-k)\), so that
\(l+v=2Q-2\). The shifted ordinate, after subtraction of a
scalar baseline, is
\[
 \frac{S s_p(2n)-R s_p(l)-S s_p(v)}{p-1}
       -R t_l-S t_v,\qquad
 t_l=v_p(l),\quad t_v=v_p(v),
\]
with \(t=q\) when the corresponding argument is zero.
Indeed each factorial is now shifted to \(Q+l\) or \(Q+v\),
and removing its last factor costs the indicated valuation.
Because \(s_p(2n)=q(p-1)\), the digit bound applied to each
of \(l/p^{t_l},v/p^{t_v}\) gives a shifted ordinate at least
\(-Rq\), with equality at the two endpoints. In the interior,
if neither is divisible by \(p\), equality needs both digit
sums maximal and hence \(l=v=Q-1\).
If an argument, say \(l\), is divisible by \(p^t\), \(t>0\),
the other is a \(p\)-unit. Equality would require the other
argument to be at least \(Q-1\), so \(l\le Q-1\).
The bound for \(l/p^t\) would force \(l=Q-p^t\). But then
\(v=Q+p^t-2\) has digit sum \(t(p-1)<q(p-1)\).
The same argument applies with the variables exchanged.
The lower polygon thus has only slope \(s_0\), with support
on its face at \(0,n/2,n\).

## Gap-two residual equations and prime exclusions

Write \(u_p(b!)\) for the \(p\)-unit of a factorial reduced
modulo \(p\). The identity
\[
          u_p(b!)=(-1)^{v_p(b!)}\prod_j b_j!
\]
for the digits \(b_j\) follows by successively removing full
blocks of \(p\)-units. In the lower polynomial these units
at \(b=n,2n,3n\) are respectively
\((-1)^{P_q},-1,(-1)^{P_q}\). The middle binomial unit is
\((-1)^{n/2}\). Introducing
\[
     y=(-1)^{D_0+R+S}z^{D_0},\qquad x=y^{P_q}
\]
on the residue face, the lower equation is
\[
                           x^2+a x+1=0,\qquad a=(-1)^{R+S}.
 \tag{1}
\]
Here \(z\) is the variable after the common slope scaling
\(t=\pi^{-(R-S)}z\), where \(\pi^{D_0}=p\).

For the upper face at \(k_r\), the factorial arguments are
\(3Q-p^r-1\) and \(Q+p^r+1\).
Their factorial valuations relative to \(r=q\) change by
\(P_q-P_r\) and \(-(P_q-P_r)\), respectively.
In the first factorial the digit factorial product ratio
relative to \(r=q\) is \(-2\) for \(r<q\).
In the second the ratio is \(1/2\) for \(0<r<q\) and \(1\)
at \(r=0\). The upper binomial unit ratio, after the
alternating sign, is \(1\) for \(0<r<q\) and \(1/2\) at
\(r=0\). These facts follow directly from the digits: for
\(r<q\) the first argument has top digit \(2\) and a single
lower digit \(p-2\) among digits \(p-1\), while the second has
top digit \(1\) and either two lower digits \(1\) or a lower
digit \(2\). At \(r=q\) the first has top digit \(1\) and
lower digits \(p-1\), and the second is \(2Q+1\).
Using \(k_r-k_q=D_0(P_q-P_r)\), the normalized upper face is
\[
       1+C\sum_{r=1}^{q-1}y^{P_q-P_r}+H y^{P_q},\qquad
             C=(-1)^R2^{R-S},\quad H=(-1)^R2^{R-1}.
 \tag{2}
\]

At a common residue root, \(y\ne0\). Let
\(Z=\sum_{r=0}^{q-1}y^{-P_r}\) and \(D=H/C=2^{S-1}\).
Equation (2) gives
\[
         Z=1-D-x^{-1}/C,\qquad Z^p=y(Z-1+x^{-1}).
 \tag{3}
\]
The second identity uses \(P_{r+1}=pP_r+1\).
The multiplier \(Z-1+x^{-1}\) cannot vanish at a common
root, since (3) would then imply \(Z=0\) and \(x=1\);
both roots in (1) are primitive third or sixth roots
for \(p>3\). Thus \(y\) lies in \(\mathbb F_p(x)\).
If \(p\equiv5\bmod6\) and \(q\) is even, then
\(y^{P_q}\) is a power of the norm to \(\mathbb F_p\),
in conflict with \(x\notin\mathbb F_p\).
If \(p\equiv1\bmod6\), then \(y\in\mathbb F_p\) and
\(x=y^q\). The latter equals one when \(p-1\mid q\).
The reduction therefore proves both asserted infinite
families.

For the finite prime exclusions at fixed \(q,R,S\), work in
the quadratic number field defined by (1), with
\(x^{-1}=t_0\). The two possibilities in (3) are the reductions of
\[
 y_+=\frac{C(1-D)-t_0}{-H+(C-1)t_0},\qquad
 y_-=\frac{C(1-D)-x}{-H+(C-1)t_0},
\]
for \(p\equiv1\bmod6\) and \(p\equiv5\bmod6\), respectively.
Both complex absolute values are less than one, since the
difference of the squared denominator and numerator moduli is
\((2+a)C(CD-1)>0\).
Thus \(y_+^q-x\ne0\), and for odd \(q\),
\(y_-\operatorname{Norm}(y_-)^{(q-1)/2}-x\ne0\).
The product of the numerators of their rational field norms,
together with the norm of the denominator, excludes only
finitely many primes. For even \(q\) omit the expression
involving \(y_-\). For each remaining prime the relevant
modular equality cannot hold. Since slope \(s_0\) is the
only possible common slope, this proves coprimality.

