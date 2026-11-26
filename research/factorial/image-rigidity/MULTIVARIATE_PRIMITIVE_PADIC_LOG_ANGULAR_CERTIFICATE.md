# Exact logarithmic coefficients for native prime angular separation

Fix an ordinary factorial binomial and adjacent orders \(n,N=n+1\)
with \(N\ge3\). Write \(F_N(z)=\sum_{k=0}^N a_kz^k\).
The following test uses the upper coefficients and a prime \(p>N+1\).
As in the [native prime chord test](MULTIVARIATE_NATIVE_PRIME_CHORD_ANGULAR_CERTIFICATE.md#a-condition-on-the-upper-newton-polygon),
assume the upper Newton polygon has a strictly separated endpoint
chord of denominator \(N\). More explicitly, with \(v_p(p)=1\)
put
\[
 A=v_p(a_0),\quad s=v_p(a_0)-v_p(a_N),\quad
 \nu_k=N(v_p(a_k)-A)+sk\quad(1\le k<N),
\]
and assume \(\gcd(s,N)=1\) and \(\nu_k>0\) for all \(k\).
Set \(\pi=p^{1/N}\). Then
\[
 p^{-A} F_N(\pi^s x)=c_0+c_N x^N+
                     \sum_{k=1}^{N-1}b_k\pi^{\nu_k}x^k
\]
with \(c_0,c_N,b_k\) rational \(p\)-adic units.
There is no hypothesis on the lower Newton polygon or on whether
its roots are units.

**Theorem.** For any \(1\le r<N\) coprime to \(N\), the scalar series
\(M_r\) specified below is in \(\mathbb Q_p\) and can be computed to
any prescribed \(p\)-adic precision by rational arithmetic.
If any \(M_r\ne0\), the ordinary pair has disjoint root arguments
under independent positive rescalings. In particular, the
condition holds if for some \(k_*\) coprime to \(N\) a divisor
\(g>1\) of \(N\) satisfies
\[
        g\mid k \quad\hbox{whenever } \nu_k<\nu_{k_*}.          \tag{1}
\]
This permits primitive angular information arbitrarily deep
above the first nonprimitive interior term. The test concerns
the original rows and uses no lower moment coefficient tests.
If a certifying prime exists for every \(3\le N<J\), the
[absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md)
and positivity of the quadratic moment give separation in
every degree.

## The exact upper logarithmic coefficient

Let \(\lambda_i\) be the distinct unit roots of
\(c_0+c_Nx^N\) in an unramified extension of \(\mathbb Q_p\).
Since \(p>N\), these have distinct residues, and each lifts
to a unique root \(x_i\) of the scaled polynomial near
\(\lambda_i\). Write \(x_i=\lambda_i\exp W_i\) for the principal
logarithm \(W_i=\log(x_i/\lambda_i)\). The smallest interior
exponent \(\beta=\min_k \nu_k\) is positive, and \(p>N+1\)
ensures \(\beta/N>1/(p-1)\). Thus the exponential and
logarithm used here converge and are inverse at this scale.

For formal variables \(t_k\), let \(U=\sum_k t_k(1+U)^{k/N}\)
be its unique formal solution with zero constant term.
Put \(W=N^{-1}\log(1+U)\). For the multiindex
\(\ell=(\ell_1,\ldots,\ell_{N-1})\), use
\[
       L=\sum_k\ell_k,\qquad K_\ell=\sum_k k\ell_k,\qquad
       w_\ell=\sum_k\nu_k\ell_k .
\]
When \(L>0\), the exact coefficient is
\[
 [t^\ell]W=\frac1{NL}\binom{L}{\ell_1,\ldots,\ell_{N-1}}
                       \binom{K_\ell/N-1}{L-1}.                \tag{2}
\]
To see this, insert a variable \(z\) multiplying the
right side of the equation for \(U\). The coefficient of
\(z^L\) in the logarithm is \(1/(NL)\) times the coefficient
of \(u^{L-1}\) in
\((1+u)^{-1}(\sum_k t_k(1+u)^{k/N})^L\), by formal
Lagrange inversion. Expanding the multinomial gives (2).
The coefficient valuation is at least \(-v_p(L)\), since
\(K_\ell/N-1\in\mathbb Z_p\) and its binomial coefficients
are integral. This last fact follows by continuity from
nonnegative integer arguments to the binomial polynomial.

To obtain \(W_i\), substitute
\(t_k=(b_k/c_0)\lambda_i^k\pi^{\nu_k}\).
Indeed the scaled polynomial equation reads
\(e^{NW_i}-1=\sum_k t_k e^{kW_i}\).
The formal manipulation evaluates to the Hensel solution.
For example \(U(t)\) has integral coefficients by its
defining recursion and converges for the indicated
positive term valuations. The recursion gives
\(v_p(U)\ge\beta/N\), and (2) converges after substitution
because \(w_\ell/N-v_p(L)\) tends to infinity with \(L\).

For a residue \(r\in\{1,\ldots,N-1\}\), let
\(k_r\in\{1,\ldots,N-1\}\) be congruent to \(s^{-1}r\bmod N\).
Since \(w_\ell\equiv sK_\ell\bmod N\), the coordinate of
\(W_i\) in the \(\pi^r\) direction is exactly
\(\lambda_i^{k_r}M_r\), where
\[
 M_r=\sum_{\substack{L>0\\w_\ell\equiv r\pmod N}}
   \frac1{NL}\binom{L}{\ell_1,\ldots,\ell_{N-1}}
      \binom{K_\ell/N-1}{L-1}
      \left(-\frac{c_0}{c_N}\right)^{(K_\ell-k_r)/N}
      p^{(w_\ell-r)/N}\prod_k(b_k/c_0)^{\ell_k}.
                                                               \tag{3}
\]
The power-basis coordinate lies over the unramified field.
Factoring \(\lambda_i^{k_r}\) out of all terms leaves
the rational scalar series shown, since
\(\lambda_i^N=-c_0/c_N\).
The terms of (3) have valuations at least
\[
             \frac{L\beta-r}{N}-\frac{L-1}{p-1}.
\]
This uses \(v_p(L)\le (L-1)/(p-1)\), obtained from
\(p^{v_p(L)}-1\ge(p-1)v_p(L)\).
The bound increases unboundedly with \(L\), so it gives
an explicit finite truncation for computing (3) modulo
each \(p^b\). A nonzero value is detectable at finite precision.

## Lower ratios and the subgroup test

For distinct upper roots their ratio has residue
\(\omega=\lambda_i/\lambda_j\). The logarithm of the ratio
divided by \(\omega\) is \(W_i-W_j\). Its coordinate at
a primitive residue \(r\), meaning \(\gcd(r,N)=1\),
equals \((\lambda_i^{k_r}-\lambda_j^{k_r})M_r\).
The factor preceding \(M_r\) is a unit because \(k_r\)
is coprime to \(N\) and the binomial root residues are distinct.

Any two lower roots \(\alpha,\alpha'\) have local degrees at most
\(n<p\). After adjoining a common finite unramified extension
\(K\) they lie in \(K(p^{1/e})\) and \(K(p^{1/e'})\) for their
ramification indices \(e,e'\le n\), respectively. Indeed, for
a uniformizer in a tame extension, Hensel lifting an \(e\)-th
root of its unit quotient with \(p\) exhibits a root of
\(x^e-p\) in the unramified base change. Its degree is exactly
\(e\), so it generates that base change. This is the
[single-row tame basis comparison](MULTIVARIATE_SINGLE_MOVABLE_ROW_PADIC_PROGRESSION.md#the-lower-ratios-cannot-have-that-displacement).

If \(\alpha/\alpha'\) were equal to the upper ratio, the lower
roots would have the same valuation \(v\).
The power \(p^v\) can be chosen in the intersection of the
two tame subfields. Divide by this same power in both fields,
then choose unramified unit lifts of the residues with quotient
\(\omega\). The log of the resulting ratio is a difference of
two logs of principal units in the respective tame subfields.
In the basis of the combined tame extension of degree
\(\operatorname{lcm}(N,e,e')\) over \(K\), its support has
fractional parts only from \((1/e)\mathbb Z\) and
\((1/e')\mathbb Z\). Neither subgroup includes a coordinate
with reduced denominator \(N\). This contradicts any nonzero
primitive coordinate just found. Hence upper and lower
off-diagonal quotients are disjoint. Equal nonreal arguments
would give a common quotient of roots by their conjugates.
A common real argument is excluded by strict positivity of the
even ordinary moment on the real line.

For (1), set \(r=s k_*\bmod N\), with \(1\le r<N\).
In every multiindex of (3), the integer \(K_\ell\)
is congruent to \(k_*\bmod N\), so it is not divisible by \(g\).
At least one term in that multiindex comes from an index
not divisible by \(g\), and its \(\nu_k\) is at least
\(\nu_{k_*}\) by (1). For \(L\ge2\), this gives
\(w_\ell\ge\nu_{k_*}+(L-1)\beta\), and
\((L-1)\beta>Nv_p(L)\).
The latter follows from \(p-1>N\) if \(v_p(L)>0\),
and from \(L\ge2\) otherwise. All nonlinear terms in
(3) therefore have strictly larger valuation than the
single linear term from \(k_*\). This proves \(M_r\ne0\).
The argument uses the residue group and imposes no bound
of the form \(\nu_{k_*}<2\beta\).

## A subfield obstruction independent of the log expansion

There is also an exact condition that only uses the tame field
degrees. This part requires \(p>N\) and the same primitive
upper chord. Pass to a common unramified field \(K\) containing
the roots of the residual binomial and all prime-to-\(p\) roots
of unity needed below. The upper roots lie in the cyclic extension
\(L=K(p^{1/N})\) by the binomial rescaling and Hensel lifting.
Each original upper root generates \(L\): its valuation \(s/N\)
has reduced denominator \(N\), so its degree over \(K\) is at
least \(N\) and at most \([L:K]=N\).

Let \(\sigma\) generate \(\operatorname{Gal}(L/K)\), let
\(\alpha\) be an upper root, and put \(q=\sigma^j(\alpha)/\alpha\)
for \(1\le j<N\). Write \(r=[K(q):K]\), \(d=\gcd(N,j)\).
Then
\[
                          \operatorname{lcm}(r,d)=N.           \tag{4}
\]
To prove it, \(\sigma^r(q)=q\), so
\(\theta=\sigma^r(\alpha)/\alpha\) is fixed by \(\sigma^j\).
Indeed the two automorphisms commute, and
\(\sigma^j(\theta)=\sigma^r(q\alpha)/(q\alpha)=\theta\).
Let \(t=\operatorname{lcm}(r,d)\). The ratio
\(\sigma^t(\alpha)/\alpha\) is a product of
\(\sigma^{ar}(\theta)\), \(0\le a<t/r\), hence belongs to
the fixed field of \(\sigma^j\). It is fixed by \(\sigma^t\)
because \(d\mid t\). Applying successive powers of \(\sigma^t\)
shows that this ratio is an \((N/t)\)-th root of unity.
If \(t<N\), the ratio would be a nontrivial root of unity
\(\xi\in K\), since \(\alpha\) generates \(L\).
But then \(F_N(\xi z)\) and \(F_N(z)\) share the root \(\alpha\).
They have the same degree and \(F_N\) is irreducible over
\(K\), so they differ by a scalar. Comparing their nonzero
constant and linear coefficients forces \(\xi=1\).
This proves (4).

The residue of the unit ratio \(q\) has order \(N/d\), since
\(\sigma\) multiplies \(p^{1/N}\) by a primitive \(N\)-th
root of unity and \(\gcd(s,N)=1\).
Thus for every prime \(\ell\mid N/d\), the full power
\(\ell^{v_\ell(N)}\) divides \(r\).
A lower quotient from ramification indices \(e,e'\)
would lie in an extension of degree \(\operatorname{lcm}(e,e')\)
over the same unramified field. Hence a primitive upper chord
alone suffices for angular separation if no lower root
ramification index is divisible by any full prime-power
divisor \(\ell^{v_\ell(N)}\) of \(N\).
This last condition is automatic when \(N\) is a prime power,
because every lower ramification index is smaller than \(N\).
It also holds when all lower roots are unramified.
The logarithmic test above can apply when this field-degree
condition does not hold.

## An ordinary example with a later primitive coordinate

At \(N=6\), take the ordinary rows
\[
 (9,10),(9,10),(5,6),(1,2),(2,1),(6,5),(10,9).
\]
At \(p=11\) the scalar-free coefficient valuations from
\(k=0\) through \(6\) are \((3,5,4,5,4,5,4)\).
To check these, each row has gap one. The three ascending
bases \(9,5,1\) cross a factor divisible by \(11\) at offsets
\(1,3,5\), respectively, each with valuation one.
There are two copies of the first upward step.
The decreasing copies have the reflected downward steps.
No upper binomial factor is divisible by \(11\).
Thus \(s=-1\) and the interior exponents are
\((11,4,9,2,7)\).
The nearest interior exponent is \(2\), at \(k=4\),
and no primitive index has exponent below \(2\beta=4\).
But \(k_*=5\), \(\nu_5=7\), has only even indices
below it in the exponent ordering. Taking \(g=2\)
in (1) gives the exact separation at orders five and six
for the displayed original profile. Formula (3) also gives
\(M_1\equiv11\pmod{121}\): the tail with \(L\ge9\) has
valuation above two by the stated truncation bound.

## A primitive chord alone can leave a local quotient collision

For comparison, use \(p=13\), \(N=6\), and
\(\pi=13^{1/6}\). Let \(\zeta\in\mathbb Q_{13}\) be a
primitive sixth root of unity. The six numbers
\[
          -\zeta^j\pi(1+\zeta^{2j}\pi^2)(1+\zeta^{3j}\pi^3),
                         \qquad 0\le j<6,
\]
are the conjugate roots of a monic degree-six polynomial in
\(\mathbb Z[z]\). For example it is the resultant of
\(t^6-13\) and \(z+t(1+t^2)(1+t^3)\).
All six roots have valuation \(1/6\), so the endpoint valuation
difference is \(-1\) and every intermediate coefficient lies
strictly above the chord by integrality of its valuation.
The polynomial is
\[
 z^6+78z^5+2496z^4+42328z^3+402480z^2+2044224z+4402944 .
\]

The upper quotient for \(j=1\) and \(j=0\) is
\(\zeta A B\), where
\[
 A=\frac{1-\pi^3}{1+\pi^3},\qquad
 B=\frac{1+\zeta^2\pi^2}{1+\pi^2}.
\]
It is also the quotient of \(\zeta A\) and \(B^{-1}\).
Those elements lie in quadratic and cubic extensions of
\(\mathbb Q_{13}\), respectively. The product of their two
local minimal polynomials, padded to degree five with monic
linear factors if necessary, gives a lower polynomial with
that exact common quotient. The logarithm of the upper ratio
divided by \(\zeta\) is \(\log A+\log B\), supported in the
two smaller tame subfields. This is a collision of local
quotients for a relaxed polynomial pair. The lower
coefficients are local, and this construction does not
give an ordinary factorial common root or common ray.
