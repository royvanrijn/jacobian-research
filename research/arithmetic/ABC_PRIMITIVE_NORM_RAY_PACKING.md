# Short real ratios and high powers in principal ideal norms

The angular packing for primitive quadratic norms extends to
higher-degree fields. The comparison uses the norm of a minor
in the normal closure and does not need a size bound for a
chosen unit or a fundamental domain.

## A close-ray count uniform in the field

Fix integers \(d\ge2\) and \(k\ge2\).
Let \(K/\mathbb Q\) have degree \(d\), let
\(\mathfrak a\) be a fixed nonzero fractional ideal of \(K\),
and consider elements \(\alpha\in K^\times\) satisfying
\[
        (\alpha)=\mathfrak a\mathfrak b^k,\qquad
        \mathfrak b\subseteq\mathcal O_K,\qquad
        N\mathfrak b\le X^\rho .
\]
For every embedding \(\sigma:K\hookrightarrow\mathbb C\),
fix a nonzero complex number \(\eta_\sigma\). Suppose that
each \(\alpha\) has a positive real scale \(r_\alpha\) with
\[
     \sigma(\alpha)
       =r_\alpha\eta_\sigma\bigl(1+O(X^{-\delta})\bigr)
                   \quad\text{for all }\sigma .
 \tag{NR1}
\]
The implied constant is fixed, uniformly in \(K\), \(\alpha\)
and \(\eta_\sigma\). If
\[
                              d\delta>2\rho ,
 \tag{NR2}
\]
then the number of distinct \(\mathbb Q^\times\)-rays of
\(\alpha\) is \(O_{d,k}(|\mathrm{Cl}(K)[k]|)\), for
sufficiently large \(X\) in terms of the margin in (NR2)
and the constant of (NR1).

The classes of \(\mathfrak b\) belong to a coset of
\(\mathrm{Cl}(K)[k]\) when any class occurs. In each
occurring class select an actual pair
\((\alpha_0,\mathfrak b_0)\).
Write \(\mathfrak b=\gamma\mathfrak b_0\) so that
\(\gamma\in\mathfrak b_0^{-1}\), and
\(\alpha/\alpha_0=u\gamma^k\) with \(u\) a unit.
Use one fixed representative \(u\) of its class modulo
\(k\)-th powers of units, absorbing the remaining unit
into \(\gamma\). There are \(O_{d,k}(1)\) unit classes
by Dirichlet's unit theorem and the fixed-degree bound
on the number of roots of unity. No estimate on the
absolute values of a unit representative enters.

For each \(\sigma\), fix one complex \(k\)-th root of
\(\sigma(u)^{-1}\). By (NR1), \(\sigma(\gamma)\) equals
that root, multiplied by a \(k\)-th root of unity and
\((r_\alpha/r_{\alpha_0})^{1/k}(1+O_k(X^{-\delta}))\).
Refine into the at most \(k^d\) possible branch choices.
If \(\gamma,\gamma'\) occur in one class and one such
unit and branch record, all complex conjugates of
\(\gamma/\gamma'\) have the same positive real scale
up to relative error \(O_k(X^{-\delta})\).

Let \(E\) be the normal closure of \(K\) and let
\(D=[E:\mathbb Q]\). For embeddings \(\sigma_i,\sigma_j\)
consider the difference
\[
        \Delta_{ij}
       =\sigma_i(\gamma)\sigma_j(\gamma')
                     -\sigma_i(\gamma')\sigma_j(\gamma).
\]
This lies in the product of the two conjugates of
\(\mathfrak b_0^{-1}\) extended to \(E\). If it is
nonzero, ideal inclusion gives
\[
                |N_{E/\mathbb Q}\Delta_{ij}|
                   \ge (N\mathfrak b_0)^{-2D/d}.
\]
On the other hand the common scale of the conjugates
gives
\[
        |\Delta_{ij}|\ll_k
             X^{-\delta}|\sigma_i(\gamma)\sigma_j(\gamma')|.
\]
The same comparison holds for every conjugate of this
minor. In the product over \(\operatorname{Gal}(E/\mathbb Q)\),
each embedding occurs exactly \(D/d\) times in each
position by transitivity. Consequently
\[
 |N_{E/\mathbb Q}\Delta_{ij}|
      \ll_{d,k} X^{-D\delta}
              \left(\frac{N\mathfrak b\,N\mathfrak b'}
                          {(N\mathfrak b_0)^2}\right)^{D/d}.
\]
This contradicts the lower bound for large \(X\) by
(NR2). Thus every minor vanishes. All embeddings of
\(\gamma/\gamma'\) then agree, so this quotient is
rational. Hence \(\alpha/\alpha'\) is rational too.
Each class, unit and branch record supplies one ray,
proving the claim. The fractional ideal denominator in
the norm lower bound cancels the representative ideal's
norm in the upper bound.

## Application to signed weighted powers

Fix \(q,k\ge2\) and \(0\le h<1\) with
\(k>2/(1-h)\). Choose a fixed sign \(s\in\{1,-1\}\).
Among positive primitive triples with largest value
\(c\le X\), suppose two letters are
\(A x^q\) and \(B y^q\), where \(A x^q\le c^h\).
The third is \(T=|A x^q+s B y^q|=C R^k\), and the
signed binomial \(A U^q+s B V^q\) is irreducible over
\(\mathbb Q\). Assume
\[
        \operatorname{rad}(AB)\le X^{f\lambda},
                    \qquad
        \operatorname{rad}(C)\le X^{g\lambda}.
\]
For every fixed budget and \(\varepsilon>0\) the count is
\[
      O_{q,k,h,\lambda,\varepsilon}
                \!\left(X^{((q+1)f/2+g)\lambda+\varepsilon}\right).
 \tag{NR3}
\]
More precisely, after fixing \(A,B,C\), the cost is
\(X^{o(1)}|\mathrm{Cl}(K)[k]|\), for the pure field
defined below. This permits replacing the class-number
loss whenever there is a stronger uniform torsion bound.
The constants are uniform in the radical caps \(f,g\).

Fix a dyadic shell \(c\asymp X\); the larger of the
two weighted powers has \(B y^q\asymp X\).
Enumerate the coefficients at a total cost
\(X^{(f+g)\lambda+\varepsilon}\) by the fixed-support
count. Take \(\theta\) to be a root of
\(Z^q+s A^{q-1}B\), which is irreducible by scaling
the given binomial. Set \(K=\mathbb Q(\theta)\) and
\(\alpha=A x-\theta y\). Its absolute norm is
\(A^{q-1}T\).

For \(p\nmid qAB\) dividing \(T\), the norm exponent
comes from exactly one degree-one prime ideal.
Indeed \(x,y\) are both units modulo \(p\), and the
equation for \(\alpha\) selects the single root
\(Ax/y\) of the monic polynomial modulo \(p\).
The polynomial discriminant is prime to \(p\),
so its distinct factors describe the prime ideals in
the full integer ring. Record which prime lies over
each such \(p\mid C\), at a cost \(X^{o(1)}\) by the
uniform bound on the number of prime factors. For
\(p\mid AB\), primitivity gives \(p\nmid T\); the
whole norm exponent of \(A^{q-1}\) on that support
is fixed. Its allocation is likewise subpower, since
the number of ideals of a fixed norm is at most a
fixed power of the divisor function. For primes
dividing the fixed \(q\) one can record the whole
ideal exponent vector in \(X^{o(1)}\) ways, as all
exponents are \(O_q(\log X)\).
These choices produce one fixed \(\mathfrak a\) and
a remaining \(\mathfrak b^k\) with
\(N\mathfrak b\le X^{1/k}\), up to an absolute dyadic
constant which does not affect the strict inequality (NR2).

Every conjugate of \(\theta\) has the same modulus,
and
\[
      \left|\frac{A x}{\sigma(\theta)y}\right|
         = \left(\frac{A x^q}{B y^q}\right)^{1/q}
         \ll_{q,h} X^{-(1-h)/q}.
\]
Thus (NR1) holds with \(r_\alpha=y\),
\(\eta_\sigma=-\sigma(\theta)\), and
\(\delta=(1-h)/q\). Since \(d=q\), (NR2) is exactly
\(k>2/(1-h)\). One rational ray gives one primitive
pair \(x,y\), because \(1,\theta\) are independent
and \(\gcd(x,y)=1\). It therefore gives one triple
for fixed coefficients and sign.

The usual ideal class count gives
\[
 |\mathrm{Cl}(K)|\ll_{q,\varepsilon}
          |D_K|^{1/2}X^\varepsilon
          \ll_{q,\varepsilon}
          \operatorname{rad}(AB)^{(q-1)/2} X^\varepsilon .
\]
For the first inequality use Minkowski's bound and
the elementary fixed-degree estimate for ideals of
prescribed norm. For the second, ramification lies
above \(qAB\). At primes \(p>q\) ramification is
tame and the discriminant exponent is at most
\(q-1\). For \(p\le q\) the discriminant exponent
is bounded in terms of \(q\): over the maximal
unramified local subfield the derivative of an
Eisenstein uniformizer polynomial gives different
exponent at most \(e-1+e v_p(e)\).
Sum these different exponents with local residue
degrees. This proves the class number estimate and
then (NR3), by summing the coefficient records and
height shells.
[Local different bound, Serre, Chapter III, §6](https://www.sas.rochester.edu/mth/sites/doug-ravenel/otherpapers/Serre-Local.pdf)

## Prime-power towers remove the class-number cost

Assume \(q=p^a\) and \(k=p^b\) for the same prime \(p\),
with \(a,b\) fixed and \(k>2/(1-h)\).
Then the count in (NR3) improves to
\[
                O_{q,k,h,\lambda,\varepsilon}
                          (X^{(f+g)\lambda+\varepsilon}).
 \tag{NR4}
\]
This includes irreducible moving coefficient fields, with
no restriction on whether they have real places.

First adjoin the \(p\)-th roots of unity. The extension
\(K(\mu_p)/K\) has degree dividing \(p-1\).
Extension of ideals is injective on \(p\)-primary class
torsion because the composition with the class norm is
raising to that degree. Over \(K_0=\mathbb Q(\mu_p)\),
put \(K_i=K_0(\theta^{p^{a-i}})\) for \(0\le i\le a\);
\(K_a=K(\mu_p)\). Irreducibility over \(\mathbb Q\)
gives \([\mathbb Q(\theta^{p^{a-i}}):\mathbb Q]=p^i\):
the degree is at most \(p^i\) and adjoining \(\theta\)
from that field has degree at most \(p^{a-i}\).
Since \([K_0:\mathbb Q]\) is prime to \(p\), it follows
that \([K_i:K_0]=p^i\). Each step \(K_i/K_{i-1}\)
is cyclic of degree \(p\), since it adjoins a \(p\)-th
root with all its conjugates in the extension.
The ramified finite primes all lie above \(qAB\), as
seen from the generator polynomials and the cyclotomic
base. The number of ramified finite primes per step
is \(O_q(1+\omega(AB))\).

For a cyclic degree-\(p\) step \(L'/K'\) with \(t\)
ramified finite primes and generator \(\sigma\), let
\(V=\mathrm{Cl}(L')[p]\), an \(\mathbb F_p\)-vector
space. The operator \(\sigma-1\) is nilpotent of index
at most \(p\). Its Jordan blocks, or the filtration
by its kernels, give \(\dim V\le p\dim V^\sigma\).
We bound the invariant classes explicitly.
For an invariant class represented by \(I\), write
\(I=(u)\sigma I\). The product of the conjugate
generators gives \(N_{L'/K'}u\in\mathcal O_{K'}^\times\).
Its class modulo norms of units is a homomorphism on
the invariant classes: changing the ideal representative
multiplies \(u\) by an element of the form \(v/\sigma v\),
and changing a generator multiplies by a unit.
The image size is bounded in terms of the degree,
since the norms contain all \(p\)-th powers of
\(\mathcal O_{K'}^\times\).
In the kernel adjust \(u\) by a unit to have norm one.
Cyclic Hilbert 90 writes it as \(v/\sigma v\), and
then \(I/(v)\) is invariant as an ideal. Invariant
ideals are extended \(K'\)-ideals times products
of ramified prime ideals with exponents smaller
than their ramification indices. Thus the subgroup
of classes represented by invariant ideals has
quotient of size at most \(p^t\) modulo the
extended \(\mathrm{Cl}(K')\). Its \(p\)-torsion
has size at most \(p^t|\mathrm{Cl}(K')[p]|\),
because the \(p\)-torsion size of an image of a
finite abelian group is at most that of the group.
Combining these bounds gives
\[
       |\mathrm{Cl}(L')[p]|
                  \le (O_q(1)\,p^t\,|\mathrm{Cl}(K')[p]|)^p .
\]
Starting with the fixed cyclotomic class group and
iterating at most \(a\) times gives
\(|\mathrm{Cl}(K_a)[p]|\le
\exp(O_q(1+\omega(AB)))=X^{o(1)}\) uniformly.
The same holds for \(p^b\)-torsion by successively
taking the kernels of multiplication by \(p\),
and then for \(K\) by the injection.
The fixed-coefficient ray count now proves (NR4).

For \(q=2\) the real norm case with \(A=B=1\) is
reducible and can also be handled directly. Write the
small opposite square as \(x^2\), the tall square
as \(y^2\), and the norm value \(T=y^2-x^2=C R^k\).
The factors \(y-x,y+x\) have common divisor at
most two. At fixed \(C\), split its prime powers
across the two factors, and record their valuations
at two, at subpower cost. Then
\(y\pm x=C_\pm z_\pm^k\) with fixed \(C_\pm\).
Since \(x/y=O(X^{-(1-h)/2})\), two ratios
\(z_+/z_-\) in this record have relative distance
\(O_k(X^{-(1-h)/2})\). Their cross determinant
is \(O(X^{1/k-(1-h)/2})\), since both root
heights are \(O(X^{1/(2k)})\). It is zero for
large \(X\), and the ratio recovers one primitive
pair \(x,y\).

## Consequence for minimal target repair

The short letter above is a condition on the original
quotient. The norm value \(T\) may be the primitive
additive target \(c\), or the other tall summand in
the real quadratic difference case. Impose the
canonical sector lower bound
\(\operatorname{rad}(c)\ge c^{(1/3-\xi)\lambda}\) in
either case. Then \(H_M(c)\le Y\) gives
\(c\le Y^{1/(M(1/3-\xi)\lambda)}\), and the count in
(NR4) yields repaired exponent
\[
                      \frac{f+g}{M(1/3-\xi)}+\varepsilon .
\]
It is strictly below \(1/M\) when \(f+g<1/3-\xi\).
The ideal powers and real closeness are from the
primitive equation. Raising the target valuations by
a common multiplier does not create these hypotheses.
[Minimal target-repair conversion](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)
[Quadratic angular packing and primitive norm data](ABC_MINIMAL_REPAIR_CYCLOTOMIC_SPLIT.md#angular-packing-of-primitive-square-kernels-and-near-power-targets)


## Full embedding determinant for binomial rays

For the binomial argument of (NR3)–(NR4) there is a stronger
comparison using all \(q\) embeddings. The shortness
condition on the target order can be weakened to
\[
                         k(1-h)(q-1)>2 .                       \tag{NR5}
\]
For either sign one may instead omit shortness and spend
a numerical chart exponent \(2/(kq(q-1))\). More precisely,
in the prime-power tower case the count without a height gap is
\[
            O_{q,k,\lambda,\varepsilon}
                (X^{(f+g)\lambda+2/(kq(q-1))+\varepsilon}).
 \tag{NR6}
\]
In the general fixed-order case the same argument uses
\((q+1)f/2+g\) in place of \(f+g\).
The input orders may depend on the fixed budget, with
implied constants then depending on those orders. This
asserts no simultaneous limit in order and height.

To prove this, write
\[
       t=Ax/(|\theta|y),\qquad
       \omega_i=\sigma_i(\theta)/|\theta|,
\]
so that the distinct \(\omega_i\) are fixed roots of \(1\)
or \(-1\). Consider a real interval of length at most
\(X^{-\delta}\) in a fixed compact set a fixed positive
distance from all the \(\omega_i\).
Keep the same fixed residual ideal, actual reference ideal,
unit coset and branch records as in the close-ray proof.
For each candidate \(\gamma_j\) its embedding vector has
coordinates
\[
         \sigma_i(\gamma_j)=a_i\,s_j f_i(t_j),
             \qquad f_i(t)=(t-\omega_i)^{1/k},
             \qquad s_j=y_j^{1/k}.
\]
Here the nonzero constants \(a_i\) are fixed within the
record and the \(f_i\) have fixed analytic branches on the
interval. These equations follow by taking the chosen root
of \(\sigma_i(\alpha/\alpha_0)/\sigma_i(u)\), using
\(\sigma_i(\alpha)=y|\theta|(t-\omega_i)\).
All derivatives through order \(q\) and the reciprocal
absolute values of \(f_i\) are bounded in terms of \(q,k\)
and the fixed distance.

Take \(q\) distinct ratios \(t_j\) in one interval.
Newton divided differences give
\[
 \left|\det_{i,j} f_i(t_j)\right|
             \ll_{q,k}\prod_{j<l}|t_j-t_l|
             \ll_{q,k} X^{-\delta q(q-1)/2}.
 \tag{NR7}
\]
Moreover this determinant is nonzero for all such
distinct nodes, once \(X\) is large. Indeed the divided
difference matrix converges uniformly, for nodes in
a shrinking interval about \(t_0\), to the Taylor
matrix with columns
\(\binom{1/k}{j-1}(t_0-\omega_i)^{1/k-j+1}\).
This last matrix has a nonzero determinant by
the Vandermonde in the distinct \((t_0-\omega_i)^{-1}\),
after removing nonzero row and column factors.
Its absolute value is bounded away from zero uniformly
on the indicated compact set. This argument also
handles nodes whose mutual separation is much smaller
than the interval width, since the estimates concern
divided differences.

By multiplying the absolute embedding products in
each column, the factored row and column scales satisfy
\[
 \left(\prod_i|a_i|\right)\left(\prod_j s_j\right)
   \asymp_{q,k}
      \left(\prod_j |N_{K/\mathbb Q}\gamma_j|\right)^{1/q}
   \ \le\ \frac{X^{1/k}}{N\mathfrak b_0},
\]
up to fixed dyadic constants. The \(\gamma_j\) are all
in the fractional ideal lattice \(\mathfrak b_0^{-1}\).
A nonzero full embedding determinant has absolute value
at least \(\sqrt{|D_K|}/N\mathfrak b_0\), by the
discriminant formula and the integral index in this
fractional ideal. By (NR7) the same determinant is
\(o(1/N\mathfrak b_0)\) if
\(\delta q(q-1)/2>1/k\), a contradiction. Hence at most
\(q-1\) distinct ratios occur per class, unit and branch
record in a chart. One ratio gives one primitive pair.

When \(Ax^q\le c^h\), the parameter \(t\) lies within
\(O(X^{-(1-h)/q})\) of zero, giving (NR5) and leaving
the coefficient count unchanged.
When shortness is omitted in the unsigned sum, interchange
the summands so \(t\in[0,1]\). The roots of \(-1\)
have a fixed positive distance from this interval.
Partition into intervals of length \(X^{-\delta}\)
for any \(\delta>2/(kq(q-1))\), pay
\(O(X^\delta)\) for the intervals and absorb the
arbitrarily small strict margin in \(\varepsilon\).
This gives (NR6).

In the difference case interchange the binomial terms so that
\(t\in[0,1)\). Away from the root \(1\) use the same compact
charts. Near \(1\) put \(t=1-\epsilon_0 r\) in a dyadic
interval with \(1\le r\le2\). The number of these scale
intervals is \(O_q(\log X)\), since the nonzero integral norm
value gives \(1-t\gg_q X^{-1}\). The root \(1\) and its
embedding are unique. Its normalized function is
\(r^{1/k}\) after extracting the factor \(\epsilon_0^{1/k}\)
into the constant row scale. All other rows are analytic
functions of \(\epsilon_0 r\), with their arguments bounded
away from zero. Their first \(q-1\) Taylor coefficients form
an invertible matrix by the same Vandermonde computation.
Take constant row combinations using its inverse and divide
the rows by the consecutive powers
\(1,\epsilon_0,\ldots,\epsilon_0^{q-2}\). The resulting
\(q-1\) rows are \(1,r,\ldots,r^{q-2}\) plus
\(O_{q,k}(\epsilon_0)\) in the \(C^{q-1}\) norm on
\([1,2]\). Alongside \(r^{1/k}\) these rows have a uniformly
nonzero divided-difference determinant for small
\(\epsilon_0\): its limit in the polynomial rows is a
nonzero constant times the divided difference of order
\(q-1\) of \(r^{1/k}\). The latter is bounded away from zero
by its derivative of the same order on \([1,2]\).
Use the compact charts for \(\epsilon_0\) bounded away
from zero. Thus in small \(\epsilon_0\) bins the analogue
of (NR7), using intervals of length \(X^{-\delta}\) in
\(r\), gains the further factor
\(\epsilon_0^{(q-1)(q-2)/2}\). The factors \(r^{1/k}\)
and the other untransformed normalized rows have absolute
values bounded away from zero. Hence the count of rows'
scales via \(N\gamma_j\) still applies before the row
combinations, including the extracted small first-row
factor. This proves (NR6) also for differences.

If the binomial norm is short while the two powered
terms are comparable and tall, write
\(|T|=X^{h_T+o(1)}\) in a dyadic height shell with
\(h_T<1\). Then \(\epsilon_0=X^{-(1-h_T)+o(1)}\)
and the ideal norm cap is \(N\mathfrak b\le X^{h_T/k+o(1)}\).
In degree \(q\ge3\) all the relative charts can therefore
be omitted when
\[
              h_T/k < (1-h_T)(q-1)(q-2)/2 .
\]
This is separate from the short-coordinate condition
(NR5), since here the binomial output is the short letter.

For example, with \(q\) a dyadic integer comparable to
\(\lambda^{-\alpha}\), \(\alpha>1/2\), and \(k=2\),
the archimedean chart cost is \(o(\lambda)\).
The third letter needs only a square representation with
a small squarefree coefficient radical, while the first
two have a common exact power order in this example.
If \(f+g<1/3-\xi\) by a fixed margin, this gives an
unsigned minimally repaired saving for all sufficiently
small fixed budgets at the lower target radical cap,
without requiring one summand to have smaller numerical
height. The primitive three-letter radical caps may all
be imposed. Even in this case, irreducibility of the
weighted binomial and the actual valuation conditions
are essential to the stated argument.
