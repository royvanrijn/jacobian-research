# Weighted quadratic norm descent for odd power targets

Two weighted square summands and an odd exact power target admit
a cyclotomic factor count when the ideal allocation has a
nontrivial principal power below the target order. The resulting
count uses the two summand radicals after paying for their parity
kernels. It only needs a bounded normalized cap on the target
radical. In the Gaussian case the principal power is automatic.

## Statement and primitive ideal allocation

Fix an odd integer \(q\ge3\) and fixed nonnegative caps
\(f\le \ell_{ab}\), \(\ell_c\). Count positive primitive triples
\(a+b=c\) in \(X/2<c\le X\) with
\[
 a=A x^2,\quad b=B y^2,\quad c=z^q,\quad z\text{ odd},\qquad
 AB\le X^{fL},\quad
 \operatorname{rad}(ab)\le X^{\ell_{ab}L},\quad
 \operatorname{rad}(c)\le X^{\ell_c L},
\]
where \(A,B\) are positive squarefree integers. Work at fixed
small budget \(L>0\), then take the height limit. Set \(D=AB\)
and \(K=\mathbb Q(\sqrt{-D})\). Primitivity makes \(D\) squarefree.
Use the full ring of integers of \(K\). With
\(\alpha=A x+y\sqrt{-D}\) and with \(\mathfrak a\) the product
of the ramified prime ideals over \(A\), one has
\[
      (\alpha)=\mathfrak a\mathfrak b^{\,q},\qquad
      \mathfrak a^2=(A),\qquad N\mathfrak b=z,                 \tag{WG1}
\]
for an integral ideal \(\mathfrak b\) coprime to \(2D\).
Indeed \(N\alpha=A z^q\). A prime over an odd \(p\nmid D\)
can divide both \(\alpha\) and \(\bar\alpha\) only if \(p\)
divides both \(x\) and \(y\), by the sum and difference of the
two elements. Thus a prime \(p\mid z\) is split and only one
of its two ideals divides \(\alpha\), at depth \(q v_p(z)\).
The primes dividing \(A\) occur at norm exponent one and are
ramified. This also handles \(2\mid A\), since then \(z\) is
odd and \(D\) is even squarefree. If \(2\nmid A\), the norm
is odd. No other primes occur.

Let \(j\) be the order of \([\mathfrak a\mathfrak b]\) in
\(\mathrm{Cl}(K)\), so \(j\mid q\) by (WG1). Write \(k=q/j\).
Assume \(k>1\).
For every fixed admissible \(k\), every \(\eta>0\) and all
sufficiently small fixed \(L>0\), the count in this sector is
\[
 \ll_{q,\ell_{ab},\ell_c,f,\eta,L,\varepsilon}
 X^{\left(f+\frac{\ell_{ab}-f}{2\tau(k)}+\eta\right)L+\varepsilon}.
                                                               \tag{WG2}
\]
One may sum over all \(k>1\) and use the worst case
\(\tau(k)\ge2\). The assertion includes the principal-allocation
case \(j=1\), for which \(k=q\). In particular it applies if
the \(q\)-torsion of the class group is trivial.

## Principal power and factor budgets

Since \(j\) is odd, \(\mathfrak a^j=(A^{(j-1)/2})\mathfrak a\).
The ideal \(\mathfrak a\mathfrak b^j\) is principal. Choose a
generator \(\gamma\) and adjust it by a unit so that
\[
             \gamma^k=A^{(k-1)/2}\alpha,\qquad
             N\gamma=A z^j .
\]
The equality of ideals follows from \(jk=q\) and (WG1).
The norm identity follows directly from \(N\mathfrak a=A\).
Raising to odd \(k\) is onto on the unit group except
possibly for \(D=3\) and \(3\mid k\). In that case
\(\mathcal O_K/(2)=\mathbb F_4\). The element
\(\alpha=A x+y\sqrt{-3}\) reduces to \(1\), because it
reduces into \(\mathbb F_2\) and has odd norm.
Also \(\gamma\) has odd norm and \(\gamma^k\) reduces to
\(1\). The unit comparing \(\gamma^k\) with
\(A^{(k-1)/2}\alpha\) therefore reduces to \(1\), so is
\(\pm1\). These units are \(k\)-th powers.
Write
\(\gamma=(A r+s\sqrt{-D})/2\) with integers \(r,s\).
The trace of \(\gamma\) is divisible by \(A\), since it is
in each ramified prime over \(A\) and \(A\) is squarefree.
The imaginary coordinate \(s\) is integral in the full ring
of integers of a quadratic field. Both \(r,s\) are nonzero,
since an odd power of a real or pure imaginary \(\gamma\)
could not give both nonzero original coordinates. We obtain
\[
                       A r^2+B s^2=4z^j.                      \tag{WG3}
\]

The gcd of the three terms in (WG3) divides \(4\).
An odd prime dividing both terms on the left cannot divide
\(D\), since \(z\) is coprime to \(D\). It would divide both
\(r,s\), so \(\gamma\) and hence \(\gamma^k\) would be divisible
by that rational prime in \(\mathcal O_K\). Since the prime
does not divide \(2D\), the equality for \(\gamma^k\) then forces
it to divide both \(x,y\). This contradicts primitivity.
The exponent at \(2\) is at most two since \(z\) is odd.
A derived primitive triple determines \((r,s,z)\) up to the
bounded choices of common divisor and signs at fixed \(A,B,j\).
Its complete height is comparable to \(X^{1/k}\) on a dyadic
original target shell \(c\asymp X\).

Set \(U=\gamma\), \(V=\bar\gamma\).
The two expressions \(U^k+V^k\) and \(U^k-V^k\) separate
the two coordinates of \(\alpha\), up to fixed factors
supported on \(2D\) and the factor \(\sqrt{-D}\) on the
imaginary side. As \(k\) is odd, the two linear factors
are \(U+V=A r\) and \(U-V=s\sqrt{-D}\). Together the
expressions have \(m=2(\tau(k)-1)\) further factors
\(\Phi_e(U,V)\), for \(e\mid 2k\), \(e\ge3\). These
values are rational integers because the forms are symmetric
and \(U,V\) are conjugate algebraic integers.
They are nonzero by the two nonzero original coordinates.

Let \(P=2Dq\). At primes \(p\nmid P\) the evaluated factor
radicals, including the two linear factors with radicals
from \(r,s\), are disjoint. If both \(U,V\) vanish modulo a
prime ideal over such \(p\), then \(p\mid r,s\), already
excluded. If one vanishes, none of the cyclotomic factors
vanishes. Otherwise the ratio \(U/V\) is a unit and is a
root of at most one of the factors of \(T^{2k}-1\),
which is separable at \(p\).
Each retained prime belongs to the radical of the corresponding
original \(x\) or \(y\), by the exact factor identities.
Consequently if the outside-\(P\) cap exponents in base
\(X^L\) are \(\beta_r,\beta_s,\beta_e\), their sum is at
most \(\ell_{ab}-f_0\) up to arbitrarily small bin errors,
where \(f_0=\log_X(AB)/L\).

## Fixed values and the two smaller dictionaries

Fix \(A,B\), the signature and dyadic caps. Each nonlinear
factor gives, after multiplying by a bounded power of two,
an integral homogeneous binary form \(F_e(r,s)\) of degree
at least two. It is squarefree as a form since the linear
change from \((r,s)\) to \((U,V)\) is invertible over \(K\).
Its coefficients have polynomial height in \(X\), uniformly
at fixed \(q\). At fixed nonzero value it gives \(X^{o_X(1)}\)
parameter pairs uniformly. If it has an irreducible factor
of degree at least two, enumerate that factor's value by
the divisor bound and apply the bounded-degree ideal norm
and unit argument. If every factor is linear, two independent
linear values recover the pair, again at divisor cost.
All relevant coefficient heights remain polynomial by the
fixed-degree factor bound. Fixing the radical outside \(P\)
then gives \(X^{o_X(1)}\) supported values on the union of
that radical and \(P\). Thus a conditional normalized count
exponent \(d_0\) satisfies \(d_0\le\beta_e\) for each of the
\(m\) factors, up to the small exponent allowances.
[Uniform norm count for moving binary forms](ABC_MINIMAL_REPAIR_CYCLOTOMIC_SPLIT.md#direct-target-bounds-for-reducible-weighted-powers)

For the derived primitive triple in (WG3), apply the
two-cheapest-letter dictionary bound (GD3).
Split each exact valuation root into its part on \(P\) and its
outside part. The parts on \(P\) have uniformly subpower
dictionaries, by fixed-support enumeration. Keep them as
cost-zero singletons in the anatomic box rather than moving
their large powers into the initial coefficients.
This splitting leaves the numerical root cap sum unchanged,
up to the usual height allowances. The outside costs for the
two summands are at most \(\beta_r+\beta_s\) together.
The raw radical cap on the entire derived triple is
\(X^{O_{\ell_{ab},\ell_c}(L)}\), since the newly introduced
\(\operatorname{rad}(D)\) and the retained primes all come from
the original budget, apart from fixed primes. At derived
height \(X^{1/k}\) use a fixed enlarged relative budget
\(O_{q,\ell_{ab},\ell_c}(L)\). Its height exponent cancels
when expressing the dictionary costs in base \(X^L\).
The cap on \(z\) is used here. The (GD3) cost is no larger than
half the sum of the two specified summand costs, regardless
of which two costs are actually smallest. Thus
\[
     2d_0\le\beta_r+\beta_s
       \le \ell_{ab}-f_0-m d_0,\qquad
     d_0\le \frac{\ell_{ab}-f_0}{2\tau(k)} .
\]
For the dyadic and radical bins, fix all outer allowances
before the height limit and all mesh widths independently
of \(X\). The root decomposition at each derived budget has
fixed truncation with arbitrarily small normalized loss.
[Two-cheapest-letter dictionary bound (GD3)](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md#group-dictionary-count)
· [Uniform fixed-support enumeration](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#a-radical-count-contradiction-to-the-claimed-smooth-number-lower-bound)

There are at most \(X^{f_0L+o_X(1)}\) kernel choices in a
small bin for \(f_0\), with the bin-width error in the exponent,
since the squarefree product and its prime allocation determine
\(A,B\). The bound \(f_0+(\ell_{ab}-f_0)/(2\tau(k))\) is increasing
in \(f_0\). Sum over bins and absorb the losses in \(\eta\)
to prove (WG2).

## Transfer and residual ideal classes

For cumulative bounds the absolute caps of the shell theorem
can be imposed relative to \(c\) instead of \(X\), then summed
over dyadic shells. At the balanced summand cap \(\ell_{ab}=2/3\), (WG2) gives
a strict saving for every \(k\ge3\) when \(f<2/9\),
with any fixed margin. In the central least-target-repair
band with joint summand cap \(2/3+\xi\) and target radical
at least \(c^{(1/3-\xi)L}\), the repaired exponent is at most
\[
 \frac{f+(2/3+\xi-f)/(2\tau(k))+\eta}{M(1/3-\xi)}
    +\varepsilon ,
\]
if \(AB\le c^{fL}\). The kernel and ideal decompositions
refer to the primitive quotient, and its least repair
has height \(H_M(c)\ge\operatorname{rad}(c)^M\).
In the exact square case \(A=B=1\), \(\mathrm{Cl}(K)\) is
trivial and \(z\) is automatically odd, since two odd
squares sum to an integer of 2-adic valuation one.
At \(q=3\) the normalized original exponent is \(1/6\)
under balanced caps, and the repaired coefficient relative
to \(1/M\) approaches \(1/2\).
[One-target minimal repair](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)

An ideal class of exact order \(q\) only supplies \(k=1\)
through this descent, leaving no nonlinear factors to pay.
For even \(z\) the factorization at \(2\) may supply a
non-real twist. The current argument does not control those
fibres or show that the relevant ideal orders are generically
proper divisors of \(q\).

## Projective recovery with arbitrary ideal classes

There is a separate count when the parity kernel of a summand
already carries enough numerical weight. It allows arbitrary
ideal classes and permits a small defect on the target power.
Fix any integer \(q\ge3\), \(f,g\ge0\) and \(\delta>0\).
In a primitive shell \(X/2<c\le X\), write
\[
       a=A x^2,\quad b=B y^2,\quad c=C z^q,
       \qquad AB\le X^{fL},\quad
       \operatorname{rad}(C)\le X^{gL},
\]
with squarefree \(A,B\), and take \(C\) to be the canonical
\(q\)-power defect, so \(v_p(C)<q\).
Let \(h=\min(\log_X a,\log_X b)\). The initially recorded
prime supports on \(a,b,c\) are those of \(A,B,C\), respectively.
Optionally add primes on the original letters with combined
squarefree product at most \(X^{s_0L}\), for a fixed \(s_0\ge0\),
and write \(w_a,w_b,w_c\) for the logarithmic numerical weights
in base \(X\) of all full prime powers on the resulting supports.
Assume
\[
               w_a+w_b+(2/q)w_c\ \ge\ h-1+2/q+\delta .
\]
For every \(\eta>0\) and every sufficiently small fixed \(L>0\),
the count is
\[
        O_{q,f,g,s_0,\delta,\eta,L,\varepsilon}
             \left(X^{(\chi_q f+g+s_0+\eta)L+\varepsilon}\right),
                                                               \tag{WG4}
\]
where \(\chi_q\) may be any established exponent such that
\(\sum_{D\le T,\ D\ \mathrm{squarefree}}|\mathrm{Cl}(\mathbb Q(\sqrt{-D}))[q]|
\ll_{q,\varepsilon} T^{\chi_q+\varepsilon}\).
One may always use \(\chi_q=3/2\). In particular
\(\chi_3=1\) by Davenport–Heilbronn, and the same bound
holds for any fixed power of two and for three times a
fixed power of two by quadratic genus theory. For a fixed
prime \(q\ge5\), Heath-Brown and Pierce give
\(\chi_q=3/2-3/(2q+2)\).
This bound does not require an additional cap on
\(\operatorname{rad}(abc)\). The sector gives a strict saving
from exponent \(L/3\) if \(\chi_q f+g+s_0<1/3\)
by a fixed margin.
[Davenport–Heilbronn, Theorem 3](https://www.math.utoronto.ca/~ila/DH2.pdf)
· [Heath-Brown–Pierce, Theorem 1.1](https://arxiv.org/pdf/1409.3177)

To prove it, put \(D=AB\), \(K=\mathbb Q(\sqrt{-D})\), and
\(\alpha=A x+y\sqrt{-D}\). As in (WG1), the prime ideals of
\(\alpha\) over the odd primes of \(c\) lie on one split branch
for each prime, and the ramified part over \(A\) has norm \(A\).
At odd primes use the exact powers of \(C\) and \(z\) on their
chosen branches. For even \(c\), the original summands are odd,
so \(A,B,x,y\) are all odd. If \(D\equiv1\bmod4\), then
\(v_2(c)=1\) and \(2\) ramifies in \(K\). If \(D\equiv3\bmod8\),
then \(v_2(c)=2\) and \(2\) is inert. In both cases put the
entire contribution above \(2\) into the defect ideal.
If \(D\equiv7\bmod8\), then \(2\) splits and \(v_2(c)\ge3\).
The residues of \(\alpha/2=(A x-y)/2+
y(1+\sqrt{-D})/2\) at the two primes above \(2\) are different.
Its two valuations therefore give valuations of \(\alpha\)
equal to \(1\) and \(v_2(c)-1\). If \(e=v_2(C)\), put a
valuation \(v_2(z)\) in \(\mathfrak b\) on the larger branch.
Put valuations \(1,e-1\) in the defect ideal across the two
branches, in that order. This may give a fractional defect
ideal with denominator at most \(2\), and its norm at \(2\)
remains \(2^e\).

The exact representation for either parity is
\[
                   (\alpha)=\mathfrak a\,\mathfrak c\,
                                  \mathfrak b^q,
              \qquad N\mathfrak c=C,\quad N\mathfrak b=z .
\]
Fix \(A,B,C\) and the split branch for each odd prime dividing
\(C\). At \(2\) there are at most two allocation choices,
including when a fractional defect remains after removing
the \(q\)-th power. This has cost \(X^{(f+g)L+o_X(1)}\). Indeed the kernel product
and its prime allocation determine \(A,B\). Then
\(C\le X^{(q-1)gL}\) and the already chosen radical has
uniformly subpower many supported values with subpower many
prime-branch choices.

Fix also the ideal class of \(\mathfrak b\). The feasible
classes solve \( [\mathfrak b]^q=[\mathfrak a\mathfrak c]^{-1}\),
so they form a coset of the \(q\)-torsion if nonempty.
In its inverse class choose a
fixed integral ideal \(\mathfrak t\) of norm \(m\ll\sqrt D\),
using the same Minkowski bound. There is a generator \(\gamma\)
of \(\mathfrak b\mathfrak t\), and we have
\[
          \alpha=\beta\gamma^q,\qquad
          (\beta)=\mathfrak a\mathfrak c\mathfrak t^{-q},
                    \qquad N\gamma=m z .
\]
For a nonempty record, fix a generator of the principal
fractional ideal on the right as \(\beta\), up to the finitely
many imaginary quadratic units. Its norm is \(AC/m^q\).
As \(m\in\mathfrak t\), the element \(2m^q\beta\) is integral.
Writing \(\beta\) in the standard quadratic basis, its absolute
complex modulus and its denominator bound its two rational
coordinate heights by \((DC)^{O_q(1)}\). No principal generator
of \(\mathfrak b\) is asserted.

Write \(\gamma=(r+s\sqrt{-D})/2\). The real and imaginary
coordinate forms for \(x,y\), respectively, are the rational
homogeneous forms
\[
 x(r,s)=\frac{\beta\gamma^q+\bar\beta\bar\gamma^q}{2A},
 \qquad
 y(r,s)=\frac{\beta\gamma^q-\bar\beta\bar\gamma^q}{2\sqrt{-D}},
 \qquad Z(r,s)=\frac{r^2+D s^2}{4m}.
\]
The three forms \(A x(r,s)^2,B y(r,s)^2,C Z(r,s)^q\)
have equal degree \(2q\), satisfy the additive identity by
the norm formula, and give the primitive original triple at
the projective parameter \([r:s]\). They have no common
geometric zero, because \(\gamma,\bar\gamma\) are independent
linear forms over \(\mathbb C\) and \(\beta\ne0\).
After clearing denominators by a common \(2q\)-th power we obtain
integral forms of height \((DC)^{O_q(1)}\). Their two summand
polynomials are squares of degree-\(q\) forms times \(A,B\).
Each underlying degree-\(q\) form has \(q\) distinct geometric
projective roots, since it is an invertible homogeneous
change of a nondegenerate binomial of degree \(q\). The
target polynomial is \(C\) times the \(q\)-th power of
a quadratic with two distinct roots.

Take a reduced primitive integral parameter \((r_0,s_0)\)
and \(H=\max(|r_0|,|s_0|)\). The resultant identities give a
\((DC)^{O_q(1)}\) upper bound for the gcd of the three integral
values. They also give upper and lower bounds for the maximum
of the values in terms of \(H^{2q}\) with the same polynomial
coefficient losses. Hence
\[
                    H=X^{1/(2q)+O_q((f+g)L)+o_X(1)}.
\]
Partition the shorter summand into numerical height bins
of any fixed small width, denoting their common exponent by
\(h\) up to the bin error. In an affine chart whose denominator
has magnitude \(H\), the underlying degree-\(q\) summand form
divided by \(H^q\) has size at most
\(X^{-(1-h)/2+O_q((f+g)L)+o_X(1)}\).
If the gap is larger than the coefficient losses, root
separation puts the chart coordinate in a real interval of
that width, with the same \(O_q((f+g)L)\) allowance in its
exponent. For the shallow gap the trivial chart suffices
with an enlarged error of that order. After recording one of
the boundedly many charts and roots, cross determinants of
two reduced parameters are at most
\[
         X^{1/q-(1-h)/2+O_q((f+g)L)+o_X(1)}
\]
with the numerical bin error also allowed.
[Projective height, root isolation and residue records](ABC_RATIONAL_PROJECTIVE_CARRIERS.md#height-local-branch-weights-and-the-real-correction)

Now record all exact complete prime powers on the chosen
supports. Enumerate the extra primes and their allocation at
cost \(X^{s_0L+o_X(1)}\). For all specified supports the exact
powers cost \(X^{o_X(1)}\) by fixed-support enumeration.
If \(p^e\Vert T\) is recorded on a summand \(T\), the
underlying degree-\(q\) form has valuation at least
\((e-v_p(D_T))/2\), with \(D_T=A\) or \(B\). If it is recorded
on the target, the underlying quadratic has valuation at least
\((e-v_p(C))/q\). In fact each integral evaluation at the
reduced parameter is a common positive integer multiple of
the primitive row.

For each chart of the reduced parameter use the derivative
resultant of the underlying form, and record a unit denominator
chart and a Hensel root branch. Fix a bad-prime integer \(J\)
of height \((DC)^{O_q(1)}\) containing \(2DC\), the leading
coefficients and the nonzero derivative resultants in both charts.
At valuation exceeding \(2v_p(J)\) the derivative has
valuation at most \(v_p(J)\) and Hensel gives precision
with loss at most \(v_p(J)\) from the value's valuation.
At lower valuations simply discard the precision.
The total loss from bad primes, coefficient factors and
integer rounding is \(O_q(\log J+s_0 L\log X)\), since
all recorded primes divide \(D C\) or the extra prime product.
There are boundedly many chart and root records per prime,
so their total cost is subpower. The combined modulus
dividing every cross determinant of the record has size
\[
    X^{w_a/2+w_b/2+w_c/q-O_q((f+g+s_0)L)-o_X(1)}.
\]
Here different complete letters use disjoint primes by primitivity.

By the assumed strict numerical margin, the modulus
exceeds the archimedean determinant cap for sufficiently
small fixed \(L\) and sufficiently fine fixed numerical bins.
Thus each record has one reduced projective parameter and
one primitive triple. Sum over \(D\le X^{fL}\) and its
subpower many allocations weighted by \(|\mathrm{Cl}(K)[q]|\).
Counting the remaining defect and support records gives (WG4).
The estimate with \(\chi_q=3/2\) follows from Minkowski's
bound and uniform ideal-divisor counting. For a fixed
\(2\)-power order the genus-theory \(2\)-primary torsion
factor is subpower, so it does not change the \(3\)-torsion
average when the order is \(3\) times a \(2\)-power.
[Short supported records and fixed-support enumeration](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#a-radical-count-contradiction-to-the-claimed-smooth-number-lower-bound)

In particular at \(q=3\) this uses the average class
number input for all kernel choices and needs no restriction
on the individual class of \(\mathfrak b\).
The theorem still requires a fixed numerical margin
from the recorded powers.

### Automatic purchases in a diffuse summand pair

Suppose in addition that
\(\operatorname{rad}(ab)\le X^{(2/3+\xi)L}\) and that
every prime of \(ab\) outside \(AB\) has
\(\log_X p\le \rho L\). For small enough fixed \(\rho,\xi>0\),
(WG4) implies a strict uniform saving below \(L/3\) whenever
\[
                    (q\chi_q-1)f+qg<(q-2)/3                 \tag{WG5}
\]
by a fixed margin. In particular the condition is
\(2f+3g<1/3\) for cube targets and their cubic defects.

Indeed bin the actual kernel mass near \(f_0=\log_X(AB)/L\).
Spend at most \(s_0=1/3-\chi_q f_0-g-\epsilon_0\) in that bin
on the additional prime supports of \(a,b\). Choose the
bin allowance and \(\epsilon_0>0\) small compared with the
fixed margin in (WG5). Greedily select primes in decreasing
order of \(v_p(ab)\), stopping before the budget is
exceeded. The missed last prime costs at most \(\rho\)
in normalized radical. Unless all primes fit, this selects
at least the fraction
\[
                  t=\frac{s_0-\rho}{2/3+\xi-f_0}
\]
of the total numerical weight outside \(AB\). This follows
by comparing the average numerical weight per unit radical
of the initial segment with the average over the whole list.
If all primes fit, take \(t=1\). Including the already paid
kernel powers gives \(w_a+w_b\ge t(\log_X a+\log_X b)\).
The displayed budget ratio is decreasing in \(f_0\) in this
range. By (WG5) it exceeds \(1/q\) by a fixed margin at all
\(f_0\le f\). Since the taller summand has exponent
\(1-o_X(1)\) and
\[
                 \frac{h-1+2/q}{1+h}\le\frac1q ,
\]
the numerical threshold in (WG4) holds with fixed slack
at large height. The total paid exponent is less than
\(1/3\) by a fixed margin after the bin and exponent losses.
This argument selects only actual prime powers of the original
summands and counts each primitive quotient once if the
least-repair cutoff is imposed as an additional constraint.
