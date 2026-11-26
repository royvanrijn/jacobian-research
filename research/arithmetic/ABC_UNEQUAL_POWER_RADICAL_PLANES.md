# Unequal powers on a plane of radical coordinates

When the radicals of three actual letters satisfy a plane relation,
two different pure valuation orders give an exact integer quotient.
After fixing that quotient, the original additive equation is
a polynomial in one root. Removing the factors that account for
the third radical leaves a nonconstant factor on fixed prime
support. This gives a subpower count even for adjacent small orders,
with the third letter's actual valuation spectrum arbitrary.

## The exact quotient from the third support

Fix \(H_0>0\), fixed integers \(p>q\ge2\), and a prime pool
\(\mathcal P\) whose prime product is at most \(X^{H_0}\).
All constants below are uniform in the data of this polynomial
height. Fix nonzero integers \(a,b,\ell_1,\ell_2,\ell_3\) of
absolute value at most \(X^{H_0}\). Consider positive integers
\(x,y,z\le X^{H_0}\) such that
\[
       \ell_1x+\ell_2y+\ell_3z=0,\qquad
       W=a x^p+b y^q\ne0,\qquad |W|\le X^{H_0},
       \qquad
       \operatorname{rad}_{\overline{\mathcal P}}(W)=z .
\]
Assume \(\gcd(y,z)=1\). In particular these conditions hold for
an original primitive row with two weighted exact powers and
with \(z\) the remaining complete letter's radical outside the
recorded pool. The additive signs can be included in \(a,b\).
The first two root variables may be required to be squarefree
when they arise as residual radical coordinates.

Enlarge \(\mathcal P\) by the primes of \(ab\ell_1\ell_2\).
If this removes some primes from \(z\), record their product
at height-subpower cost and absorb it in \(\ell_3\). Thus \(z\)
again denotes the exact radical outside the enlarged pool.
All statements after this operation are uniform on the new
record. Put \(d=p-q\). Let \(L(t)\) be the primitive integral
polynomial proportional over \(\mathbb Q\) to
\[
                   a(-\ell_2)^p t^d+b\ell_1^p .
\]
Both coefficients are nonzero. Modulo \(z\) we have
\(\ell_1x=-\ell_2y\). Multiplying the additive equation by
\(\ell_1^p\) then gives
\[
        y^q\bigl(a(-\ell_2)^p y^d+b\ell_1^p\bigr)=0
                                \pmod z.
\]
The removed coefficient content has all its primes in the
enlarged pool, and \(y\) is a unit modulo \(z\). It follows
that
\[
                             k=L(y)/z\in\mathbb Z .           \tag{UP1}
\]
This integer is a quotient by the actual third radical.
There is no bound on the valuations of \(W\) on that support.
[Polynomial-height supported-value counts](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

## Uniform recovery from one quotient

**UP2.** For each fixed \(k\), there are
\(O_{H_0,p,q,\varepsilon}(X^\varepsilon)\) triples in
the stated record. The constants are uniform in \(k\),
the coefficients and the prescribed pool. If a solution
occurs then \(k\) has polynomial height, with exponent
bounded in terms of \(H_0,p,q\). The same uniformity holds
if the pool and coefficients have a larger fixed
polynomial height after the initial stripping.

First let \(k\ne0\), put \(t=y\), and write
\[
  z=\frac{L(t)}k,\qquad
  x=f(t)=-\frac{\ell_2}{\ell_1}t
         -\frac{\ell_3}{k\ell_1}L(t),\qquad
  F(t)=a f(t)^p+b t^q.
\]
At each complex root \(w\) of \(L\), \(w\ne0\) and
\(f(w)=-(\ell_2/\ell_1)w\ne0\). The definition of \(L\)
gives \(F(w)=0\). The roots of \(L\) are simple, since it
is a nonzero binomial with constant term nonzero in
characteristic zero.

The total multiplicity of the roots of \(L\) in \(F\)
is at most
\[
                         d+\deg\bigl(pt f'(t)-qf(t)\bigr).
\]
In fact near each root \(w\) the quotient \(F(t)/t^q\)
is regular. Its derivative is
\[
          a f(t)^{p-1}t^{-q-1}\bigl(pt f'(t)-qf(t)\bigr).
\]
The visible factors in front of the bracket are nonzero
at \(w\), so each unit of multiplicity beyond the first
is charged to a zero of the bracket. The bracket is
nonzero since \(f(0)\ne0\).

If \(d\ge2\), then \(\deg f=d\) and \(\deg F=pd>2d\).
The multiplicities accounted for by \(L\) are at most
\(2d\), so \(F\) has another factor. For \(d=1\), if
\(f\) is nonconstant, \(\deg F=p\ge3\) and the same
bound is two. If \(f\) is constant, the bracket has
degree zero and \(\deg F=q\ge2\), so again an extra
factor remains.

Clear the fixed rational denominators in \(F\), and
remove from the resulting integer polynomial all
primitive irreducible factors associated with a root
of \(L\), to their full multiplicity. The remaining
primitive factor product \(G\) is nonconstant and
coprime to \(L\). Its coefficients and its
nonzero resultant with \(L\) have polynomial height
uniformly at the fixed degrees. For the coefficient
claim, the norm of a product of two real polynomials
of bounded degrees is bounded below by a positive
constant times the product of their coefficient norms.
This follows by restricting both coefficient norms to
one and using compactness and lack of zero divisors.
Apply it to the integral factors after removing content.
The resultant height follows from its determinant.

At any actual \(t\), a prime dividing \(G(t)\) outside
the prescribed pool and the denominator support divides
the nonzero \(W=F(t)\). Hence it divides \(z\), and thus
\(L(t)\). Such a prime divides the fixed resultant.
Consequently the nonzero integer \(G(t)\) is supported
on a fixed polynomial-height pool after the resultant
is included. It is of polynomial height. Uniform
fixed-support enumeration gives \(X^{o_X(1)}\) possible
values, and a nonconstant polynomial has at most
its degree many preimages of each value. This recovers
\(y\), then \(z,x\), proving the count in the nonzero
quotient branch.

For \(k=0\), first fix \(y\) among the at most \(d\)
integer roots of \(L\). In the plane, \(x\) is an
affine function of \(z\) with nonzero slope and
nonzero constant term. The polynomial
\(a x(z)^p+b y^q\) vanishes at \(z=0\) and has nonzero
derivative there. After factoring out \(z\) it has
a nonconstant factor of degree \(p-1\) whose value
at zero is nonzero. Extra prime divisors outside the
original pool of this factor would have to divide
\(z\), and hence its fixed constant term. Supported-value
enumeration as above counts \(z\) at subpower cost.
This proves (UP2) in both branches without requiring
any distribution of valuations on the third letter.

## Count with a height or radical price on the quotient

Let \(n=X^L\) at a fixed \(L\). Suppose the coefficients
of the primitive \(L\) have height at most \(n^\theta\),
and \(y\le n^{r_y}\), \(z\ge n^{r_z}\). Equation (UP1)
gives
\[
 |k|\le 2 n^{\theta+d r_y-r_z}.
\]
Thus the row count in this record is
\[
                 O\!\left(
                   n^{(\theta+d r_y-r_z)_+}X^\varepsilon
                  \right).                                   \tag{UP3}
\]
Alternatively, if the nonzero \(k\)'s have outside
radical at most \(n^b\) relative to a further prescribed
polynomial-height pool, the bound is \(O(n^b X^\varepsilon)\).
Enumerate the outside radical, then its complete value
and the part supported on the pool at subpower cost.
The zero case carries no quotient cost.

One may sum over a dictionary of original coefficient
parts, plane normals and pools by paying its cardinality.
The supported parts stripped from \(z\) cost only an
inner height-subpower factor. In the height bound this
uses the lower cap \(r_z\) for the radical after that
stripping and the primitive height of the resulting
quotient polynomial. This distinction is relevant
when many primes on the third letter lie in the
original record.

In particular, for adjacent fixed orders \(p=q+1\),
a lower third-radical exponent at least the second
root exponent leaves only the cost \(\theta\).
A dictionary of size \(n^s\) with \(s+\theta<1/3\)
by a fixed margin gives a strict small-budget sector
bound if the stated original caps and plane relation
hold in that dictionary. The theorem works for orders
\(2,3\), where the high-order product hypothesis of
(TG1) is unavailable. For larger fixed orders the
product gap \(d r_y-r_z\) specifies the sector for the
height version. The quotient version instead leaves
only the actual radical price of the affine numerator
divided by the third support. Neither version presumes
that such planes or cheap quotients occur universally.
[Two high valuation grids without the radical-plane relation (TG1)](ABC_TWO_VALUATION_GRIDS.md#count-and-range)

## A deeper quotient forced by the original third valuations

For adjacent orders \(p=q+1\), primes of order at least \(J\)
on the original third letter give an exact quotient of depth
\(J\). Fix \(J\ge2\). After the coefficient-pool stripping in
(UP1), record the product of primes on \(z\) whose valuation
on \(W\) is below \(J\). Move these primes into the pool and
absorb their product into \(\ell_3\). The remaining \(z\)
is squarefree and every prime on it has valuation at least
\(J\) on \(W\). This recording costs at most \(n^s\) if the
removed prime product is at most \(n^s\).

Put
\[
 P(Y,Z)=a(-\ell_2Y-\ell_3 Z)^p+b\ell_1^pY^q,
 \qquad  P(y,z)=\ell_1^p W.
\]
At \(Z=0\) the polynomial is
\(Y^q(C_0Y+D_0)\), with \(C_0=a(-\ell_2)^p\),
\(D_0=b\ell_1^p\), and simple nonzero root
\(y_0=-D_0/C_0\). Let \(Q_J(Z)\) be the polynomial
of degree less than \(J\) with constant coefficient
\(y_0\) obtained by solving \(P(Q_J(Z),Z)=0\bmod Z^J\)
successively in the coefficients. The recursion divides
at each step by \(C_0y_0^q\), which is nonzero. The
coefficients have polynomial height for fixed \(p,q,J\),
with denominators supported on \(ab\ell_1\ell_2\).
Choose a positive common denominator \(D_J\) with
this support, and put \(Q_J^*=D_JQ_J\) in
\(\mathbb Z[Z]\).

**UP4.** Each actual row yields the integer
\[
             k_J=\frac{D_Jy-Q_J^*(z)}{z^J}.                  \tag{UP4}
\]
For each fixed \(k_J\) the original row count is
\(O_{H_0,p,q,J,\varepsilon}(X^\varepsilon)\) uniformly
at polynomial height in the recorded data. If
\(D_J,\|Q_J^*\|_\infty\le n^\Theta\),
\(y\le n^{r_y}\), and \(z\ge n^{r_z}\), the count
including the record of the removed primes is
\[
           O\!\left(n^{s+\max\{0,\Theta+r_y-Jr_z,\Theta-r_z\}}
                     X^\varepsilon\right).                 \tag{UP5}
\]
There is also a bound \(O(n^{s+b}X^\varepsilon)\) when
the nonzero \(k_J\)'s have outside radical at most
\(n^b\) relative to an additional prescribed
polynomial-height prime pool.

Indeed, for a prime \(v\mid z\), both \(y\) and \(y_0\)
are units at \(v\), and \(P(y,0)=0\bmod v\) implies
\(y=y_0\bmod v\). The coefficients of \(Q_J\) are
\(v\)-integral and \(D_J\) is a unit there. Both
\(P(y,z)\) and \(P(Q_J(z),z)\) vanish modulo \(v^J\).
Their difference factors as \(y-Q_J(z)\) times a
\(v\)-adic unit, since the other factor reduces to
\(P_Y(y_0,0)=C_0y_0^q\). Thus \(v^J\mid
D_Jy-Q_J^*(z)\). Coprimality of the prime powers gives
(UP4).

Now fix \(k_J\). Write \(t=z\), set
\[
 v(t)=Q_J(t)+\frac{k_J}{D_J}t^J,\qquad
 u(t)=-\frac{\ell_2v(t)+\ell_3t}{\ell_1},\qquad
 W_0(t)=a u(t)^p+b v(t)^q.
\]
The constant terms of \(u,v\) are nonzero. The
polynomial \(W_0\) vanishes at zero. Near zero,
differentiating \(W_0(t)/v(t)^q\) shows that the
multiplicity at zero is at most one plus the degree of
\[
                      p u'(t)v(t)-q u(t)v'(t).
\]
If \(\deg v=e\ge2\), then \(\deg u=e\); more
generally, whenever both have common positive degree
\(e\), the derivative numerator has nonzero leading
coefficient of degree \(2e-1\), and \(\deg W_0=pe>2e\).
If one is constant, the other has degree one, since
\(\ell_3\ne0\) and \(u,v\) obey the plane relation.
The derivative numerator is then a nonzero constant
and \(\deg W_0\ge2\). These cases exhaust the degrees.
After removing the factor at zero to full multiplicity,
a nonconstant factor with nonzero value at zero
remains. Clear denominators and apply the same
supported-value argument as in (UP2), now with the
cleared denominator and the constant term of that factor
in the pool. A prime
outside the pool on this factor would divide the true
\(W\), hence \(z=t\), hence the constant term. This
is impossible. Uniformity follows from fixed degrees
and polynomial coefficient height. For any occurring
quotient the displayed formula (UP4) bounds its
height polynomially, so the new coefficients have
the required height.

For \(z\ge1\), the magnitude of (UP4) is at most a
constant depending on \(J\) times
\[
               n^\Theta\bigl(y/z^J+1/z\bigr).
\]
Counting the integer quotients gives (UP5). For the
radical version use uniform supported enumeration
on the quotient and add the one zero choice.

There is a direct coefficient bound with the cancellations in the
input ratios retained. Write \(w_0=\ell_3/(\ell_2y_0)\) and write
\(\mathcal H\) for reduced rational height. After scaling \(Y=y_0v\)
and \(Z=w/w_0\), the equation at the simple root becomes
\((v+w)^p=v^{p-1}\), with \(v(0)=1\). Put \(h=(v+w)/v\).
Then \(h=1+wh^p\) and \(v=h^{-p}\). The coefficients of \(h\)
and \(v\) are integral by recursion on the degree. Thus
\(Q_J(Z)\) is the truncation of \(y_0 v(w_0Z)\).
Writing the two rational inputs in lowest terms yields a common
denominator for which
\[
   \max(D_J,\|Q_J^*\|_\infty)
     \ll_{p,J}\mathcal H(y_0)\mathcal H(w_0)^{J-1}.
\]
The denominator support in this bound agrees with the pool already
used in the local argument.

At depth two, the improved price can be read without
a recursion. The coefficient of \(Z\) in \(Q_2\) is
\(-p\ell_3/\ell_2\). For unit-sized plane coefficients
and no removed primes, the denominator and coefficient
height of \(Q_2\) after clearing denominators are
bounded by a fixed multiple of the primitive
\(L\)-height. When the second root and third radical
have comparable normalized exponent \(r>0\), this
replaces the quotient height price \(\theta\) in
(UP3) by \((\theta-r)_+\), up to the comparability
allowance. More generally (UP5) pays for low
third-letter valuations and for the actual size of the
depth-\(J\) reconstruction polynomial. This mechanism
uses divisibility of the original third letter at the
recorded powers, while the two free letters may have
small and different residual orders.

### Large adaptive orders

If the two complete weighted powers also have magnitude at
most \(X^{H_0}\), one can allow orders growing with height
without a new error on any large-order branch. Fix an inner
exponent allowance. If both orders exceed a sufficiently
large fixed cutoff, enumerate the two roots directly at
subpower cost for that allowance. If \(p\) exceeds the
cutoff and \(q\) does not, fix \(x\) and the first complete
value at the same cost. For \(x>1\) there are only
\(O_{H_0}(\log X)\) possible \(p\)'s, while \(x=1\)
makes that complete value independent of \(p\).
The plane makes \(y\) an affine function of \(z\) with
nonzero slope and constant term. The original expression
\(a x^p+b y(z)^q\) has degree \(q\ge2\). If it vanishes
at zero, its derivative there is nonzero, so removing
the factor at zero still leaves a nonconstant factor.
Otherwise retain the whole polynomial. The same
fixed-support argument counts \(z\) after including the
cleared denominator and the nonzero constant term in the
pool. This proves a uniform height-subpower bound on all
such large-order branches, with cutoffs set using the
chosen inner allowance. The quantitative quotient prices
are needed only for the remaining fixed orders.

Browning–Verzobio also analyze polynomial restrictions of
weighted powers in their generalized Fermat slices. Their
degree comparisons bound when a sum can itself be a
perfect power on a curve. Here the third letter's radical
is the polynomial parameter after reconstruction, and the
nonconstant remaining factor is controlled by supported
integer enumeration. This leaves the third valuations
free above the required recorded depth.
[Browning–Verzobio, Lemmas 2.7 and 2.15](https://arxiv.org/html/2608.24512v1#S2.SS3)
