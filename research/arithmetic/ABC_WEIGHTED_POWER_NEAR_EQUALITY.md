# Uniform near-equality counts for weighted powers

Fix \(B>0\), \(0<\theta<h\), and positive integers \(d_1,d_2\). Put
\(s=(d_1+1)(d_2+1)\) and \(D=s(s-1)/2\). Let \(u,v,w>0\) satisfy
\[
 w(s-1)>d_1u+d_2v,\qquad
 \max\{3\max(u,v),\,(s-1)w\}<\frac{h-\theta}{8D}.
 \tag{WP0}
\]
The second condition is used only to make a fixed power of the physical
precision dominate all auxiliary root-grid scales. Constants in the
following estimate may depend on these fixed parameters.

## Moving two-base power fibres

Uniformly for positive real \(a,b,T\) and positive integers
\(j,k\le(\log X)^B\) with \(j>d_1\) or \(k>d_2\), the number of positive
integer pairs satisfying
\[
 x\le X^u,\quad y\le X^v,\quad
 \min(a x^j,b y^k)\ge X^h,\quad
 |a x^j+b y^k-T|\le X^\theta
\]
is \(O(X^w(\log X)^{O_{d_1,d_2,B}(1)})\). In particular, for all
fixed \(C,\eta>0\) there is an integer \(d=d(C,\eta)\) such that for all
sufficiently small fixed \(L>0\), two bases at most \(X^{CL}\) with
\(\max(j,k)>d\) give
\[
 \#\{(x,y)\}\ll_{B,C,h,\theta,\eta,L,\varepsilon}
                  X^{\eta L+\varepsilon}.                     \tag{WP1}
\]
For this statement take \(d_1=d_2=d\), \(u=v=CL\) and \(w=\eta L\).
The required smallness of \(L\) depends on \(C,\eta,h-\theta\).
The factor polynomial in \(\log X\) absorbs summation over all
orders in the specified range. If the coefficients are positive
integers and the two summands are at most \(X\), every order for a
base greater than one is at most \(\log_2 X\). When a base equals one,
its summand is fixed and spacing gives at most one value of the
other base for each remaining order. Thus arbitrary orders give
the same bound for distinct summand pairs in that integral setting.

### Normalization and spacing

Write \(\tau=h-\theta\), \(m=\max(u,v)\) and
\(J=\lceil(\log X)^B\rceil\). Split by which summand is larger.
We describe the orientation \(b y^k\ge a x^j\). The opposite orientation
exchanges \(u,v\), \(d_1,d_2\) and the letters. Since \(T\gg X^h\), the
relative discrepancy is \(O(X^{-\tau})\).
If \(a x^j/T<X^{-2m}\), then \(b y^k/T\) is within
\(X^{-2m}+O(X^{-\tau})\) of one. Two distinct integers
\(y\le X^v\) with \(b y^k\asymp T\) give two values separated by
\(\gg T X^{-v}\), by monotonicity and the mean-value theorem
or the binomial expansion. Thus this part has at most one \(y\).
For that \(y\), consecutive eligible \(a x^j\ge X^h\) are separated by
\(\gg X^{h-u}\), so it has at most one \(x\).

On the remaining part split the bases into geometric boxes of ratio
\(1+1/J\). There are \(O(J\log X)\) boxes per base. Write \(p=x/U\),
\(q=y/V\) for the positive real left endpoints \(U,V\), and let
\(\alpha=aU^j/T\), \(\beta=bV^k/T\). On a nonempty box
\[
 \beta\gg1,\quad 1\le p,q\le 1+1/J,\quad
 t=\alpha p^j\ge X^{-2m},\quad 1-t\gg1
\]
at its eligible grid points. Indeed the larger summand occupies at
least a fixed fraction of \(T\), and the \(j\)-th and \(k\)-th powers of
normalized coordinates are each between one and \(e\). Divide the
displayed near-equality by \(T\) and use the positive branch
\(f(p)=((1-\alpha p^j)/\beta)^{1/k}\).
Then \(|q^k-f(p)^k|\ll X^{-\tau}\) at every eligible point.
Since \(q^k\ge1\), the derivative of the \(k\)-th power between
\(q\) and \(f(p)\) is bounded below uniformly even for moving \(k\).
Thus \(|q-f(p)|\ll X^{-\tau}\). Monotonicity
gives \(1-t\gg1\), \(t\ge X^{-2m}\), and \(f(p)\asymp1\) throughout
the interval spanned by their \(p\)-coordinates. On this interval all
fixed-order derivatives of \(f\) have bounds polynomial in \(J\).
This follows by differentiating \(\log f=(\log(1-\alpha p^j)-\log\beta)/k\)
and using \(t\ll1\), \(1-t\gg1\), and \(p\asymp1\).

### Small auxiliary relation

Partition the interval into cells of length \(X^{-w}\). In one cell
form an \(s\)-row evaluation determinant of the monomials
\(p^i q^\ell\), \(0\le i\le d_1\), \(0\le\ell\le d_2\).
Multiplying the column \((i,\ell)\) by \(U^i V^\ell\) gives integer
entries even when \(U,V\) are real. The product of these factors is
at most \(X^{s(d_1u+d_2v)/2}\).
The analytic rows with \(q=f(p)\) admit Taylor expansions in \(p\)
through degree \(s-2\), with remainder
\(O(X^{-(s-1)w}J^{O_s(1)})\).
The original rows have the same expansion with that remainder,
by \(|q-f(p)|\ll X^{-\tau}\) and (WP0).
In a nonzero multilinear Taylor term the chosen Taylor orders
must be distinct. If it uses \(r\) remainder rows, the total order
is at least \((s-r)(s-r-1)/2+r(s-1)\ge D\).
The determinant is therefore \(O(X^{-wD}J^{O_s(1)})\).
The strict first inequality in (WP0) forces the cleared integer
determinant to vanish for all sufficiently large \(X\).
Thus some nonzero polynomial \(F(p,q)\) of bidegree at most
\((d_1,d_2)\) vanishes at all the original grid points of the cell.
We may choose real coefficients normalized to have maximum modulus
one. Put \(H=F(p,f(p))\). At every original grid point of the cell,
\[
                         |H(p)|\ll_{d_1,d_2} X^{-\tau}.        \tag{WP2}
\]
No bound for the heights of the coefficients of \(F\) is needed.

### Uniform contact bound

On this graph use the differential operator
\(\mathcal D=p(1-t)\,d/dp\) with \(t=\alpha p^j\).
Its iterates on \(p^i f^\ell\) are \(p^i f^\ell P_{i,\ell,n}(t)\),
where \(P_{i,\ell,0}=1\) and
\[
 P_{i,\ell,n+1}=(i(1-t)-\ell jt/k)P_{i,\ell,n}
                       +jt(1-t)P'_{i,\ell,n}.
\]
Here \(\deg P_{i,\ell,n}\le n\) and \(k^nP_{i,\ell,n}\) has integral
coefficients. The determinant \(W(t)\) of this polynomial array
for \(n=0,\dots,s-1\) has degree at most \(D\) and \(k^DW(t)\)
is integral. It is nonzero. To see this, the plane polynomial
\(q^k-(1-\alpha p^j)/\beta\) is irreducible over \(\mathbf C(p)\)
by Eisenstein at any simple zero of \(1-\alpha p^j\), and also
irreducible in \(\mathbf C[p,q]\). Its bidegree is \((j,k)\).
No nonzero polynomial of bidegree \((d_1,d_2)\) can contain its
positive analytic branch, because one of \(j,k\) is greater than
the corresponding degree. Consequently the displayed monomials on
the graph are linearly independent real analytic functions.
Their ordinary Wronskian is not identically zero, as follows by
choosing linear combinations of distinct Taylor leading orders at
an interior point. Passing to \(\mathcal D\) multiplies that Wronskian
by a nonvanishing factor and factoring the monomial from each column
leaves \(W(t)\).

For the points where \(|W(t)|<X^{-\tau/4}\), factor \(W\) over
\(\mathbf C\) and use the lower bound \(k^{-D}\) on its leading
coefficient. If \(\deg W>0\), each such \(t\) is within
\(X^{-\tau/(4D)}J^{O_s(1)}\) of one of at most \(D\) complex roots.
The case of constant \(W\) has no such points for large \(X\).
Because \(dt/dp=jt/p\gg X^{-2m}\), each of these exceptional real
intervals contains at most one grid coordinate \(p\), by (WP0)
and the spacing \(1/U\ge X^{-u}\).
For each \(p\) there is at most one \(q\) in the whole box, by
the original power spacing on the larger summand.

At every remaining candidate at least one of the quantities
\(\mathcal D^n H\), \(0\le n<s\), has modulus
\(\gg_s X^{-\tau/4}J^{-O_s(1)}\). This follows by inverting the
Wronskian matrix, whose entries and cofactors are polynomially bounded
in \(J\), with determinant of modulus \(\gg_s |W(t)|\).
The coefficient vector has maximum modulus one.
Conversion to ordinary derivatives, using the polynomial-in-\(J\)
bounds for \(p(1-t)\) and its fixed derivatives, implies
\(\max_{n<s}|H^{(n)}(p)|>X^{-\tau/2}\) for large \(X\).

For each \(n<s\), the superlevel sets
\(H^{(n)}(p)>X^{-\tau/2}\) and \(H^{(n)}(p)<-X^{-\tau/2}\)
have \(O_s(J^2+1)\) interval components. Indeed \(H^{(n)}\)
is a rational expression in \(p,f(p)\), with denominators cleared
by powers of \(f\), since \(f'=-\alpha j p^{j-1}/(k\beta f^{k-1})\).
The resulting polynomial equations for the two boundary levels
have degrees \(O_s(J+1)\). If a boundary polynomial vanishes
identically on the graph, the derivative is identically that level.
Otherwise Bézout with the irreducible plane polynomial above
bounds the number of level crossings by \(O_s(J^2+1)\).
For \(n=0\), no point in either superlevel set satisfies (WP2).
For \(n\ge1\), any \(n+1\) grid coordinates in one component
give an \(n\)-th divided difference of magnitude
\(O_s(X^{-\tau+nu})\) by (WP2) and spacing. The mean-value identity
for divided differences contradicts the lower derivative bound
on the interval between them, by (WP0). Thus each cell has
\(J^{O_s(1)}\) grid points. Summation over cells and boxes proves (WP1).

## Both signs and all omitted-letter orientations

The same bound (WP1), with a possibly smaller budget threshold, holds for
\[
        |a x^j-b y^k-T|\le X^\theta,\qquad
        \min(a x^j,b y^k,T)\ge X^h.                            \tag{WP4}
\]
Use \(d_1=d_2=d\) and put \(m=CL\), \(w=\eta L\),
\(\tau=h-\theta\), and \(E=\tau/(32D)\). Take \(L\) smaller if needed
so that \(3m+E<\tau/(4D)\), \(2dm<E\), and \((s-1)w<E/2\).
Write \(U=a x^j\), \(V=b y^k\). If \(V<U X^{-2m}\), spacing
first determines \(x\) from \(U=T(1+O(X^{-2m}))\)
and then determines \(y\). The spacing of the second step uses
\(V\ge X^h\).

Use the same geometric boxes for the other pairs, now with \(p\)
the normalized \(y\), \(q\) the normalized \(x\),
and write \(U=S q^j\) and \(V=B_0 p^k\) for positive \(S,B_0\)
fixed within a box. In the part \(T/U\ge X^{-E}\) use the exact graph
\(f(p)=((T+B_0p^k)/S)^{1/j}\).
Here \(S\gg X^h\). Since \(q^j\ge1\), the \(O(X^{-\tau})\)
error in \(j\)-th powers gives the same bound on \(|q-f(p)|\)
by the uniform spacing of the real root map near \(q^j\).
On the interval spanned by the candidates put
\[
     t=\frac{B_0p^k}{T+B_0p^k},\qquad
     \mathcal D=p\,d/dp,\qquad
     \mathcal Dt=kt(1-t),\qquad
     \mathcal D\log f=(k/j)t .
\]
Both \(t\gg X^{-2m}\) and \(1-t\gg X^{-E}\) hold throughout
this interval by monotonicity and the endpoint bounds.
The logarithmic derivative identity gives the required
polynomial-in-\(J\) derivative bounds for \(f\).
Apply the same auxiliary determinant, with the same cell width.
In the contact calculation the recurrence becomes
\[
 P_{i,\ell,n+1}=(i+\ell kt/j)P_{i,\ell,n}
                         +kt(1-t)P'_{i,\ell,n}.
\]
Its determinant is nonzero by the same simple-zero argument, now
for \(q^j=(T+B_0p^k)/S\). Its degree is still at most \(D\),
with leading coefficient at least \(j^{-D}\) in absolute value.
For the determinant cutoff \(X^{-\tau/4}\), the exceptional intervals in
\(t\) have the same power width as before. Their inverse images have
length \(O(X^{2m+E-\tau/(4D)}J^{O_s(1)})\), smaller than grid
spacing. The derivative argument for the remaining points is unchanged.

In the part \(T/U<X^{-E}\), one has \(U/V=1+O(X^{-E})\).
Put \(g=\gcd(j,k)\), \(j'=j/g\), \(k'=k/g\).
If both \(j',k'\le d\), the reduced rational \(x^{j'}/y^{k'}\)
lies in a relative \(O(X^{-E})\) strip about \((b/a)^{1/g}\).
A nonzero cross determinant would have magnitude
\(O(X^{-E+(j'+k')m})<1\), so the ratio is unique.
It fixes \(K=U/V\). Any solution has \(K>1\) for large \(X\),
since \(T\gg X^\theta\). The relation
\(|(K-1)b y^k-T|\le X^\theta\) with \(T\ge X^h\)
then determines \(y\) by spacing, and hence determines \(x\).
This argument does not require coprime bases.

If \(\max(j',k')>d\), the comparison graph in the box is the
monomial \(f_0(p)=(B_0/S)^{1/j}p^{k/j}\) and the error is
\(O(X^{-E})\). The error in powers gives the same bound in
roots as above, since \(1\le q^j\le e\). Its fixed-order derivatives are polynomially
bounded in \(J\) on the candidate interval, so the same Taylor
determinant applies with this precision. For the monomials
\(p^i f_0^\ell\), the operator \(p\,d/dp\) has eigenvalues
\(i+\ell k/j\). These are pairwise distinct on the chosen
rectangle because \(\max(j',k')>d\), and are separated by at
least \(1/j\). The Wronskian after removing the bounded positive
monomial factors is a Vandermonde determinant of magnitude
\(\gg J^{-D}\). Inverting it gives a polynomial-in-\(J\) lower
bound on some ordinary derivative of the auxiliary \(F(p,f_0(p))\).
Use levels \(\pm X^{-E/2}\) for the divided-difference argument.
The defining polynomial of the positive graph can be taken as
\(q^{j'}-(B_0/S)^{1/g}p^{k'}\), irreducible because
\((j',k')=1\). For instance, an algebraic root over
\(\mathbf C(p)\) has valuation \(k'/j'\) at \(p=0\), so
the ramification index and hence its degree are at least \(j'\).
The derivative level sets have the same polynomial-in-\(J\)
component bound by Bézout. This proves (WP4).

Consequently the two-base estimate applies whether the fixed
complete letter is a summand or the target, as long as all three
physical heights exceed \(X^h\). It uses no congruence depth
or common-multiplier choices.

## What the bound does on actual target dictionaries

Write \(N=X^L\) and let the fixed target shell have height comparable
to \(X\). Suppose candidate summand pairs have the above
factorizations with \(x,y\le N^C\) and one order greater than \(d\).
Impose \(\min(a x^j,b y^k)\ge X^h\) with \(h>\theta\).
If for each complete target the coefficient records available
to the filtered test number at most \(N^{f+o_X(1)}\), then the real
cutoff leaves at most \(N^{f+\eta+o_X(1)}\) candidate support pairs
per target radical. Here the cost of exact target representatives
on a fixed radical is subpower. In particular, the record hypothesis
holds if the two coefficient radicals have product at most \(N^f\)
and the coefficients are integers at most \(X\). Choose their radical
product and assign each prime to its membership pattern in the two
coefficients, then use the fixed-support estimate. In the primitive
case the pattern has only two choices per prime. Without primitivity,
enumerate the pair of radicals as divisors of their product, still at
subpower cost. Thus overlap of supports does not require primitivity
for this bound. The polylogarithmically many adaptive orders add a subpower factor.

The independent prime-block test now has mean false tail
\[
 O\!\left(N\min\{1,N^{f+\eta-\beta-u_0+o_X(1)}\}^{\,r}\right),
 \tag{WP3}
\]
for target radical cap \(N\), threshold \(N^{u_0}\), and \(r\)
independent prime blocks of total exponent \(\beta\) each.
Apply the candidate-by-candidate discrepancy divisibility estimate
and the same argument as (SI3). In balanced boxes this replaces
the free-letter support exponent near one by the paid coefficient
exponent \(f\), whenever \(f<1\). If \(f+\eta<u_0\), even the real
filter alone precludes popularity at this threshold in the sector.
The count concerns the original summands and target. For exact primitive rows, (SH1)
already gives a uniform saving when the shorter side has numerical
height below a fixed absolute power. Thus the tall assumption used
here is compatible with the remaining saturation sector. In the
one-target minimally repaired problem it therefore counts primitive
quotients, without including the extra choices of common multiplier.
[Short-interval target test (SI3)](ABC_SHORT_INTERVAL_RADICAL_FILTER.md#consequences-for-real-filtered-popularity-tests)
· [Minimal repair transfer](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)
· [Strict saving for short summands (SH1)](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-strict-saving-when-the-short-term-is-small)

The minimum summand height matters even for a tall target. Set
\(a=1\), \(y=k=1\), and \(b=T\) with \(T\asymp X\). Fix \(j>d\)
and take \(L\) small enough that \(CLj<\theta\). Every \(x\le X^{CL}\)
then passes the near-equality filter. This violates (WP1) with
\(\eta<C\) and sufficiently small \(\varepsilon\) if the lower height
bound on \(a x^j\) is removed.
At fixed \(T\), geometry of the actual coefficient records matters.
Choosing exact power decompositions does not by itself show that their
effective record count is less than the elementary support count
on a diffuse profile.

Huxley develops determinant bounds for rational and integer points
near smooth curves. Adiceam and Marmon use approximate determinants
for small values of nonsingular ternary forms, with separate
treatment of neighbourhoods of high contact. The first determinant
step here follows this general method. The coefficient-uniform
terminal step uses the displayed binomial Wronskian recurrence and
a bound on its degree independent of the moving orders. The
coefficient-uniform exact plane-curve estimate already used in
the fine-grid bound takes equality as its input. This result
treats a moving real interval of numerically large width.
[Huxley, The Rational Points Close to a Curve](https://www.numdam.org/item/ASNSP_1994_4_21_3_357_0/)
· [Adiceam–Marmon, Rational points close to non-singular algebraic curves](https://arxiv.org/abs/2204.11607v2)
· [Browning–Verzobio, Lemma 2.2](https://arxiv.org/html/2608.24512v1#S2.SS1)

## Scope of fixed-target literature

Pastén's Mordell estimate gives a square–cube gap of order at least
\((\log X)^2/(\log\log X)^8\) in its height parameter. This has no
power-sized exclusion window. His cubic Thue–Mahler bound retains
a cubic dependence on the largest support prime in the displayed
upper bound for logarithmic height. This does not exclude solutions
at the heights of the balanced diffuse boxes, where support primes
can be polynomial in the original height.
[Pastén, Theorems 1.2 and 1.8](https://arxiv.org/html/2608.23559v2#S1)

Ziegler's \((1,1,c)\) case starts at at least \(|S|+1\) ordered
solutions on the same support set \(S\). The logarithmic height
cutoff depends on the largest support prime, the regulator and
factorials in \(|S|\). The popular targets in the balanced endpoint
may have many different pairs of free supports and need not meet
that hypothesis on any one small support set.
[Ziegler, Theorem 1](https://arxiv.org/html/2609.15487v1#S2)

Fouvry and Waldschmidt count exact represented values for families of
binary fewnomials under the coefficient restriction
\(A^*(F)\le\exp(\eta\,\deg(F)/\log\deg(F))\). At a fixed order
this coefficient bound is independent of the height limit. It does
not control all the coefficient records of (WP3), which can have
positive numerical height exponent even at fixed budget and order.
[Fouvry–Waldschmidt, Theorem 1.4](https://webusers.imj-prg.fr/~michel.waldschmidt/articles/pdf/Fewnomials-EF-MW.pdf)
