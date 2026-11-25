# Recovering the kernel gcd by an elementary real-quadratic norm count

Let an actual selected simple bipartite graph have squarefree positive
left kernels \(b\le K_1\), right kernels \(d\le K_2\), and one fixed
positive representation per edge,
\[
 \kappa_1 a^2b^3+\kappa_2 c^2d^3=n,\qquad
 \gcd(\kappa_1 a^2b^3,\kappa_2 c^2d^3)=1,\quad
 a\le R_1,\ c\le R_2.                                      \tag{1}
\]
The coefficients are fixed positive integers during each graph count.
Assume \(n\ge2\), all four caps lie between1 and \(n\), and
\[
                         2(R_1R_2)^2<n.                    \tag{2}
\]
The path count below does not need (2); it enters only the existing
actual codegree bound used to pass from paths to cycles.

At a shared-right path \(b_1-d-b_2\), with \(b_1\ne b_2\), put
\[
 \alpha=\gcd(a_1,a_2),\quad a_i=\alpha r_i,\qquad
 g=\gcd(b_1,b_2),\quad b_i=g u_i,
\]
\[
 L=\frac{c_1^2-c_2^2}{\alpha^2g^3},\qquad
 e'=\operatorname{sf}(|L|),\qquad
 {\cal V}=\frac{\alpha^3g^{3/2}}{\sqrt{e'}}.                 \tag{3}
\]
Here sf is the squarefree kernel: \(|L|=e'v^2\) uniquely with
\(e'\) squarefree and \(v\) a positive integer. Square roots in
counting weights are positive real numbers.

The actual coprimality and squarefree conditions imply \(L\) is a
nonzero integer; a proof is included below.

## 1. A uniform elementary norm lemma

**Lemma 1.** Let \(D>1\) be squarefree, let \(M\ne0\) be an integer,
and let \(H\ge1\). The number of integer solutions
\[
                  x^2-Dy^2=M,\qquad |x|,|y|\le H             \tag{4}
\]
is
\[
 O\!\left(\tau(|M|)\,2^{\omega(|M|)}
              [1+\log(H(1+\sqrt D))]\right),                \tag{5}
\]
with an absolute implied constant. In particular it is
\(O_\epsilon(n^\epsilon)\), uniformly, when \(|M|,D,H\) are bounded
by fixed powers of \(n\). The constant in this consequence can depend
on the powers and on \(\epsilon\), but not on \(D,M\).

No ideal factorization, class number estimate, fundamental-unit
existence theorem, or Dirichlet unit theorem is needed.

Put \(s=\gcd(x,y)>0\). Then \(s^2\mid M\). Write
\(x=sx',y=sy'\), \(M'=M/s^2\), so \(\gcd(x',y')=1\).
A prime dividing \(y'\) cannot divide \(M'=x'^2-Dy'^2\);
thus \(\gcd(y',M')=1\).

Let \(m=|M'|\). For \(m>1\), each primitive solution determines
\[
                  r=x'/y'\pmod m,\qquad r^2=D\pmod m.       \tag{6}
\]
For \(m=1\), use one residue class; all divisibility statements
below are then automatic.

There are at most \(2^{\omega(m)+1}\) possible classes in (6).
For an odd prime \(p\nmid D\), unit square roots have at most two
classes modulo any \(p^a\): the factors of a difference of squares
show uniqueness of each sign, or one can use the invertible
derivative \(2r\). At \(p=2\nmid D\) there are at most four classes:
for two odd roots the minimum of the valuations of their sum and
difference is one, giving the two signs modulo \(2^{a-1}\), each
with at most two lifts modulo \(2^a\). Small exponents give no more
classes.

If \(p\mid D\) and \(p\mid M'\), primitiveness gives \(p\mid x'\)
and \(p\nmid y'\). Since \(D\) is squarefree,
\(v_p(x'^2-Dy'^2)=1\), so \(v_p(M')=1\); there is only the zero
root modulo \(p\). This argument includes \(p=2\).
The Chinese remainder theorem proves the claimed root bound.
No bad prime or high target valuation is omitted.

Take two primitive solutions \((x',y')\), \((X,Y)\) in the SAME
class (6). Their algebraic numbers have the same norm \(M'\), and
\[
 \frac{x'+y'\sqrt D}{X+Y\sqrt D}
   =\frac{(x'X-Dy'Y)+(y'X-x'Y)\sqrt D}{M'}
                                                   \in\mathbb Z[\sqrt D].
                                                                  \tag{7}
\]
Both numerator coefficients are divisible by \(m\), by (6), and
hence by the signed integer \(M'\). The quotient has norm one.
This is an integral norm-one element, not a fractional-order
assertion.

Put \(B=H(1+\sqrt D)\). Both real embeddings of each
\(x'+y'\sqrt D\) have absolute value at most \(B\), and their
product has absolute value \(m\). Its first embedding therefore
has absolute value in \([m/B,B]\). If two such absolute values
differ, their ratio, after replacing the algebraic quotient by
its negative if needed, is a positive real number \(u>1\).
Equation (7) makes \(u\) an algebraic integer of norm one, with
conjugate \(u^{-1}\). Its integer trace satisfies
\[
                u+u^{-1}\ge3,\qquad
                u\ge(3+\sqrt5)/2.                           \tag{8}
\]
If the quotient is rational it is \(1\) or \(-1\), so it cannot
have absolute value greater than one. Equal absolute values have
multiplicity at most two, since a real field embedding is injective.
Consequently the values in each class have a uniform multiplicative
spacing and there are \(O(1+\log B)\) of them. This uses neither
a chosen fundamental unit nor a regulator lower bound.

There are at most \(\tau(|M|)\) choices of \(s\), and
\(2^{\omega(m)}\le2^{\omega(|M|)}\). Summing the class bound proves
(5). If \(m>B^2\), no solution exists; this empty case is harmless.
The usual elementary divisor bound absorbs both divisor factors and
the logarithm in the polynomial-height consequence.

## 2. The actual cubic and norm equations

Actual within-edge coprimality gives
\[
 \gcd(\alpha,\kappa_2d)=\gcd(g,\kappa_2d)
       =\gcd(\kappa_1,\kappa_2d)=1.
\]
Subtract the two edge equations. First \(\alpha^2\) divides
\(c_1^2-c_2^2\); then \(g^3\) divides the quotient. If \(c_1=c_2\),
prime-valuation parity in \(a_1^2b_1^3=a_2^2b_2^3\) would force
the two squarefree kernels \(b_1,b_2\) to coincide. Thus \(L\ne0\).
The reduced kernels satisfy
\[
 \kappa_1r_1^2u_1^3-\kappa_1r_2^2u_2^3
                                 +\kappa_2Ld^3=0.            \tag{9}
\]
Here \(u_1,u_2\) are coprime squarefree positive integers,
\(u_1\ne u_2\), and each is coprime to \(d\).
In particular \((u_1,u_2,d)\) is primitive and all three signed
cubic coefficients are nonzero.

The sole deep counting input is Heath-Brown,
*Counting Square-full Solutions to x+y=z*, arXiv2601.07817v1,
Theorem1.3: for fixed nonzero signed integer coefficients,
primitive diagonal cubic zeros of maximum height \(H\) number
\(O_\epsilon(H^\epsilon)\), uniformly in those coefficients.
Its exact scope, rather than a global square-full triple estimate,
is used here. For (9), \(H\le n\).

Now fix ONLY the coefficient data \((r_1,r_2,L)\), and fix one
reduced kernel solution \((u_1,u_2,d)\) of (9) capable of producing
an actual path. Because \(\gcd(\kappa_1,\kappa_2d)=1\), (9)
implies \(\kappa_1\mid L\).
Set
\[
 A_i=\kappa_1r_i^2u_i^3,\qquad X=\alpha^2g^3.
\]
The two actual edge equations are
\(A_iX+\kappa_2c_i^2d^3=n\). Eliminating \(X\) gives
\[
                  A_2c_1^2-A_1c_2^2=nL.                    \tag{10}
\]
Consequently the integers
\[
 x=r_2u_2^2c_1,\quad y=r_1u_1c_2,\quad
 D=u_1u_2,\quad M=u_2nL/\kappa_1
\]
satisfy
\[
                         x^2-Dy^2=M.                        \tag{11}
\]
The integer \(D>1\) is squarefree because the reduced kernels are
distinct, coprime and squarefree. The integer \(M\) is nonzero.
The height bounds are uniform:
\[
               D\le n^2,\quad |M|\le n^4,\quad
                       |x|,|y|\le n^4.                     \tag{12}
\]
Indeed \(r_i,u_i,c_i\le n\) and \(|L|<R_2^2\le n^2\);
division by the positive \(\kappa_1\) can only reduce the bound.

Lemma1 gives \(O_\epsilon(n^\epsilon)\) possible pairs \(c_1,c_2\),
because \(x,y\) recover them at most once. The value
\[
                  X=\frac{n-\kappa_2d^3c_1^2}
                            {\kappa_1r_1^2u_1^3}            \tag{13}
\]
is then fixed. Discard nonpositive or nonintegral values.
A positive integer has at most one expression
\(X=\alpha^2g^3\) with \(\alpha>0\) integral and \(g\) squarefree:
the parity of every prime valuation fixes \(g\), and then fixes
\(\alpha\). Valuation one gives no expression.
Thus no independent \(g\), \(\alpha\), root sign, or
representation multiplicity is added.

Recovered \(b_i=gu_i\) determine at most one selected ordered actual
path. All further restrictions, such as dyadic ranges, squarefreeness
of \(gu_i\), the prescribed gcd, and agreement with the selected
edge labels, only discard candidates.

**Proposition 2.** Each fixed coefficient tuple
\((r_1,r_2,L)\) has \(O_\epsilon(n^\epsilon)\) actual selected
shared-right path completions. The two subpower counts are the
primitive cubic theorem for the reduced kernels and Lemma1 for
the roots. This is uniform in the fixed positive coefficient pair.

## 3. Coefficient frequency without a g sum

For real \(A,G,E\ge1\), restrict
\[
       \alpha\in[A,2A),\qquad g\in[G,2G),\qquad e'\le E.
\]
Actual positivity gives
\[
                 0<|L|<\frac{R_2^2}{A^2G^3}=L_0.
\]
The number of possible nonzero signed integers \(L\) with squarefree
kernel at most \(E\) is at most
\[
 2\sum_{\substack{1\le t\le E\\t\ {\rm squarefree}}}
                \left\lfloor\sqrt{L_0/t}\right\rfloor
                       \le4\sqrt{L_0E}.                    \tag{14}
\]
If \(L_0<1\), there are no lifts. Enlarging \(E\) beyond its actual
range is harmless. The variables \(r_i\le R_1/A\) give at most
\((R_1/A)^2\) ordered pairs; if \(A>R_1\) the path set is empty.

Use Proposition2 AFTER this coefficient count. It gives
\[
 \#\{\hbox{selected paths in this range}\}
      \ll_\epsilon
       \frac{n^\epsilon R_1^2R_2\sqrt E}{A^3G^{3/2}}.          \tag{15}
\]
In particular there is NO sum over \(g\) in (14)--(15): it is
recovered after the two completion bounds. This is the extra
factor \(G^{-1}\) compared with counting every integer \(g\)
as a coefficient index.

## 4. The stronger weighted path and cycle theorem

**Theorem 3.** Under the stated actual graph assumptions,
\[
          \sum_{\text{shared-right paths}}{\cal V}
                        \ll_\epsilon n^\epsilon R_1^2R_2.   \tag{16}
\]
Under (2) also,
\[
 \sum_{C_6}\sum_{\text{shared-right paths in }C_6}{\cal V}
                        \ll_\epsilon n^\epsilon K_1R_1^2R_2,\tag{17}
\]
and for EVERY real \(Y>0\),
\[
 \#\{C_6:\hbox{some shared-right path has }{\cal V}\ge Y\}
                        \ll_\epsilon
                           n^\epsilon K_1R_1^2R_2/Y.        \tag{18}
\]
There are side-swapped versions.

Partition \(\alpha,g,e'\) into dyadic ranges. They are bounded by
\(n,n,n^2\), so there are \(O(\log(n)^3)\) range triples.
On a block with \(e'\in[E,2E)\), the maximum weight is
\(O(A^3G^{3/2}/\sqrt E)\). Multiplying (15), with upper bound
\(2E\), cancels all three dyadic variables and proves (16) after
absorbing logarithms with a smaller initial epsilon.

The exact general codegree theorem for the actual selected graph
under (2) gives \(O_\epsilon(n^\epsilon)\) common neighbors for
each same-side pair. A selected path has
\(O_\epsilon(K_1n^\epsilon)\) cycle continuations: choose the
third left kernel, then its two common right neighbors.
This proves (17). Positivity gives (18) by Markov, including all
small positive thresholds. This uses the existing general
rectangular codegree proof, not an extension of a balanced-only
statement.

The primitive cubic theorem and the actual codegree theorem are
separate inputs. The new norm-recovery and coefficient-frequency
arguments are supplied here. A bounded completion count by itself
would not bound path frequency; (14) is essential.

## 5. Comparisons and counted sectors

Write
\[
 {\cal W}=\alpha^3\sqrt{g/e'}
\]
for the preceding kernel-gcd-refined weight. Then
\[
                         {\cal V}=g{\cal W}\ge{\cal W}.      \tag{19}
\]
For the earlier actual square quotient
\(e_0=|\ell|/\gamma^2\), where
\(\ell=g^3L\) and \(\gamma=\gcd(c_1,c_2)\), put
\(W_0=\alpha^3/\sqrt{e_0}\).
Actual coprimality gives
\(\gcd(\alpha,\gamma)=\gcd(g,\gamma)=1\), so
\(\gamma^2\mid L\). Since \(|L|=e'v^2\) with \(e'\) squarefree,
\(\gamma\mid v\). Hence
\[
                  {\cal V}/W_0=g^3(v/\gamma)\ge g^3.         \tag{20}
\]
All old path/cycle moments and every positive-threshold tail
follow termwise, without taking them as premises.

The root cap gives
\[
 \alpha g^{3/2}\sqrt{e'}\,v<R_2,\qquad
 \alpha^4g^3<R_2{\cal V},\qquad
                    {\cal V}<\alpha^2R_2.                 \tag{21}
\]
In particular \({\cal V}>1/R_2\), and a sector
\({\cal V}\ge Y\) requires
\(\alpha>\sqrt{Y/R_2}\). A threshold merely below one does not
automatically count every cycle.

For the coefficient-shaped critical rectangle,
\[
 K_1R_1^2R_2
  =16n^{131/153}/(\kappa_1^{19/9}\kappa_2^{14/9}).
\]
The general codegree condition holds for all sufficiently large
\(n\), uniformly, since
\(2(R_1R_2)^2\le32n^{44/51}<n\); the bounded initial range is
absorbed uniformly. Therefore the critical budget
\(n^{26/51+\epsilon}\) holds on
\[
 {\cal V}\ge
      n^{53/153}/(\kappa_1^{19/9}\kappa_2^{14/9}).            \tag{22}
\]
For coefficient one, membership requires
\(\alpha>2^{-1/2}n^{10/153}\).

In the balanced maximum box \(R_1=R_2=K_1=T\), the tail is
\(O_\epsilon(n^\epsilon T^4/Y)\), giving the \(T^2n^\epsilon\)
cycle budget on \({\cal V}\ge T^2\). This sector requires
\(\alpha>\sqrt T\) and contains both preceding weight sectors.
These are counted sectors and necessary size conditions, not
actual strict-gain constructions or an exhaustion of the residual.

No full six-cycle hypothesis, all-target prescribed-sum saving,
or unrestricted aggregate over varying coefficient pairs is
asserted. The squarefree kernels, distinct same-side vertices,
actual coprimality and selected representations remain essential.


## 6. The coefficient-divisibility gain

Write the fixed positive coefficient uniquely as
\[
 \kappa_1=h_\kappa^2 k_\kappa,\qquad k_\kappa\ {\rm squarefree},
 \qquad \chi_1=h_\kappa k_\kappa.                            \tag{23}
\]
The exact divisibility \(\kappa_1\mid L\), already established in
Section2, improves the coefficient frequency.

**Corollary 4.** The right sides of the path frequency(15), path
moment(16), cycle moment(17), and every positive-threshold tail(18)
may each be divided by \(\chi_1\), with constants still uniform in
the fixed positive coefficient pair. Equivalently the weight
\(\chi_1{\cal V}\) has the former weighted bounds.
If \(\kappa_1\) is squarefree, the gain is \(1/\kappa_1\).
There is a side-swapped statement with the corresponding coefficient.

For \(L\ne0\), write \(|L|=e'v^2\), where \(e'\) is squarefree.
Prime-by-prime valuation comparison gives
\[
       \kappa_1\mid e'v^2
 \quad\Longleftrightarrow\quad
       \frac{h_\kappa k_\kappa}{\gcd(k_\kappa,e')}\mid v.      \tag{24}
\]
This also covers primes dividing both \(h_\kappa\) and \(k_\kappa\):
at a prime with \(v_p(\kappa_1)=2r+\eta\), \(\eta=0\) or1,
the minimum required exponent of \(v\) is \(r\) when
\(\eta=0\) or \(p\mid e'\), and \(r+1\) otherwise.

Consequently the number of signed lifts with \(0<|L|<L_0\),
\(\operatorname{sf}(|L|)\le E\), and \(\kappa_1\mid L\) is at most
\[
 \frac{2\sqrt{L_0}}{\chi_1}
       \sum_{\substack{e'\le E\\e'\ {\rm squarefree}}}
                  \frac{\gcd(k_\kappa,e')}{\sqrt{e'}}
 \le
       \frac{4\tau(k_\kappa)\sqrt{L_0E}}{\chi_1}.             \tag{25}
\]
Indeed \(\gcd(k,e)=\sum_{a\mid k,\ a\mid e}\varphi(a)\).
For each \(a\mid k\), dropping squarefreeness gives
\[
 \sum_{\substack{e\le E\\a\mid e}}\frac1{\sqrt e}
       \le\frac{2\sqrt E}{a},
 \qquad
 \sum_{a\mid k}\frac{\varphi(a)}a\le\tau(k).
\]
Empty ranges and real \(E\ge1\) obey the same inequalities.

Count only these necessary coefficient data BEFORE the primitive
cubic and norm completions. Replacing (14) by (25) proves the claimed
gain in (15). A nonempty graph has \(\kappa_1\le n\), so
\(\tau(k_\kappa)\) is absorbed by starting with a smaller epsilon.
If the graph is empty all conclusions are immediate.
The dyadic cancellation and the unchanged codegree/Markov arguments
give the gains in (16)--(18), proving the corollary.

In the critical shaped rectangle the counted sector can therefore
be enlarged to
\[
 {\cal V}\ge
 \frac{n^{53/153}}
 {\chi_1\,\kappa_1^{19/9}\kappa_2^{14/9}}.                    \tag{26}
\]
For squarefree \(\kappa_1\), its denominator is
\(\kappa_1^{28/9}\kappa_2^{14/9}\).
This is still a conditional membership sector furnished by an
unconditional weighted count. It does not show that every cycle
has a path in the sector, or permit an unrestricted coefficient
aggregate.


## 7. Classical context and the weaker refined-weight consequence

Uniform Pell counting has a classical antecedent. T. Estermann,
[Einige Saetze ueber quadratfreie Zahlen](https://doi.org/10.1007/BF01455836),
Math. Ann.105(1931),653--662, Hilfssaetze2--3, gives uniform bounds
for positive integer coefficients a,b and positive integer X, counting
positive solutions of \(ax^2-bz^2=l\), \(ax^2\le X\), with
\(ab\) nonsquare and nonzero signed \(l\); the final statements have
no coprimality hypothesis. The original
[positive-norm statement](https://gdz.sub.uni-goettingen.de/content/PPN235181684_0105/1000/0/00000664.jpg)
and [negative-norm statement](https://gdz.sub.uni-goettingen.de/content/PPN235181684_0105/1000/0/00000666.jpg)
are on printed660 and662. This is context only: Lemma1 above supplies
the needed count locally by primitive residue classes and trace
spacing. No claim of priority for uniform Pell counting is made.

The [preceding kernel-gcd-refined cycle bound](ERDOS_943_KERNEL_GCD_REFINED_SQUARE_QUOTIENT.md),
Theorem1, also follows by a real alternative route here.
Its weight is \({\cal W}\), and (19) gives
\({\cal W}\le{\cal V}\). Thus its weighted path and cycle sums
are termwise bounded by (16)--(17). For every \(Y>0\), a cycle
with some \({\cal W}\ge Y\) has some \({\cal V}\ge Y\), so (18)
gives the preceding tail with the same prefactor.
The old theorem is a consequence, not a premise of this argument.
Since \(\chi_1\ge1\), Corollary4 also implies the cycle bounds without the coefficient gain in Theorem3; the coefficient-gain branch can furnish
them independently of the direct dyadic argument without the coefficient gain.
