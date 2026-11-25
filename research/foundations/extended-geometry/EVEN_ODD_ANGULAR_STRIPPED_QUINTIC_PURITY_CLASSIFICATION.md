# Even height and odd angular parity in the stripped quintic space

Let \(K\) be any field of characteristic zero. Work in
\[
 A_K=K[x,y,t]/(t^2+2xy-1),\qquad
 \Lambda(F)=\operatorname{CT}_x I_t
 F\left(x,\frac{1-t^2}{2x},t\right),
\]
where \(I_t(t^{2r})=1/(2r+1)\) and \(I_t(t^{2r+1})=0\).
A phase is **pure** if \(\Lambda(f^m)=0\) for every integer \(m\geq1\).
It is eventually pure if these equalities hold for all sufficiently
large \(m\). It is safe against every ambient marker if, for every
\(g\in A_K\), \(\Lambda(gf^m)=0\) for all sufficiently large \(m\).

Consider exactly the eight-dimensional space
\[
\begin{split}
 f={}&a_1y+a_3y^3+D y^5+b_1yt^2+C y^3t^2+Ayt^4\\
    &\quad +x t^2(q_0+B t^2)
       =P(y,t)+xt^2(q_0+B t^2).                         \tag{1}
\end{split}
\]
Thus \(P\) is odd in \(y\), even in \(t\), and has total degree at most5.
This is a subspace of the even14-dimensional stripped quintic space,
not that entire space.

The normal monomials \(y^jt^r\) (\(j\geq0\)) and \(x^it^r\)
(\(i\geq1\)) are linearly independent. Their filtration degrees are the
displayed total degrees: their highest parts are the same normal
basis modulo the homogeneous relation \(2xy+t^2\), using \(xy\) as
leading monomial. Consequently (1) has least degree at most3 when its
four degree5 coefficients \(A,B,C,D\) vanish, and otherwise has least
degree5. There is no separate degree2 or degree4 stratum in (1).
Its angular antipode sends \(f\) to \(-f\), so every odd moment is zero.

## 1. A local cubic boundary lemma

**Lemma 1.** Over every characteristic-zero field, if
\[
 h=a y+c y^3+b yt^2+q xt^2,\qquad q\ne0,                \tag{2}
\]
and its moments of orders2,4,6 all vanish, then \(a=b=c=0\).
In particular every pure member of (2) is positive one-sided.

For nonnegative \(r,s\), the native trace is
\[
 \Lambda(x^ry^rt^{2s})
 =\frac{2^r r!(2s)!(r+s)!}{s!(2r+2s+1)!}.              \tag{3}
\]
It is zero for unequal \(x,y\) exponents or odd \(t\) exponent.
Formula (3) follows by substituting \(xy=(1-t^2)/2\), expanding the
finite polynomial, and integrating; equivalently it is the beta
integral identity, an identity between rational numbers.

The entire finite expansion needed below uses no stored moments.
Divide (2) by \(q\), and temporarily write its coefficients as \(a,b,c\),
so its positive coefficient is1. At order \(2k\), select \(j\) copies of
\(cy^3\), \(n_b\) copies of \(byt^2\), \(n_a=k-2j-n_b\) copies of \(ay\),
and \(k+j\) copies of \(xt^2\). Angular balance forces exactly this.
Hence, for this normalized phase,
\[
 \Lambda(h^{2k})=
 \sum_{\substack{0\leq j\leq\lfloor k/2\rfloor\\
                 0\leq n_b\leq k-2j}}
 \frac{(2k)!a^{k-2j-n_b}b^{n_b}c^j}
 {(k+j)!(k-2j-n_b)!n_b!j!}\,
 \Lambda(x^{k+j}y^{k+j}t^{2(k+j+n_b)}).                \tag{4}
\]
Use (3) with \(k=1,2,3\). In the original unnormalized coefficients,
\[
 M_2=\frac{2q}{105}(7a+3b).
\]
Thus \(a=-3b/7\). On this relation,
\[
 M_4=\frac8{3003}cq^3+\frac{1088}{735735}b^2q^2,
\]
forcing \(c=-136b^2/(245q)\). On both relations,
\[
 M_6=-\frac{2304}{554498945}q^3b^3.                    \tag{5}
\]
Characteristic zero and \(q\ne0\) imply \(b=0\), then \(a=c=0\).
This proves the lemma. It uses no winding-one classification: \(y^3\)
is not in the space \(U(t)+xS(t)+yR(t)\).

## 2. Precisely scoped imported quintic inputs

The following established inputs are used only in their displayed domains.

**N (Newton necessity).** An even-in-\(t\), pure stripped
phase of actual degree5 with both positive and negative angular terms
has \(B\ne0\) and at least two of \(A,C,D\) nonzero. This is the even
consequence of Theorem2 in *Stripped quintic Newton face reduction*,
the node SNF-EVEN-TWO-NEGATIVE-TOPCOEFFICIENTS-NECESSARY.
It includes the positive-height2, degree-\(P\)-at-most4, and sole
negative degree5 leader exclusions. It is not a theorem about all
phases of degree at most4.

**BC14 (finite obstruction).** Every member of (1) with \(BC\ne0\)
has at least one nonzero moment among orders2,4,...,14.
This is Corollary2 of *Sharp finite-moment obstruction for
odd-angular quintics*. Its scope is pure, not eventual; it does not
include the \(B=0\) or \(C=0\) boundary.

The complementary finite certificate covers the remaining boundary.

**C0-12 (finite obstruction).** Over every
characteristic-zero field, for every choice of five parameters
\(a_1,a_3,b_1,q,\lambda\), the phase
\[
 a_1y+a_3y^3+y^5+b_1yt^2+\lambda yt^4+xt^2(q+t^2)       \tag{6}
\]
has at least one nonzero moment among orders2,4,...,12.
This is Theorem1 of *Sharp finite-twelve obstruction for the C=0
odd-angular quintic*. Its complete chart proof covers every parameter,
including all pivot-zero and \(\lambda=0\) fibers. It is a finite
pure obstruction, not an eventual-purity theorem.

## 3. The remaining normalization has degree six

**Lemma 2.** No phase (1) with \(C=0\)
and \(BD\ne0\) is pure.

Adjoin a sixth root \(\rho\) of \(D/B\), and set \(\mu=1/(\rho B)\).
The actual quadric automorphism
\[
 x\longmapsto\rho x,\qquad y\longmapsto\rho^{-1}y,
 \qquad t\longmapsto t
\]
preserves \(\Lambda\) by angular constant terms. Apply it to \(f\),
then multiply by \(\mu\). The coefficients of \(xt^4\) and \(y^5\)
both become1, since
\[
 \mu\rho B=1,\qquad \mu\rho^{-5}D=D/(\rho^6B)=1.
\]
The result is (6). Its normalized parameters are
\[
 \left(\frac{a_1}{\rho^2B},\frac{a_3}{\rho^4B},
       \frac{b_1}{\rho^2B},\frac{q_0}{B},
       \frac{A}{\rho^2B}\right)
\]
in the order \((a_1,a_3,b_1,q,\lambda)\).
Every transformed moment is \(\mu^m\Lambda(f^m)\).
C0-12 makes one nonzero. Injectivity of the field extension makes
the corresponding original moment nonzero. No sixth root in \(K\),
ordering, or complex embedding is assumed.

## 4. Exact pure classification and ambient-marker safety

**Theorem 3.** For (1) over any
characteristic-zero field, the following are equivalent:
\[
 \begin{array}{ll}
 \text{(i)}& f\text{ is pure};\\
 \text{(ii)}& P=0\ \text{or}\ (q_0=B=0).
 \end{array}                                             \tag{7}
\]
Every pure phase in this space is safe against every fixed marker in
the FULL algebra \(A_K\). More precisely, for every nonzero \(g\in A_K\),
\[
 \Lambda(gf^m)=0\qquad(m>\deg g).                        \tag{8}
\]

Suppose \(P\ne0\) and \((q_0,B)\ne(0,0)\).
If \(A=B=C=D=0\), the phase has degree at most3. Then \(q_0\ne0\),
and Lemma1 excludes purity.

Otherwise the phase has actual degree5. N gives \(B\ne0\) and
at least two nonzero coefficients among \(A,C,D\). If \(C\ne0\),
BC14 gives a contradiction for every value of \(D\). If \(C=0\),
the two-leader condition forces both \(A,D\ne0\), so Lemma2 gives
a contradiction. This exhausts mixed phases and proves
(i) implies(ii).

For the converse and marker statement, the nonzero positive part
has angular support \(\{1\}\) after Laurent substitution.
A nonzero \(P\) has angular support contained in \(\{-1,-3,-5\}\).
A one-sided \(f^m\) therefore has all weights at least \(m\), or
all weights at most \(-m\). A marker of least sphere degree \(d\)
has a Laurent representative with weights in \([-d,d]\), as is seen
by substituting each monomial of a degree-\(d\) representative.
For \(m>d\), no term of \(gf^m\) has weight zero, proving (8).
Taking \(g=1\) proves all positive moments vanish. The zero phase
and zero marker are immediate.

The necessity proof uses all-order purity. It does not give a uniform14
nonvanishing cutoff on the whole mixed eight-dimensional space:
N is an analytic necessary-condition theorem without a supplied bounded prefix.
Nor does it classify eventual purity; a finite early nonzero moment
does not exclude eventual vanishing. The safety direction here is
pure implies safety, not a converse from eventual safety to purity.

## 5. Exact dependencies and scope

The local cubic lemma, normalization, reduction and full-marker
support argument are supplied above. The three imported results are
N, BC14 and C0-12 with exactly their displayed domains. In
particular the cubic boundary is handled locally, without extending
the winding-one theorem to a \(y^3\) term. The necessity proof uses
neither an ordinary-prime assertion for a generic sixfold nor a global
purity-recognition cutoff. The Newton necessary-condition theorem
and the two exact finite chart theorems are separately scoped inputs;
the balanced-Prym D0 obstruction is not a premise.

This concerns eight coefficients, not the other six even stripped
coefficients, all24 stripped coefficients, the36-dimensional full
degree-at-most5 sphere space, or arbitrary primitive quintics.
It gives no global degree5 existence or nonexistence result.
