# A conditional Laurent-coefficient invariance lemma

This note independently proves a scalar Laurent-series step relevant to
Lemma 6.3 of Zheglov's
[arXiv:2410.06959v5](https://arxiv.org/html/2410.06959v5). It repairs the
raw/compressed exponent convention and proves invariance of all the
distinguished coefficients under explicit hypotheses. It does **not**
assert that the earlier head-order and correspondence constructions
supply those hypotheses for every hypothetical DC pair. The main DC1
claim remains unverified in the
[primary-source audit](ZHEGLOV_DC1_PRIMARY_SOURCE_AUDIT.md).

The exact [cached v5 PDF](sources/zheglov-2410.06959v5.pdf), dated
19 January 2026, has 78 pages and SHA256
`897913ef95e011332e86e348803c6c9da5dbd196db17a04f8f951110e9b095e3`.
The relevant primary locations are Remark 5.2, Section 6.2 Step 3,
and Lemma 6.3 on printed pages 70--71.

## Exact hypotheses

Let \(k\) be a characteristic-zero field and let \(t\) be a commuting
Laurent variable. The order of a nonzero element of
\(k((t^{-1}))\) is its highest exponent. Fix integers
\[
 p\geq2,\quad d_0\geq2,\quad 1\leq l_0<d_0,\quad
 \gcd(d_0,l_0)=1,\quad I\geq1,
\]
where \(d_0\mid p\). Fix nonnegative integers \(q_0,\ldots,q_{I-1}\),
each divisible by \(d_0\), and constants
\(\epsilon_0,\ldots,\epsilon_{I-1}\in k^\times\), such that
\[
                 q_{j+1}<p q_j\qquad(0\leq j<I-1).       \tag{1}
\]
Define fixed polynomials in commuting variables by
\[
 H_0(X,Y)=Y,\qquad
 H_{j+1}(X,Y)=H_j(X,Y)^p-\epsilon_j^p X^{q_j}.             \tag{2}
\]
The source application has \(d_0=d_2\), \(l_0=l_2\), and uses the
fixed common chopping exponent \(p\). Its usual monic initial head has
\(\epsilon_0=1\); the lemma does not require that normalization.

For an integer \(N>1\), put
\[
 D=d_0N+l_0,\quad \sigma=D/d_0,\quad p_N=\sigma p,
 \quad q_{N,j}=\sigma q_j\quad(j<I),\quad P_N=t^{p_N}.
                                                               \tag{3}
\]
In particular \(D>N+1\), and all exponents in (3) are integral.
Let \(Y_N\in k((t^{-1}))\) satisfy **all** the following conditions:
\[
 H_j(P_N,Y_N)=\epsilon_jt^{q_{N,j}}+\text{lower terms}
                         \qquad(0\leq j<I),              \tag{4}
\]
\[
 \begin{split}
 L_{N,k}&=(p-1)\sum_{j<k}q_{N,j},\qquad L_{N,0}=0,\\
 q_{N,I}&:=\operatorname{ord}H_I(P_N,Y_N)
                  =L_{N,I}-p_N+N+1,\\
 [t^{-p_N}]Y_N&=0.                                      \tag{5}
 \end{split}
\]
The order and integrality hypotheses imply exactly one of the alternatives
\[
 \delta_N:=p q_{N,I-1}-q_{N,I}>N+1,                      \tag{6a}
\]
or
\[
 d_0=2,\ l_0=1,\ C=1,\ \delta_N=N,\qquad
 C=\frac{p q_{I-1}-(p-1)\sum_{j<I}q_j+p}{d_0}.
                                                               \tag{6b}
\]
In the exceptional case (6b), additionally require
\([t^{q_{N,I}-1}]H_I(P_N,Y_N)=0\). This zero coefficient is the only
additional hypothesis here; it is not inferred from the head orders.
There is no separate gap assumption, and the boundary
\(\delta_N=N+1\) cannot occur.

To prove the dichotomy, (3) and (5) give
\(\delta_N=DC-N-1\), where \(C\) is an integer because \(p\)
and every \(q_j\) are divisible by \(d_0\). The shared heads
imply \(\delta_N>0\), as also shown by (13)--(15) below, so
\(C\geq1\). For \(C\geq2\),
\(\delta_N\geq2D-N-1\geq3N+1>N+1\). If \(C=1\) and
\(d_0\geq3\), then \(\delta_N\geq2N>N+1\), since \(N>1\).
The sole remaining case is \(d_0=2,l_0=1,C=1\), giving
\(\delta_N=N\).

## The invariant coefficients and their common head root

For each \(0\leq k\leq I\), the coefficient
\[
                 C_{N,k}=[t^{L_{N,k}-p_N}]H_k(P_N,Y_N)   \tag{7}
\]
depends only on the fixed data (1)--(2), independently of \(N\) and
of the remaining coefficients of \(Y_N\) satisfying (4)--(5) and
the exceptional zero-coefficient condition just stated.
This is a written all-\(N\) conclusion, not an extrapolation from the
finite controls below.

Construct Laurent series \(A_j(z)\) backwards by
\[
 A_{I-1}(z)=\epsilon_{I-1}z^{q_{I-1}/d_0},\qquad
 A_j(z)=\epsilon_jz^{q_j/d_0}
 \left(1+\frac{A_{j+1}(z)}
                  {\epsilon_j^p z^{p q_j/d_0}}\right)^{1/p}.
                                                               \tag{8}
\]
The unit binomial expansions exist by (1) in characteristic zero.
Put \(W=A_0\), \(A_I=0\), and
\[
 \lambda_k=p^k\prod_{j<k}\epsilon_j^{p-1},\quad
 \lambda_0=1,\quad
 \ell_k=(p-1)\sum_{j<k}q_j/d_0.                          \tag{9}
\]
The exact invariant is
\[
 C_{N,k}=[z^{\ell_k-p/d_0}]A_k(z)
                    -\lambda_k[z^{-p/d_0}]W(z).        \tag{10}
\]
For \(k=0\), this is zero, as required by (5). For the final index,
\[
                    C_{N,I}=-\lambda_I[z^{-p/d_0}]W(z). \tag{11}
\]
These formulas require no integer-polynomial assertion about the
\(\epsilon_j\). The binomial construction gives fixed rational
functions of these constants; the stronger source wording about
integer polynomial coefficients is not imported here.

### Proof: the linear term

Write \(W_N(t)=W(t^D)\). Construction (8) gives
\[
 H_j(P_N,W_N)=A_j(t^D),\quad j\leq I,\qquad
                         H_I(P_N,W_N)=0.                \tag{12}
\]
Let \(E=Y_N-W_N\). Factor the difference of the two values of
\(H_I\) using (2). Each factor at stage \(j\) has leading term
\(p\epsilon_j^{p-1}t^{(p-1)q_{N,j}}\), because both values have
the same head in (4) and (12). Consequently
\[
                  \operatorname{ord}E=\gamma_N=-p_N+N+1.
                                                               \tag{13}
\]
The derivative
\(J_{N,k}=\partial_YH_k(P_N,W_N)\) has leading term
\(\lambda_kt^{L_{N,k}}\). Every exponent in its support belongs
to \(D\mathbb Z\), since (12) has that property. Its next possible
exponent is at most \(L_{N,k}-D\). Because \(D>N+1\), the
coefficient of \(t^{L_{N,k}-p_N}\) in \(J_{N,k}E\) is exactly
\[
 \lambda_k[t^{-p_N}]E
                   =-\lambda_k[z^{-p/d_0}]W.             \tag{14}
\]
Lower derivative coefficients would require an exponent of \(E\)
strictly greater than (13). The baseline term in (12) contributes
the first term of (10).

### Proof: all nonlinear terms

Expand as a polynomial in an auxiliary commuting variable \(U\):
\[
 H_j(P_N,W_N+U)-H_j(P_N,W_N)=\sum_m K_{j,m}U^m.
\]
Give \(U\) the order \(\gamma_N\), and set
\[
                  g_j=q_{N,j}-\gamma_N-L_{N,j}\quad(j<I).
\]
Factoring differences as above shows that \(g_j>0\). Moreover
\[
 g_j-g_{j+1}=p q_{N,j}-q_{N,j+1}>0,\qquad
                         g_{I-1}=\delta_N.              \tag{15}
\]
Induction through the powers in (2) gives, for \(j\geq1,m\geq1\),
\[
 \operatorname{ord}K_{j,m}+m\gamma_N
       \leq\gamma_N+L_{N,j}-(m-1)g_{j-1}.                \tag{16}
\]
For the first step this is the binomial expansion of
\((W_N+U)^p\). At a later step, a product of \(r\) variations
with total \(U\)-degree \(m\) has order at most
\(p q_{N,j}-r g_j-(m-r)g_{j-1}\). Since
\(g_{j-1}>g_j\), this is at most
\(p q_{N,j}-m g_j\), the next instance of (16).

Under (6a), (15)--(16) put all terms with \(m\geq2\), for every
prefix \(k\geq1\), strictly below \(L_{N,k}-p_N\). This proves
(10) in the nonresonant case.

Under (6b), each prefix \(k<I\) has
\(g_{k-1}\geq N+D>N+1\), since every earlier difference in
(15) is a positive multiple of \(D\). Its nonlinear terms are
therefore again too low. At the final stage, terms with \(m\geq3\)
have order at most \(q_{N,I}-2N<L_{N,I}-p_N\), since \(N>1\).
The quadratic term has order at most
\(r_N+1\), where \(r_N=L_{N,I}-p_N\).

Every coefficient \(K_{I,2}\) has exponents in \(D\mathbb Z\).
Its possible top exponent
\[
              r_N+1-2\gamma_N=L_{N,I}+p_N-D
\]
also belongs to \(D\mathbb Z\), because \(D=2N+1\). Hence
the quadratic contribution at \(r_N\) can only use the coefficient
of \(E^2\) at \(2\gamma_N-1\), which is
\(2E_{\gamma_N}E_{\gamma_N-1}\).

The linear contribution at \(q_{N,I}-1\) is
\(\lambda_IE_{\gamma_N-1}\); derivative coefficients below the
head have a gap at least \(D\). All nonlinear terms have order at
most \(q_{N,I}-N<q_{N,I}-1\). The zero coefficient in (6b) thus
forces \(E_{\gamma_N-1}=0\). The dangerous quadratic contribution
vanishes. Equations (14) and (12) now prove (10) also in this case.

## Three precise normalization controls

These elementary statements explain how the trace and missing-head
conditions can follow **if** the relevant earlier constructions exist.
They do not establish existence of those constructions.

**Puiseux trace.** Suppose a monic polynomial in \(Y\), of degree
\(p_N\) and with coefficients in \(k[P_N]\), has exactly the
branches \(Y_N(\zeta t)\), counting multiplicity, as \(\zeta\)
runs over the \(p_N\)-th roots of unity after constant extension.
Their sum is polynomial in \(P_N\) by Vieta. The coefficient of
\(P_N^{-1}\) in the sum is
\(p_N[t^{-p_N}]Y_N\), so it is zero. This is the trace input in
Remark 5.3. The asserted rank-one commuting-curve construction that
would supply this branch hypothesis is not audited by this lemma.

**Cyclic-vector scalarization.** In the skew Laurent ring with
coefficient ring \(k^s\), let \(\tau\) cyclically permute the
components and \(D a=\tau(a)D\). Suppose
\[
 L=\epsilon D^q+0D^{q-1}+\text{lower},\qquad
 S=1+s_1D^{-1}+\text{lower},\qquad \epsilon\in k^\times.
\]
The coefficient at \(D^{q-1}\) in \(S^{-1}LS\) is
\(\epsilon(\tau^q(s_1)-s_1)\). Its component sum is zero.
If this coefficient is scalar, all \(s\) entries are equal, so
characteristic zero forces the scalar to be zero. Thus a claimed
monic auxiliary scalarizing conjugation preserves the missing
head-minus-one coefficient. The cyclic-vector model and existence
of that conjugation remain separate hypotheses of this application.

**Filtered Schur and tail steps.** In an algebra with an additive
order filtration, conjugation by \(1+O(D^{-N})\) changes an element
of order \(q\) only in orders at most \(q-N\), hence preserves its
head-minus-one coefficient for \(N>1\). If the initial second
operator is changed by a tail of order at most \(-p_N\), telescoping
the powers in (2) bounds the change in its \(k\)-th chopped operator
by \(L_{N,k}-p_N\). The heads are unchanged: at every earlier
stage this bound lies below the head by \(g_j+N+1>0\).
At the final stage the bound is \(q_{N,I}-N-1\), so the missing
head-minus-one coefficient is also preserved. This statement requires
the shared heads (4), order-submultiplicativity, and the same fixed
first operator; it is not a claim that an arbitrary completion or
conjugation satisfies them.

## The repaired indexing and the source distinguished exponents

The telescoping identity
\[
 \sum_{j<k}(p q_{N,j}-q_{N,j+1})
                 =L_{N,k}+q_{N,0}-q_{N,k}
\]
shows that the source's distinguished exponent \(r(k)\) is exactly
\(L_{N,k}-p_N\). Formula (10) therefore treats all its prefixes,
including the exceptional final index.

The source's raw Laurent indices scale by
\(\sigma=(d_2N+l_2)/d_2\), whereas its displayed coefficient formula
uses \(d_2N+l_2\). The consistent compressed formula is
\[
                  c_{N,Dj}=c_{d_0j}.                     \tag{17}
\]
For example, if an original series has the same head chain and first
deviation from its common head root at exponent \(-p+1\), (17)
holds for integers \(j\geq1-p/d_0\). Both sides then equal the
coefficient of \(z^j\) in \(W\). Equivalently, on the relevant raw
indices divisible by \(d_0\), use
\(c_{N,\sigma i}=c_i\), with \(i>-p+1\).
Both the coefficient index and any threshold must follow the chosen
convention. The proof above bypasses the ambiguous displayed formula.

## Exact finite controls and their limits

The [independent Fraction checker](check_conditional_laurent_invariance.py)
and [retained receipt](conditional_laurent_invariance_controls.json)
give two explicit families at \(N=3,5,9\). They are commutative
Laurent coefficient controls, **not** DC-pairs or ordinary differential
operator constructions. The all-\(N\) result is proved above.

For \(p=4,q_0=6,d_0=2,l_0=1,q_1=22\), the original head
distances \(2\) and \(7\) sum to \(p+q_0-1=9\). The compressed
root is \(W=z^3(1+z^{-1})^{1/4}\). Its coefficient at \(z^{-2}\)
is \(231/8192\), and \(\lambda_I=16\). The final invariant is
the nonzero number \(-231/512\). The prefix invariant is
\(-231/2048\). The controls truncate this root immediately before
its \(t^{-p_N}\) coefficient and add \(t^{\gamma_N}\). They have
zero trace coefficient and common final leading coefficient 16.

For the resonant data \(q_1=16\), the original distances are
\(8\) and \(1\). Use
\[
 Y_N=t^{3D}+\tfrac14t^{-D}
                     +t^{\gamma_N}+v t^{\gamma_N-1}.
\]
Its final head-minus-one coefficient is exactly \(16v\), while
its distinguished coefficient is exactly \(192v\). The required
missing coefficient therefore forces \(v=0\) and removes the
resonance. The controls retain nonzero \(v\) as a deliberate check
that omitting this hypothesis changes the conclusion.

The main theorem of the primary paper is not adopted. The unresolved
application inputs include production of the shared polynomial chain,
all scaled head orders and the final-order equality (5), the exceptional
missing-head normalization, the requisite Schur/auxiliary conjugations,
and the rank-one commuting curve with the stated branch/trace property.
The gap dichotomy itself follows from the retained lemma's order and
integrality hypotheses. No DC1 resolution, refutation,
formal verification, external review or novelty priority is claimed.
