# Dyadic isolation of split real-cubic Frobenius fields

Let \(A/K\) be a threefold with full geometric rational endomorphism
algebra a totally real cubic field \(F\), acting over \(K\), and let
\(L/\mathbb Q\) be the normal closure of \(F\). Work at good absolute
degree-one places \(v\) of odd characteristic \(p\), split completely
in \(F\), with exactly one failed coefficient factor. Omit the finite
set needed for integral real multiplication and a polarization.
For the relative trace \(a_v\) put
\[
 D_{v,i}=\sigma_i(a_v)^2-4p,\qquad
 d_v=-\frac{\prod_i D_{v,i}}{p},\qquad
 M_v=L(\sqrt{D_{v,1}},\sqrt{D_{v,2}},\sqrt{D_{v,3}}).
\]
Here \(i\) labels the three embeddings of \(F\) into \(L\).
The [singleton inertia calculation](CM_NEWTON_FROBENIUS_RAMIFICATION_AND_DENOMINATORS.md#inertia-and-the-frobenius-weyl-group)
shows that \(M_v\) is the full Frobenius splitting field, and that
every singleton sign change is realized by inertia at some place
above \(p\). The Weil bounds and the valuations at \(p\) give
\(0<d_v\le64p^2\) and \(p\nmid d_v\).

## Isolating the residue characteristics

**Proposition (Dyadic Frobenius-field isolation).** In every
interval \(X/2<p\le X\), choose at most one of the indicated places
above each rational prime, and call their number \(n\). At least
\(\lceil n/5\rceil\) of their fields are mutually linearly disjoint
over \(L\) for all sufficiently large \(X\), with every selected
characteristic unramified in all the other selected fields.
After discarding a set of places of density zero in the split
degree-one class, one can instead take at least \(\lceil n/3\rceil\)
of the remaining fields. In this reduced set an edge of
cross-ramification \(p\to q\) has \(v_q(d_v)=1\) and
\(\left(\frac p q\right)=1\). Thus two characteristics congruent to
three modulo four cannot ramify in each other's fields.

An odd prime \(q\ne p\) splitting completely in \(F\) can
ramify in \(M_v/L\) only if \(q\mid d_v\). In fact \(q\) also
splits in \(L\), and a ramified quadratic generator at a prime
over \(q\) needs odd valuation there. Put a directed edge
\(p\to q\) whenever \(q\mid d_v\) for two
distinct characteristics in the interval. There are at most
two outgoing edges for sufficiently large \(X\) because three
prime factors greater than \(X/2\) would have product exceeding
\(64X^2\). The underlying graph has at most \(2k\) edges
on any induced subgraph of \(k\) vertices. Repeatedly deleting
a vertex of degree at most four and coloring in reverse order
gives five colors. The largest color class gives the first bound.
The singleton inertia at each selected characteristic provides
three independent valuation tests on its own generators, all
vanishing on the other fields' generators. Thus the joint
extension of \(L\) has degree \(8\) to the power of the number
selected.

If \(d_v\) is divisible by a product of two rational primes in
\((X/2,X]\), counted with multiplicity, the remaining integer is
at most \(256\). Such \(d_v\) is prime to all fixed rational
primes greater than \(256\) once \(X\) exceeds twice those primes.
For a sufficiently large split auxiliary prime \(\ell\),
the condition \(\ell\mid\prod_i D_{v,i}\) has probability
\(\gg1/\ell\) in the full common-determinant residual image.
For example the determinant is square with probability \(1/2\)
and then a single factor has either of two trace values making
its discriminant zero, each of probability \(\asymp1/\ell\).
The corresponding events are independent for any fixed set of
sufficiently large auxiliary primes by product image and
Chebotarev. The sum of reciprocals of auxiliary primes split
in \(L\) diverges. Taking larger fixed sets shows that the
two-large-prime-factor locus has density zero, including squares
of one such prime. The remaining graph has outdegree at most one
and the same coloring argument now uses three colors.
[Effective and qualitative auxiliary product distributions](CM_NEWTON_FROBENIUS_DISCRIMINANT_SIEVE.md#effective-auxiliary-distributions-on-the-identity-class).

The dyadic exclusion applies in natural density, since it holds
on every sufficiently large interval \((X/2,X]\). Equivalently
one can discard the locus where \(d_v\) has two prime factors
at least \(p/2\) and at most \(2p\), counted with multiplicity.
On the remaining set an edge gives \(v_q(d_v)=1\), so some
quadratic generator is indeed ramified above \(q\).
The [quadratic norm identity](CM_NEWTON_SPLIT_COFACTOR_NORM_AND_TRACE_COUNT.md#the-cofactor-as-a-short-quadratic-norm)
gives \(d_v=V^2-pU^2\) with rational integers \(U,V\).
If \(q\mid U\), then \(q\mid V\) and \(q^2\mid d_v\).
Hence \(q\nmid U\), and \(p\equiv(V/U)^2\bmod q\).
Quadratic reciprocity proves the assertion about two-cycles.

For comparison, Burungale–Hida–Lai bound the number of places
with a *prescribed* Frobenius field, uniformly in that field under
GRH and their connectedness hypothesis. The argument here uses
the residue-characteristic ramification of this specific Newton
pattern to separate fields within growing groups of reductions.
[Burungale–Hida–Lai, Theorem 1.3](https://arxiv.org/pdf/2402.07935).

## Adjoining the normalizing square roots

**Proposition (Exact Kummer rank for isolated characteristics).**
Let \(\mathcal P\) be a finite set of the indicated places with
distinct characteristics. Suppose \(p_w\nmid d_v\) for distinct
\(v,w\in\mathcal P\). Write \(M_{\mathcal P}\) for the compositum
of the \(M_v\) and put \(C_{\mathcal P}=L(\sqrt{-d_v}:v\in\mathcal P)\).
Then
\[
 \begin{split}
 M_{\mathcal P}(\sqrt{p_v}:v\in\mathcal P)&=M_{\mathcal P}C_{\mathcal P},\\
 M_{\mathcal P}\cap C_{\mathcal P}&=L,\qquad
 [M_{\mathcal P}C_{\mathcal P}:L]=8^{|\mathcal P|}[C_{\mathcal P}:L].
 \end{split}
\]
The field \(C_{\mathcal P}\) is precisely the largest subextension
unramified at all places over the selected characteristics.
In particular the augmented fields \(M_v(\sqrt{p_v})\) are
mutually linearly disjoint over \(L\) exactly when the classes
of \(-d_v\) are linearly independent in \(L^\times/L^{\times2}\).

The product relation \(\prod_i D_{v,i}=-p_vd_v\) gives the
compositum equality. At a place of \(L\) above \(p_v\), each
cofactor class has zero valuation and each \(D_{w,i}\) for
\(w\ne v\) has zero valuation. The valuation vectors of \(D_{v,i}\)
run over the three singleton vectors as the place varies. Thus
a squareclass relation between all these generators and the
cofactor classes has all \(D\)-exponents zero. This proves
disjointness and the degree formula. It also proves the maximal
unramified statement since these are elementary abelian
quadratic extensions and all selected residue characteristics
are odd. Each augmented field has degree \(16\) over \(L\):
its additional cofactor \(-d_v\) is totally negative in \(L\),
so is not a square.

For the level-\(97\) reductions at \(1583\) and \(1709\), the
raw cofactors are \(7\cdot1567\cdot9689\) and
\(7\cdot19\cdot2184571\). Neither characteristic divides the
other cofactor. The primes \(1567\) and \(19\) give independent
squareclass tests on the two cofactors over \(\mathbb Q\).
They remain independent over \(L\), which in this example is
cyclic cubic. Thus the augmented compositum has degree \(256\)
over \(L\).
[Exact level-\(97\) cofactor calculations](CM_NEWTON_SPLIT_COFACTOR_NORM_AND_TRACE_COUNT.md#exact-checks).

## Growing disjoint subfamilies of augmented fields

**Theorem (Dyadic normalized Frobenius-field independence).**
Fix \(\delta>0\). For every sufficiently large \(X\) such that
at least \(\delta X/\log X\) places in \(X/2<p\le X\) have the
split one-failure pattern, one can choose
\(\gg_{A,K,\delta}\log\log X\) of their augmented fields
\(M_v(\sqrt{p_v})\) mutually linearly disjoint over \(L\).
Under GRH for the Dedekind zeta functions of \(LK(A[m])\) and
their subfields for all positive integers \(m\), the lower bound
improves to \(\gg_{A,K,\delta,\alpha}(\log X)^\alpha\) for every
\(1<\alpha<76/75\).

First use the dyadic graph to obtain an isolated subfamily of
at least \(c_K\delta X/\log X\) places with distinct characteristics
and \(p_w\nmid d_v\) for every pair of different places.
The five-color bound suffices after keeping at most one place
above each rational prime. In this subfamily write
\(c_v=\operatorname{sf}(d_v)\) for the positive squarefree part. The
[unconditional cofactor sieve](CM_NEWTON_UNCONDITIONAL_COFACTOR_SIEVE.md#a-squareclass-test-that-includes-zero-discriminants)
bounds the size of each prescribed cofactor fibre by
\[
 O_{A,K,\theta}\!\left(
 \frac{X}{\log X}\frac{\log\log X}{(\log X)^\theta}\right),
 \qquad 0<\theta<1/18.
\]
The local test \(U_\ell(c)\) is one on the designated fibre
including zero discriminants, and its second moment, including
the exceptional quadratic character term, is uniform in \(c\).
The small cofactor restriction in that sieve is used only when
summing over the possible classes. Similarly the
[GRH cofactor sieve](CM_NEWTON_FROBENIUS_DISCRIMINANT_SIEVE.md#a-square-sieve-for-the-cofactor)
gives a uniform bound \(O_{A,K,\theta}(X^{1-\theta})\) per fibre
for \(0<\theta<1/38\), even when \(c\) reaches the full Weil
size \(64X^2\). Indeed for each place only \(O_\theta(1)\)
auxiliary primes \(\ell\asymp X^\theta\) divide
\(c\prod_i D_{v,i}\), a nonzero integer of size \(O(X^5)\).
The loss in the square-sieve sum is therefore uniform in \(c,v\).

The unconditional fibre bound leaves
\(\gg_{A,K,\theta,\delta}(\log X)^\theta/\log\log X\) distinct
cofactor classes in the isolated subfamily. Over \(\mathbb Q\)
their span has rank at least the base-two logarithm of this
number. Restriction to \(L^\times/L^{\times2}\) lowers rank by
at most a fixed constant since its kernel consists of the
squareclasses defining quadratic subfields of \(L\). The exact
Kummer rank formula then gives \(\gg_{A,K,\delta}\log\log X\)
independent augmented fields.

Under GRH the number of distinct rational classes \(-c_v\)
is \(\gg_{A,K,\theta,\delta}X^\theta/\log X\).
Let \(r\) be their rank. A row basis of their sign and prime
valuation parity vectors has \(r\) pivot coordinates whose
projection is injective on their span. Apart from a possible
sign coordinate the pivots are distinct primes \(q_1,\ldots,q_j\)
with \(j\le r\), and the product of the selected primes in a
projected vector divides \(c_v\le64X^2\). For \(0<s<1\), the
number of such projected vectors is at most
\[
 2(64X^2)^s\prod_{i=1}^j(1+q_i^{-s})
                 \le X^{2s}\exp(O_s(1+r^{1-s})).
\]
The first bound weights each admissible subset by
\(((64X^2)/\prod q_i)^s\) and accounts for the sign by a factor
of two. The second orders the pivot primes and uses \(q_i\ge i+1\).
If \(r\le(\log X)^\alpha\) with \(\alpha(1-s)<1\), this gives
at most \(X^{2s+o(1)}\) distinct classes. For
\(1<\alpha<76/75\), choose \(\theta<1/38\) and \(s\) with
\(1-1/\alpha<s<\theta/2\). This contradicts the lower bound
on their number. Subtracting the fixed rank loss over \(L\)
gives the asserted size of the augmented subfamily. This uses a
positive fraction of one-failure places as input and gives no
upper bound on that fraction.
