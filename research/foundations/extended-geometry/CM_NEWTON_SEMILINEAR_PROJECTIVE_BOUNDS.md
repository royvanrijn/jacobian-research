# Projective trace bounds when Galois permutes the CM action

Let \(A/K\) be a polarized abelian variety over a number field of
dimension \(2e\). Suppose it is not potentially CM and that
\(\operatorname{End}^0(A_{\bar K})\) contains a CM field \(E\)
of degree \(2e\) on which Rosati acts by conjugation. Assume
\(G_K\) preserves \(E\) as a subalgebra. It may act nontrivially
on \(E\). Let \(F=E^+\), and use the
[connected projective trace field \(F_0\)](CM_NEWTON_PROJECTIVE_TRACE_FIELD.md)
and balanced-pair number \(b\) after a finite extension defining
the action. The projective classes are the fibres of restriction
of the embeddings of \(F\) to \(F_0\). Put \(d=[F:F_0]\),
so each class has size \(d\) and has a constant balanced or definite
signature type.

Let \(v\) be a good absolute degree-one prime of \(K\) of residue
characteristic \(p\), outside ramification in the field defining
the action and in \(E\), and away from the coefficient order
index. Let \(\gamma\in\operatorname{Aut}(E)\) be the permutation
of the field action induced by Frobenius, and \(h\) its order.
Choose the coefficient labels at a place above \(v\). Let
\(\theta\) be a Frobenius of \(p\) in the normal closure of \(E\)
for that choice. The \(\theta\)-orbits on \(E\)-embeddings
are their groups at places \(\lambda\) of \(E\) over \(p\).
Write \(f_\lambda=[E_\lambda:\mathbb Q_p]\) and
\(g_\lambda=\gcd(h,f_\lambda)\).
The split case below has \(f_\lambda=1\) and trivial \(\theta\).

**Theorem.** Suppose \(b>0\). Outside a set of natural density zero
over \(K\), at every such prime at which the oriented Hodge
multiplicity is constant within each coefficient place, the
nonordinary balanced projective classes form a union of orbits
under \(\gamma\) and \(\theta\). Give each balanced class \(C\)
the cost
\[
 w_C=\max_{\sigma\ {\rm over}\ C}
          \left\lceil\frac{h}{\max(2,g_{\lambda(\sigma)})}\right\rceil .
\]
Here \(\sigma\) ranges over the two oriented \(E\)-embeddings
of each conjugate pair in \(C\). The failing classes satisfy
\[
 \sum_{C\ {\rm failing}} w_C
       \ \le\
       \left\lfloor\frac{h(b/d)-1}{2}\right\rfloor .
\]
Each failed class costs \(2d\) in the \(p\)-rank. If each joint
orbit of balanced projective classes has cost exceeding this
bound, the prime is ordinary.

In particular, at split primes the number \(r\) of failed balanced
pairs satisfies
\[
 r\ \le\
 d\left\lfloor\frac{h(b/d)-1}{2\lceil h/2\rceil}\right\rfloor .
\]
A transitive Frobenius permutation on the balanced projective
classes gives ordinary reduction there. The same holds at
nonsplit primes satisfying the multiplicity condition and
\(f_\lambda\le2\) on the balanced classes, if the joint action
of \(\gamma\) and \(\theta\) is transitive. In the nonsplit assertion it
suffices that \(g_\lambda\le2\) on the balanced classes.
For \(h=1\) the cost
is one at all coefficient places and the bound recovers the
connected projective count, including nonsplit places. If \(b=0\),
homogeneity of the Hodge multiplicities in the coefficient
places gives ordinary slopes there directly.

## Monodromy on the projective classes

Let \(G\) denote algebraic monodromy for one auxiliary prime,
and work over its algebraically closed coefficient field.
The identity component \(G^\circ\) centralizes \(E\) and acts
on the rank-two spaces \(V_\sigma\). By the
[rank-two irreducibility argument](CM_NEWTON_PROJECTIVE_TRACE_FIELD.md#density-over-the-field-of-definition-of-the-action)
each such action is irreducible and its projective image is
\(\mathrm{PGL}_2\). The linked representations over each class
can be identified projectively: in a pair of rank-two factors
the connected projective image is a product or a graph, and
equality of the projective trace rules out the product.
Conjugate embeddings are identified by the polarization.
Choose a projective line \(W_C\) for each class with all these
identifications. The derived group of \(G^\circ\) surjects
projectively onto \(\prod_C\mathrm{PGL}(W_C)\). A connected
subdirect product of these simple groups is a product of
diagonals, by the simple Lie algebra decomposition. Any
diagonal would identify projective traces across distinct
classes.

An element \(g\in G\) permutes the \(E\)-embeddings and the
projective classes. For \(\sigma\) in \(C\) its projective
block map \(\mathbb P(V_\sigma)\to\mathbb P(V_{gC})\), in the
chosen identifications, intertwines the action of \(u\in G^\circ\)
on \(W_C\) with the action of \(gug^{-1}\) on \(W_{gC}\).
These representations are the same for all \(\sigma\) in \(C\).
The centralizer of \(\mathrm{PGL}_2\) in itself is trivial,
so the block map depends only on \(C\). Thus \(G\) acts by
permutations and projective block maps on these class lines.

The order \(h\) of the action on \(E\) is constant on each
component of \(G\). On a component of order \(h\), define
\[
 H(g)=\prod_{\sigma:E\hookrightarrow\overline{\mathbb Q}_\ell}
   \frac{\operatorname{tr}(g^h\mid V_\sigma)^2}
        {\det(g^h\mid V_\sigma)}.
\]
This regular class function is nonconstant on every component,
with its appropriate exponent \(h\). To see this, for a fixed \(g\)
write the projective block product around each cycle of class
lines. Since \(G^\circ\) has the full independent projective
product, multiplication of \(g\) by a suitable element of
\(G^\circ\) varies each return block arbitrarily, independently
across cycles. For a cycle of length \(a\), \(g^h\) acts
projectively in that cycle by conjugates of the \((h/a)\)-th
power of the block product. A block of eigenvalue ratio \(t\)
has projective trace \(2+t^{h/a}+t^{-h/a}\) after this power.
Choose the other cycles with nonzero projective traces and
vary \(t\). The resulting product of positive integer powers
of these functions is nonconstant.

## Strict norm bound on the Frobenius power

At a good degree-one prime with Frobenius \(g\), the power
\(\pi_v^h\) commutes with \(E\) and is the \(p^h\)-Frobenius of
the base change to \(\mathbb F_{p^h}\). Let \(\tau,\delta\in E\)
denote its relative trace and determinant. These are independent
of cohomology theory, by degree polynomial comparison with
\(a\pi_v^h\) and \(a\pi_v^{2h}\) for \(a\) in the coefficient
order, followed by Newton identities. All traces are integral.
Each complex determinant has absolute value \(p^h\), so
\(N_E(\delta)=p^{2eh}\) by CM conjugation. Hence
\[
       H(g)=\frac{N_E(\tau)^2}{p^{2eh}}\in\mathbb Q,\qquad
       0\le H(g)\le 4^{2e}.
\]
Nonnegative \(p\)-valuation would make this rational number
an integer, because its denominator is a power of \(p\).
The finitely many integral levels contain no component.
Algebraic Chebotarev discards them on a density-zero set,
including the trace-zero level. Thus
\(v_p(N_E(\tau))<eh\) at all remaining degree-one primes.
This use of algebraic Chebotarev is on the original field:
raising Frobenius to the component exponent is inside the
class function and does not change the primes being counted.
[Sawin, Theorem 2.1](https://www.numdam.org/item/10.1016/j.crma.2016.01.025.pdf).

At a coefficient place \(\lambda\), the \(p\)-divisible summand
has height \(2f_\lambda\). Its two relative slopes, with possible
equality, are repeated \(f_\lambda\) times on the underlying
\(p\)-divisible group. Indeed, after splitting \(E_\lambda\) in
the unramified coefficient extension of the isocrystal, the
semilinear Frobenius cycles its \(f_\lambda\) idempotents, so the
Newton data agree on them. Over \(\mathbb F_{p^h}\) the \(h\)-fold
crystalline Frobenius acts linearly and equals \(\pi_v^h\) on
the rational crystalline module over its base field. Its scalar
extension preserves the coefficient and slope decompositions.
Thus each relative block of \(\pi_v^h\) realizes \(h\) times
the two normalized slopes. The endpoint dimension on this summand
is \(f_\lambda m_\lambda\), where \(m_\lambda\in\{0,1,2\}\)
is the common oriented Hodge multiplicity in that place.
The determinant valuation in the relative component is
\(h m_\lambda\). All local trace valuations are integers since
\(E_\lambda/\mathbb Q_p\) is unramified.

On a definite component the slopes are both zero or both one.
On a balanced component they sum to one. If balanced and ordinary,
the trace is a unit because the root valuations are distinct,
\(0\) and \(h\). If balanced and not ordinary, the lower
normalized slope is at least \(1/\max(2,g_\lambda)\).
For equal slopes this follows from their sum. For distinct
slopes the smallest slope \(s\) has \(f_\lambda s\) integral
by the first Newton vertex at abscissa \(f_\lambda\). Also
\(h s=v_p(\sigma(\tau))\) is integral by the distinct root
valuations. Hence \(g_\lambda s\) is integral. Thus the trace
valuation in the nonordinary case is at least
\(\lceil h/\max(2,g_\lambda)\rceil\). The integral breakpoint
is on the underlying height-\(2f_\lambda\) component.

For each balanced class set \(u_C=v_p(\sigma(\tau))\), using
any \(\sigma\) in that class. It is independent of the choice
because the projective trace of \(\pi_v^h\) is constant there,
and the determinant valuation is \(h\). For a definite class
set \(u_C=v_p(\sigma(\tau))-hm_\lambda/2\), nonnegative by
the slopes. Again twice this value is the projective trace
valuation, which is constant. Under the identification of the
conjugate embeddings the normalized projective ratio is
\(T=\tau^2/\delta\in F\), with
\(N_F(T)=p^{-eh}N_E(\tau)\), by Rosati.
Consequently the strict norm inequality is
\[
                  -hb+2d\sum_C u_C <0.
\]
For a failing balanced class we have \(u_C\ge w_C\).
This proves the cost bound. The Frobenius endomorphism \(\pi_v\)
is a quasi-isogeny commuting with \(\pi_v^h\), and permutes
the coefficient pieces by \(\gamma\), so it preserves the
Newton pattern and the balanced type. The arithmetic action
\(\theta\) also preserves the slopes, by the isocrystal
cycling just used. Both act on the projective classes,
the former by the monodromy calculation and the latter by
restriction of field embeddings to \(F_0\).
This proves the orbit constraint and the \(p\)-rank assertion.
\(\square\)

The conditions give a finite Chebotarev criterion on any
combined splitting and endomorphism Frobenius types on which
the Hodge multiplicity condition and the cost test hold.
This extends the [cyclic RM generator-class argument](CM_NEWTON_DISCONNECTED_RM_ORDINARITY.md)
to geometric CM actions preserved as algebras, without requiring
a single cyclic permutation of all the oriented embeddings.
Cantoral Farfán and coauthors use a componentwise normalized
Frobenius test for simple signatures with their own
endomorphism and monodromy hypotheses. Here the projective
classes and their cycle costs arise before any density
calculation. In the independent split identity class,
\(h=1\) and both actions are trivial, so with three balanced
pairs the single-failure pattern still fits the bound.
[Cantoral Farfán et al., Proposition 4.8 and Theorems 4.18–4.19](https://arxiv.org/html/2508.11174v1).

## Averaged Hodge types and a single balanced pair

One can remove the homogeneity restriction by comparing to the
averaged Hodge polygon instead. The resulting conclusion is
\(\mu\)-ordinarity at density one over \(K\) when there is
exactly one balanced pair, without requiring the coefficient
action or the full geometric endomorphism algebra to be defined
over \(K\).

For a place \(\lambda\mid p\) write \(b_{\lambda,i}\) for the
number of its \(f_\lambda\) coefficient embeddings of oriented
Hodge multiplicity \(i\), for \(i=0,1,2\). With the convention
that \(m_\sigma\) is the Hodge endpoint, the relative
\(\mu\)-ordinary slope pair is the average of the sorted
Hodge pairs:
\[
 L_\lambda=\frac{b_{\lambda,2}}{f_\lambda},\qquad
 U_\lambda=\frac{b_{\lambda,1}+b_{\lambda,2}}{f_\lambda}.
\]
This is the usual \(\mu\)-ordinary polygon relative to the
geometric \(E\)-action defined over the residue extension of
degree \(h\). The formula also shows directly that the choice of
dual orientation just interchanges the slopes with their
complements.
[Cantoral Farfán et al., Proposition 3.6 and note 12](https://arxiv.org/html/2508.11174v1).

The actual \(p\)-divisible polygon lies above this reference.
One can see the needed rank-two comparison directly on its
Dieudonné lattice. After splitting the coefficient places,
semilinear Frobenius carries each height-two lattice to the next
one with elementary divisors in \(\{0,1\}\). Their sum is the
oriented Hodge multiplicity. In each of the \(b_{\lambda,2}\)
steps with sum two, the image lies in \(p\) times the target.
One complete cycle is therefore divisible by
\(p^{b_{\lambda,2}}\), giving smallest normalized slope
\(s_\lambda\ge L_\lambda\). The determinant has normalized
valuation \(L_\lambda+U_\lambda\), from the sum of the Hodge
multiplicities. If \(b_{\lambda,1}=0\), the two reference
slopes coincide, and equality of polygons follows.

For a conjugate projective class \(C\) put
\(\beta_C=b_{\lambda,1}/f_\lambda\) for any coefficient place
induced by an oriented embedding above it. This is well defined:
the oriented \(\theta\)-orbit above \(C\) projects onto the
\(\theta\)-orbit of \(C\), visiting each class there equally
often. The original balanced status is constant in each
projective class, so \(\beta_C\) is just the balanced fraction
in that orbit. In particular \(d\sum_C\beta_C=b\).

The arithmetic semilinear action on the Hodge space of the
special fibre permutes its labels by the composite of
\(\theta\) with precomposition by \(\gamma^{-1}\), using
arithmetic Frobenius for \(\gamma\). These two permutations
commute. This follows by applying arithmetic Frobenius
to the equation \(e x=\sigma(e)x\) for an eigenvector, using
its conjugation of the endomorphism \(e\) on the left and
its scalar Frobenius on the right. Thus \(\gamma\) preserves
the averaged multiplicities within the coefficient places.
Also the \(h\)-th power of this action preserves original
multiplicities along \(\theta^h\)-orbits. Those suborbits
have length \(f_\lambda/\gcd(f_\lambda,h)\). In particular
\(hL_\lambda\) and \(hU_\lambda\) are integers.

**Proposition.** Keep the hypotheses on \(A,E\) of the theorem,
and use any of the allowed good degree-one primes, including
nonhomogeneous Hodge multiplicities at its coefficient places.
At the primes outside the same density-zero exception, the
projective classes on which the \(\mu\)-ordinary polygon fails
form a union of joint \(\gamma,\theta\)-orbits. They have
\(\beta_C>0\). Give each such class the cost
\[
 W_C=\max_{\sigma\ {\rm over}\ C}
     \left\lceil h\min\!\left(
       \frac{1}{g_{\lambda(\sigma)}},\frac{\beta_C}{2}\right)\right\rceil .
\]
The total cost of the failing classes is at most
\(\lfloor(h(b/d)-1)/2\rfloor\) when \(b>0\). In particular
\(b=1\) implies \(\mu\)-ordinary reduction on a set of natural
density one over \(K\).

**Proof.** Keep the \(T=\tau^2/\delta\) from the strict norm
argument, for the original prime's \(h\)-th Frobenius power.
For every \(\sigma\) above \(C\),
\[
              u_C=\frac{v_p(\sigma(T))+h\beta_C}{2}
                  =v_p(\sigma(\tau))-hL_{\lambda(\sigma)} .
\]
The identity follows from the determinant valuation
\(h(L_\lambda+U_\lambda)\). This is a nonnegative integer
because the trace has valuation at least \(h s_\lambda\).
The first formula makes it constant on the class. The
projective ratios of Frobenius are constant under \(\gamma\)
by conjugation of the relative blocks, and under \(\theta\)
because embeddings in one unramified coefficient place give
the same valuation. The same invariance of \(\beta_C\) was
proved above. When \(\beta_C>0\), equality \(u_C=0\) holds
exactly for the reference slopes: the trace valuation is
exactly that of the smaller root when the slopes are
distinct. If the polygon fails, its smallest slope is at
least \(L_\lambda+\min(1/g_\lambda,\beta_C/2)\).
For distinct slopes use the two integrality constraints
\(f_\lambda s_\lambda\in\mathbb Z\) and
\(h s_\lambda=v_p(\sigma(\tau))\in\mathbb Z\).
Together with \(f_\lambda L_\lambda,hL_\lambda\in\mathbb Z\)
they give an increment of at least \(1/g_\lambda\).
For equal slopes use their sum.
Hence any failing class has \(u_C\ge W_C\).
The strict norm inequality is again
\[
                  -hb+2d\sum_C u_C<0 ,
\]
which proves the cost constraint.

If \(b=1\), then \(d=1\). The only \(\theta\)-orbit with positive
\(\beta_C\) contains the unique balanced pair. Let \(a\) be its
length on conjugate pairs. Each oriented embedding orbit above
it has \(f_\lambda=a\) or \(2a\), since the passage from
oriented embeddings to pairs has fibres of size two.
Consequently \(b_{\lambda,1}=f_\lambda/a\le2\), and each
class in this orbit has cost at least \(h/(2a)\).
Failure in one class would force failure in all \(a\)
classes and cost at least \(h/2\). This exceeds the strict
budget. All other places have equal reference slopes and
their polygons agree with them. \(\square\)

This proposition specializes to the homogeneous bound since at
its balanced places \(b_{\lambda,1}=f_\lambda\). In the general
case the rank loss from a failed positive-\(\beta_C\) pair is
exactly the number of zero slopes in its two reference components
(counted once per pair of oriented embeddings). It is one for
each of \(b_{\lambda,0},b_{\lambda,2}\) that equals zero.
Both smaller slopes increase on failure. Thus the cost bound
can also be used for a weighted \(p\)-rank estimate, including
places where the reference itself has positive fractional slopes.

The exclusion of potentially CM varieties is necessary even
for the averaged conclusion. Set \(B:y^2=x^3-1\) and
\(A=B^2/\mathbb Q\), with \(\mathbb Q(i)\) acting by the
rotation of the two copies and Rosati for the product
polarization. It gives one balanced conjugate pair.
The order-three automorphism of \(B\) provides its CM
over \(\overline{\mathbb Q}\). At \(p\equiv5\bmod12\) the
acting \(\mathbb Q(i)\) splits, while \(x\mapsto x^3\)
permutes \(\mathbb F_p\) and gives \(\#B(\mathbb F_p)=p+1\).
Thus the reference is ordinary but the two copies are
supersingular. These primes have density \(1/4\).
