# Canonical quadratic norm elements and the cost of moving fields

Every primitive additive triple gives an imaginary quadratic norm
element with no extra coefficient factor in its norm. At a deep
target divisor its ordinary ray residue is a fixed multiple of the
absolute summand residue. Thus the absolute ray count applies to
arbitrary rows in a packet once their fields are specified. There
are primitive quotients showing that the field count cannot simply
be omitted, even when the complete target is fixed.

## The quadratic norm element and its actual packet residue

Let \(a+b=c\) be positive and primitive. Write \(ab=D z^2\) with
\(D\) positive squarefree and \(z>0\). In \(K=\mathbb Q(\sqrt{-D})\)
choose the square root with positive imaginary part and put
\[
                     \alpha=(a-b)+2z\sqrt{-D}.
\]
Then
\[
             N_{K/\mathbb Q}(\alpha)=c^2,\qquad
             \operatorname{Tr}_{K/\mathbb Q}(\alpha)=2(a-b).
                                                               \tag{CQ1}
\]
Both complex absolute values are \(c\). In particular there is
no archimedean log variable in the normalized ray group, and
\((N\alpha,\operatorname{Tr}\alpha)\) recovers the ordered
primitive row. In the notation \(a=A u^2,\ b=B v^2\) with
squarefree \(A,B\), primitivity gives \(D=AB\), and
\(\alpha=(A u+v\sqrt{-AB})^2/A\). This is the usual square
norm equation for the primitive pair.
[Uniform quadratic norm recovery](ERDOS_943_PELL_RECOVERED_KERNEL_GCD_WEIGHT.md#1-a-uniform-elementary-norm-lemma)

Every odd \(p\mid c\) is prime to \(Dz\) and splits in \(K\).
Indeed \(-D\equiv(a/z)^2\pmod p\) since \(b\equiv-a\pmod p\).
The trace is a unit at \(p\), since it is congruent to \(4a\)
there. Thus at exactly one prime over \(p\), \(\alpha\) is a
unit, while at the other it has valuation \(2v_p(c)\) by its
norm. On the unit branch its conjugate vanishes to that latter
depth. In the degree-one completion the precise congruence is
\[
                   \alpha\equiv 4a-2c\pmod{p^{2v_p(c)}}.
                                                               \tag{CQ2}
\]
In particular, on that branch \(\alpha\equiv4a\pmod{p^j}\)
whenever \(p^j\mid c\). This residue does not involve the
values of the two original squarefree coefficients.

Fix a prescribed odd divisor \(M\mid c\) with \(M\ge X^\nu\),
a prescribed positive integer \(P\) of polynomial height, and a
dictionary of at most \(R\) values of \(a\bmod M\).
Assume \(c\le X\), the radical of \(c\) outside \(2PM\) is
at most \(X^{\ell L}\), and there are at most \(S\) fields
\(K\) in the family, each of discriminant
\(\ll X^{gL}\). Fix \(\nu>0\), the polynomial height exponent
for \(P\), \(g,\ell\), and \(\eta>0\) before choosing \(L\).
For all sufficiently small fixed \(L>0\) and every
\(\varepsilon>0\), (AT1) gives
\[
              \#\{(a,b,c)\}
                 \ll S R^{1/2}X^{(\ell/2+\eta)L+\varepsilon}.
                                                               \tag{CQ3}
\]
For each fixed field record the unit prime over each \(p\mid M\)
at cost \(2^{\omega(M)}=X^{o_X(1)}\). The ideal modulus is the
product of its \(j\)-th powers for \(p^j\parallel M\) and has
norm \(M\). Equation (CQ2) fixes at most \(R\) ordinary ray
residues modulo this ideal. Outside the rational pool \(2PM\),
the radical ideal of \(\alpha\) has norm equal to the radical
of \(c\), because it uses exactly one degree-one prime over
every occurring odd prime. Apply (AT1), sum over the field
list, and use the exact reconstruction in (CQ1). No separate
dictionary of \(A,B\) is needed. The original summand radical
caps can be imposed, and they supply a uniform \(g\) through
\(D\le\operatorname{rad}(ab)\) and \(|\operatorname{Disc}(K)|\le4D\).
[Absolute residue lists and norm moments (AT1)](ABC_ABSOLUTE_TARGET_RAY_MOMENTS.md#algebraic-integers-with-a-ray-residue-dictionary)

If the complete target dictionary has at most \(T\) entries,
one may replace \(R^{1/2}X^{(\ell/2+\eta)L}\) in (CQ3) by its
minimum with \(T\). Indeed, fixing \(c\) and \(K\) fixes the
norm \(c^2\), with \(\tau(c^2)^2\) possible ideals and at most
six generators of each, up to the required norm. This is the
classical ideal-count bound for positive definite binary
quadratic representations.
[Green–Soundararajan, §3.2 and Lemma 3.1](https://www.cambridge.org/core/journals/journal-of-the-institute-of-mathematics-of-jussieu/article/covering-integers-by-x2-dy2/D767FC91FA688562345709324C5FD592)

The same estimate applies to a neighbourhood defined by a
joint square defect. Fix a first row, a prescribed
polynomial-height prime pool \(P_1\), and a bound \(B_1\ge1\)
on the squarefree defect of \(a_1b_1a_2b_2\) outside \(P_1\).
The next field is determined by that defect and a subset of
the primes of \(P_1\), since its square class is the product
of the first field's positive squarefree parameter and the
defect in \(\mathbb Q^\times/\mathbb Q^{\times2}\).
There are at most \(B_1 X^{o_X(1)}\) next fields.
Apply (CQ3) with this local field count and the same
actual residue dictionary. At a fixed target no coefficient
or residue record beyond the joint defect is needed.

## A fixed-packet family over distinct fields

Fix integers \(h\ge e\ge2\), take a large odd squarefree
integer \(Q\), and put
\[
                  c=M=Q^e,\qquad X=Y=Q^h,\qquad L=1/h.
\]
For each positive integer \(x\) prime to \(Q\) with
\[
 \cos(\pi/(3e))\sqrt Q<x<\cos(\pi/(6e))\sqrt Q ,
\]
write \(Q-x^2=D_x z_x^2\) with \(D_x\) positive squarefree.
Set
\[
 (x+z_x\sqrt{-D_x})^e=u_x+v_x\sqrt{-D_x},
 \qquad a_x=u_x^2,\qquad b_x=D_x v_x^2 .
\]
The base element has norm \(Q\) and its argument is between
\(\pi/(6e)\) and \(\pi/(3e)\). The powered element therefore
has real and imaginary parts between positive fixed fractions
of \(Q^{e/2}\), and \(a_x+b_x=c\). At each prime dividing \(Q\)
the base element uses exactly one of its two split primes,
because \(x,z_x,D_x\) are units there and the rational norm is
squarefree. This remains true after taking the power. A
rational prime dividing both \(u_x,v_x\) would divide \(Q\)
and then occur on both split primes, which is impossible.
Also no prime divides both \(u_x,D_x\), since it would divide
\(c\) and \(D_x\), while \((Q,D_x)=1\). This verifies
primitivity.

The chosen interval contains \(Q^{1/2-o(1)}\) integers prime
to \(Q\). Möbius inversion gives a main term equal to the
interval length times \(\varphi(Q)/Q\), with error
\(O(\tau(Q))\). Here \(\varphi(Q)/Q\ge 2^{-\omega(Q)}\).
For any one \(D\), the number of representations
\(x^2+D z^2=Q\) is at most \(6\cdot 2^{\omega(Q)}\),
by counting norm-\(Q\) ideals in its ring of integers and
at most six unit multiples of each generator. Thus select
one row for each of
\[
                               N=Q^{1/2-o(1)}
\]
distinct \(D_x\). These give distinct rows because
\(\operatorname{sf}_2(a_xb_x)=D_x\). Their \(a_x\bmod M\)
are distinct, as \(0<a_x<c=M\) and \(c\) is fixed.
The same is true of their unit residues \(4a_x\bmod M\).
Each field has discriminant at most \(4Q\), and all the
canonical norm elements have radical ideal norm one
outside \(M\).

This is also a family of distinct minimally repaired primitive
quotients. The least \(h\)-full target height is
\(H_h(c)=Q^h=Y\); for \(h>e\) the common multiplier
\(Q^{h-e}\) is mandatory and unique. One may take \(Q\)
to be the product of the odd primes up to a growing
cutoff. In that case every packet prime has size
\(Y^{o(1)}\), since \(\log Q\) grows faster than the
logarithm of the largest prime factor. For example,
Bertrand's postulate applied to any fixed number of
disjoint dyadic intervals below the cutoff suffices
for this last assertion. A fixed choice of \(e/h>0\)
up to bounded rounding keeps \(M\ge Y^\nu\) for a
uniform positive \(\nu\) as \(h\) grows.

Consequently a field-uniform variant of (CQ3) with
\(S^\gamma\) in place of \(S\), retaining the same \(R^{1/2}\)
term, is false for every fixed \(\gamma<1/2\) on these data
alone. Here \(S=R=N\), \(\ell=0\) and \(X^L=Q\).
The proposed estimate would force
\(N^{1/2-\gamma}\ll Q^{\eta+h\varepsilon}\).
Choose \(\eta<(1/2-\gamma)/2\), take sufficiently large
fixed \(h\) with \(e/h\) bounded below, then take
\(\varepsilon>0\) sufficiently small at that budget.
The estimate contradicts the lower bound as \(Q\to\infty\).
This obstruction omits the two original summand radical
caps. Any bound that removes the field cost has to use
additional information from the summands.

The normalized norm elements \(\alpha_x/c\) in this family
are multiplicatively independent. In the compositum of the
distinct quadratic fields,
each such element transforms under Galois conjugation by
its own quadratic character, with conjugation acting by
inversion. Apply each automorphism to a multiplicative
relation and multiply the resulting relations weighted by
one of the quadratic characters. Character orthogonality
leaves a power of only the element belonging to that
character. It is not a root of unity since at any prime
above an odd rational prime dividing \(Q\) its valuation
is \(\pm e\) in its own field. Thus the original exponent
of that element must vanish. The field cost persists even
for such multiplicatively independent rows.


## Exact class separation at every order in a deep packet

Fix a complete letter \(T\) of a primitive triple with all
letters at most \(X\), and write \(U,V\) for the other two.
For every \(\varepsilon>0\), each prescribed \(q\)-power class
of \(UV\) or \(U/V\) contains \(O_\varepsilon(X^\varepsilon)\)
such rows, uniformly for every integer \(q\ge2\) and \(T\).
This does not require the radical caps.

To prove it write \(U=A u^q,\ V=B v^q\) with canonical
\(q\)-free coefficients \(A,B\). Since \((U,V)=1\), the primes
of \(AB\) are precisely the primes of the prescribed power
class. For the product, assigning each prime to its letter
determines the reduced exponent from the class. For the
ratio, the two assignments give complementary reduced
exponents. On any occurring class the support is contained
in one \(UV\le X^2\), so the assignments number
\(2^{\omega(UV)}=X^{o_X(1)}\) uniformly in \(q\). At fixed
\(A,B\), for all \(q\) above a threshold depending on
\(\varepsilon\), count the two roots by \(X^{2/q}\).

It remains to treat a finite list of orders. For each order,
one has a separable binary equation
\[
                         A u^q\pm B v^q=\pm T ,
\]
with the signs depending on the fixed role. If the binomial
has a nonlinear irreducible factor \(F\), choose a primitive
integral factor. Its exact nonzero value divides the right
side, giving \(X^{o_X(1)}\) possibilities. For each value,
the fixed-degree norm-element argument in (CS4) counts
\(X^{o_X(1)}\) inputs, uniformly in the coefficients.
Specifically, the algebraic integer \(a_0 u-\theta v\) in the
root field of \(F\), where \(a_0\) is the leading coefficient
of \(F(T,1)\) and \(\theta\) is \(a_0\) times a root of that
polynomial, has fixed nonzero integral norm \(a_0^{\deg F-1}F(u,v)\).
The coefficients, embeddings and norm have polynomial height
at each fixed order. Count principal ideals of that norm,
then generators in the polynomial embedding box by
fixed-degree unit log packing. Irrationality of \(\theta\)
recovers both inputs from the element. If all factors are
linear, fix the signed values of two distinct factors. The
value of each divides \(T\) and the two forms determine the
inputs. Factor coefficient bounds are uniform in fixed
degree by the same multiplication argument as in (CS4).
[Uniform fixed-value norm count (CS4)](ABC_MINIMAL_REPAIR_CYCLOTOMIC_SPLIT.md#direct-target-bounds-for-reducible-weighted-powers)

More generally, in any family with at most \(T_*\) complete
values of the fixed role, prescribe an anchor row
\((U_0,V_0)\) and a prime pool supported on a positive integer
at most \(X^C\), for fixed \(C\). For \(B_*\ge1\), the number of
query rows for which some integer \(q\ge2\) satisfies
\[
      \prod_{\substack{p\notin\mathcal P\\
           q\nmid v_p((U/V)/(U_0/V_0))}}p\le B_*
\]
is at most \(O_{C,\varepsilon}(T_* B_* X^\varepsilon)\).
The same statement holds with the rational ratio quotient
replaced by the product \(UVU_0V_0\).
For each order enumerate the displayed squarefree support at
cost \(B_*\). The query's canonical coefficients \(A,B\) are
supported on its union with the prime pool and the anchor support.
Their exact values are counted at height-subpower cost
uniformly, since their heights are at most \(X\) on a fixed
support of polynomial prime product. In particular one does not
pay for all possible reduced exponents up to \(q\) separately.
Then use the preceding fixed-letter count, or its proof with
\(A,B\) already fixed. Orders above \(4\log_2 X\) give the
same defect support as the first integer above this threshold,
since every nonzero valuation has absolute value at most
\(4\log_2 X\). The union over the orders therefore costs only
another height subpower. The pool may depend on the anchor
and on the height. For the ratio comparison, multiplying
either original row by a common integer leaves its odds
unchanged, including under any least target repair.

The same query bound holds for defect comparisons of
\(U^r V^s/(U_0^{r_0}V_0^{s_0})\), allowing the orders and
all four integer weights to vary, provided the query
weights \(r,s\) are units modulo \(q\). The prime support
test off the anchor and the pool still detects every query
coefficient prime. Enumerate the support and the exact
coefficients as before. For \(q\le\log_2 X\) sum the orders,
without summing over any weights. For \(q>\log_2 X\),
both query roots are one, so the coefficient values
already determine the query for all such orders at once.
[Uniform supported-value enumeration](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

Apply this to a surviving packet of (AP2), with \(n=X^L\)
and \(c_0=1/6\). The final thinning in its proof retains at
least \(n^{2c_0-s-\eta/2}\) marked edges per packet and the
target dictionary has size below \(2n^{c_0-s}\). Thus each
such packet requires
\[
                           \gg_\varepsilon
                                 n^{c_0-\eta/2} X^{-\varepsilon}
                                                               \tag{CQ4}
\]
distinct \(q\)-power classes of the complete summand product
or odds for every \(q\ge2\). Moreover for every fixed anchor
the fraction of query edges in that packet with an outside
joint defect at most \(n^t\) at any order is
\[
                        O_\varepsilon(
                              n^{t-c_0+\eta/2}X^\varepsilon).
                                                               \tag{CQ5}
\]
Inside each occupied residue block the corresponding
bound is \(O_\varepsilon(n^{t-(c_0-s)+4\eta}X^\varepsilon)\),
using its target dictionary upper bound and edge lower
bound. In particular a greedy selection in one packet gives
\(\gg_\varepsilon n^{c_0-t-\eta/2}X^{-\varepsilon}\) rows
with pairwise defect greater than \(n^t\) at every order,
when the same pool is used for all queries and anchors.
These inequalities hold with all three original
radical caps used to obtain (AP2). In particular, at late
budgets cheap joint power classes, even with adaptively chosen
orders, cannot cover a positive fraction of the pairs in a
near-full packet. The proof supplies a dispersion bound on
the actual classes; controlling how those dispersed classes
interact with the other original radical constraints remains
open.
[Composite packet count and final thinning (AP2)](ABC_AGGREGATE_DEEP_PACKET_INVERSE.md#composite-packet-inverse)
