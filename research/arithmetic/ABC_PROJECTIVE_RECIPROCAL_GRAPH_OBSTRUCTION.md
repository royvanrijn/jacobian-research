# Projective reciprocal cancellations in dense additive graphs

Every eight-edge trade supplied by common-target rectangle
supersaturation has a nonzero projective reciprocal defect. More
generally, there are dense primitive integer graphs at the endpoint
height and dictionary scales in which the projective coordinate has
unique six-fold sums, even though the complete letters have small
higher sumsets. The two claims distinguish label counting from the
first moments of labels along their incidence lines.

## The defect of a common-target rectangle trade

Write \(F(z)=z-1/z+1/(1-z)\). On the two edges in each
cell \(i,j\in\{0,1\}\) of the common-target trade (BT2), use the
opposite signs \((-1)^{i+j}\) and \(-(-1)^{i+j}\), taking the
first sign when \(h\) is on the first summand. The two summand pairs
and their common target are
\[
 (a+it+h,b+js), \quad (a+it,b+js+h),\qquad c_{ij}=a+b+h+it+js .
\]
Here \(h,t,s\) are nonzero integers. All complete values are positive,
and each of the three dictionaries contains four distinct values.
Let \(P_A,P_B,P_C\) be the products of the four distinct values
in the three roles. Put
\[
                    A_*=2a+h+t,\qquad B_*=2b+h+s .
\]
Then the signed projective defect satisfies
\[
 \boxed{\quad
 \sum_{\text{eight edges }e}\epsilon_e F(a_e/c_e)
    =hst\left(\frac{A_*+B_*}{P_C}
                   -\frac{A_*}{P_A}-\frac{B_*}{P_B}\right).
 \quad}                                                        \tag{FD1}
\]
In particular its sign is \(-\operatorname{sgn}(hst)\). If all
twelve labels lie in \([\sigma X,X]\) with fixed \(\sigma>0\),
its magnitude is comparable to \(|hst|/X^3\), with constants
depending only on \(\sigma\).

Indeed \(F(a/c)=a/c+a/b-b/a\) for \(a+b=c\).
For
\[
 R(x;u,v)=\frac1x-\frac1{x+u}-\frac1{x+v}+\frac1{x+u+v}
       =\frac{uv(2x+u+v)}{x(x+u)(x+v)(x+u+v)},
\]
the \(a/c\) part contributes \(h R(a+b+h;t,s)\).
For the \(a/b\) part, write the cell difference with \(a_i=a+it\)
and \(b_j=b+js\) as
\[
 a_i\left(\frac1{b_j}-\frac1{b_j+h}\right)+\frac h{b_j}.
\]
The term \(h/b_j\) cancels across \(i\); the alternating sum
of \(a_i\) is \(-t\). The contribution is \(-t R(b;h,s)\).
The \(-b/a\) part gives \(-s R(a;h,t)\) by the same calculation
with the roles exchanged. This proves (FD1). In one of the two
four-edge matchings each target is strictly larger than its
matched first summand and second summand. Hence
\(P_C>P_A,P_B\). The summand numerators \(A_*,B_*\) are positive
sums of opposite corner values, proving the strict sign. In the
stated box every matched quotient \(c_e/a_e\) and \(c_e/b_e\) is
at least \(1+\sigma\), so the bracket in (FD1) has magnitude
bounded above and below by positive constants times \(X^{-3}\).
[Common-target rectangle trades (BT1–BT2)](ABC_BALANCED_LATIN_TRADE_OBSTRUCTION.md#eight-edge-supersaturation-for-actual-additive-values)

## Linear identities stable under two independent translations

Let \(S\) be a finite set of nonnegative integers. For each
\((x,y)\in S^2\), consider the rational function of independent
variables \(U,V\) obtained by using summands \(U+x,V+y\).
With coefficients \(w_{xy}\in\mathbb Q\), define
\[
 \Phi_w(U,V)=\sum_{x,y}w_{xy}
 \left(\frac{U+x}{U+V+x+y}+\frac{U+x}{V+y}-\frac{V+y}{U+x}\right).
\]
The following conditions are equivalent:

- \(\Phi_w\) is constant.
- Every row \(x\) has \(\sum_y w_{xy}=\sum_y y w_{xy}=0\),
  every column \(y\) has \(\sum_x w_{xy}=\sum_x x w_{xy}=0\),
  and every diagonal \(z=x+y\) has
  \(\sum_{x+y=z} w_{xy}=\sum_{x+y=z}x w_{xy}=0\).
- The coefficient polynomial \(W(P,Q)=\sum_{x,y} w_{xy}P^xQ^y\)
  is divisible by \((P-1)^2(Q-1)^2(P-Q)^2\).

In these cases \(\Phi_w=0\). Call this a first-moment
cancellation. The first condition forces zero residues on each
distinct linear pole \(U+x=0\), \(V+y=0\), \(U+V+z=0\).
Their numerator polynomials are respectively
\(-\sum_y w_{xy}(V+y)\), \(\sum_x w_{xy}(U+x)\), and
\(\sum_{x+y=z}w_{xy}(U+x)\). These vanish exactly under the
second condition. After clearing the product of the distinct
linear poles, the numerator has total degree at most the
degree of this product. Vanishing residues give divisibility
by each factor and hence a constant rational function.
On \(U=V\to+\infty\) its limit is \(F(1/2)\sum_{x,y}w_{xy}=0\).

For the polynomial equivalence, a row's two conditions are
\(W(P,1)=\partial_QW(P,1)=0\). Similarly the column conditions
fix double vanishing at \(P=1\). The diagonal conditions give
\(W(P,P)=\partial_PW(P,P)=0\). Viewing \(W\) as a polynomial
in its first variable with the second as parameter gives
double divisibility by \(P-Q\). The three linear factors
are distinct, so their squared product divides \(W\).
This is the first-moment version of the line-sum divisibility
in algebraic discrete tomography.
[Hajdu–Tijdeman, §4, divisibility lemma for line sums](https://math.unideb.hu/sites/default/files/inline-files/disctom8.pdf)

A nonzero first-moment cancellation has at least \(13\)
active edges. Every occupied line in any of the three
directions has at least three nonzero entries, since with
two distinct positions the sum and first moment determine
both weights. The convex hull of the support has dimension
two, since a collinear support has an isolated point on a
line of another direction. Its six supporting lines at
the two extrema of each of \(x,y,x+y\) are nondegenerate
sides, as each contains at least three entries. The support
therefore has at least six convex vertices and six
nonvertex side points. A nonvertex point on the maximum
\(y\) side has \(x\)-coordinate strictly between the extrema
of \(x\). Its vertical line intersects the boundary at
only two points. It must contain at least one strictly
interior support point, giving the bound \(13\).
When the nonzero coefficients are all \(1\) or \(-1\),
every occupied line has at least four entries because its
sum is zero. Now the boundary has at least eighteen support
points. The maximum \(y\) side has two nonvertex points
with distinct intermediate \(x\)-coordinates. Each of their
vertical lines contains at least two strictly interior
support points, giving the bound \(22\).

The degree condition along a label also explains a
terminal-pair obstruction. In a prescribed first-moment
cancellation pattern every one of the labels of the last
two edges already occurs on the preceding edges, since
each occupied line has degree at least three. The preceding
edges determine the exact final edges when the pattern is
fixed. Recovery of a final pair from its projective
coordinates gives no further loss for this pattern.

## A primitive graph at the endpoint height scale

Fix an integer \(K\ge1\), a finite alphabet of nonzero rational
weights, and a real \(0<L<3/K\). There are arbitrarily large
\(X\) and primitive graphs of ordered positive triples \(a<b\),
\(a+b=c\), \(X/2<c\le X\), with these properties:

- They have \(X^{L/3+o(1)}\) edges. Their three complete-letter
  dictionaries have sizes \(X^{L/6+o(1)}\), and a central
  complete-target fibre also has size \(X^{L/6+o(1)}\).
- In each dictionary every fixed-order sumset has cardinality
  \(X^{L/6+o(1)}\).
- Every zero relation in \(F(a/c)\) with at most \(K\) distinct
  rows and coefficients from the alphabet is a first-moment
  cancellation on the unshifted labels.

To see this, take integers \(d,H\to\infty\) with
\(\log d=o(\log H)\), put \(N=H^d\) and
\[
 Q=4H\lceil dH^{6/L-1}\rceil,\qquad
 M=Q^d,\qquad
 S=\left\{\sum_{j=0}^{d-1}u_jQ^j:0\le u_j<H\right\}.
\]
Then \(|S|=N\), \(S\subset[0,M)\) and
\(\log M=(6/L+o(1))\log N\). The raw first and second summand
dictionaries are \(U+S,V+S\), with integer shifts chosen from
\([4M,5M)\) and \([7M,8M)\). The targets lie in
\(U+V+S+S\) and the whole graph lies in the claimed shell with
\(X=15M\). The interval separation gives \(a<b\).
The raw dictionaries have fixed-order sumsets of size at most
\((2kH)^d=N^{1+o(1)}\) for each fixed \(k\), by carry-free
addition. Retaining primitive edges or pruning inactive
labels only decreases these upper bounds.

For any particular non-first-moment relation among at most
\(K\) candidate edges, clearing the distinct linear
denominators gives a nonzero polynomial in \(U,V\) of
degree at most \(3K\). A nonzero degree-\(D\) bivariate
polynomial vanishes on at most \(2DM\) pairs from the shift
square. To see this, choose a nonzero coefficient of a
power of \(V\); it is a degree-at-most-\(D\) polynomial
in \(U\) and has at most \(D\) zeros. For each other \(U\),
there are at most \(D\) zeros as a function of \(V\).
There are \(O(N^{2K})\) candidate relations at the finite
alphabet. Thus the measure of shifts that admit any
accidental relation is \(O(N^{2K}/M)=o(1)\).

It remains to choose good shifts with plenty of primitive
edges. For a fixed pair \(x,y\), the probability that the
two shifted values share a prime is at most
\[
 \sum_{p\le 9M}(p^{-1}+M^{-1})^2
 \ \le\ \frac14+\sum_{\substack{j\ge3\\j\ {\rm odd}}}j^{-2}+o(1)
 \ =\ \frac{\pi^2}{8}-\frac34+o(1)<\frac12 .
\]
The bound follows by independence of the shifts and by
counting their residues in intervals of length \(M\).
The central target \(U+V+\max S\) has exactly \(N\)
candidates, given by pairing each digit vector with its
digitwise complement. Apply Markov's inequality to the
number of nonprimitive candidates in the whole graph
and at the central target. The two probabilities of
having fewer than \(N^2/100\) or \(N/100\) primitive
edges are together at most
\(2(\pi^2/8-3/4+o(1))/(99/100)<1\).
Thus there are good shifts outside the accidental
relation set. Each complete dictionary of active
values has size at least \(N/100\), because one
complete value occurs on at most \(N\) edges. This
also gives the matching lower bound on the sizes
of the fixed-order sumsets.
[Carry-free integer model and primitive-edge selection](ABC_CRITICAL_TRIANGLE_ENTROPY_MODEL.md#exact-addition-height-and-primitivity)

In particular take \(K=20\) and the alphabet
\(\{-10,\ldots,-1,1,\ldots,10\}\). Two multisets of
six projective values with equal sum would give a
nonzero relation with at most twelve active edges
after cancellations. This is impossible. The same
conclusion holds for two subsets of at most ten
distinct projective values: the resulting signed
coefficients are \(1,-1\), so every nonzero relation
needs at least twenty-two active edges. The projective
values are distinct on these primitive ordered rows
because \(F\) is injective on \((0,1)\) and the reduced
ratio determines the row. Thus the model gives scalar
sum uniqueness alongside ordinary additive compression
of the complete dictionaries. The construction supplies
no total radical restriction on a row. Relation forcing
for the capped families would need to use the prime
supports or numerical weights inside the complete
values, beyond their graph and height exponents.

## The denominator of a forced trade in a balanced capped graph

For every \(\eta>0\) there are \(\kappa_0,\tau_0,L_0>0\) with
the following property. For fixed \(0<L<L_0\) put \(n=X^L\) and
let a family of primitive triples in \(X/2<c\le X\) have at least
\(n^{1/3-\tau_0}\) rows, three complete dictionaries of cardinality
at most \(n^{1/6+\kappa_0}\), and individual radical bounds
\[
                n^{1/3-\kappa_0}\le \operatorname{rad}(W)
                                  \le n^{1/3+\kappa_0}.
\]
For large enough \(X\) it contains at least \(n^{5/6-\eta}\) distinct
reduced denominator radicals of the scalar defects (FD1) from its
faithful trades (BT2). Each of those denominator radicals is at
least \(n^{4-\eta}\).

Here a faithful trade has twelve distinct role labels as in (BT2).
Treat their integer values as the twelve linear forms on the five
parameters \(a,b,h,t,s\). This parameter space has rank five, since
the labels \(a,a+h,a+t,b,b+s\) form a basis. No two label forms
are proportional: across roles their leading coefficients in \(a,b\)
are in different directions, and within a role they have the same
leading coefficients but distinct shifts.

For a label form \(W'\), the two original edges incident with it in
the trade span distinct two-dimensional spaces of label forms,
intersecting exactly in the line spanned by \(W'\). To see this for
an \(A\)-label, the opposite \(B\)-labels in its two edges differ
by \(s\), which is independent of the \(A\)-label and either
\(B\)-label. For a \(B\)-label use the change by \(t\).
For a \(C\)-label the incident \(A\)-labels differ by \(h\),
which is independent of that \(C\)-label and either \(A\)-label.
The two companions in any single edge span its plane together
with \(W'\), since that edge satisfies the additive equation.

Fix two distinct labels \(W,W'\). In one of the planes through
\(W'\), the form \(W\) lies outside the plane. Choose a companion
label \(T\) in this edge. Then \(W,W',T\) are independent. Extend
them to a basis with two more label forms drawn from the trade.
Enumerate the four basis labels other than \(W'\). The count is
at most \(n^{4(1/6+\kappa_0)}\). The original edge containing
\(W'\) and \(T\) has the latter value fixed. If the common radical
of \(W,W'\) is at least \(n^\xi\), the outside radical cap for
\(W'\) with pool supported on \(W\) is at most
\(n^{1/3+\kappa_0-\xi}\). Apply the uniform fixed-letter bound
(FL1) to that edge. The number of possibilities for \(W'\) is
\(O(n^{1/6+\kappa_0/2-\xi/2+\gamma})\), for any prescribed
\(\gamma>0\) at sufficiently small fixed \(L\). The five basis
values determine the entire trade. For each label pair the count
of violations is therefore
\[
                O(n^{5/6+(9/2)\kappa_0-\xi/2+\gamma}).
                                                               \tag{FD2}
\]
The bound uses the full original edge after fixing \(T\); it does
not need a separate cap on the varying shift.
[Fixed-letter bound with a prescribed prime pool](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

A prime is private on a label if it divides no other label of the
trade. Such a prime of an \(A\)-label cannot divide \(h\) or \(t\),
because each is a same-role difference from that label. Its possible
intersection with the shift support \(\operatorname{rad}(s)\)
can also be removed at the cost (FD2). Fix the two \(B\)-labels of
its incident edges, whose difference is \(s\). Their span does not
contain the \(A\)-label, and either one is a fixed companion in
one of those edges. Extend to a basis as before, fix the other
four basis values, and apply (FL1) with pool supported on \(s\).
For a \(B\)-label use the \(A\)-pair with difference \(t\).
For a \(C\)-label the two same-role differences are \(t,s\);
use the two \(A\)-labels in its incident edges, with difference
\(h\). Their span excludes the \(C\)-label. This establishes (FD2)
as well for each of the twelve remaining possible intersections
of a private label support with \(hts\).

Choose \(\xi>0\) small in terms of \(\eta\), and then choose
\(\gamma,\kappa_0,\tau_0\) sufficiently smaller than \(\xi\).
The count (BT1) supplies \(\gg n^{5/6-10\tau_0-15\kappa_0}\)
faithful trades because the ratio of its edge count to the
\(7/4\) power of the maximal dictionary size tends to infinity.
Deleting the finitely many classes counted by (FD2) leaves the
same lower bound up to a constant. On every remaining trade each
label has private radical outside \(hts\) of mass at least
\(1/3-\kappa_0-12\xi\). At a prime in this private radical,
the corresponding term of (FD1) has the unique negative pole
in the bracket. Its numerator is a \(p\)-unit: the role's
numerator \(A_*,B_*\), or \(A_*+B_*\) is the sum of two opposite
corner labels and can be written with either corner pair,
including the single label divisible by \(p\). Also \(hst\) is
a \(p\)-unit here. Consequently the reduced defect denominator
retains the full private power of every such prime. Summing
over the labels gives the lower bound
\(n^{4-12\kappa_0-144\xi}\) for its radical.

For multiplicity, fix one denominator radical. The exact
denominator has polynomial height in \(X\), by (FD1), so it has
only height-subpower many possibilities on that fixed support.
Assign its prime powers to the twelve label slots or an unused
slot. This has height-subpower many choices and allows the
assignment of all the retained private powers to their labels.
Now on each of the five basis labels \(a,a+h,a+t,b,b+s\),
the total raw mass missing from its assigned powers is at most
\(12\xi\). Specify each missing support and the full powers on it
at cost \(n^{60\xi+o_X(1)}\) across the five labels. These values
determine the entire trade. Choosing the tolerances sufficiently
small in terms of \(\eta\) gives the asserted count of distinct
denominator radicals.
