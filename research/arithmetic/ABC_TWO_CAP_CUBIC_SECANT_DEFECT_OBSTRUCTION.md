# Two radical caps and the cubic defects of projective secants

Two balanced letter caps do not bound the corrected cubic defects
of typical tuples. The obstruction persists when the two tracked
letters are exact cubes, each has a diffuse degree law, and the
target multiplicity is prescribed. Within an affordable three-cap
fibre there is also a lower bound on the price of the corrected
classes on typical triples.

## Two exact-cube dictionaries

Use \(n=X^L\), taking \(X\to\infty\) at each fixed \(L>0\) before
letting \(L\to0\). For a nonzero integer \(z\) and prime pool \(P\)
put
\[
 d_3(z;P)=\prod_{\substack{p\notin P\\3\nmid v_p(z)}}p.
\]
A pool may be prescribed at all heights, with prime product at most
\(X^C\) for a fixed \(C\). All subsequent height counts are uniform
over this choice. The constructed family can depend on the pool.

**CD1.** Fix \(0\le u\le1/2\). There are families of
\(n^{1/3+o_X(1)}\) distinct positive primitive triples \(a+b=c\le X\)
with \(\operatorname{rad}(a),\operatorname{rad}(c)\le n^{1/3}\).
Both \(a\) and \(c\) are exact cubes, and each occurring \(c\) has
exactly \(\lfloor n^{u/3}\rfloor\) representations. The degree laws
of \(a,c\) tend to the diffuse measure of density \(1/6\) on
\([2,4]\). Put \(S_{ij}=a_i c_j-a_j c_i\), up to sign, on
distinct rows. There are \(B_L\to\infty\) and \(\varepsilon_L\to0\)
such that outside at most \(n^{5/6+\varepsilon_L+o_X(1)}\) ordered
distinct triples of rows, the targets are distinct and
\[
                       d_3(S_{12}S_{13}S_{23};P)>n^{B_L}.
                                                               \tag{CD1}
\]
In particular both secant-corrected cubic classes on \(A,C\)
in (PM5) are the class of \(S_{12}S_{13}S_{23}\), with the
same conclusion if the opposite sign convention for the
secants is used. The construction makes no small radical
assumption on \(b\).

Take
\[
 t=\left\lfloor\tfrac12\log_2(1/L)\right\rfloor,\qquad
 E=3\cdot 2^t,\qquad r=t-3,\qquad m=\lfloor L^{-1/4}\rfloor .
\]
For small \(L\) all the parameters used below are positive.
Set
\[
 k_i=\left\lfloor\frac{2+(2i-1)/m}{EL}\right\rfloor,\qquad
 H_L=\frac{1}{3m}\sum_{i=1}^m EL k_i,\qquad
 R=n^{1/(3m)},\qquad Q=\lfloor R/(\log R)^2\rfloor .
\]
For each of the \(m\) slots choose a disjoint palette of \(Q\)
primes in fixed subintervals of \([R/8,R/4]\) for the \(a\)-dictionary,
and a second set of disjoint palettes in \([R/2,R]\) for the
\(c\)-dictionary. The prime number theorem provides the palettes
at fixed \(L\) and large height. Within either dictionary let a
word be \(\prod_{i=1}^m p_i^{E k_i}\) with one prime from each
palette. The dictionaries \(\mathcal A_0,\mathcal C_0\) each have
\(Q^m=n^{1/3+o_X(1)}\) words. Their radicals are at most \(n^{1/3}\).
Moreover \(H_L\le1\), \(H_L\to1\), and every word has height
\(\asymp_L X^{H_L}\). Each \(c>2a\) for any \(a,c\) in the
respective dictionaries. This gives actual positive primitive
rows \((a,c-a,c)\) on the full Cartesian product, since the
two prime alphabets are disjoint. The normalized degree measure
of each tracked letter is
\[
                   \frac{1}{3m}\sum_{i=1}^m\delta_{EL k_i}
                                 +o_X(1)
\]
in total variation at fixed \(L\). This converges to the claimed
diffuse measure. This choice of a power palette follows the
two-cap construction in (DO1).
[Diffuse power dictionaries](ABC_DIFFUSE_TWO_CAP_PRIMITIVE_SECANT_OBSTRUCTION.md#diffuse-power-dictionaries)

## Small defects on individual secants

Among ordered pairs of distinct Cartesian rows, those with
\(d_3(S_{ij};P)\le n^F\) number at most
\[
                          n^{2/3+F/(2r)+o_X(1)}                \tag{CD2}
\]
for any \(F\) fixed at the current budget. To see this, write
\(a_i c_j=W^E\), \(a_j c_i=Z^E\), and put
\(g=\gcd(W,Z), x=W/g,y=Z/g\). Division by \(g^E\) preserves
the cubic defect since \(3\mid E\). The positive factors
\[
                       x^{2^j}+y^{2^j}\qquad(3\le j\le t-1)
\]
all divide \(|x^E-y^E|\). Also they are full factors of
the difference, away from the primes dividing \(E/2^t=3\)
and \(2\). In detail, they are \(\Phi_{2^{j+1}}(x,y)\)
in the homogeneous cyclotomic factorization. If a prime
\(p\nmid E\) divides \(\Phi_{2^{j+1}}(x,y)\), the ratio
\(x/y\) has exact order \(2^{j+1}\bmod p\). The factors
dividing \(x^E-y^E\) are pairwise coprime at such a
prime, since \(T^E-1\bmod p\) has simple roots. Hence the
valuation of the full difference equals that of the
selected factor there. Enlarge the pool by \(\{2,3\}\).
The outside cube-defect supports of the selected factors
are disjoint and each contributes its defect to the full
difference. At least one factor therefore has cube
defect at most \(n^{F/r}\).

For this factor enumerate its cube-free coefficient at cost
\(n^{F/r+o_X(1)}\), pricing its prescribed pool at subpower
cost. Apply (BT1) to the squarefree homogeneous form
\(x^{2^j}+y^{2^j}\), of degree \(2^j\ge8\), with prime order
three. The form and its splitting field with \(\zeta_3\)
adjoined are fixed at the chosen \(L\). Thus the term involving
the class group in (BT1) is constant in the inner height
limit. The cube-free coefficient has polynomial height,
and its number of distinct primes is \(o(\log X)\), so
(BT1) gives \(X^{o_X(1)}\) possible primitive pairs \((x,y)\)
per coefficient. We do not need constants uniform as the
degree grows in the subsequent budget limit.
[Primitive power values of binary forms (BT1)](ABC_BINARY_POWER_TWISTS.md#primitive-power-and-powerful-values-of-binary-forms)
· [Uniform Mordell–Lang input to BT1](https://arxiv.org/html/2105.15085v4)

It remains to restore the Cartesian rows from \((x,y)\).
The disjoint \(a\)- and \(c\)-palettes distinguish their
prime factors even after taking the cross-products. At
a slot where the two \(a\)-words or two \(c\)-words choose
different primes, the surviving prime in each of \(x,y\)
recovers the respective choice. At a slot where they
agree, \(g\) cancels the whole prime power in that slot,
and there are at most \(Q\) choices to restore it.
If \(s\) of the \(2m\) slots of the two letters agree,
there are at most \(Q^s\) lifts. The trivial count of
original pairs with this many shared slots is at most
\(\binom{2m}{s}Q^{4m-s}\). Compare the two bounds to get
\[
       n^{o_X(1)}
       \min\{n^{4/3-s/(3m)},\,n^{F/r+s/(3m)}\}
                           \le n^{2/3+F/(2r)+o_X(1)}.
\]
There are only finitely many slots and factor choices at
fixed budget. This proves (CD2).

## Overlap between secants and selection of a sparse family

For \(n<d\le X^2\), any unit class modulo \(d\) contains
at most \(n^{1/6+o_X(1)}\) values of \(\mathcal A_0\),
uniformly in \(d\) and the class. For completeness, let
\(e_i=E k_i\), and for each slot sum
\(S_i(\chi)=\sum_p\chi(p^{e_i})\) over the unit primes of
that palette. The normalized character mean of
\(|S_i|^{2m}\) is at most \(Q^m X^{o_X(1)}\). Fix one
ordered \(m\)-tuple of primes. The product of the other
lies in a coset of the kernel of the \(e_i\)-th power map,
of size \(X^{o_X(1)}\) uniformly by the structure of the
unit groups at prime powers and the fixed value of \(e_i\).
Its ordinary product is at most \(R^m=n^{1/3}<d\)
and has at most \(m!\) ordered factorizations within the
palette. Hölder bounds the character mean of
\(|\prod_i S_i|^2\) by \(Q^m X^{o_X(1)}\). Orthogonality
and Cauchy–Schwarz give the claimed unit-class count.
[The same moment for diffuse power dictionaries](ABC_DIFFUSE_TWO_CAP_PRIMITIVE_SECANT_OBSTRUCTION.md#thin-dictionaries-with-small-unit-classes)

Fix two distinct anchor rows \(i,j\). Suppose
\(\operatorname{rad}\gcd(|S_{ij}|,|S_{ik}|)>n^2\).
Remove the primes dividing \(a_i c_i\), of radical at
most \(n^{2/3}\). There is a squarefree divisor \(d>n\)
of \(|S_{ij}|\) coprime to \(a_i c_i\) which also divides
\(|S_{ik}|\). There are divisor-many possibilities for \(d\).
The congruence gives \(a_k\equiv a_i c_i^{-1}c_k\bmod d\).
If \(p\mid c_k,d\), this would force \(p\mid a_k\)
against primitivity. Thus for each prescribed \(c_k\)
this is a unit class for \(a_k\). It follows that there
are at most \(n^{1/2+o_X(1)}\) possible Cartesian query rows,
since \(|\mathcal C_0|=n^{1/3+o_X(1)}\). Fixing the two
anchors costs at most \(n^{4/3+o_X(1)}\), so the number
of Cartesian triples failing any of the three pairwise
secant overlap caps is at most \(n^{11/6+o_X(1)}\).

Select targets independently with probability \(n^{-u/3}\).
For each selected target independently select a uniform
subset of exactly \(J=\lfloor n^{u/3}\rfloor\) members of
\(\mathcal A_0\). Then the row count is
\(n^{1/3+o_X(1)}\) with probability tending to one.
A prescribed tuple of \(k\) rows with distinct targets has
selection probability \(n^{-k/3+o_X(1)}\). Take
\(B_L=\sqrt r\) and \(F=B_L+5\) in (CD2).
The expected number of selected bad secant pairs with
distinct targets is at most \(n^{F/(2r)+o_X(1)}\).
The expected number of overlap-failing triples with
distinct targets is at most \(n^{5/6+o_X(1)}\). Markov's
inequality gives both bounds with an exponent allowance
\(r^{-1/4}\) simultaneously with the row-count event.
Triples with repeated targets are at most \(O(N^2 J)\)
where \(N\) is the selected row count, hence at most
\(n^{5/6+o_X(1)}\) by \(u\le1/2\). Triples containing a
bad secant pair have exponent at most
\(1/3+F/(2r)+r^{-1/4}+o_X(1)\).

On each remaining triple any one secant has cubic defect
greater than \(n^{B_L+5}\) outside \(P\). Removing its
primes that divide another secant loses at most \(4\)
in normalized radical mass. Its private defect therefore
survives in the product with mass greater than \(B_L\).
This proves (CD1), taking \(\varepsilon_L\) a fixed
multiple of \(r^{-1/4}\).

The omitted actual summand has large cubic defect as
well in this construction. In the Cartesian difference
\(c-a\) the two bases are coprime and there are no shared
slots to restore. The same factorization and (BT1)
therefore give at most \(n^{F/r+o_X(1)}\) Cartesian rows
with \(d_3(c-a;P)\le n^F\). Their expected selected
count is at most \(n^{F/r-1/3+o_X(1)}\), which tends
to zero for small fixed \(L\). Thus the construction also
allows every selected row to have \(d_3(b;P)>n^F\).
The third original radical cap excludes these examples.
[Secant-corrected cubic records (PM5)](ABC_PROJECTIVE_SHAPE_SECANT_CLASS_COUNTS.md#the-cubic-case-with-no-shape-record)

## A generic price floor in affordable three-cap fibres

**CD3.** Work in the balanced affordable-fibre setting of (CR1)
and (DC3), with retained fibre exponent \(r_*>0\) in base \(n\).
Use a prescribed pool including the fixed original factors
and the common projective modulus after peeling. Up to
arbitrarily small fixed losses, the union of the three
corrected cubic defects in (PM5) has mass at least \(4r_*\)
on typical ordered triples of the retained fibre. Each union
of two corrected defects has mass at least \(3r_*\).
If all three two-class unions are within a small fixed
margin of \(3r_*\), then every edge and every single original
letter on that tuple has outside cubic defect mass within
a comparable margin of \(r_*/3\).

Here typical refers to deleting a power-smaller fraction at
each fixed budget. The exceptional power for secant privacy
need not be uniform in the small-budget limit. We take the
support-removal tolerances and each further thinning loss
small enough in terms of the desired errors. In particular
perform (DC3) first and choose the loss in the (CR1) thinning
smaller than the exceptional power in (DC3). Then residual
letter supports of distinct rows intersect with arbitrarily
small normalized mass on typical tuples. The (PS2) terminal
bound with height tolerance a sufficiently small multiple
of \(L\) gives the same overlap bound between different
secants and between a secant and an original letter,
off the common projective modulus. For two edges sharing
an anchor, fix the anchor and the other edge and use its
support as the prescribed pool. For a third row's letter,
fix that row and one anchor of the edge. For an endpoint
letter use the anchor's support, since a secant prime off
the fixed factors makes the two endpoints equal projectively
modulo it, and hence vanishing of a complete letter at one
endpoint implies vanishing of the same letter at the other.
These are the pool arguments used for (CR1). Dividing raw
secants by fixed factors and pairwise residual gcds affects
the outside defect masses by at most the prescribed overlap
allowance.
[Removing concentrated support (DC3)](ABC_AFFORDABLE_FIBRE_DETERMINANT_COMPLEXITY.md#removing-concentrated-support-from-a-large-fibre)
· [Private secant packets (CR1)](ABC_CROSS_RATIO_FRESH_PRIME_OBSTRUCTION.md#exclusive-prime-packets-on-bounded-tuples)

Write \(e_{ij}\) for the outside defect mass of secant \(S_{ij}\)
and \(d_{i,U}\) for that of original letter \(U\) in row \(i\).
Put \(s=\sum_{i<j}e_{ij}\) and \(g_i=\sum_U d_{i,U}\).
By the preceding overlap conclusions the three-class union mass
is \(s+\sum_i g_i\) up to arbitrarily small error. The two-class
union in roles \(U,V\) is \(s+\sum_i(d_{i,U}+d_{i,V})\) up to
the same type of error. In fact a private prime has a contribution
from exactly one original letter or exactly one secant. Its
nonzero residue modulo three survives in the appropriate
corrected class. In the secant case it appears in all the
corrected classes and is charged once.

By (CR1), \(e_{ij}\ge r_*/3\) up to the fixed allowance.
Also \(g_i\ge r_*\) typically. Indeed if the union cubic
defect of the original \(abc\) is at most \(n^{r_*-\eta}\),
one can fix the three cube-free coefficients at that cost
up to subpower, then count primitive roots of the diagonal
cubic using Heath-Brown's uniform theorem. This is a
power-smaller original-row set for every fixed \(\eta>0\),
and pairwise coprimality makes the union defect mass equal
to \(g_i\). Thus \(s+\sum_i g_i\ge4r_*\) with the stated errors.
[Uniform diagonal cubic bound, Theorem 1.3](https://arxiv.org/html/2601.07817v1#S1)

For each ordered edge \(i,j\), (SD2) with one anchor and
two selected query letters implies \(e_{ij}+d_{j,U}+d_{j,V}
\ge r_*\) typically, again allowing the arbitrarily small
error from gcd stripping and the count exponent. Sum this
over an oriented triangle to obtain the two-class floor.
Finally sum the three two-class identities. The result is
\(3s+2\sum_i g_i\ge9r_*\) up to the same errors. If all three
unions are within margin \(\kappa>0\) of \(3r_*\), the edge
and vertex lower bounds force \(e_{ij}=r_*/3+O(\kappa)\)
and \(g_i=r_*+O(\kappa)\), with the auxiliary errors chosen
small compared with \(\kappa\). In the one-anchor estimate
the query's two selected letters then have mass at least
\(2r_*/3-O(\kappa)\) for each role pair. Subtract from \(g_i\)
to get an upper bound \(r_*/3+O(\kappa)\) on each single role.
The equality for \(g_i\) gives the matching lower bounds.
[Secants as linear forms, one-anchor count (SD2)](ABC_SECANT_VALUATION_DEFECT_GRAPHS.md#secants-as-linear-forms-on-a-primitive-pair)

Thus the three-class cube record cannot be cheap on the
typical tuples of an affordable fibre. A strictly cheaper
two-class relation would also have to select from the
exceptional tuples in these counts or use a different
arithmetic record.

## Uncorrected cubic classes on genuine projective triples

**CD4.** There are constants \(c,\kappa>0\) such that for
all sufficiently large \(X\), at least \(cX^\kappa\)
distinct ordered triples of positive primitive rows of
height at most \(X\) have distinct projective directions
and satisfy
\[
  \frac{a_i/c_i}{a_j/c_j},\
  \frac{b_i/c_i}{b_j/c_j}\in\mathbb Q^{\times 3}
                                \qquad (1\le i,j\le3).
\]
Thus both uncorrected threefold quotient products are exact
cubes, even though the joint uncorrected record has a
power-sized fibre. In the shape-free cubic count all these
tuples have the same projective shape. No small original
radical or common-gcd caps are asserted here.

Consider the smooth Fermat cubic
\[
                        E: U^3+V^3=W^3
\]
with the rational flex \(O=(1:-1:0)\) as its origin.
Over \(\mathbb Q(t)\), put \(u^3=t,\ v^3=1-t\), and take
\(P=(u:v:1)\) over \(\mathbb Q(t,u,v)\). Choose the branches
real for \(0<t<1\). Adjoin \(\zeta_3\). The cube classes of
\(t,1-t\) are independent over \(\mathbb Q(\zeta_3)(t)\),
as their valuations at \(0,1\) show. The Galois group of
the extension by \(u,v\) independently scales \(u,v\)
by third roots of unity. Let \(g\) be the corresponding
diagonal automorphism of the projective cubic. It sends
\(O\) to another flex, so \(g(O)\) is a point of order
dividing three. The map \(h(Q)=g(Q)-g(O)\) is an origin
automorphism and a group homomorphism. In this case the
last assertion follows directly from the chord law,
since \(g\) preserves collinearity. It follows that
\([m]g=g[m]\) when \(m\equiv1\bmod3\).

Use \(m=1,4,7\) and write \(u_m,v_m\) for
the \(U/W,V/W\) coordinates of \([m]P\). The group-law
formulas are defined over \(\mathbb Q\). By equivariance,
the diagonal Galois action on \(u_m,v_m\) is the same
as its action on \(u,v\), respectively. Thus \(u_m/u\)
and \(v_m/v\) belong to the fixed field
\(\mathbb Q(\zeta_3)(t)\). They are also real on the
chosen branch and hence are in \(\mathbb Q(t)\). Therefore
\[
             x_m=u_m^3,\qquad
             x_m/t,\ (1-x_m)/(1-t)\in\mathbb Q(t)^{\times3}.
\]
These assertions use rational functions of \(t\) and
specialize wherever their denominators are nonzero.

As \(t\downarrow0\), \(P\) approaches the rational flex
\(T=(0:1:1)\), of exact order three. The coordinate \(U/W\)
is a real local parameter there. Since \([m]\) fixes \(T\)
and its derivative on the tangent space there is
multiplication by \(m\), one has \(u_m\sim m\,t^{1/3}\).
Thus \(x_m\sim m^3 t\). In particular the three rational
functions \(x_1,x_4,x_7\) are distinct and in \((0,1)\)
for all small enough positive \(t\). Taking \(x_m=a_m/c_m\)
in lowest terms and \(b_m=c_m-a_m\) gives three actual
positive primitive rows with the indicated relative
classes. Substitute \(t=1/j\) for positive integers \(j\)
away from a bounded set. All three denominators are
\(O(j^D)\) for some fixed \(D\), by clearing the fixed
rational-function denominators. Since \(x_1=t\), the tuples
are distinct. Letting \(j\) run up to a small fixed
multiple of \(X^{1/D}\) proves the power-sized count.

In particular the missing corrected pole class in
(PM5) takes power-many values on this family. With both
uncorrected quotient products cubes, the other two
corrected classes equal the corrected pole class. Fixing
it would give subpower completion by (PM5). This is a
cubic counterpart to the two-class obstruction for parity
quadruples.
[Shape-free corrected cubic completion (PM5)](ABC_PROJECTIVE_SHAPE_SECANT_CLASS_COUNTS.md#the-cubic-case-with-no-shape-record)
· [Parity obstruction at fixed shape (PSO1)](ABC_PROJECTIVE_SHAPE_TWO_CLASS_OBSTRUCTION.md#positive-primitive-quadruples-from-odd-multiples)
