# Exact Eisenstein norm powers under the original radical caps

Use \(n=X^L\) with height first at fixed \(L\), and
primitive positive rows \(a+b=c\in(X/2,X]\).
Write \(N=a^2+ab+b^2\). For a fixed order \(q\ge2\)
and a total original cap \(\operatorname{rad}(abc)\le n^B\),
the count with \(N=w^q\) has limiting coefficient
\[
 \frac{B}{(2a_0+1)\tau(s)},\qquad
           3q=3^{a_0}s,\quad 3\nmid s.                  \tag{A}
\]
In particular the coefficients at orders \(2,3,4\)
are \(B/6,B/5,B/9\). Uniformly over adaptive orders \(q\ge2\),
the proper-power sector has limiting coefficient \(2B/9\).
These bounds use a total original cap. The further defect
bound below uses three individual original caps.

For a limiting coefficient \(D\), the convention throughout
is a bound \(O_{B,L,\eta,\varepsilon}(n^{D+\eta}X^\varepsilon)\)
for each \(\eta>0\) and sufficiently small fixed \(L\).
In fixed-order claims, the threshold on \(L\) and constants
may depend on the fixed order. All exponent allowances are
chosen before taking \(X\to\infty\).

## EN1. Exact factorization and counting

Let \(\omega^2+\omega+1=0\), put
\(\alpha=a-b\omega\), and use the Euclidean property of
\(\mathbb Z[\omega]\) as in CP4. The ideals of \(\alpha\)
and \(\bar\alpha\) are coprime outside 3: a common prime
ideal would divide both \(a-b\omega\) and \(a-b\bar\omega\),
hence \(b\) and \(a\) at an unramified prime, contrary
to primitivity. If 3 divides \(N\), then \(a\equiv b\bmod3\)
and both are units at 3. The identity
\[
 \prod_{j=0}^2
  (\omega^j\alpha-\omega^{-j}\bar\alpha)
       = \alpha^3-\bar\alpha^3
       = \pm(\sqrt{-3})^3abc .
\]
Each difference divided by \(\sqrt{-3}=\omega-\omega^2\)
is one of the complete letters with a sign.
Since a unit and its conjugate have the same cube,
\[
              \beta^{3q}-\bar\beta^{3q}
                          =\pm(\sqrt{-3})^3abc.         \tag{B}
\]

First use (B). The \(\Phi_d(\beta,\bar\beta)\), for
\(d\mid3q\), give distinct homogeneous factors, each
fixed at the chosen \(q\). If \(3\nmid d\), the
whole factor is rational after division by a nonzero
constant where necessary (in particular \(d=1\)
is proportional to \(y\)). If \(3\mid d\), split the
roots of \(\Phi_d\) into the two classes
\(j\equiv1,-1\bmod3\). Each product of
\(\beta-\zeta_d^j\bar\beta\) has coefficients in
\(\mathbb Q(\omega)\). Conjugation replaces each factor
by \(-\zeta_d^{-j}\) times itself, so the conjugate of
the product is a constant multiple of the product.
Divide by any nonzero coefficient of that product to
give a rational binary form, and clear denominators.
The subgroup fixing \(\mathbb Q(\omega)\) is transitive on
each set of roots and on all roots when \(3\nmid d\).
Thus these rational forms are irreducible. Their
geometric root sets are disjoint, as are those for different
\(d\), since the change from \((x,y)\) to
\((\beta,\bar\beta)\) is invertible over the algebraic
closure. Thus all resultant primes and normalization
primes can be put in a fixed finite pool. Away from
that pool, primitivity of \((x,y)\) makes the supports
of the factor values pairwise disjoint, and the union
of the supports divides the original \(abc\). The
number of these factors is \(T=(2a_0+1)\tau(s)\).

There are three linear factors for odd \(q\), from
\(d=1,3\), proportional to \(x,y,x+y\). For even \(q\)
the indices \(2,6\) add \(2x+y,x-y,x+2y\). Every
three of the resulting distinct directions yield an
additive primitive triple after applying their fixed
linear relation, signs, and bounded common gcd.
Its height is comparable to \(X^{1/q}\), since any
two of the fixed independent directions control the
height of \(\beta\). It recovers \((x,y)\) at bounded
cost and its radical cap is the sum of the three
form caps up to fixed primes. Choose exact radical
bins at each fixed budget, and apply MCB2 or the
one-third bound in the relative budget \(qL\). Thus
each three linear forms have combined normalized
radical cap at least \(3d\) at row count exponent \(d\).
Average over all such triples to get \(k_{\rm lin}d\)
as the lower bound for the total linear cap.

A fixed exact value of each nonlinear factor determines
\((x,y)\) at \(X^\varepsilon\) cost. After a fixed
change of coordinates and clearing fixed denominators,
the irreducible form is a constant times the norm of
a linear expression in a fixed field. At a fixed value,
its integral ideal divides a fixed multiple of that value
and has subpower many possibilities. Every embedding
of a candidate generator has polynomial height. The
norm is a nonzero integer, giving reciprocal polynomial
lower bounds for its embeddings. For each ideal the
logarithmic unit lattice therefore gives polylogarithmically
many generators in this range. An exact generator
recovers the two coordinates since the degree is at
least two. Consequently each
nonlinear factor costs normalized radical at least
\(d\), by supported-power enumeration of its exact
polynomial-height values. Summing the linear and
nonlinear charges gives \(B\ge Td\) and proves (A).

## EN2. Uniformity over exact power orders

If \(q\) has a fixed divisor two or three, take
that divisor. Otherwise take a prime \(p\ge5\)
dividing it. The unramified \(p\)-th power case has
\(T=6\), including exactly three linear forms. Even
without applying an abc count at relative budget
\(pL\), fixing any two linear values recovers
\((x,y)\), so the sum of the three linear caps is
at least \(3d/2\). Add the three nonlinear charges
\(3d\) to get \(d\le2B/9\). This uses no smallness
assumption on \(pL\). At each fixed budget one
need apply this argument only for \(p\le P(L)\),
with a finite cutoff large enough in terms of
\(B,L\) and the desired allowance. For larger
orders, enumerate \(x,y\ll X^{1/p}\) directly
at exponent \(2/p\), which is at most the desired
fixed positive multiple of \(L\) for this cutoff.
There are \(O(\log X)\) possible orders for the
nontrivial norms. For the count up to \(X\), sum the
shell estimates with any fixed allowance in the exponent.
The constants for the finite
order list are taken at the fixed-budget height
limit.

## EN3. Retaining a small radical of norm defects

For \(q\ge2\), put
\[
 D_q(N)=\prod_{\substack{p\mid N\\q\nmid v_p(N)}}p .
\]
Under the three individual caps \(n^{1/3+o(1)}\) at the
balanced small-budget limit, a bound \(D_q(N)\le n^f\)
with fixed \(f<1/6\) eliminates saturation, uniformly
over all choices of the order \(q\). For each fixed
\(f<1/3\) there is also a fixed \(Q_*(f)\) such that
\(D_q(N)\le n^f,\ q\ge Q_*(f)\) eliminates saturation.
Concretely each statement provides positive constants
\(\delta,\kappa,L_0\), depending on \(f\), such that the
count with three caps \(n^{1/3+\kappa}\) and
\(0<L<L_0\) is \(O_{f,L,\varepsilon}
(n^{1/3-\delta}X^\varepsilon)\).

Here is the arithmetic count used in both assertions.
A prime dividing \(abc\) cannot divide \(N\), by
primitivity and the congruence obtained by setting
the corresponding letter to zero. Split the ideal
of \(\alpha\) into quotient and remainder exponents
modulo \(q\). The UFD gives
\[
                          \alpha=\delta_0\beta^q,
\]
with an exact integral coefficient \(\delta_0\)
including the unit choice. Its support is over
\(D_q(N)\), since each split rational prime selects
just one of the conjugate primes and the ramified
exponent is at most one. For a fixed defect radical
there are height-subpower choices of the orientations
and all exact coefficient powers, uniformly up to
polynomial height. Thus coefficients cost
\(n^fX^\varepsilon\). The root is primitive and
\(|\beta|\ll X^{1/q}\).

Choose a full divisor \(D\) from any one original
letter, with \(X^w\le D\) and radical at most \(n^\mu\).
Its choices cost \(n^\mu X^\varepsilon\).
The quotient \(\alpha/\bar\alpha\) is a specified
cube root of unity modulo \(D\mathbb Z[\omega]\).
For the original positions \(a,b,c\) these are
\(\omega^2,1,\omega\), respectively. Both
\(\delta_0\) and \(\beta\) are units modulo \(D\),
so \((\beta/\bar\beta)^q\) is fixed modulo \(D\).
For bounded \(q\), the number of possible unit ratio
classes is height-subpower. At unramified primes it
is bounded locally by a constant depending on \(q\)
and the field, using cyclic finite-field units
and the principal unit filtration. At the finitely
many primes dividing \(3q\), the number of \(q\)-torsion
classes stays bounded as the modulus depth grows, by
the local logarithm on a subgroup of fixed finite index.
If two roots have the same class, their cross-product
difference \(\beta\bar\beta'-\beta'\bar\beta\) belongs to
\(D\mathbb Z[\omega]\) and has absolute value
\(O(X^{2/q})\). For \(qw>2\) by a fixed margin the
difference vanishes. The roots are rationally proportional
and primitive, hence equal up to sign. This proves the
count with exponent
\[
                             f+\mu+\eta .                \tag{C}
\]
The arbitrary normalized allowance \(\eta>0\) absorbs
adaptive orders above a finite cutoff at each fixed
\(L\): enumerate \(\beta\) in a box of height
\(O(X^{1/q})\) when \(2/q\le \eta L\). Nonunit roots
require \(q=O(\log X)\). A unit root has its whole
original norm support in the coefficient, and can
be counted directly at cost \(n^fX^\varepsilon\),
without enumerating unbounded orders.

There is a second count for each fixed \(q\), without
requiring \(qw>2\). Write the signed original letters
as degree-\(q\) homogeneous integer forms in the two
coordinates of \(\beta\). Their root sets are the three
distinct sets obtained by requiring
\((\beta/\bar\beta)^q\) to equal a fixed cube root times
\(\bar\delta_0/\delta_0\). Their roots are simple and
the three sets disjoint outside primes dividing
\(3q\,N(\delta_0)\). The latter norm is automatically
coprime to the bought divisor. Omit the fixed primes
dividing \(3q\) when choosing \(D\).
Also \(|N(\delta_0)|\le n^{(q-1)f}\), so the form
coefficients satisfy the polynomial-in-\(n\) bound
at fixed \(q\).

Let \(S\) be either other role and let \(P_S\) be
its form. Of the \(q+1\) charts with first column
\((i,1)\), \(0\le i\le q\), and second column \((1,0)\),
at least one has leading \(P_S\) coefficient a unit
at each prime of \(D\) larger than \(q\). Omit the
other fixed small primes at purchase. At least one
chart retains product of assigned full prime powers
of size at least \(D^{1/(q+1)}\). For fixed \(D\) and the
coefficient, the assignment is deterministic using
the first eligible column. On the retained submodulus
the bought form has at most \(q^{\omega(D)}\)
projective root classes, by simple-root Hensel lifting
in the two local coordinate charts. The tested value
\(P_S\) is a unit by original primitivity. Apply BR1,
after fixed sign changes if needed to make the
root coordinates positive. Primitive inputs with a
zero coordinate have only bounded possibilities.
Writing \(r_S\) for the original radical cap on \(S\)
in base \(n\), the normalized bound is
\[
                             f+\mu+r_S/2+\eta .          \tag{D}
\]
The choice of \(D\) is charged once. Its retained
submodulus depends on the same records.
[Binary-form ray estimate BR1](ABC_BINARY_FORM_RAY_AMPLIFICATION.md#br1-primitive-inputs-with-prescribed-projective-residues)

To deduce the uniform claims, argue by a sequence
at balanced saturation with vanishing positive
cap tolerance and vanishing exponent loss. On a
limiting profile GE1 supplies a positive compact
nonatomic window on some original letter. It
supplies full divisors of arbitrarily small fixed
price \(\mu\) and positive fixed depth \(w\) after
refining. Use the window away from the finitely
many bad primes at bounded \(q\). For \(f<1/6\),
any sequence of bounded orders has a constant-order
subsequence to which (D) applies with positive
slack. Any divergent order subsequence satisfies
\(qw>2\) with slack and (C) applies. For \(f<1/3\)
and orders required to be at least \(Q_*\), failure
for every \(Q_*\) would give a divergent order
sequence and the same contradiction. In this
argument one can first pick the window away from
any fixed finite prime support in the bounded
order case: its atom-free mass can be retained
outside the boundedly many large-price primes,
and small-price primes with fixed rational
names have vanishing total raw price. Counting
orders at each fixed \(L\) uses the finite cutoff
in (C) before the height limit.
[GE1 and the positive nonatomic window](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees)

## EN4. The cyclic pair count

For a fixed original row the primitive cyclic determinant
output recovers a second row at subpower cost. The
unreduced determinant pair is the image of the two
primitive coordinates of row two under an integer matrix
of determinant \(N_1\). Its content \(g\) divides \(N_1\)
by the adjugate and has subpower choices. The bounded
signed-output choices and exact \(g\) then recover the
second row. The primitive output norm is \(N_1 N_2/g^2\).

If the first-row catalogue costs \(n^{d_0}\), the
primitive output has height in \([X^\nu,X^B]\) for
fixed positive \(\nu,B\), its radical cap is \(n^s\)
and its norm is an exact fixed \(q\)-th power, the
pair exponent is at most \(d_0+s/T+\eta\), with \(T\)
from (A). Use output height shells to apply EN1.
For adaptive exact power output the corresponding
bound is \(d_0+2s/9+\eta\). Without a lower height
floor one can use \(s/(T-k_{\rm lin}/2)\) in the fixed
order bound, since any two independent linear values
recover a root pair. The criterion for beating the
original full pair exponent is that the stated
pair price, including any additional records,
sits strictly below \(2/3\). An endpoint extraction
forcing this norm condition remains open.
[Cyclic determinant transfer CP4](ABC_CYCLIC_PROJECTIVE_PAIR_TRANSFER.md)
