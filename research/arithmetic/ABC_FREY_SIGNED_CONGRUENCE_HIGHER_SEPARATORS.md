# Higher-precision separators for signed Frey congruences

Use the [single signed newspace](ABC_FREY_COEFFICIENT_PRIME_OLD_ANNIHILATOR.md#one-classical-signed-space-and-one-signed-toric-lattice)
of the Frey form \(f\) at squarefree conductor \(R\), with congruence
index \(r_f\). A sharp bound for this index is equivalent to rational
abc. The [first-order partition bound](ABC_FREY_SIGNED_CONGRUENCE_TANGENT_PARTITIONS.md#local-depths-and-independent-groups)
controls independent first differences from \(f\). A group of
competitors with repeated initial directions can still be separated
at the largest depth in the group. The additional data are the
first differences between competitors and, when those are dependent,
the saturated lattice of their differences.

## Exact linear separators at higher precision

Fix a coefficient prime \(\ell\). Work over the ring of integers
\(\mathcal O\) in a finite splitting field, with uniformizer \(\pi\),
fraction field \(K\), and ramification index \(e\). Let \(T\) be the
signed Hecke order after scalar extension and let \(\lambda_0\) be
evaluation at \(f\). Fix operators \(x_1,\ldots,x_b\in\ker\lambda_0\).
For clarity one may take \(x_j=T_j-a_j(f)\) through the
[short signed coefficient bound (OA8)](ABC_FREY_COEFFICIENT_PRIME_OLD_ANNIHILATOR.md#shorter-coefficient-vectors-for-a-fixed-sign-character).
Those operators span the augmentation ideal additively. A smaller
selection can give an upper bound with fewer columns.

For a nonempty group \(S\) of competitors put
\[
 z_i=(\lambda_i(x_j))_{j=1}^{b}\in\mathcal O^b,\qquad
 m_i=\min_{1\le j\le b}v_\pi(z_{ij}),\qquad M=\max_{i\in S}m_i .
\]
Assume every \(z_i\ne0\), since a zero row prevents a separator
using these operators. Let \(W_S\) be the saturation in
\(\mathcal O^b\) of the span of all \(z_i-z_{i'}\) with \(i,i'\in S\).
For any \(a\in S\) let \(h_S\) be the divisibility exponent of
the class of \(z_a\) in the free module \(\mathcal O^b/W_S\).
Use exponent \(+\infty\) if the class is zero. Equivalently,
\[
     h_S=\sup\{h\ge0:z_a\in W_S+\pi^h\mathcal O^b\}.
 \tag{HS1}
\]
This is the least exponent of a nonzero value on \(f\) of an
integral linear operator in \(1,x_1,\ldots,x_b\) which vanishes
on \(S\). Indeed its nonconstant part must have coefficient vector
orthogonal to \(W_S\), so the possible constant terms are exactly
the evaluations of all integral linear forms on the class of \(z_a\).
Since \(\mathcal O^b/W_S\) is free, these form the ideal
\((\pi^{h_S})\). In particular \(h_S\ge M\). For the full
additive Hecke generators this is the exact idempotent exponent
in the projection to the factors of \(S\cup\{f\}\).

There is a normalized form using a smaller matrix. Make \(B\) the
matrix with rows \(z_i/\pi^{m_i}\), put \(w=(\pi^{M-m_i})_{i\in S}\),
and write \(\Delta_q(B)\) for the minimum valuation of its
nonzero \(q\times q\) minors, where \(q=\operatorname{rank}_K B\).
If \(w\notin\operatorname{im}_K B\), there is no separator of
this form. Otherwise
\[
              h_S=M+\Delta_q(B)-\Delta_q([B\mid w]).
 \tag{HS2}
\]
In fact a value \(\pi^{M+t}\) is possible exactly when
\(\pi^t w\in B\mathcal O^b\). Adding \(w\) as a column gives
a cyclic enlarged lattice quotient, whose length is the
difference of the two determinantal valuations. In particular
the largest elementary divisor of \(B\) bounds the increment
above \(M\). A dual certificate for the increment is
\[
 h_S-M=
     \max_{y\in K^S,\ y^{\mathsf T}B\in\mathcal O^b}
           \max(0,-v_\pi(y^{\mathsf T}w)).
 \tag{HS3}
\]
This follows by duality on the saturated \(K\)-span of the
columns. For the full additive operators, a witness \(y\) gives
a lower bound for the genuine simultaneous index and can use a
small group \(S\). Restriction to the factors of \(f\) and \(S\)
makes that local exponent a lower bound for the entire signed
exponent.

The same calculation works with any finite list of monomials
of positive degree in the \(x_j\) as additional columns.
In a projected algebra of rank \(r=|S|+1\), the monomials
through degree \(r-1\) in algebra generators span the entire
order additively. Their reductions span a filtration of the
rank-\(r\) special fiber. If the span stops growing, it is stable
under multiplication by the generators and is already the whole
algebra. It starts in dimension one, so degree \(r-1\) suffices
by Nakayama. The arithmetic here uses the valuation of the
separator's value as well as its available degrees.
[Kreuzer–Long–Robbiano, Proposition 3.7 for the separator degree over a field](sources/kreuzer-long-robbiano-cayley-bacharach-1804.09469.pdf).

## Using the mutual differences between competitors

Choose \(a\in S\) with \(m_a=M\). For any tree connecting the
competitors in \(S\), normalize the difference on each edge
\(ij\) by the least valuation of a nonzero entry of \(z_i-z_j\).
If these \(|S|-1\) normalized edge rows have independent
reductions, they generate the saturated lattice \(W_S\).
If \(z_a/\pi^M\bmod\pi\) is independent of those edge rows as
well, then
\[
                             h_S=M .
 \tag{HS4}
\]
The edge vectors span the rational difference space because
the graph is a tree. Independence mod \(\pi\) makes their span
saturated. Formula (HS1) now proves (HS4). Thus all competitors
can have the same first direction from \(f\), provided the
first nonzero mutual differences resolve the collision
transversely.

For two competitors with equal depth \(m\), put
\(A=z_1/\pi^m,\ B'=z_2/\pi^m\). Suppose they are independent
over \(K\), write \(c=\min_j v_\pi(A_j-B'_j)\), and let
\(d\) be the minimum valuation of the coefficients of
\(A\wedge B'\). Then
\[
                             h_{\{1,2\}}=m+d-c .
 \tag{HS5}
\]
This is (HS2) with target \((1,1)^{\mathsf T}\) and \(d\ge c\).
The increment vanishes precisely when the first nonzero
normalized mutual difference is independent of \(A\bmod\pi\).
If the original depths differ and \(A,B'\) again denote the
separately normalized rows, the extra exponent is \(d\) instead.
One coordinate of the target is then \(1\), the other is
in \(\pi\mathcal O\), and the augmented matrix has a unit
minor of size two.

For \(\ell\ge5\) and full additive generators, the first
nonzero difference on an edge is also a leading adjoint
class. Indeed the two corresponding Galois specializations
are congruent over the jointly integral representation up to
the mutual depth. Dividing their difference by that power
of \(\pi\) gives a cocycle at the common residual representation.
The trace detection in [(TD7)](ABC_FREY_SIGNED_CONGRUENCE_TANGENT_PARTITIONS.md#adjoint-interpretation-at-the-large-image-primes)
identifies the linear relations of these classes with those of
their normalized Hecke differences. This makes (HS4) a rank test
using mutual deformation classes, including when the classes
obtained by comparing all branches to \(f\) are parallel.

## Combining projections and bounding the missing contribution

For disjoint groups \(S,U\), write \(T_{S,U}\) for the
projected order on \(f\) and both groups. Let \(K_S,K_U\)
be the kernels of projection onto \(f\) and the respective
group, and put \(C=T_{S,U}/(K_S+K_U)\). The exact pullback
\[
 T_{S,U}\simeq
   (T_{S,U}/K_S)\times_C(T_{S,U}/K_U)
\]
follows by taking the difference of any lifts of a compatible
pair. Use full additive operators and put
\(m=\max(h_S,h_U)\). The copies of \(\pi^m e_f\) in the two
projected orders have images in \(C\) with difference \(\eta\).
Both images have the same evaluation at \(f\), and after
inverting \(\pi\) the ring \(C\) is just \(K\). Thus \(\eta\)
is torsion, and the exact simultaneous exponent is
\[
             h_{S\cup U}=m+
                 \min\{j\ge0:\pi^j\eta=0\}.
 \tag{HS6}
\]
This identifies which part of the overlap of the two projections
can contribute beyond the maximum of their exponents.

For any partition into groups, multiplying linear separators
bounds the signed index as
\[
 e\,v_\ell(r_f)\le \sum_{S} h_S ,
 \tag{HS7}
\]
where larger lists of monomials may be used for some groups.
Using (HS5) or (HS4) can improve the first-order independent
partition without relying on the full-length separator.
Conversely (HS3) for a group or (HS6) for two groups gives a
lower certificate, which is needed for a disproof based on
accumulation. To prove abc by this route one needs a uniform
bound \((2+\varepsilon)\log R+O_\varepsilon(1)\) on the sum
of suitably chosen local upper bounds weighted by
\(\log\ell/e\). The geometric conditions above quantify
when higher-order separation helps. They do not currently
control the frequency of groups with a large increment.

## Obstruction at every fixed precision

Fix any precision \(k\ge1\). There are abstract split orders
with the same reduction modulo \(\pi^k\) as augmented algebras,
the same full pairwise depths, and the same normalized
character differences modulo \(\pi^k\) for every pair.
Nevertheless their simultaneous exponents differ by an
arbitrarily large amount.

Take \(a,b\ge k\), \(n\) distinct nonzero residue units \(u_i\),
and \(a+k\le H\le na\). Over a DVR with at least \(n+1\) residue
elements use
\[
 \begin{split}
 x&=(0,\pi^a u_1,\ldots,\pi^a u_n),\qquad
 R_H=\mathcal O[x]+\pi^H\mathcal O^{n+1},\\
 I_H&=\ker(\lambda_0:R_H\to\mathcal O),\qquad
 T_H=\mathcal O\cdot 1+\pi^b I_H .
 \end{split}
\]
The idempotent exponent for \(R_H\) is \(H\), by the
[Vandermonde argument](ABC_FREY_SIGNED_CONGRUENCE_TANGENT_PARTITIONS.md#limits-of-pairwise-depth-and-the-remaining-estimate).
The exponent for \(T_H\) is then \(H+b\), since its augmentation
ideal is obtained by multiplying \(I_H\) by \(\pi^b\).
Every distinct pair has depth \(a+b\), seen on the generator
\(\pi^b x\). The augmented algebra \(T_H/\pi^k T_H\) is the
square-zero extension of \(\mathcal O/\pi^k\) by a free module
of rank \(n\), since \(b\ge k\).

Write \(y=P(x)+\pi^H z\in I_H\) with \(P(0)=0\), which is
possible by subtracting its constant term in \(\pi^H\mathcal O\).
For \(t=\pi^b y\), normalized difference on the pair \(i,j\)
satisfies
\[
 \frac{\lambda_i(t)-\lambda_j(t)}{\pi^{a+b}}
                \equiv (u_i-u_j)P'(0)\pmod{\pi^k},
 \qquad u_0=0.
 \tag{HS8}
\]
Terms of polynomial degree at least two and the
\(\pi^H z\) term vanish after normalization by the inequalities
on \(a,H\). The functional \(y\mapsto P'(0)\bmod\pi^k\)
is well-defined, as (HS8) with \(j=0\) and an \(i\ne0\)
with unit \(u_i\) shows. It is primitive since it takes value
\(1\) on \(x\). It can therefore be made the same first
coordinate on all the augmented reductions. This identifies
all the normalized pair maps together with the reduced
algebras at precision \(k\). Varying \(H\) leaves those data
fixed, while \(H+b\) ranges from \(a+b+k\) to \(na+b\).
The group of matrices over \(T_H\) with determinant in the
constant units has full residual image and generates these
orders by traces, using determinant-one companion matrices
with arbitrary prescribed trace. The construction is in
abstract trace orders and makes no modular realization
assertion.

For a finite check over \(\mathbb Z_7\), take \(n=4\),
\(a=b=k=2\), and \(u_i=i\). PARI Smith normal forms of the
columns \(1,\pi^b x^j\ (1\le j\le n)\) and
\(\pi^{b+H}e_i\ (1\le i\le n)\), before and after adjoining
\(e_0\), give exponents \(6,7,8,9,10\) as \(H=4,5,6,7,8\).
Differences of every pair of coordinate rows have valuation
four throughout.

## A short higher-precision certificate at level 442

In the [level 442 packet](ABC_FREY_SIGNED_CONGRUENCE_TANGENT_PARTITIONS.md#a-signed-modular-packet-at-level-442)
take only the two operators \(T_5,T_9\). Their eigenvalues
on \(f\) and the two competitors are \((4,1),(-2,1),(2,-3)\).
The normalized difference rows at \(2\) are \((-3,0)\)
and \((-1,-2)\). Their first reductions coincide. The
mutual difference divided by \(2\) is \((-1,1)\), transverse
to either initial direction. Equation (HS5) therefore gives
cost one at \(2\) using just these two operators. Globally the
exact relation is
\[
                        1+T_5+T_9=6 e_f
\]
on this packet. At \(3\) one of the other forms agrees with
\(f\) modulo \(3\), and at \(2\) they both agree with it
modulo \(2\). Thus \(6\) is the exact congruence index.
This is a certificate in a genuine signed modular space,
outside the fixed Frey input family.
