# Additive patterns in a fixed affordable divisor fibre

A short relative additive pattern on the targets forces the same
affine pattern on both summands in a fixed divisor fibre. Counts for
the three normalized patterns then determine the original triples.
The counted records use the actual prime valuations of all three
letters.

## Relative patterns and exact affine lifting

Fix \(X\ge2\), coprime positive integers \(F,G\), and signs
\(\sigma,\tau\in\{-1,1\}\). Consider ordered \(m\)-tuples of positive
primitive triples
\[
 U_i=F u_i,\qquad V_i=G v_i,\qquad
 T_i=\sigma U_i+\tau V_i>0,\qquad
 (U_i,V_i)=1,\qquad U_i,V_i,T_i\le X.
\]
The target values are distinct. Their relative additive pattern is
the unique tuple of distinct integers \(\boldsymbol\ell=(\ell_i)\)
such that
\[
 T_i=T_1+d\ell_i,\qquad d>0,\qquad
 \ell_1=0,\qquad \gcd_i|\ell_i|=1 .
\]
Here \(d=\gcd_i|T_i-T_1|\). An arbitrary fixed rational relative
pattern can be put in this form by subtracting its first coordinate,
clearing denominators and dividing by the gcd. Its integral
height is \(H=\max_i|\ell_i|\).

Assume \(v_i\le V_0\) and
\[
                              2 H V_0<F .                    \tag{AFP1}
\]
The same argument applies with \(U,F\) and \(V,G\) exchanged.
For example, an affordable pair of full divisors with
\(F G\ge X^{h+\eta}\) and \(V_i\le X^h\) satisfies (AFP1)
when \(H\le X^{\eta-\delta}\), for any fixed \(\delta>0\)
and sufficiently large \(X\).

Reducing the additive equation modulo \(F\) gives
\(\tau G(v_i-v_1)\equiv d\ell_i\pmod F\). Therefore
\[
       F\mid \ell_j(v_i-v_1)-\ell_i(v_j-v_1).
\]
The absolute value is at most \(2HV_0\). By (AFP1) it is zero,
so \(v_i\) is an affine function of \(\ell_i\).
The signed relation makes \(u_i\) affine as well. For each
\(K\in\{U,V,T\}\), write
\[
         K_i=a_K+b_K\ell_i,\qquad
         g_K=\gcd_i K_i=\gcd(a_K,b_K),\qquad
         k_i=\frac{K_i}{g_K}=A_K+B_K\ell_i.                   \tag{AFP2}
\]
The intercept is integral because \(\ell_1=0\), and the slope is
integral by Bézout and \(\gcd_i|\ell_i|=1\). The same identity
gives the equality of the two gcds in (AFP2). The normalized
vector \((k_i)\) is positive and primitive. A constant vector
has all \(k_i=1\). Otherwise \((A_K,B_K)\) is a coprime pair
with \(A_K>0,\ B_K\ne0\), and both entries have height at most
\(X\).

## A count by the joint defect record

Fix \(\boldsymbol\ell\), \(F,G\) and a prime pool \(S_0\) with
\(\prod_{p\in S_0}p\le X^C\), for fixed \(C\).
Let \(S\) also include the primes dividing the nonzero
differences \(\ell_i-\ell_j\). This pool has polynomial height
uniformly for fixed \(m,C\), since \(H\le X\) on any occurring
pattern. For an integer \(b>0\), write
\(\operatorname{rad}_q(b)=\prod_{p:\,q\nmid v_p(b)}p\), and
let \(\operatorname{sf}_q(b)\) be its \(q\)-th-power-free part
with exponents in \(\{0,\ldots,q-1\}\).

Consider either of the following cases, with \(B\ge1\).

- For any fixed even \(m\ge4\), put
  \(D_K=\operatorname{sf}_2(\prod_i K_i)\) on each letter.
- For \(m=3\) and any fixed \(q\ge3\), put
  \(D_K=\operatorname{sf}_q(\prod_i k_i)\) on each letter,
  using the normalized entries in (AFP2).

In both cases impose the joint cap
\[
             \prod_{\substack{p\mid D_U D_V D_T\\p\notin S}}p
                                   \le B .                   \tag{AFP3}
\]
Then the ordered tuple count satisfying (AFP1) is
\[
                       O_{m,q,C,\varepsilon}(B X^\varepsilon),
                                                               \tag{AFP4}
\]
with \(q=2\) in the even case. This charges the union of the
defect supports across the letters. The original radical caps
may all be imposed in addition. The three common contents
\(g_U,g_V,g_T\) are allowed to vary without a separate cost.

There are at most \(B X^{o_X(1)}\) joint records
\((D_U,D_V,D_T)\). Choose the outside product of primes in
(AFP3), then assign the exponents for each occupied prime
on the three letters. For fixed \(q\), this costs at most a
fixed constant to the power of the number of primes involved.
The parts on \(S\) cost subpower by the same estimate. All
records have polynomial height, so their numbers of distinct
prime divisors are \(O_{m,q,C}(\log X/\log\log X)\).

In the even case
\(\operatorname{sf}_2(\prod_i K_i)=\operatorname{sf}_2(\prod_i k_i)\),
since \(g_K^m\) is a square. A nonconstant normalized template
therefore gives a rational point on
\[
                  D_K z^2=\prod_{i=1}^m (t+\ell_i),
                  \qquad t=A_K/B_K .
\]
The product is a square times \(D_K\) after dividing by
\(B_K^m\). The rational coordinates have polynomial height.
The forms \(A+\ell_i B\) have distinct projective zeroes.
For \(m=4\) the four-form elliptic argument with rational
\(2\)-torsion bounds the completions by \(X^{o_X(1)}\).
For even \(m\ge6\) the genus-\((m-2)/2\) argument with
rational Jacobian \(2\)-torsion gives the same count.
The point's \(t\) determines the primitive pair
\((A_K,B_K)\) with \(A_K>0\).
[The four-form and higher-genus bounds](ABC_MULTI_TARGET_MIXED_KERNEL_COUNTS.md#the-mixed-kernel-count)

For \(m=3\), two distinct entries of a primitive template
can share primes only in \(S\). In fact
\(k_i-k_j=B_K(\ell_i-\ell_j)\) and
\(\gcd(k_i,B_K)=1\). Off \(S\), the \(q\)-defect supports of
the three entries are disjoint. Thus their individual
defects can be recovered as extra records at subpower
cost once the joint record in (AFP3) is fixed. More
precisely, choose \(q_0\mid q\) to be an odd prime
dividing \(q\), or \(q_0=4\) if \(q\) is a power of two.
The \(q_0\)-th-power-free coefficients \(C_i\) of the
three normalized entries have outside support within the
\(q\)-defect support of their product, and their allocation
and exponents cost subpower. Write \(k_i=C_i z_i^{q_0}\).
The exact linear relation
\[
 (\ell_2-\ell_3)k_1+(\ell_3-\ell_1)k_2
                         +(\ell_1-\ell_2)k_3=0
\]
gives a diagonal equation in the three positive roots
with fixed nonzero signed coefficients of polynomial
height. The roots are primitive because the \(k_i\) are.

For \(q_0=3\), the uniform diagonal cubic bound gives
\(X^{o_X(1)}\) completions. For \(q_0=4\), the signed
diagonal quartic maps to the bounded-height elliptic
curve with rational \(2\)-torsion used in (PT6).
For an odd \(q_0=p\ge5\), the product
\(k_1k_2k_3^{p-2}\) is a fixed coefficient times a
\(p\)-th power. Writing the third of the three linear
forms as a nonzero combination of the other two gives
the cyclic twist (PT2) and its subpower bound (PT3).
In all cases two normalized entries recover
\((A_K,B_K)\). If one chooses a fixed odd prime \(p\)
at the outset, one may alternatively charge
\(\operatorname{rad}_p(K_1K_2K_3^{p-2})\) on each
original letter in (AFP3). The weights sum to \(p\),
so the contents drop out, and the three weights are
units modulo \(p\), so off \(S\) the same defect is detected.
[Uniform cyclic, cubic and quartic completions](ABC_THREE_TARGET_PRIME_ORDER_STARS.md#the-unweighted-product-for-any-fixed-order-at-least-three)

After the three templates have been fixed, the contents
are a positive integer solution of
\[
       \sigma g_U (u_i^{*})+\tau g_V(v_i^{*})=g_T(t_i^{*})
                              \quad(1\le i\le m),
\]
where the starred vectors are the three normalized templates.
The columns have rank at most two because they are affine
in \(\ell_i\). If they had rank one, the positive primitive
vectors would coincide. Each common entry would divide
\(\gcd(U_i,V_i)=1\), so all entries would be one and the
targets would be constant. The rank is therefore two.
The kernel is one-dimensional and there is at most one
positive integer vector of contents with
\(\gcd(g_U,g_V)=1\), as required by the original primitivity.
This recovers all original letters and proves (AFP4).

## Long progressions with small-prime defects bound the whole triple height

There is also a vanishing result in large height shells for long
progressions, allowing a number of missing indices. Fix \(k\)
above Pratt's absolute threshold, and
\(2\le q\le(\log\log k)^{1/5}\). Let
\(I\subseteq\{0,\ldots,k-1\}\) satisfy
\[
                   |I|\ge k-\frac{0.2k}{\log k}.
\]
Suppose the actual target values in an affordable pair record
have the form \(T_i=t_0+i d\) for \(i\in I\), with \(d>0\).
Use the same primitive signed triples and the same upper cap
\(V_i/G\le V_0\) as above, and assume \(2kV_0<F\).
On each letter \(K\), put \(g_K=\gcd_{i\in I}K_i\). Suppose
\[
    p\mid\operatorname{rad}_q\!\left(\prod_{i\in I}(K_i/g_K)\right)
                            \quad\Longrightarrow\quad p\le k
\]
for all three letters, where the product is of the normalized
entries \(K_i/g_K\).
Then every complete value in each of the observed triples is
\[
                          O(k^2(\log k)^{9/2}).              \tag{AFP5}
\]
If \(|I|\ge k-k/(\log k)^3\), the bound improves to \(O(k^2)\).
Both implied constants are absolute. For example a fibre with
maximum height comparable to \(X\) contains no such pattern if
\(k^2(\log k)^{9/2}=o(X)\) and the lift condition holds.

Here the defect hypothesis is on a product in each letter.
When \(q\mid |I|\), the three normalized products can be replaced
by the complete products because each content appears to the
power \(|I|\). Another sufficient condition, with no restriction
on \(|I|\bmod q\), is that all complete entries have
valuations divisible by \(q\) outside the primes at most \(k\).
In that case division by the gcd preserves the condition
prime by prime.

Indeed, \(I\) contains adjacent indices because it is more than
half the full interval by a fixed margin. Hence the differences
\(I-I\) have gcd one. Subtract any index in \(I\) and apply the
affine lifting argument, with \(H<k\). Translate the intercept
back to index zero. The three normalized vectors have the form
\[
                      A_K+i B_K,\qquad i\in I,\qquad
                         \gcd(A_K,B_K)=1 .
\]
Their values on \(I\) are positive. If \(B_K=0\), they equal one.
Otherwise, reflect \(i\mapsto k-1-i\) when necessary to make
the step positive. The product of the observed entries is an
exact \(q\)-th power times an integer supported on primes at
most \(k\). If the resulting initial term is zero, primitivity
already forces the sign-adjusted step to equal one. Otherwise
Pratt's Theorem 1.5 bounds the positive step by
\(3(\log k)^{3/2}\) and the absolute initial term by \(k\) times
one plus that step. With at most \(k/(\log k)^3\) missing
indices the step is at most two. Undoing a reflection adds at
most \(k\) times the step to the bound for the intercept.
Thus, writing \(D=1+3(\log k)^{3/2}\), all three templates
satisfy \(|B_K|\le D,\ |A_K|\ll kD\). One can take \(D=3\)
under the stronger missing-index bound.
[Pratt, Theorems 1.4–1.5](https://arxiv.org/html/2611.33008v1#S1)

The signed \(2\)-by-\(3\) matrix of their intercepts
and slopes has rank two, by the same primitivity and
nonconstant target argument used above. Its three signed
maximal minors are integers of size \(O(kD^2)\) and form a
nonzero vector spanning the kernel. Dividing by their common
gcd gives the unique primitive vector up to sign, so the
actual contents are all \(O(kD^2)\). The normalized values
are \(O(kD)\) and the complete values are \(O(k^2D^3)\).
This gives (AFP5) and the improved bound. In particular the
height bound has no dependence on the common contents
before the argument. Pratt bounds a single primitive
progression. The application here uses the simultaneous
lift and rank-two reconstruction to bound the heights of
all three original letters.

## Two patterned letters and a full divisor on the third

The normalized templates of any two letters in a lifted pattern
are linearly independent. If two were proportional, they would
coincide by positivity and primitivity, and the additive relation
would make the third proportional to the same vector. Primitivity
within each triple would force every common normalized entry to
be one, contrary to the distinct targets.

Suppose the small-prime condition for (AFP5) holds on \(V,T\),
and retain an arbitrary \(U\). The modulus \(F\) divides both
affine coefficients of \(T-\tau V\), by the gcd-one condition on
the selected index differences. The contents \(g_T,g_V\) are units
modulo \(F\). Eliminating them gives
\[
                       F\mid A_V B_T-A_T B_V\ne0.            \tag{AFP6}
\]
The Pratt bounds on just \(V,T\) now imply \(F\ll k(\log k)^3\),
or \(F\ll k\) with at most \(k/(\log k)^3\) missing indices.
In an affordable record with \(FG\ge X^{h+\eta}\), \(V_i\le X^h\)
and fixed \(\eta>0\), we have \(F\ge X^\eta\), so this pattern
is impossible when \(k(\log k)^3=o(X^\eta)\) and the lift holds.
Interchanging the free letters interchanges \(F,G\).

For an arbitrary prime pool, use (AFP4) with only the defects
on \(V,T\) charged, and let \(B\) bound their joint outside
support. For complete caps \(V_{\max},T_{\max}\) the count per
shape is
\[
             O_{m,q,C,\varepsilon}\!\left(
                  B X^\varepsilon
                         \left(1+\frac{V_{\max}T_{\max}}{F^2G}\right)
                                \right).                    \tag{AFP7}
\]
Enumerate just the two normalized templates by the argument
for (AFP4). Let \(M_V,M_T\) be their maximum entries. Since the
relative indices include zero and another integer, each
intercept and absolute slope is bounded by its maximum entry.
Thus (AFP6) gives \(M_VM_T\ge F/2\). Write
\(r=g_T,\ s=g_V/G\). They are coprime, satisfy \(r\le T_{\max}/M_T\)
and \(s\le V_{\max}/(GM_V)\), and \(T_i\equiv\tau V_i\bmod F\)
at any fixed index prescribes a unit class for \(r/s\bmod F\).
All normalized entries used here are units modulo \(F\) in any
nonempty record. Distinct pairs have cross determinant divisible
by \(F\). On scaling the two coordinates by their upper bounds,
the vectors have norm at most \(\sqrt2\) and angular spacing at
least \(F\) divided by twice the product of those bounds. This gives
\(O(1+V_{\max}T_{\max}/(FGM_VM_T))\) pairs. Each recovers \(V,T\)
and hence \(U\). Use \(M_VM_T\ge F/2\) and sum over the two template
records to prove (AFP7).

## Shape dictionaries and covered target occupancy

The relative shapes can range over any dictionary of
admissible \(\boldsymbol\ell\). If \(B_{\boldsymbol\ell}\)
is the permitted joint cap for that shape, the ordered
tuple count in the fibre is at most
\[
               O_{\varepsilon}\!\left(
                          X^\varepsilon\sum_{\boldsymbol\ell}
                                       B_{\boldsymbol\ell}\right)
\]
with the fixed theorem parameters understood. More generally
one can sum the actual cardinalities of the joint outside
defect dictionaries, by the same proof. Every distinct
target covered by at least one selected tuple costs at most
\(m\) times this bound in aggregate. In an affordable record
the target determines the primitive pair: (AFP1) gives
\(V_0<F\), so reduction modulo \(F\) fixes \(v\) from
\(T\) and the relation fixes \(u\).

For example, in base \(n=X^L\), if the admissible divisor
pair records have normalized count exponent at most \(t\),
each fibre uses at most \(n^s\) shapes, each joint defect cap is
at most \(n^f\), and the selected tuples cover the targets
under study, the covered original count has exponent at
most \(t+s+f\) up to arbitrary fixed allowances.
A strict margin below \(1/3\) in this sum gives a saving.
The hypothesis concerns actual tuple coverage and the actual
dictionaries. The bound does not give such coverage, or
a small shape or defect dictionary, for unrestricted fibres.

The lattice step is the multi-point version of the rank
drop used for short affine carriers. The common contents
and the offsets here are recovered across all realizations
of a relative shape. They need not be specified as an already
fixed small-parameter carrier in that argument.
[Short affine carriers (AF1)](ABC_AFFORDABLE_FIBRE_AFFINE_RANK.md#sparsity-for-finitely-many-small-integer-parameters)
Related results of Bajpai, Bennett and Chan treat existence
of powerful arithmetic progressions and conditional bounds
on their step and common gcd.
[Bajpai–Bennett–Chan, Theorems 1.1–1.2](https://arxiv.org/pdf/2302.03113)

The lift hypothesis is substantive. For \(F=G=1\) and
positive \(g\) divisible by \(5\), take the four-tuples
\[
 U=(1,1,2g+1,2g+1),\qquad
 V=(2g-1,3g-1,2g-1,3g-1),\qquad
 T=g(2,3,4,5).
\]
All four triples are primitive. The last two gcds divide
\(2\) and \(5\), respectively, and are one by parity and
\(5\mid g\). The relative pattern is \((0,1,2,3)\),
the two summand products are squares, and the target product
has squarefree kernel \(30\). Thus a fixed shape and a
fixed joint defect allow \(\gg X\) actual tuples without
(AFP1). This example makes no claim under a small original
radical budget.
