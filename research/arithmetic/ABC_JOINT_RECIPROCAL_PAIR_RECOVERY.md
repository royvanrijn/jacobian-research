# Pair recovery in joint abc reciprocal relations

For a primitive positive row \(x=(a,b,c)\) with \(a+b=c\), put
\[
 {\cal R}(x)=(1/a,1/b,1/c),\qquad
 {\cal J}(x)=(c/a,c/b,a/c).
\]
The three components of each vector record the three complete prime
supports separately. A distinguished pair in a fixed weighted equation
for either vector can be recovered after the other rows have been fixed.
Combining this recovery with effective radical counting gives a strict
saving for every fixed order on distinct active rows. The same argument
applies to the two moving coordinates in a complete-letter fibre.

## Recovering the final two rows

Let \(u,v\) be positive integers of polynomial height. For fixed
nonzero coefficients \(\alpha,\beta\) of bounded polynomial height
and a fixed rational \(t\), clearing the coefficient denominators
writes
\[
           \frac{\alpha'}u+\frac{\beta'}v=\frac{A}{B},
           \qquad \alpha',\beta',A,B\in\mathbb Z,\quad B>0.
\]
When \(A\ne0\) one has
\[
            (Au-\alpha'B)(Av-\beta'B)=\alpha'\beta'B^2 .
            \tag{JR1}
\]
If the equation has a solution, all entries have bounded polynomial
height with exponents determined by the original bounds. The divisor
estimate gives \(O_\varepsilon(X^\varepsilon)\) possible ordered pairs
\((u,v)\), uniformly in the coefficients and target. When \(t=0\)
the only possible relation is \(v=-(\beta/\alpha)u\).
This is the signed weighted form of the elementary two-unit-fraction
divisor count.
[Browning and Elsholtz, Theorem 1 and §2](https://www.math.tugraz.at/~elsholtz/WWW/papers/papers33FINAL2013.pdf)

Fix all but two terms in a weighted relation for \({\cal R}\). If two
of the three remainder coordinates vanish, the corresponding two
letters of the two remaining rows are proportional with the same
factor \(-\beta/\alpha\). The additive equation then makes the whole
rows proportional. Primitivity forces identical rows. Thus, for
distinct rows, at least two remainder coordinates are nonzero. The
divisor count in those two coordinates determines all letters of both
rows at height-subpower cost. On a fixed complete \(T\)-fibre it
suffices to use the two moving reciprocal coordinates. Both remainders
cannot vanish: the signed relation between the two moving values has
nonzero constant term \(T\) up to sign, which forces the proportionality
factor to be one. A nonzero remainder in one moving coordinate fixes both
rows at height-subpower cost.

For \({\cal J}\) the recovery is purely algebraic. Write \(u=a_i/c_i\)
and \(v=a_j/c_j\) for the two remaining ratios in \((0,1)\).
Two of the remainder coordinates give
\[
             \alpha u+\beta v=s,\qquad
             \frac{\alpha}{u}+\frac{\beta}{v}=t.
\]
Eliminating \(v\) yields
\[
             tu(s-\alpha u)=\alpha s+(\beta^2-\alpha^2)u .
             \tag{JR2}
\]
For fixed remainders the polynomial is identically zero only if
\(s=t=0\) and \(\beta^2=\alpha^2\). Positivity of \(u,v\) and
\(s=0\) then force \(\beta=-\alpha\) and \(u=v\).
For distinct rows this case is excluded. There are therefore
at most two choices for the last ordered pair: a primitive positive
row is determined by its reduced ratio \(a/c\).

One can also retain just the two projective components
\((F(z),z)\) on the oriented family \(0<z=a/c\le1/2\), where
\[
                    F(z)=z-\frac1z+\frac1{1-z}.
\]
Indeed the two remainders give \(v=(s-\alpha u)/\beta\) and a
fixed value for \(\alpha F(u)+\beta F(v)\). The only finite poles
of \(F\) are \(0,1\), both simple. The resulting rational function
of \(u\) can be constant only if the affine map
\(u\mapsto(s-\alpha u)/\beta\) permutes those poles. This forces
\(v=u\) or \(v=1-u\), neither of which gives a distinct pair
with both ratios at most \(1/2\). Otherwise clearing the four
linear denominators gives a nonzero polynomial of degree at
most five. Thus this two-component variant also gives bounded
pair recovery. Use unlabelled accounting for the original
radical in (JR3), since \(F\) alone has all three prime supports
in its reduced denominator. For example, a bound on the
size of a dictionary for the weighted \(z\)-sum allows one to
sum (JR3) over that dictionary while keeping the scalar
\(F\)-target fixed. On a fixed \(c\)-fibre the auxiliary
equation is a weighted sum of complete first summands.
[Denominator and pole structure of \(F\)](ABC_THREE_CAP_RATIO_SUMSETS.md#a-single-ratio-and-an-observable-separable-by-letters)

## Sequential pricing and a saving at every fixed order

Fix \(B>0\), an integer \(k\ge3\), a coefficient height exponent \(K\),
and \(\eta>0\). Let \(n=X^L\) with \(L>0\) sufficiently small and fixed.
Consider distinct primitive rows \(x_i=(a_i,b_i,c_i)\), \(1\le i\le k\),
in \(X/2<c_i\le X\), with \(\operatorname{rad}(a_ib_ic_i)\le n^B\).
Fix nonzero rational weights of heights at most \(X^K\) and a
rational target vector. Then, for \(J={\cal R}\), \(J={\cal J}\), or the oriented variant \((F(z),z)\),
\[
 \#\left\{(x_1,\ldots,x_k):\sum_{i=1}^k w_i J(x_i)=t\right\}
   \ll_{B,k,K,\eta,L,\varepsilon}
           X^{\{kB/6-B/[3(k-1)]+\eta\}L+\varepsilon}.
           \tag{JR3}
\]
The bound is uniform in the weights and in the target. At \(k=2\)
pair recovery gives \(O_\varepsilon(X^\varepsilon)\), and at \(k=1\)
there is at most one row.

On a fixed complete \(T\)-fibre let the moving letters be \(U,V\). Give
each moving role an arbitrary prescribed prime pool \(P_U^*,P_V^*\)
supported on integers of bounded polynomial height. Suppose the
outside radical caps in base \(n\) are \(\ell_U,\ell_V\), respectively,
write \(B=\ell_U+\ell_V>0\), and assume \(\operatorname{rad}(T)\le n^C\)
for fixed \(C\). For \(J=(1/U,1/V)\), the bound becomes
\[
 \#\left\{(x_1,\ldots,x_k):\sum_{i=1}^k w_i J(x_i)=t\right\}
   \ll X^{\{kB/8-B/[4(k-1)]+\eta\}L+\varepsilon}.
           \tag{JR4}
\]
The constant is uniform in \(T\), the pools and the targets, under
the stated bounds. The \(k=1,2\) conclusions remain valid. In the
balanced primitive target fibre where the two caps are both \(1/3\)
and the pools are empty, the exponent in base \(n\) is
\(k/12-1/[6(k-1)]\) up to the allowances.

For the proof, let a labelled occurrence \((W,p)\) represent a prime
\(p\) in its role \(W\). Use this accounting for \({\cal R}\), since each
component has the primes of its own letter as its reduced denominator
support. For \({\cal J}\), the supports of the components correspond
in order to \(a,b,c\) by primitivity. For the fixed fibre, count only
moving-role occurrences outside the initial prescribed pool in that
role. The total normalized mass \(b_i\) of initially charged
occurrences on a row is at most \(B\).

Enlarge the pool in each role by the supports of all weight numerators
and denominators, and of the target denominator in that role. These
pools have polynomial height, since the target in a nonempty solution
set is a sum of a fixed number of rationals of polynomial height. A
labelled prime outside its enlarged pool cannot occur just once among
the rows, by the unique negative valuation test in the corresponding
component. Write \(m_\pi\) for the multiplicity of a labelled prime
\(\pi=(W,p)\), and use weight \(m(\pi)=\log_n p\). Define
\[
 S=\sum_{\pi\ {\rm in\ the\ added\ pools}}m_\pi m(\pi)
    +\sum_{\substack{\pi\ {\rm outside\ the\ enlarged\ pools}\\m_\pi\ge2}}
                       (m_\pi-2)m(\pi),\qquad
 D=kB-\sum_i b_i .
\]
The first sum concerns only the initially charged occurrences.
For each pair \(i,j\) let \(q_{ij}\) be the mass of outside-pool
labelled primes occurring precisely on those two rows and let
\(q=\max_{i<j}q_{ij}\). Pool incidence and multiplicity at least three
have total incidence at most \(3S\). With \(N_k=\binom{k}{2}\) this gives
\[
                    kB\le D+3S+2N_kq.
                    \tag{JR5}
\]
Both \(D\) and \(S\) are nonnegative.

Choose a maximizing pair and place its two members last in the exposure
order. There are only finitely many pair choices. Charge each labelled
prime outside the enlarged pools on its first occurrence. The total
charge is \((\sum_i b_i-S)/2\). At least \(q\) is charged in the last two
places. For the first \(k-2\) rows count sequentially by their new charges.
The allowed effective radical count charges \(v/3+\eta'\) in base \(n\)
for global new mass \(v\). In the global case no prime is recorded in two
roles within a single row. Taking the union of the earlier primes over
all roles only decreases the effective mass for this upper bound.
In the fixed fibre, the two role charges \(v_U,v_V\) cost
\(\min(v_U,v_V)/2+\eta'\), which is at most
\((v_U+v_V)/4+\eta'\). Here the recorded pools may be kept separate on the
two roles by (FL1). These bounds are uniform after the earlier rows are
fixed. Bin all first-occurrence costs on a fixed finite mesh; the
losses from rounding and the individual allowances can be made smaller
than \(\eta\), with the mesh chosen before the height limit.
[Effective radical pricing (EB4)](ABC_EFFECTIVE_RADICAL_BUDGET.md#counting-by-a-sublogarithmic-prime-partition-sum)
· [Fixed-letter outside radical count (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

The last two rows have only height-subpower cost by pair recovery.
Thus the global exponent in base \(n\) is at most
\(kB/6-(D+S+2q)/6\) before the allowances, and the fixed-fibre
exponent is at most \(kB/8-(D+S+2q)/8\). Since \(N_k\ge3\),
(JR5) gives \(D+S+2q\ge kB/N_k=2B/(k-1)\), proving both bounds.
The proof only counts joint equations. A scalar weighted relation in
one predetermined linear combination can allow cancellation between
the different components.

## Variable targets of small denominator

There is a version in which the target is not fixed. On the original
box in (JR3) or (JR4), let the target denominators in reduced form be at
most \(Q_j\) in component \(j\), with total normalized denominator radical mass at most \(r\) in
the sense specified below. Once the first \(k-2\) rows and the final weight
indices are known, suppose possible target values lie in an interval
of length \(\Delta_j\) in component \(j\). A fixed containing interval
for the last two contributions supplies such an interval by translation.
The height bounds may be enforced through finitely many dyadic boxes
at height-subpower cost. The corresponding estimates are
\[
 \begin{split}
 X^{\{kB/6-B/[3(k-1)]+(1+1/N_k)r/6+\eta\}L+\varepsilon}
      \prod_j(1+Q_j^2\Delta_j) &\quad\text{globally},\\
 X^{\{kB/8-B/[4(k-1)]+(1+1/N_k)r/8+\eta\}L+\varepsilon}
      \prod_j(1+Q_j^2\Delta_j) &\quad\text{on a fixed fibre}.
 \end{split}                                                     \tag{JR6}
\]
The displayed bounds suppress constants of the same scope as before.
Only polynomial-height target data are needed, which is automatic on
a nonempty solution set.

For this version the pools contain the weight supports and the
prescribed initial supports. Do not add the unknown target
denominators. Let \(M_1\) be the total outside-pool mass of labelled
primes having multiplicity one. For each component, a singleton prime
divides its target denominator, so \(M_1\le r\). For labelled prime
accounting use
\[
                    r\ \ge\ \sum_j\log_n\operatorname{rad}(q_j)
\]
where \(q_j\) is the reduced target denominator in component \(j\).
Equivalently this is the logarithm of the product of the individual
radicals, counted with their component labels. The quantity \(S\)
keeps its pool and higher-incidence terms. The total new charge is now
\((\sum_i b_i-S+M_1)/2\), and
\(kB\le D+3S+M_1+2N_kq\).
Thus \(D+S-M_1+2q\ge kB/N_k-(1+1/N_k)M_1\). Count the first
\(k-2\) rows as before. Reduced fractions of denominator at most \(Q_j\)
are spaced by at least \(Q_j^{-2}\). The shown factors bound the number
of target vectors after the first rows have been fixed, and pair
recovery then applies to each.

For the complete reciprocal coordinates and uniformly bounded
absolute weights there is no extra exponent when the targets are
arbitrary integers. In that case \(r=0\), \(Q_j=1\), and the final
two contributions lie in intervals of bounded length. On a box with
all recorded complete letters at least \(X^h\), bounded absolute
weights, and \(Q_j\le n^\tau\) for fixed \(h,\tau>0\), the interval
cost is bounded at sufficiently small \(L\). Only the singleton
penalty in (JR6) remains. In the unlabelled version of the accounting,
used below for repaired reciprocal coordinates, one can instead use
the radical of the product of the target denominators.

## Diagonals and least repairs

For a homogeneous even moment on a set \(\mathcal E\) of rows, first
merge indices using equal rows. A nontrivial \(m\)-against-\(m\)
equality leaves \(k\ge2\) distinct active rows and \(s\) cancelled
groups with \(k+2s\le2m\). The signed weights after merging are fixed
integers of size at most \(m\) in each of finitely many equality
patterns. Each cancelled group costs at most
\(n^{B/3+\eta'}X^{\varepsilon'}\) globally by effective radical
counting, or \(n^{B/4+\eta'}X^{\varepsilon'}\) on the fixed fibre
by (FL1). Thus for every fixed \(m\ge2\), the nondiagonal
count for either global vector is
\[
                    \ll X^{\{mB/3-B/[3(2m-1)]+\eta\}L+\varepsilon}.
                    \tag{JR7}
\]
On the fixed fibre with two moving reciprocals it is at most
\(X^{\{mB/4-B/[4(2m-1)]+\eta\}L+\varepsilon}\).
The \(k=2\) case uses the stronger subpower completion bound.

The diagonal count is \(m!|\mathcal E|^m+O_m(|\mathcal E|^{m-1})\).
Consequently on a sufficiently near-saturating family the
\(2m\)-row moment is asymptotic to the diagonal count with a
power-saving error. This also applies to a scalar encoding at
any prescribed finite moment order. For instance, use
\[
             J_1(x)+M J_2(x)+M^2 J_3(x)
\]
with integer \(M\) larger than twice the uniform bound for the
integer coefficients obtained by multiplying a \(2m\)-term
difference of sums by all its denominator products. This bound
is a fixed power of \(X\), so \(M\) can have bounded polynomial
height. Reducing the resulting identity modulo \(M\), and then
dividing by \(M\) twice, shows that each of the three coordinate
differences separately vanishes.

For a minimally repaired row \((A,B',C)=(Za,Zb,Zc)\), the
projective vector \({\cal J}\) is unchanged. The full repaired
reciprocal vector \((1/A,1/B',1/C)\) also satisfies the global
bounds when \(Z\) has bounded polynomial height and its prime
support is in \(\operatorname{rad}(abc)\). Use unlabelled
prime accounting on the original radical \(abc\), since a
prime can occur in more than one repaired role. A prime with
a singleton original occurrence would still give a unique
negative valuation in at least one of the corresponding
reciprocal equations. The proof of (JR5) and the new-charge
bound then use the original prime support union, with the same
global cost. Pair recovery gives repaired complete values, and
division by their common divisor recovers the primitive quotient.
Two zero remainder coordinates would make the two repaired
rows proportional and thus give the same primitive quotient.
For least target repair \(Z=Z_M(c)\) its required height
bound follows from bounded \(ML\) under the original cap
\(\operatorname{rad}(c)\le n^B\). These statements count
primitive quotient rows with their prescribed repair.
[Least target repair](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)

The fixed-fibre bound applies in particular to the original
primitive representations of a fixed \(r\)-full target \(T\).
For \(p\)-full and \(q\)-full summands one may take
\(L=1/p+1/q\), \(\ell_U=(1/p)/L\) and
\(\ell_V=(1/q)/L\), with \(1/r=O(L)\) to keep the target cap
uniform. A fixed original target also has fixed least target
repair multiplier, so its two repaired reciprocals differ
from the original ones by a single constant factor. The
same argument applies to a fixed-gcd fibre of nonprimitive
representations after dividing by the recorded gcd and
checking the primitive radical caps separately.
[Sums of three powerful numbers in primitive coordinates](https://arxiv.org/html/2608.24512v1)

The divisor factorization in (JR1) is standard. The point of the
bounds is that a terminating pair can be put at the end of the
prime exposure order and pays no new radical cost there. This
uses the joint reciprocal constraints, while scalar denominator
energy estimates use total height or the cardinality of rational
dictionaries. Directly forcing enough nondiagonal joint relations
from the positive-height three-letter graph remains open.
[Bloom and Kuperberg, Theorem 3](https://arxiv.org/html/2312.09021v2)
· [Single projective coordinate and matching residual](ABC_THREE_CAP_RATIO_SUMSETS.md#near-maximal-relations-pair-the-repeated-prime-supports)


## Powerful representations of one target with arbitrary gcd

Fix \(C,K,\eta>0\) and \(k\ge3\). Let \(p,q\ge2\) be fixed with
\(L=1/p+1/q\) sufficiently small, put \(X=N\), and suppose
\(\operatorname{rad}(N)\le n^C\) for \(n=N^L\). Let
\(\mathcal A\) consist of positive representations \(A+B=N\)
with \(A\) \(p\)-full and \(B\) \(q\)-full. No coprimality is
required. For fixed nonzero rational weights of polynomial
height at most \(N^K\), the number of distinct-row solutions to
\[
 \sum_{i=1}^k w_i(1/A_i,1/B_i)=t
\]
is at most
\[
      O_{C,K,k,\eta,p,q,\varepsilon}
        \left(N^{\{k/8-1/[4(k-1)]+\eta\}L+\varepsilon}\right).
        \tag{JR8}
\]
In particular this applies to \(r\)-full targets if \(1/r\le CL\).
The pair count at order two is subpower even without the fullness
conditions.

Fix the tuple of gcds \(g_i=\gcd(A_i,B_i)\mid N\) at divisor cost.
For row \(i\), write \(W_i=N/g_i\), \(u_i=A_i/g_i\) and
\(v_i=B_i/g_i\). Outside the primes of \(N\), the \(u_i,v_i\)
parts are still \(p\)-full and \(q\)-full. Hence their outside
radicals are at most \(W_i^{1/p},W_i^{1/q}\). The prime support
of \(N\) is put in the initial pool on both roles. Away from it
the prime supports of the complete reciprocals are exactly those
of the primitive moving letters.

First assume every \(W_i\ge N^{1/8}\). When pricing a row with
assigned new masses \(d_U,d_V\) in base \(n\), apply (FL1) at
height \(W_i\) with the fixed budget \(8L\) and caps
\[
 \frac{d_U\log N}{8\log W_i},
 \qquad
 \frac{d_V\log N}{8\log W_i}.
\]
These are bounded uniformly. The fixed letter \(W_i\) satisfies
\(\operatorname{rad}(W_i)\le\operatorname{rad}(N)\le W_i^{8CL}\).
The resulting cost in base \(n\) is
\(\min(d_U,d_V)/2\) up to the arbitrarily small normalized
allowance. The prime union proof of (JR4) now applies with
\(B=1\), using the complete reciprocal equations to test for
singleton primes. The last pair is recovered in its complete
values by (JR1) and \(A+B=N\), so its gcds and primitive target
values need not agree.

If some \(W_i<N^{1/8}\), expose those rows first. For each such
row the \(N\)-supported parts of both primitive summands cost
height-subpower to fix. Counting the smaller of the two outside
radical dictionaries and using the fixed \(W_i\) then costs at
most
\[
                   N^{\min(1/p,1/q)/8+\varepsilon'}
                       \le n^{1/16}N^{\varepsilon'}
\]
per row. The residual target for the remaining equations is now
fixed. If at least three tall-quotient rows remain, apply their
bound just proved. Its gap relative to their diagonal exponent
is at least \(1/[4(k-1)]\), and each exposed short-quotient row
used less than \(1/8\) in normalized exponent. If at most two
remain, their count is subpower by pair or single-coordinate
recovery. In this case even the bound using \(n^{k/16}\) is
within the exponent in (JR8) for every \(k\ge3\). Sum over the
finitely many short-quotient position sets and over the gcds.
This proves the stated uniform estimate.

For \(2m\)-row energy the same cancellation argument as in
(JR7) yields a nondiagonal bound
\[
                O\!\left(N^{\{m/4-1/[4(2m-1)]+\eta\}L+\varepsilon}\right).
                \tag{JR9}
\]
Indeed the one-row count at a fixed original target costs
\(n^{1/4+\eta'}N^{\varepsilon'}\) by the tall and short
quotient estimates above, so each cancelled group can be
charged at that rate. Both (JR8) and (JR9) retain the original
reciprocal coordinates and the common target throughout.
