# Full-divisor separation and parity moments on mixed Euclidean triples

Let \(n=X^L\), with \(L>0\) fixed before \(X\) tends to infinity. Work with any set of positive primitive triples \(a+b=c\le X\), optionally imposing all three individual radical caps, a common shell, or a least-repair cutoff. Fix an assignment \(S,T,F\) of the complete letters to orders \(e_S=2,e_T=3,e_F=6\). Interpret the corresponding equation with the prescribed signs as
\[
 A x^2+B y^3=Cz^6.
\]
Use the signed coefficient magnitudes from the exact power-free decompositions; primitivity makes their supports disjoint. For each role let \(P_J\) be the complete full part supported on its defective primes \(e_J\nmid v_p(J)\). Put
\[
 R=\prod_J\operatorname{rad}(P_J),\qquad
 \Delta=A^3B^2C,\quad d=\operatorname{sf}_{\pm}(\Delta),\quad
 w_{\rm def}=\sum_J\frac{\log_X P_J}{e_J}.
\]
Here \(\operatorname{sf}_{\pm}\) is the signed squarefree kernel. In particular \(|d|\mid R\). Work in strata \(R\le n^\rho\), \(|d|\le n^\sigma\). In each triple allow three full divisors \(G_J\) consisting of nondefective prime powers. Let \(b\) be a raw bound
\[
 \prod_J\operatorname{rad}(G_J)\le n^b,
\]
or an actual uniform joint dictionary exponent for the record \((G_S,G_T,G_F)\) conditional on \(R\) and the full defective factor choices, with a separate allowed \(X^\varepsilon\) loss. It is enough that a witnessing assignment of \(G_J\) exists for each row. Suppose uniformly on the selected triples
\[
 w_{\rm def}+\sum_J\frac{\log_X G_J}{e_J}\ge w_0>0.
\]
For the signed squarefree \(d\), write \(h_3(d)=|\mathrm{Cl}(\mathbb{Q}(\sqrt d))[3]|\), with split kernel weight one. Use the signed moment exponent from the [original packing bound](ABC_MORDELL_LOCAL_PACKING.md), extended to \(0\le k\le1\) by the first moment and Hölder:
\[
 f(k)=
 \begin{cases}
 1& (0\le k\le1),\\
 (5k+13)/18& (1\le k\le4),\\
 11/6+(1597/5000)(k-4)& (k\ge4).
 \end{cases}
\]
For fixed real \(k\ge0\) with \(3^{-k}\le\min(1,3w_0)\), the number of selected triples is
\[
 O_{\rho,\sigma,b,w_0,k,L,\varepsilon}
       (n^{\rho+b+(f(k)-1)\sigma}X^\varepsilon).
\]
In particular \(w_0\ge1/9\) gives the count \(n^{\rho+b}X^\varepsilon\). The bound counts primitive quotients, including when the input triples are indexed by their least-repaired targets.

## MP1. Packing with additional full divisors

Fix the squarefree \(R\). Each prime of \(R\) has at most three role assignments and at most five possible nonzero valuation residues on its role. Fix the corresponding full defective parts \(P_J\) and powers on \(2,3\) at uniform height-subpower cost, by fixed-support factor enumeration with the full factors at most \(X\). This fixes coefficients \(A,B,C\) including their signs. Fix one witnessing \((G_S,G_T,G_F)\). Raw enumeration of its radical and then exact exponents and colors costs \(n^bX^{o_X(1)}\). A smaller available conditional dictionary can replace it. Every original complete letter in the indicated record is left unchanged.

Use the rational point
\[
 U=-ABy/z^2,\qquad V=A^2 Bx/z^3
 \quad\text{on }E_\Delta:V^2=U^3+\Delta,
 \qquad t=U^3/\Delta=-By^3/(Cz^6).
\]
The other signed ratio is \(t+1=Ax^2/(Cz^6)\). The primitive triple is recovered from \(t\) and the fixed assignment of roles and signs, and \(h(t)\le\log X\). Scale by a choice of sixth root of \(\Delta\) to \(E_1:v^2=u^3+1\), so \(t=u^3\) and \(\widehat h=h(t)/6+O(1)\), uniformly. This uses the fixed-model canonical height comparison and functoriality under the isomorphism. The poles of \(t\) and its vanishing orders at the branch points are precisely the associated torsion discs: on \(S,T,F\) the parameters are \(v,u,u/v\) respectively, with valuations \(i/2,i/3,i/6\) when \(v_p(J)=i\) and \(p\nmid6\). The calculation of the orders uses the exact rational \(t,t+1\) and primitivity. It applies equally at defective and nondefective primes. Over the chosen algebraic closure of each local field we use the extension of \(v_p\) normalized by \(v_p(p)=1\). There are at most three choices of center at each prime; partition the point set by them, at cost \(3^{\omega(R\prod\operatorname{rad}G_J)}=X^{o_X(1)}\). At each locked prime \(p\nmid6\), translation by the selected torsion center is integral because \(E_1\) has good reduction for \(p\nmid6\). The difference of two points in the same local part lies in the infinity disc with parameter valuation at least \(i/e_J\), by the integral formal group law. These statements apply if only a lower bound \(i\) is fixed, though we have fixed the full valuation powers. The pole order six of \(t\) there gives \(v_p(t(P-Q))\le-6i/e_J\) off the torsion diagonal. Here \(P-Q\) is computed on \(E_\Delta(\mathbb{Q})\); the rational function \(U^3/\Delta\) gives a rational number, so all these denominator contributions on distinct rational primes add. The finitely many rational torsion points can be factored out at uniform height-subpower cost, or simply use Northcott after scaling to \(E_1\) in degree at most six. At \(2\) and \(3\) use potential good reduction of the fixed model.
At \(2\) set \(u=4^{1/3}x,\ v=2y+1\), giving the good model
\(y^2+y=x^3\). At \(3\), take \(\pi^4=3\) and set
\(u=\pi^2x-1,\ v=\pi^3y\), giving
\(y^2=x^3-\pi^2x^2+x\), smooth modulo the maximal ideal.
The roots of \(u^3+1\) and the points \(v=\pm1\) are a fixed finite
set. If the relevant original valuation \(i\) is sufficiently large,
the identities \(u^3=t,\ v^2=t+1\) place the point near one of the
same torsion centers, or near infinity in the sixth role. For
example the factors of \(u^3+1\) have pairwise differences of fixed
valuation, so a highly divisible product forces one root to be
strictly closer than the others. At each fixed center the original
parameter \(v,u,u/v\), respectively, has a simple zero. Its ratio
with a translated good-model formal parameter is a local analytic
unit with constant valuation in a sufficiently small disc. Thus the
formal depth is \(i/e_J+O(1)\), uniformly in the twist and the points.
Take differences in the integral formal group and compare the
parameter at infinity in the same way. This adds
\(6i/e_J-O(1)\) to the denominator order of \(t(P-Q)\). Bounded \(i\)
can be ignored at these two primes at a total \(O(1)\) log loss.
The exact powers of \(2,3\) can be fixed at polylogarithmic cost and
used as full locks at zero normalized price. Thus distinct
Mordell–Weil vectors in a fixed part have
\[
 \widehat h(P-Q)\ge w_0\log X-O(1).
\]

Apply the [sphere cap proof](ABC_MORDELL_LOCAL_PACKING.md) with radius squared \(\frac16\log X+O(1)\) and this separation. The rational Mordell–Weil rank satisfies the same [CM bound](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#the-3-class-group-in-the-rank-bound)
\[
 3^r\ll_C C^{1+\omega(R)} h_3(d)^2.
\]
It depends only on the twist \(\Delta\), and the good-prime divisors used for extra separation do not enter that rank cost. The elementary cap bound, with any fixed extra angular allowance chosen after \(\varepsilon\), gives \(O_{\varepsilon,k,w_0}(X^\varepsilon h_3(d)^k)\) points per fixed part. This uses \(3^{-k}\le\min(1,3w_0)\), including equality as in the source proof. All \(X^{o_X(1)}\) charges can be incorporated within a separate inner allowance at each fixed \(L\).

For each signed \(d\), the defect radicals \(R\) are multiples of \(|d|\); their number is at most \(n^\rho/|d|\). After summing the auxiliary lock costs uniformly within each coefficient choice, use the signed moment theorem \(\sum_{|d|\le Z}h_3(d)^k\ll_{k,\varepsilon} Z^{f(k)+\varepsilon}\) and dyadic summation. This gives \(n^{\rho+b+(f(k)-1)\sigma}X^\varepsilon\). If bounding minimally repaired quotients separately in height \(Y\), maintain any actual lower balanced target radical \(R_c\ge c^{L(1/3-\xi)}\). Then \(H_M(c)\le Y\) implies \(c\le Y^{1/(ML(1/3-\xi))}\), just as in the standard Euclidean transfer. If the primitive shell exponent is \(E\), sum the shells at fixed \(L,M\) to get \(O(Y^{E/(M(1/3-\xi))+\varepsilon})\), choosing the inner height allowance after \(L,M\). Each primitive quotient supplies one least repair.

## MP2. All subcritical mixed defects on a broad high-degree band

Fix \(0<\eta<1/3\). Let \(0<\tau\ll\eta\) be fixed. In a balanced height shell and at any fixed budget \(L\), suppose:
- the exact mixed defective raw total \(r=\log_n R\) is at most \(1/3-\eta\);
- at most \(\tau\) of its raw mass comes from primes \(2,3\) or from scaled degree \(Lv_p(J)<2-\tau\);
- the square role has total raw mass at least \(1/3-\tau\). After deleting \(2,3\) and degrees below \(2-\tau\), at most a further \(\tau\) of its nondefective raw mass is lost. In the retained nondefective set every prime has normalized logarithmic price at most \(\tau\).

Use fixed raw bins \([r_-,r_+]\) of width at most \(\tau\). Then
\[
 w_{\rm def}\ge \frac{2-\tau}{6}\max(0,r_--\tau)
 \quad\text{and the retained square-role nondefective raw mass is at least }
 \frac13-r_+-2\tau.
\]
Put \(a=2-\tau\) and buy square-role separation by a lock of raw at least
\[
 b_*=\max\left(0,\frac{2}{9a}-\frac{\max(0,r_--\tau)}3\right).
\]
This is affordable: for \(\tau=0\), the difference between available raw mass \(1/3-r\) and the required raw \(1/9-r/3\) is \(2(1/3-r)/3\ge2\eta/3\). Removing the bin and degree tolerances retains a positive margin by taking \(\tau\) sufficiently small in terms of \(\eta\). Greedily select whole nondefective prime powers in the retained square-role set until their raw mass first reaches \(b_*\); the overshoot is at most \(\tau\). A lock of raw \(b\le b_*+\tau\) has weighted numerical contribution at least \(ab_*/2\). Hence the total local separation is at least \(1/9\). The theorem above bounds the normalized exponent in the bin by
\[
 r_++b_*+\tau\le \frac19+\frac23r_++O(\tau).
\]
Summing the finitely many bins gives a fixed saving. For instance, after choosing \(\tau\) sufficiently small, the normalized exponent is at most \(1/3-\eta/2\). The proof works uniformly without any upper bound on individual scaled degrees, provided the stated raw portions and the maximum prime-price condition hold. It has no condition on the two complete dictionary prices or the balanced multiplicity. All original radical caps and any original least-repair cutoff remain available.

The high-degree criterion applies at every original dictionary price. The method only requires a mixed signature whose total defect has subcritical raw price, or a correspondingly cheap actual defect catalogue. It cannot extract such a signature from arbitrary exact valuations. For example, the positive interval \([2,4]\) can be approximated by every exact degree congruent to \(1\bmod 6\). For that formal profile no nontrivial role under any one of the three Euclidean signatures is nondefective. This formal choice rules out extracting such a signature from the weak measures alone. The occurrence of genuine balanced triples with these valuation data remains open.



## MP3. Selection by exact degree without a band hypothesis

Let \(w\) be the full defect weight of a chosen orientation. On one
role \(J\), let \(h_J=\log_X J\), \(a_J=\log_n\operatorname{rad}(J)\),
\(r_J=\log_n\operatorname{rad}(P_J)\), and
\(w_J=\log_X(P_J)/e_J\). Suppose each prime on \(J\) outside the
defect and the two free primes has normalized price at most \(\tau\).
Set \(t_k=3^{-k-1}\). If \(w<t_k\le h_J/e_J\), order those
nondefective primes by decreasing \(Lv_p(J)/e_J\), their ratio of
separation to price. First compute the fractional bound on all nondefective primes with
\(w\) the defect weight alone. Buying the largest ratios up to
the threshold requires at most
\[
 b\le (a_J-r_J)\frac{t_k-w}{h_J/e_J-w_J}+\tau
     \le (a_J-r_J)e_Jt_k/h_J+\tau
\]
The average of a top fractional segment is at least the full
nondefective average. Including the free powers separately can
only reduce the price of this feasible purchase. One nonfree
prime accounts for the rounding. The second inequality
uses \(w\ge w_J\) and \(t_k\le h_J/e_J\).
If the defect weight already reaches the threshold, use no purchase.

For example, a balanced chart with a tall square role
\(h_S=1-o(1)\), cap \(a_S\le1/3+o(1)\), and vanishing maximum
nonfree prime price has a strict saving whenever its mixed defect
price is below \(7/27\) by a fixed margin. At \(k=1\) the purchase
cost is at most \(2/27+o(1)\), independent of the valuation profile.
“Tending to zero” here means that the shell and prime-price
tolerances are chosen below fixed fractions of the saving margin
before the inner height limit. More generally, selecting the
top ratios jointly from all three roles offers an exact
fractional selection bound with one such rounding loss.

## MP4. Joint parity moments over the three CM orientations

For each complete role \(J\) let \(O_J=\operatorname{sf}_2(J)\)
and let \(\sigma_J\) cap \(\log_n O_J\). Let \(\tau_J\) cap the
raw mass on primes with positive even valuation not divisible by
three. Fix the three \(O_J\) and their full odd-valuation parts
at subpower cost after fixing the kernels. For each choice of cube
role \(T\), pick one orientation \((S,T,F)\) of the other roles
as square and sixth. In that orientation, the new defective
primes are precisely the even nonmultiples of three on \(T,F\).
Their full powers cost at most \(n^{\tau_T+\tau_F}X^\varepsilon\).
The fixed parity parts on all roles can contribute to separation,
including on the cube role if its valuations are multiples of
three. Let \(b_T\) cap the conditional dictionary price of
additional full powers and suppose the combined fixed powers
provide separation at least \(3^{-k_T-1}\) for a fixed \(k_T\ge0\).
Write \(q_T=\tau_T+\tau_F+b_T\). The signed squarefree kernel
of the twist is \(\epsilon_T O_S O_F\), where \(\epsilon_T\)
is fixed by the orientation. The kernels on each primitive row
are pairwise coprime.

At each fixed parity triple, MP1's pointwise count before averaging
gives \(X^\varepsilon n^{q_T}h_3(\epsilon_T O_S O_F)^{k_T}\) for
each \(T\) separately. The count is therefore at most the geometric
mean of these three bounds. Sum over the parity triples. If
\(F_{ij}\) is supported on the allowed kernel pairs of roles \(i,j\),
ordinary matrix Cauchy–Schwarz gives
\[
  \sum_{x,y,z} F_{12}(x,y)F_{13}(x,z)F_{23}(y,z)
       \le \prod_{i<j}\left(\sum_{x,y}F_{ij}(x,y)^2\right)^{1/2}.
\]
Indeed the left side is the Frobenius inner product of the
\(13\)-matrix and the product of the \(12\)- and \(23\)-matrices.
The Frobenius norm of the product is at most the product of the
two Frobenius norms, by entrywise Cauchy–Schwarz. Here \(F_{SF}\)
is \(h_3(\epsilon_T O_S O_F)^{k_T/3}\) on the allowed pairs.
Every fixed squarefree product \(O_SO_F\) has subpower many
factorizations, and it is at most \(n^{\sigma_S+\sigma_F}\).
The signed class moments now give the exponent
\[
 \frac13\sum_T q_T+
 \frac12\sum_{(S,T,F)}(\sigma_S+\sigma_F)f(2k_T/3).
\]
In particular, if all three cyclic orientations have separation at
least \(3^{-5/2}\), the first moment suffices. For cyclic roles
the exponent simplifies to
\[
             \sum_J\sigma_J+\frac23\sum_J\tau_J
                                +\frac13\sum_T b_T .
\]
A single orientation needs \(1/9\) separation to remove the class
moment cost in this manner. Averaging separate complete joint
records before summing the kernels gives the lower threshold.

For unequal prices, let \(n^{a_i}\) and \(n^{a_{ij}}\) cap the
actual one- and two-role kernel support sizes. Choose geometric
mean weights \(\lambda_T\ge0\), \(\sum_T\lambda_T=1\), and edge
weights \(\vartheta_T\ge0\) with each pair summing to at most one
and \(\vartheta_T>0\) when \(\lambda_T k_T>0\). In normalized
counting measure the pair-function inequality with exponents
\(\vartheta_T\) follows from ordinary Hölder. The feasible polytope
has vertices zero, the three unit vectors, and
\((1/2,1/2,1/2)\). The respective bounds are the constant bound,
one-edge integration, and the matrix inequality above. Interpolate
the vertex bounds on the product of the same nonnegative functions
raised to the chosen weights. With \(s_T=\lambda_T k_T/\vartheta_T\)
and \(p_T\ge\max(1,s_T)\) for each active edge, another Hölder
application on its support gives the exponent
\[
 \sum_i a_i+\sum_T\lambda_Tq_T+
 \sum_{(S,T,F)}\vartheta_T\left[
       (1-s_T/p_T)a_{SF}
       +(s_T/p_T)f(p_T)(\sigma_S+\sigma_F)-a_S-a_F\right].
\]
Terms with \(\vartheta_T=0\) are omitted. One may take projections
containing all kernel records used by the selection. This permits
unequal signed moments and actual joint kernel support prices.

For an exact-profile comparison, on each original role take odd
valuation \(1\bmod 12\) at raw mass \(.075\), even valuation
\(2\bmod 12\) at mass \(s\in[.05,13/240)\), both of scaled degree
in \([.65,.66]\), and valuation \(6\bmod 12\) on the remainder.
Take the last scaled degrees in \([4.3,4.5]\) with means
chosen for total height one. These residues permit diffuse weak
profiles along decreasing budgets. Every cyclic orientation has
defect separation at least \(.65(.075+s/2)\ge .065>3^{-5/2}\),
so the joint exponent is \(\rho=.225+2s<1/3\). The
\((3,3,3)\) defect is at least \(.375\), and the two fourth roles
alone in a \((2,4,4)\) factorization use total defect mass \(2/3\).
For a single orientation, even spending the full remaining
budget \(b\le1/3-\rho\le1/120\), the separation is at most
\(.66(.075+13/480)+2.25/120=.086125<3^{-11/5}\).
Thus \(k>6/5\) is needed and the single-orientation
penalty \((f(k)-1)(.15)\) exceeds \(1/120\). This comparison
concerns the stated raw-price bounds on valuation data. It makes
no occurrence claim for arithmetic families.

## MP5. The limit of marginal class moments

Consider three formal kernel alphabets of size \(N\), with a marked
subset of size \(\lfloor N^{1-a/2}\rfloor\) on each alphabet, for
fixed \(0<a<5/9\). Give each pair a weight \(H_{ij}\) comparable to
\(N^a\) on the product of the two marked subsets and weight one
elsewhere. Its first moment is \(O(N^2)\), and for every fixed
\(p\ge1\) its \(p\)-th moment is \(O(N^{2+a(p-1)})\), within the
three separate bounds \(N^{2f(p)}\). But the sum of
\(\min_{i<j}H_{ij}^k\) over triples is at least a constant times
\(N^{3+a(k-3/2)}\) when \(k>3/2\), just from the three marked
subsets. Independent pair moments alone therefore cannot carry
the class-cost removal past the symmetric threshold of MP4.
These are formal weights. Genuine joint class-group bounds could
preclude their concentration.

The three signed squareclasses in the CM orientations have product
minus a square. Their quadratic fields lie in a multiquadratic
compositum of degree at most eight. The first-moment theorem for
\(2\)-extensions of Lemke Oliver, Wang and Wood requires a
transposition in the degree-action; the regular multiquadratic
action of degree four or eight does not have one. Their separate
discriminant-ordered bound for groups without transpositions uses
a coarser size parameter than the required joint kernel
catalogue. It does not here supply a three-edge moment at the
price \(\sum_J\sigma_J\).
[Average 3-torsion in 2-extensions, Theorems 1.1 and 6.5](https://www.cambridge.org/core/journals/forum-of-mathematics-pi/article/average-size-of-3torsion-in-class-groups-of-2extensions/6AABE169B5BA8192FE036455DE4F806F)

## Remaining valuation extraction

A sufficient unrestricted advance would extract from every
hypothetical saturation sequence a subfamily of exponent
\(1/3-\varepsilon_0\) on which MP1 or MP4 returns exponent
at most \(1/3-\Delta\) with \(\Delta>\varepsilon_0\),
after all chart, shell and tolerance costs. The high-band and
parity conditions concern actual residues and full divisors,
even when the weak measures are identical. Formal profiles with
all valuations \(1\bmod6\) have complete mixed defects and
full parity kernels, so neither extraction follows from weak
measures alone. The numerical and radical constraints on actual
primitive rows remain part of this obligation.

The point repulsion inputs are [Helfgott–Venkatesh](https://arxiv.org/abs/math/0405180)
and the fixed-model [local packing proof](ABC_MORDELL_LOCAL_PACKING.md).
The signed moment inputs are [Heath-Brown–Pierce, Theorem 1.2
and Corollary 1.4](https://arxiv.org/pdf/1409.3177) and
[Chan–Koymans, Theorem 1.1](https://arxiv.org/html/2505.00611v1#S1).
