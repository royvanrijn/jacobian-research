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

## MP6. A ramification estimate for asymmetric parity kernels

Let \(A\le B\le C\) be dyadic lower endpoints for positive
pairwise-coprime squarefree kernels \(x,y,z\), respectively.
Write \(h_{xy}=h_3(\epsilon_{xy}xy)\) for fixed signs on the
three edges and use the analogous notation for the other two.
For fixed \(0\le k\le2\), put \(\kappa=1597/5000\).
Then
\[
 \sum_{x\asymp A,\ y\asymp B,\ z\asymp C}^{*}
        \min(h_{xy},h_{xz},h_{yz})^k
 \ \ll_{\varepsilon,k}\ (ABC)^{1+\varepsilon}
  \left(1+(A/B)^{1/3}(AB)^{\kappa(k-3/2)_+}\right),
\]
where the star indicates the pairwise coprimality and squarefree
conditions. The estimate uses the ramified primes of the actual
shared kernel \(x\).

For \(x\) in its dyadic range, let \(S_B(x)\) count the first
\(h_{xy}\) moment over allowed \(y\). The index-three subgroups of
a quadratic ideal class group correspond to cubic fields with
the same fundamental discriminant. In particular \(h_3(d)-1\)
is twice this cubic-field count. At every odd prime dividing \(x\),
such a cubic field is partially ramified because its fundamental
discriminant has valuation one there. Dropping other constraints
gives a cubic-field count with discriminant bound \(O(AB)\).
Apply Bhargava–Taniguchi–Thorne, Theorem 1.4 with \(t=U=1\),
partially ramified modulus equal to the part of \(x\) prime to \(6\)
and no
conditions at other primes. The two residue coefficients in
their Table 1 are \(O(x^{-1+\epsilon})\) after multiplication
over this modulus. Together with the trivial class term this gives
\[
 S_B(x)\ll_\epsilon (ABC)^\epsilon B+E_B(x),\qquad
 \sum_{x\asymp A}E_B(x)\ll_\epsilon (ABC)^\epsilon A^{4/3}B^{2/3}
\]
for nonnegative \(E_B\) bounded by a constant times their absolute
errors. Deleting \(2,3\) has multiplicity at most four. The same
holds with \(C\) in place of \(B\), including either sign of
the discriminant.

Call a center \(x\) good if both bounds give first moments at most
a fixed multiple of \((ABC)^\epsilon B\) and
\((ABC)^\epsilon C\), respectively. The bad set has size at most
\(O_\epsilon((ABC)^\epsilon A\mu)\), where
\(\mu=(A/B)^{1/3}\), and summing the two first moments over it
costs at most \((ABC)^{O(\epsilon)}AB\mu,\ AC\mu\), respectively.
If every edge weight is at least \(H\), each good \(x\) has
\(O((ABC)^\epsilon B/H)\) allowed \(y\) and the corresponding
\(C/H\) bound. Their contribution to the number of triples is
at most \(ABC(ABC)^{O(\epsilon)}H^{-2}\).
For bad centers use the three-edge matrix inequality from MP4.
The first two edge support counts are at most
\(AB\mu(ABC)^{O(\epsilon)}/H,\ AC\mu(ABC)^{O(\epsilon)}/H\).
The third is at most \(BC(ABC)^{O(\epsilon)}/H\), by the
ordinary first moment and divisor enumeration.
This gives the bound \(ABC(ABC)^{O(\epsilon)}\mu H^{-3/2}\)
on this part. Dyadic summation stops at \(H\ll(AB)^\kappa\)
by the pointwise quadratic torsion bound on the \(xy\) edge.
This proves the stated estimate.

Consequently in MP4, if all three orientations have conditional
record cost at most \(q_{\max}\) and separation at least
\(3^{-k-1}\), one can count using the minimum of the three
pointwise class bounds, at exponent
\[
 q_{\max}+\sigma_1+\sigma_2+\sigma_3+
  \max\{0,(\sigma_1-\sigma_2)/3+
              \kappa(\sigma_1+\sigma_2)(k-3/2)_+\},
\]
where the \(\sigma_i\) now describe the actual dyadic kernel
sizes in \(n\), in increasing order. A fixed exponent asymmetry
allows \(k>3/2\) with no class penalty. For \(A^3\le B\) this
includes \(k=2\), giving separation \(1/27\). There is also a
two-orientation proof of this last range. Proposition 4.7 of
the same paper gives \(S_B(x)\ll_\epsilon (ABC)^\epsilon B\)
uniformly for \(A^3\le B\), by placing the cubic-field count
in a discriminant interval enlarged by an arbitrarily small
power to meet the strict modulus cutoff. Use the same bound
for \(C\), and sum the product of the two neighbor counts.

For a comparison on exact valuation profiles, give the three
odd masses \(0.02,0.09,0.09\), and give each role even noncube
mass \(0.066\). Put these defects in scaled degrees
\([.45,.46]\), with residues \(1,2\bmod12\) respectively.
Place the remaining mass at valuation \(6\bmod12\) in degrees
\([3.8,5.4]\) with mean making each numerical height one.
Using the center of odd mass \(0.02\), both relevant orientations
have separation greater than \(1/27\), while
\(\sigma_1+\sigma_2+\sigma_3+2(0.066)=0.332\).
This gives a strict bound on the stated sector, including
diffuse weak measures. The three separation levels of the
unaugmented defects are all below \(3^{-5/2}\).
These are formal profile parameters illustrating the criteria;
the bound itself is on genuine primitive triples satisfying them.
[Bhargava–Taniguchi–Thorne, Theorem 1.4, Table 1 and Proposition 4.7](https://arxiv.org/pdf/2107.12819)

## MP7. A discriminant formulation of the balanced moment

The following implications are conditional on weak Malle upper
bounds for explicit permutation actions. They give a field-count
form of the missing balanced joint moment. Work with the signed
squareclasses in MP4, whose product is minus a square. Triangles
with any \(h_{ij}=1\) contribute at most \(ABC\) to the minimum
moment.

For the other triangles choose one cubic field of fundamental
discriminant associated with each edge. The three quadratic
characters are distinct and independent. Indeed each odd-kernel
prime occurs in precisely two edge supports, the product signed
character is that of \(\mathbb Q(i)\), and a trivial character
or an equality between two characters would leave one of the
three fields with squarefree kernel \(1\) or \(-1\) in the
pairwise-coprime setting. This has already been excluded.
The compositum of the three \(S_3\) closures has group \(G=S_3^3\).
Its sign image is full, and its kernel inside \(C_3^3\) projects
onto each coordinate since the commutator does. The signs act
on this kernel by the three distinct coordinate characters.
Their idempotents over \(\mathbb F_3\) show that each coordinate
line is in the kernel.

In \(G\), let \(H\) be the order-four subgroup of even sign
flips in the standard sign complement. It is core-free,
since no nonidentity double flip stays in \(H\) under all
conjugations by \(C_3^3\). The action on \(G/H\) therefore
has degree 54 and is faithful. A double flip acts as an
involution with six fixed cosets: its centralizer has size
\(2\cdot2\cdot6=24\), and its conjugacy class meets \(H\) in
one element. Its permutation index is 24. All other
nonidentity elements have no fixed cosets and index at least
27. At an odd prime dividing \(xyz\), the three cubic
discriminants are fundamental, with the prime appearing on
exactly two edges. Tame inertia in the compositum is a double
flip, so the degree-54 subfield discriminant has exponent
24 there. The contribution at 2 is bounded in fixed degree.
The total discriminant is therefore \(O((xyz)^{24})\); its
part outside 2 equals the 24th power of the odd part of \(xyz\).

The degree-54 subfield determines the compositum as its Galois
closure. There are boundedly many cubic subfield choices in a
fixed closure, and the three edge radicands recover the three
coprime kernels by gcd. Thus the count of ordered cubic
triples has only bounded multiplicity over the degree-54
fields. Conversely the number of choices of cubic triples
per kernel triple is \(\prod_{i<j}(h_{ij}-1)/8\), which is
comparable to \(\prod h_{ij}\) here. It suffices to prove
\(O_\varepsilon(Y^{1/24+\varepsilon})\) for just these
degree-54 fields with the specified sign relation and inertia,
ordered by discriminant \(Y\). Weak Malle for the full
permutation group predicts this bound. It would give
\[
                 \sum{}^{*}\min_{i<j} h_{ij}^k
                       \ll_\varepsilon (ABC)^{1+\varepsilon}
                    \qquad (0\le k\le3).
\]
In particular it would allow separation \(1/81\) in all three
orientations at the common conditional record cost.
[Weak Malle formulation for solvable groups](https://arxiv.org/abs/1804.11318)

For \(k\le2\) the corresponding two-edge problem has smaller
degree. Use two cubic fields with independent quadratic
characters and their Galois compositum of group \(S_3\times S_3\).
In the regular degree-36 action every nontrivial involution
has index 18, the minimum. At odd primes the inertia is an
involution, with single or double flips according to which
quadratic discriminants contain the prime. Its discriminant
is \(O((xyz)^{18})\) when using the \(xy,xz\) edges, and the
vertex factors are recovered by gcd. The weak \(Y^{1/18+\varepsilon}\)
bound in this regular action would therefore give the exponent
free \(k\le2\) moment for balanced kernels as well. Degenerate
cases with a trivial edge weight have already been removed.
Current results for a symmetric group times an abelian group,
or for the product-action wreath families, have different
groups or permutation representations from these two counts.
[Masri–Thorne–Tsai–Wang, symmetric times abelian](https://doi.org/10.1112/S0010437X20007587)
· [Inductive methods for counting number fields](https://doi.org/10.1515/crelle-2025-0074)

## MP8. A balanced sector counted by short ideal representatives

Let \(x,y,z\) run dyadically at a common scale \(N\), with
the pairwise coprime squarefree restrictions of MP6.
Suppose the \(xy\) edge is imaginary, and let \(h_M(-xy)\)
count the nontrivial classes in
\(\mathrm{Cl}(\mathbb Q(\sqrt{-xy}))[3]\) whose least integral
ideal norm is at most \(M\). For \(M\ll N\) one has, uniformly
in \(x\),
\[
               \sum_{y\asymp N}^{*} h_M(-xy)
                              \ll_\varepsilon N^\varepsilon M^{5/2}/N .
\]
Choose a minimum-norm ideal \(\mathfrak a\) in such a class.
It is primitive as an ideal, in the sense that no rational
prime ideal \(p\mathcal O_K\) divides it. If \(m=N\mathfrak a\)
and \(\mathfrak a^3=(\alpha)\), the imaginary quadratic norm
equation gives
\[
                        4m^3=t^2+xy w^2 .
\]
Here \(t,w\) are integers with \(|t|\ll M^{3/2}\) and
\(0<|w|\ll M^{3/2}/N\). The nonzero condition on \(w\) follows
because a rational \(\alpha\) with every ideal valuation
divisible by three is itself a rational cube, even at the
ramified primes where the ideal valuation has the factor two.
Then \(\mathfrak a\) would be principal. For fixed \(m,t,w\)
with \(w\ne0\), the positive \(y\) is determined and the
ideal cube, hence the ideal \(\mathfrak a\), is determined.
This gives an upper-bound count by these integer solutions.

Let \(w_0\) be the odd part of \(|w|\). For an odd prime
\(p\mid(m,w_0)\) the norm equation gives \(p\mid t\).
Thus \((\alpha)\) is divisible by \(p\mathcal O_K\).
If \(p\) is unramified this contradicts primitivity of
\(\mathfrak a\), since every prime over \(p\) then divides
\(\mathfrak a\). At a ramified prime the primitive ideal
has valuation one there, so \(v_p(m)=1\).
The two terms \(t^2,xy w^2\) have valuations of different
parity, and their sum has valuation three. Thus \(v_p(w)=1\).
For integer \(m,w\) satisfying these necessary conditions,
the congruence
\[
                   t^2\equiv4m^3\pmod{x_{\rm odd}\,w_0^2}
\]
has at most \(N^\varepsilon\gcd(m,w_0)\) residues. For unit
right-hand side any odd prime power has at most two roots.
When \(p\mid m\) but \(p\nmid w_0\), the modulus has at most
a single factor \(p\) and the root is zero. At a prime common
to \(m,w_0\), the modulus exponent is at most three and the
square-zero congruence has at most \(p\) roots. Take the product
by CRT.

Put \(V=4M^{3/2}/N\). If \(V<1\) the sum is empty. Otherwise the number
of admissible \(t\)'s is at most
\[
 N^\varepsilon\gcd(m,w_0)
                  \left(1+\frac{M^{3/2}}{N w_0^2}\right).
\]
On summing over \(m\le M\), use
\(\sum_{m\le M}\gcd(m,w_0)\le M\tau(w_0)\), by the divisor
formula for the gcd. The sum
\(\sum_{0<|w|\le V}\tau(w_0)/w_0^2\ll1+\log V\),
since each odd part repeats for at most \(O(1+\log V)\)
powers of two and the odd weighted divisor sum converges.
The unweighted sum is \(O_\varepsilon(VN^\varepsilon)\).
This proves the row bound with arbitrarily small losses.
[Norm equation method for quadratic class groups, Pierce](https://sites.math.duke.edu/~pierce/papers/Pierce_class_numbers.pdf)

On these balanced triples the ordinary first moment on the
\(xz\) edge is \(O(N^{2+\varepsilon})\) after summing over
\(x,z\), by factorization multiplicity and the signed first
moment. Hence the joint sum \(h_M(-xy)h_3(\epsilon xz)\)
is \(O_\varepsilon(N^{1+\varepsilon}M^{5/2})\).
For kernel triples with \(h_3(-xy)>1\) on which at least
an \(N^{-\delta}\) fraction of nontrivial classes have
such representatives, the actual two-edge product moment
is bounded by \(O_\varepsilon(N^{1+\delta+\varepsilon}M^{5/2})\).
The corresponding minimum moment for \(k\le2\) is obtained
by adding the trivial \(O(N^3)\) contribution where any edge
has weight one. Thus \(M\le N^{4/5-\eta}\) and
\(\delta<(5/2)\eta\) suffice for the price-free \(k=2\)
mixed packing moment on this sector, even with comparable
parity kernels. The same bound counts the MP7 degree-36
fields when their imaginary resolvent meets this fraction
condition. It refers to the fraction of ideal classes within
a resolvent and requires no identification of individual
classes with specific cubic fields.

Using the imaginary quadratic Minkowski bound \(M\ll N\)
for all classes gives an unconditional \(O(N^{7/2+\varepsilon})\)
bound for the corresponding balanced degree-36 family
with an imaginary edge. The Cauchy–Schwarz bound from the
separate quadratic second moments has exponent \(1+2f(2)=32/9\).
The \(7/2\) bound remains larger than the
\(N^{3+\varepsilon}\) field count needed in MP7.
The count of classes with long minimum ideals, and extraction
of the short-ideal sector from actual balanced abc families,
remain open.

## MP9. A degree-six relative discriminant for the two-edge family

Take two fundamental \(S_3\) cubic discriminants with squarefree
radicands \(\epsilon_1 xy,\epsilon_2 xz\) in the degree-36
family of MP7. Write \(G=S_3\times S_3\) for the group of
the compositum. Its action on the cosets of the diagonal
\(S_3\) is the degree-six action \(b\mapsto g_1 b g_2^{-1}\)
on \(b\in S_3\). This action is faithful, since an element in
its kernel has equal components which commute with all of
\(S_3\). The associated field \(M_6\) contains the quadratic
field \(K_0=\mathbb Q(\sqrt{\epsilon_1\epsilon_2 yz})\),
fixed by the kernel of the product sign. The extension
\(M_6/K_0\) is cubic with relative quadratic resolvent the
biquadratic compositum of the original two resolvents.
Indeed on a sign block the kernel of the product sign has
cubic action with one order-three diagonal line in its kernel,
and its quotient is \(S_3\).

At an odd prime of \(x\), inertia is a double transposition
in the product group. The representative \((t,t)\), \(t\)
a transposition, fixes \(1,t\) in the degree-six action and
one in each sign block. Its permutation index is two
overall, and one on each relative cubic block. At an odd
prime of \(yz\), a single transposition acts without fixed
points, with index three. It interchanges the two sign
blocks, so the extension over \(K_0\) is unramified there.
In particular
\[
  |\operatorname{Disc}(M_6)|^{(2)}=x_{\rm odd}^{\,2}
                              (y_{\rm odd}z_{\rm odd})^3,
 \qquad
 \mathfrak d(M_6/K_0)^{(2)}
                             =x_{\rm odd}\mathcal O_{K_0}.
\]
The superscript deletes powers or ideals over 2. At primes
of \(x\) this last assertion holds both for a split prime
and for an inert one, since the local action on each cubic
block has index one. This relative discriminant is a
rational principal ideal with the specified squarefree
generator. The degree-six field determines its Galois
closure, so its original pair of cubic subfields has only
bounded multiplicity.

In balanced boxes \(|\operatorname{Disc}(K_0)|\asymp N^2\)
and the relative discriminant has norm \(\asymp N^2\), but
its odd rational generator lies in a box of size \(N\).
Taking opposite signs gives an imaginary \(K_0\), available
among actual signed CM orientations. The reformulation is a
count over moving quadratic bases with this restricted ideal
condition. General relative cubic discriminant counts do not
enforce the rational generator restriction.
[Relative cubic counts and 2-extensions](https://arxiv.org/abs/2110.07712)

## MP10. Transverse norm catalogues at Minkowski height

Use MP8's imaginary fields with \(x,y\asymp N\), and choose
one minimum-norm ideal and one cube generator for each class.
At fixed \(x\) restrict the actual nonzero magnitudes \(|w|\)
in their representations \(4m^3=t^2+xyw^2\) to a catalogue
of at most \(K\) values. If \(m\le M\ll N\), then
\[
       \sum_{y\asymp N}^{*} h_{M,\mathrm{cat}}(-xy)
                              \ll_\varepsilon N^\varepsilon M K,
\]
where the left side counts only represented nontrivial classes.
Unlike MP8's unrestricted \(w\) count, this estimate can use
even a single large transverse height. The entire modulus
\(q=xw^2\), including its prime at 2, is needed.

For fixed \(m,w\), the condition \(y\asymp N\) puts \(t\) in
two intervals of length \(O(q)\). To see this take \(N\le y<2N\).
If \(4m^3\ge4Nq\), rationalize the difference of the square
roots to bound the length by \(O(Nq/m^{3/2})\). Otherwise bound
the full interval by \(O(m^{3/2})=O(Nq/m^{3/2})\).
Nonemptiness gives \(4m^3\ge Nq\ge N^2\), so \(N/m^{3/2}=O(1)\).

The odd prime-power root count for \(t^2\equiv4m^3\bmod q\)
is bounded by \(N^\varepsilon\gcd(m,w_{\rm odd})\), as in MP8.
At 2 there are uniformly boundedly many roots for the
genuine primitive-ideal data. If \(w\) is odd, the modulus
has depth at most one. If \(m\) is odd, write \(t=2t_0\) at
depth at least three and use the unit square-root bound
modulo the modulus divided by four, with the bounded lifts.
If \(m,w\) are both even and 2 is unramified, then
\(\alpha=(t+w\sqrt{-xy})/2\) is a rational scalar modulo
\(2\mathcal O_K\), using the integral basis
\(1,(1+\sqrt{-xy})/2\). Its even norm makes this scalar zero.
This contradicts primitiveness of its ideal cube at the
unramified prime, as in MP8. If 2 is ramified, primitiveness
gives \(v_2(m)=1\). A valuation \(v_2(w)\ge3\) would force
\(v_2(t^2)=5\) in the norm equation, which is impossible.
Then the modulus has depth at most five. This handles all cases.

The interval and residue bounds give
\(O_\varepsilon(N^\varepsilon\gcd(m,w_{\rm odd}))\) choices
of \(t\) per \(m,w\). Sum over \(m\le M\) using the same
gcd divisor bound as in MP8, then over the \(K\) values.
Each triple determines at most one class as before.
In particular the joint weighted count with another edge
\(h_3(\epsilon xz)\) is \(O_\varepsilon(N^{2+\varepsilon} M K)\).
For a resolvent sector where at least an \(N^{-\delta}\) fraction
of the nontrivial classes admit these actual records,
\(M K\le N^{1-\delta+o(1)}\) suffices for the price-free
\(k\le2\) minimum moment. This includes the full Minkowski
range \(M\ll N\) when the transverse catalogue and reciprocal
fraction cost are both subpower. There is no bound here forcing
such a catalogue for all balanced kernels.

## MP11. Uniform balanced rows via cubic coefficient lattices

The [fixed-divisor cubic count](ABC_FIXED_DIVISOR_CUBIC_DISCRIMINANT_COUNTS.md)
improves the whole balanced row to \(N^{4/3+\varepsilon}\)
and first gives the two-edge moment \(N^{10/3+\varepsilon}\),
uniformly including prime shared kernels. It applies the repeated-root
condition at each prime of \(x\) directly to successive minima
of the cubic coefficient lattice. In the rank-three case
the integral dual cubic has discriminant divisible by the whole
modulus, and a sufficiently small wedge forces a rational
double root. At the present pointwise class bound this
does not yet lower the symmetric threshold in MP4.

The averaged short-certificate argument improves this further
to \(N^{29/9+\varepsilon}\) on the whole balanced two-edge moment.
Nonsingular short certificates restrict the shared squarefree
kernel by their discriminants. Singular certificates restrict the
root assignments directly. Combining the two-edge estimate with
MP4 gives \(N^{3+(4/9)(k-3/2)+\varepsilon}\) for the balanced
three-edge minimum moment on \(3/2\le k\le2\). The \(k=2\) target
\(N^{3+\varepsilon}\) remains open.
[Averaged cubic certificates](ABC_FIXED_DIVISOR_CUBIC_DISCRIMINANT_COUNTS.md#an-averaged-saving-from-short-cubic-certificates)

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
