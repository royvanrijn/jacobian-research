# Spliced sums and numerical carries on shallow abc registers

The shallow Latin reduction forces many signed splices of
genuine rows, whether or not the third shallow edge exists.
Most such splices have a fresh radical far above the original
cap on the computed role. Even a full shallow Latin square
can have tall primitive lifts with all three deep locks and
maximal mixed residue and numerical-bin entropies, while all
its off-diagonal triangles fail the numerical test.
[Actual shallow partial Latin squares (SRL5)](ABC_SHALLOW_RESIDUE_REGISTER_SATURATION.md#a-three-register-partial-latin-square)

## The outside radical of a forced splice

Put \(n=X^L\), taking the height limit at each fixed \(L\).
Consider a complete-label matching of
\(N=n^{d+o(1)}\) original positive primitive rows
\(A+B=C\), with \(X/2<C\le X\) and
\(\operatorname{rad}(T)\le n^{r_T+o(1)}\) on each role.
The allowances in the count and caps can be made arbitrarily
small at sufficiently small fixed budgets. Fix any full
locks \(D_T\mid T\) throughout the matching, allowing \(D_T=1\),
with normalized radical masses \(\mu_T+o(1)\).
Let \(\mathcal P\) be the union of their prime supports and
write \(\operatorname{rad}_{\mathcal P^c}\) for the product
of the prime divisors outside this union. All stated masses
are fixed and bounded before the height limit, and the locks
have polynomial height. Retain the original eligibility
restrictions on the matching, including any least-repair
cutoff.

We can require that every prime occurs in at most one
role across the retained family. Colour the primes
independently with three colours. A row is retained
with probability \(3^{-\omega(ABC)}=X^{-o_X(1)}\), by
primitivity and polynomial height. Some colouring loses
only height-subpower many in the exponent. The matching
condition then ensures distinct ordered pairs of retained
rows give distinct ordered source pairs in any two roles.
[Colouring an actual complete matching (IM1)](ABC_OCCURRENCE_INDUCED_ADDITIVE_MATCHINGS.md#extracting-an-exact-induced-matching-under-all-three-caps)

Mark the computed role \(W\), with opposite source roles
\(U,V\). Give each row its common absolute register on
\(U,V\), with the sign chosen so that the residues agree.
It is taken modulo a fixed divisor \(M_W\mid D_W\).
Suppose the number of occupied register values is at most
\(n^{\sigma+o(1)}\) with a fixed \(\sigma<d\).
When \(W=C\), call \(A_i+B_j\) a splice output for
distinct rows \(i,j\) having the same register.
When \(W=A\) or \(W=B\), use the positive cross
difference of the opposite target and summand.
For each unordered distinct pair in a register class,
at least one of the two cross differences is positive.
In the \(W=A\) orientation their sum is \(A_i+A_j>0\),
and likewise with the roles exchanged.

**SRS1.** There are at least \(n^{2d-\sigma-o(1)}\) such
positive splices. For every fixed \(\delta>0\) with
\[
 t_* =3(2d-\sigma)-r_U-r_V+\mu_U+\mu_V-\delta>0 ,
\]
all but a power-small proportion of the splices satisfy
\[
                  \operatorname{rad}_{\mathcal P^c}
                    (\text{splice output})>n^{t_*}.
                                                               \tag{SRS1}
\]
The height and budget thresholds can depend on the fixed
caps, \(\delta\), and the margin \(d-\sigma\). In the balanced
case \(r_U=r_V=d=1/3\) at a shallow depth \(\sigma=d/2\),
the threshold is \(5/6+\mu_U+\mu_V-\delta\).
This lower threshold is on the output of a cross pair.
The original matched outputs still satisfy their original
cap \(r_W\).

Indeed Cauchy–Schwarz on the register class sizes gives at
least \(N^2/n^{\sigma+o(1)}-N\) ordered off-diagonal
pairs. Keep at least half of them if choosing the positive
difference orientation. All source pairs are coprime by
the colouring, so the new sums, or the corresponding equations
with a cross difference moved to the left, are positive
primitive triples. Their largest value is comparable to \(X\).
This is automatic for a difference splice because the original
target is one of the two operands. For a sum splice we can
first retain a fixed side \(A\) or \(B\) that exceeds \(X/4\)
on at least half of the original rows and use the restricted
matching. Each splice is divisible by \(M_W\). When all three
shallow moduli are present, a splice output also has the two
opposite residues inherited from its two operands. In particular
it belongs to the target pair class that would close the shallow
triangle. That class need not contain an original target edge.

Give the fixed pool \(\mathcal P\) zero price and use the
ordinary logarithmic price elsewhere. A polynomial-height
prime pool has \(\omega(\mathcal P)=o(\log X)\), as follows
by comparing its product with the product of the first primes.
Thus these prices satisfy (EB5), including the zero-price
choice. The two complete operand radicals outside the pool
have caps at most \(r_U-\mu_U+o(1)\) and
\(r_V-\mu_V+o(1)\) in units of \(\log n\).
A splice with outside output radical at most \(n^t\) therefore
has total effective radical budget at most
\(r_U+r_V-\mu_U-\mu_V+t+o(1)\).
The raw budget is also \(O(L)\), since the pool radical mass
is bounded and \(t\) is fixed. Apply (EB4) to the derived
primitive triples on the finitely many required height shells.
The resulting bound is
\[
 O_{L,\varepsilon}\!\left(
       n^{(r_U+r_V-\mu_U-\mu_V+t)/3+\eta} X^\varepsilon
                    \right)
\]
for any fixed \(\eta>0\) and sufficiently small fixed \(L\).
Every derived triple in the indicated orientation determines
its complete operands and hence the ordered original pair.
At \(t=t_*\) choose the chart and counting losses smaller
than the gap \(\delta/3\). This proves (SRS1).
[Effective radical budgets and prescribed pools (EB4–EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

## Near splices in both residue and numerical scales

**SRS2.** In a balanced zero-multiplicity complete matching,
suppose the tall source role used below has a fine maximum
exact-layer cap as in SRL4. Take a shallow lock on the
computed role of depth \(\sigma\), and fix
\(0\le \xi<d-\sigma\). Put \(H=Xn^{-\xi}\).
There are \(n^{2d-\sigma-\xi+o(1)}\) off-diagonal positive
splices, in the exponent sense, whose two source rows
agree on this register and whose splice output is within
\(H\) of the genuine output of one source row. In this
claim select the splices by requiring the corresponding
two tall complete values to differ by at most \(H\).
The power-small exception to the following outside
radical inequality has the same scope as SRS1:
\[
    \operatorname{rad}_{\mathcal P^c}(\text{splice output})
      > n^{\,3(2d-\sigma-\xi)-r_U-r_V+
                      \mu_U+\mu_V-\delta}.
\]
Here the displayed threshold is assumed positive and
\(\delta>0\) is fixed before taking small budgets.
In the balanced case with \(\sigma=d/2\) it is
\(5/6-3\xi+\mu_U+\mu_V-\delta\).

For a sum splice \(A_i+B_j\), select the same tall
summand on at least half the rows by pigeonhole, say
\(A\ge X/4\). Compare the splice with \(C_j\) by testing
\(|A_i-A_j|\le H\). For a difference splice compare its
two target values \(C_i,C_j\), both at least \(X/2\).
The computed difference then lies within \(H\) of
the genuine computed role on the row whose summand
was kept. Use either orientation of each unordered pair
when it makes the computed difference positive, as in
SRS1. The corresponding tall role determines the row
in a complete matching. All its log values range over
a fixed bounded interval after subtracting \(\log X\).
Two values in one log bin of width sufficiently small
times \(n^{-\xi}\) are within \(H\). Cauchy--Schwarz over
at most \(n^{\sigma+\xi+o(1)}\) joint cells gives the
lower bound after subtracting the diagonal. Conversely
a pair within \(H\) has log bins at bounded distance.
JRM2 and Cauchy--Schwarz between neighboring bin occupancies
bound the number of such pairs by
\(n^{2d-\sigma-\xi+o(1)}\), using the full mixed
gain at depth \(\sigma+\xi\) for fine layers. The
colouring and the tall-side pigeonhole lose only
subpower or constant factors. They retain the
height shell and the original caps. Apply the effective
budget count in SRS1, now against this smaller
reference population. This proves the same power
loss with the stated threshold for numerical near
splices, without requiring any loose triangle in the
original shallow array.
[Joint real and congruence collision bound (JRM2)](ABC_JOINT_REAL_RESIDUE_MOMENTS.md#joint-collision-entropy-and-conditional-full-factors)

## A full shallow Latin square with large spurious carries

The following exact arithmetic construction concerns the joint
constraint from shallow residues and short real bins. It omits
the three original radical caps. Fix positive \(w_A,w_B,w_C\)
with \(w_A+w_B+w_C<1\), and fix \(d=1/3\).
At every sufficiently small fixed budget \(L>0\), as
\(X\) grows, it supplies \(n^{d+o_X(1)}\) actual distinct
positive primitive rows of three mutually comparable
numerical heights. All three full locks \(D_T\mid T\) are
fixed on these rows with \(D_T=X^{w_T+O(L)+o_X(1)}\).
The full register on each lock is injective on the rows.
The shallow registers form an entire Latin square, and
every pair among the three shallow registers and the three
individual numerical-bin labels has entropy \(d+o_X(1)\)
in base \(n\). Nevertheless every spurious shallow triangle
has numerical defect \(\gg X n^{-d/2}\).

For the modulus choices, fix an integer \(j_0\ge1\), put
\(\mu=d/(2j_0)\), and choose distinct primes \(p_T\)
from disjoint fixed relative intervals of scale \(n^\mu\).
The prime number theorem provides these primes at late
heights. Set
\[
      M_T=p_T^{j_0},\qquad
      D_T=p_T^{\,\lceil w_T/(L\mu)\rceil},\qquad
                        D_0=D_AD_BD_C .
\]
The powers in \(D_T\) are fixed at each budget.
For sufficiently small fixed \(L\), \(M_T\mid D_T\),
\(D_T/M_T\gg n^{d/2}\), and
\(D_0=X^{\sum_T w_T+O(L)+o_X(1)}\).
Choose a prime \(q\asymp n^{d/2}\) small enough that
\(q\le\min_T\varphi(M_T)\), and take injective maps
from \(\mathbb F_q\) into the three shallow unit groups.
The row indices are \((i,j,k)\) with \(k=i+j\) in
\(\mathbb F_q\). At \(M_A\), for instance, the symbol
\(i\) gives the common unit of the two opposite signed
coordinates with opposite signs. Use the same convention
at \(M_B,M_C\). Each shallow symbol occurs on \(q\) rows.
Choose distinct full residue lifts modulo \(D_T\) on these
rows. This is possible since each unit residue modulo \(M_T\)
has \(D_T/M_T\) unit lifts. At every mark the marked complete
coordinate is prescribed to be zero modulo its full lock.

For numerical bins write \([z]_q\) for the representative
of \(z\in\mathbb F_q\) in \(\{0,\ldots,q-1\}\). Define
\[
          \kappa_A(i,j)=[i+2j]_q,\qquad
          \kappa_B(i,j)=[3i+5j]_q,\qquad
          \kappa_C(i,j)=\kappa_A(i,j)+\kappa_B(i,j).
\]
The last sum is an integer. Take \(H=\lfloor X/(100q)\rfloor\)
and tall offsets \(A_0=\lfloor X/4\rfloor,\ B_0=\lfloor X/3\rfloor\).
Use the interval from \(A_0+H\kappa_A\) to
\(A_0+H\kappa_A+H/10\) for the complete \(A\)-value, and
the corresponding \(B\)-interval with \(B_0,\kappa_B\).
The three prescribed congruences on each operand combine
by the Chinese remainder theorem modulo \(D_0\).
These intervals have \(X^{\kappa}\) available progression
steps for some \(\kappa>0\) at sufficiently small fixed
budgets, since \(H\asymp X n^{-d/2}\) and
\(\sum_T w_T<1\). Choose \(A\) in its class.
No common prime of the two operand classes can divide \(D_0\),
since at each marked prime at least one operand is a unit.
For the part \(A'\) of \(A\) outside \(D_0\), inclusion–exclusion
on the \(B\)-progression leaves its length in progression
steps times \(\varphi(A')/A'\), with error
\(O(2^{\omega(A')})\). The main term is larger than the
error for late heights, using \(\varphi(A')/A'\ge
2^{-\omega(A')}=X^{-o(1)}\). Choose such a \(B\) and take
\(C=A+B\). This realizes all residues, all numerical bins,
and primitivity. The complete values are distinct in each
role since each sees two opposite shallow symbols
and each ordered pair occurs only once in the Latin square.
The chosen offsets give \(X/2<C<X\) and \(A,B,C\asymp X\).

For \(q>7\), the six linear forms
\[
                   i,\ j,\ i+j,\ i+2j,\ 3i+5j,\ 4i+7j
\]
on uniform \((i,j)\) have pairwise nonproportional coefficient
vectors. Every pair is therefore uniform on \(\mathbb F_q^2\).
Reducing the integer \(C\)-token modulo \(q\) gives the last
form. Passing from that form to the integer token introduces
at most two choices. Passing between the tokens and actual
logarithmic bins of width comparable to \(n^{-d/2}\) introduces
only bounded ambiguity in both directions, because the values
are comparable to \(X\). All asserted pair entropy exponents
follow, for either Shannon or collision entropy.

For a spurious shallow triple \((i,j,k)\), its \(A\)-edge
uses the row \((k-j,j,k)\) and its \(B\)-edge uses
\((i,k-i,k)\), with indices in the field. Its \(C\)-edge
comes from the true row \((i,j,i+j)\). The integer token
difference for the spurious splice satisfies
\[
 \kappa_A(k-j,j)+\kappa_B(i,k-i)-\kappa_C(i,j)
                                     \equiv 6(k-i-j)\pmod q .
\]
If \(k\ne i+j\), it is a nonzero integer. The errors within
the chosen intervals can change the numerical defect
\(A_{jk}+B_{ik}-C_{ij}\) by at most \(H/5\) from its token
multiple of \(H\). Thus its absolute value is at least \(4H/5\).
For every fixed \(\theta<1\) this exceeds \(X^\theta\) at
all sufficiently small fixed budgets and large heights.
The \(q^3-q^2\) spurious triangles all pass the shallow
modular test but fail the short real test.
[Unconstrained CRT integer lifts of full locks](ABC_SHALLOW_RESIDUE_REGISTER_SATURATION.md#sharpness-and-the-lifting-requirement)

## Sparse shallow graphs with no spurious triangles

Kong, Shangguan and Wang, Theorem 1.1, give linear tripartite
hypergraphs with \(q^{2-o(1)}\) edges on parts of size at
most \(q\) and no loose triangles. If all three pair
projections of a vertex triple occur and the triple is
absent as a hyperedge, the three edges witnessing its
pairs must be distinct and form a loose triangle.
Thus a shallow graph at the full pair-count exponent
need have no spurious triangle at all.

The same integer lifting argument applies to these partial
Latin squares. Independently assign to each edge a uniform
token \(K\in\{0,\ldots,q-1\}\). A fixed shallow symbol has
degree at most \(q\), so each joint fibre of that symbol
and \(K\) is binomial of mean at most one. The binomial upper
tail and a union bound over all \(O(q^2)\) such fibres
give an assignment in which every fibre has size
\(O(\log q)\). Set both source tokens equal to \(K\)
and the target token to \(2K\), and use the preceding
full-lock and real-bin lifting. This preserves every mixed
shallow and real pair entropy at exponent \(d\).
The combined real labels can nevertheless have entropy
only \(d/2+o(1)\). Again the original radical caps are
omitted. Together the two constructions delimit the
consequences of shallow and real-bin entropy without the
caps. The bound (SRS1) gives a restriction on splicing
genuine capped rows, but supplies no reason for a
spliced output to inherit the third cap.
[Kong, Shangguan and Wang, Theorem 1.1](https://arxiv.org/html/2611.05261v2#S1.SS1)
