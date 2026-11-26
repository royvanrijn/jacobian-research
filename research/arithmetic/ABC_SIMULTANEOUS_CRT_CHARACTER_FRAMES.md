# Character frames for simultaneous prime congruences and real windows

Products of actual complete values can mix jointly modulo a fixed
number of auxiliary primes, even when their root coordinates are
correlated. A positive Pell correction removes all pure quadratic
modes jointly while moving a real logarithmic window by a controlled
amount. Taking a ratio of two independent products guarantees mass
in a short window about one.

The proof uses the weighted tuple character frame, the primitive
large sieve for composite conductors, and a separation of the
nonquadratic and quadratic conductor factors. The mixing supplies
auxiliary units. Its exponent does not pay for independent
replication of a genuine occurrence.

[One-prime character frames and windows (CFM1–CFM6)](ABC_COMPLETE_VALUE_CHARACTER_FRAMES.md)
· [Exact bilinear split with an auxiliary scale (BS1)](ABC_AUXILIARY_BILINEAR_SPLIT.md)

## Dictionary and prime progressions

Fix \(L>0\), put \(n=X^L\), and take the height limit at this
fixed budget. Let
\[
 \Omega=\prod_{i=1}^J[1,B_i]_{\mathbb Z},\qquad
 B_i=\lfloor n^{b_i}\rfloor,\quad b_i>0,\qquad
 T=\sum_i b_i,\quad b_*=\max_i b_i .
\]
Let the positive dictionary \(\mathcal D\) have at least \(n^d\)
distinct values. Choose one representation
\(W=F\prod_i x_i^{e_i}\) in \(\Omega\) for each value, with
the same positive integer \(F\) of polynomial height and
positive integral orders \(e_i\) fixed before the height limit.
Use the uniform probability law on these representations.
Write \(\Delta=T-d\ge0\). The positive exponents, the finite
lists and all moment orders below are fixed during that limit.
The constant in each \(o_X(1)\) allowance may depend on this
data, including \(L\).

Fix an integer \(m\ge1\) and \(\beta>b_*\).
Choose distinct marker primes \(\ell_1,\ldots,\ell_m\equiv1\bmod4\)
outside the primes dividing the orders. Take the \(j\)-th
auxiliary prime \(q_j\) from a fixed reduced progression in
\([n^\beta,2n^\beta]\) with
\[
 q_j\equiv3\pmod8,\qquad q_j\equiv2\pmod p
       \ \ (p\text{ odd},\ p\mid\textstyle\prod_i e_i),
 \qquad
       \left(\frac{\ell_k}{q_j}\right)=(-1)^{{\bf1}_{j=k}}
       \quad(1\le k\le m).
\]
Quadratic reciprocity makes the last requirements conditions
on reduced classes modulo the marker primes, so the
progressions exist by CRT and they are disjoint. The prime
number theorem in each fixed progression supplies
\(n^{\beta-o_X(1)}\) primes. Remove primes dividing \(F\).
There are \(O_{L,\beta}(1)\) such primes at polynomial
height. The remaining prime choices are independent across
the \(m\) labelled progressions. Put \(Q_{\boldsymbol q}=\prod_jq_j\)
and \(G_{\boldsymbol q}=(\mathbb Z/Q_{\boldsymbol q}\mathbb Z)^\times\).
The roots and the values are units in every component.
For all \(i,j\) one has \(\gcd(e_i,q_j-1)\mid2\).

## A phase moment separated by conductor type

Index characters \(\chi=\prod_j\chi_j\) of \(G_{\boldsymbol q}\)
by their local characters. Call a component *quadratic* when
it is the nonprincipal quadratic character. Let \(S\) be the
indices of the nonprincipal, nonquadratic components and
\(R\) those of the quadratic components. Write
\[
 I_{\rm odd}=\{i:e_i\ {\rm odd}\},\qquad
 s_{j,I}=\sum_{i\in I}\frac{1}{2\lceil 2j\beta/b_i\rceil},
 \qquad x_+=\max(x,0).
\]
For an integer \(h\ge1\) set
\[
 g_{j,I,\nu}(h)
     = j\beta-h\Delta-\nu j\beta(1-hs_{j,I})_+ .
\]
The nonprincipal characters under control are chosen in
one of two ways. Always include those with \(S\ne\varnothing\),
using \(I=\{1,\ldots,J\}\), \(j=|S|\), and \(\nu=2\).
One may also include the pure quadratic modes with
\(S=\varnothing\), \(R\ne\varnothing\), using \(I=I_{\rm odd}\),
\(j=|R|\), and \(\nu=1\). Let \(\delta_h\) be the minimum of
\(g_{j,I,\nu}(h)\) over the selected types and \(1\le j\le m\).
Only use the conclusion below when \(\delta_h>0\).

**CCF1.** For \(\tau\in\mathbb R\) put
\(\Phi_{\boldsymbol q,\chi}(\tau)
  =\mathbb E[\chi(W)e^{2\pi i\tau\log W}]\).
Uniformly in \(\tau\),
\[
  \mathbf E_{\boldsymbol q}
       \sum_{\chi\ {\rm under\ control}}
                 |\Phi_{\boldsymbol q,\chi}(\tau)|^{2h}
                 \ll n^{-\delta_h+o_X(1)}.
 \tag{CCF1}
\]
This includes the incomplete-moment arm with
\(hs_{j,I}<1\). A simpler sufficient depth condition is the
following. If every selected \(I\) has mass at least \(a_0>0\),
take
\[
                    h a_0\ge 2(2m\beta+b_*).
\]
All positive-part penalties then vanish and
\(\delta_h\ge\beta-h\Delta\). Thus one needs a sufficiently
small box cardinality defect. The finite formula for
\(\delta_h\) allows smaller depths when it stays positive.

### Primitive characters and the active slot frame

First let \(S\ne\varnothing\). Freeze all \(q_j\) at indices
outside \(S\). The factors at \(R\) then give a fixed bounded
tuple twist. Absorb them and the real phase into the probability
weights on \(\Omega\), which preserves their squared norm.
The varying character
\(\widetilde\chi=\prod_{j\in S}\chi_j\) is primitive modulo
\(D=\prod_{j\in S}q_j\), and \(\widetilde\chi^{e_i}\) is still
primitive for every slot \(i\). Indeed each local component has
order greater than two, and no odd divisor of any \(e_i\)
divides the corresponding \(q_j-1\). Its power therefore stays
nonprincipal. Distinct choices of these labelled primes
give distinct conductors by the disjoint progressions.
For a pure quadratic type instead let \(D=\prod_{j\in R}q_j\).
Now \(\widetilde\chi\) and its powers on odd slots have
conductor \(D\), and its even powers are principal on the units.

Here is the frame estimate for either case. Write \(j\) for the
number of varying primes and \(I\) for the corresponding set of
slots whose powers stay primitive. Set \(k_i=\lceil2j\beta/b_i\rceil\).
On \(\Omega\) with probability measure use the vectors
\(\xi_f(\boldsymbol x)=\prod_i\widetilde\chi_f(x_i)^{e_i}\)
indexed by the varying conductor–character pairs \(f\).
For a fixed anchor \(f\), the slot Gram factor
\[
 r_i(f,g)= B_i^{-1}\sum_{x\le B_i}
            \widetilde\chi_f(x)^{e_i}
            \overline{\widetilde\chi_g(x)^{e_i}}
\]
satisfies
\[
                    \sum_g |r_i(f,g)|^{2k_i}
                         \ll n^{o_X(1)}\qquad(i\in I).
\]
Indeed, taking \(k_i\)-th powers gives coefficients supported
on integers at most \(B_i^{k_i}\), with squared coefficient
sum at most \(B_i^{k_i}n^{o_X(1)}\) by the fixed-order divisor
bound. The anchor gives a bounded coefficient twist. The
primitive character large sieve applies to the varying
\(\widetilde\chi_g^{e_i}\) on conductors at most \(2^j n^{j\beta}\).
The power map has multiplicity at most \(|e_i|^j\) by CRT.
The resulting bound is
\(n^{o_X(1)}(1+n^{2j\beta}/B_i^{k_i})\).
This applies also when the varying conductor equals the
anchor conductor, since the large sieve tests the varying
powered character by itself.

For clarity, the primitive-character large sieve used here
follows from the primitive Gauss-sum identity. At modulus \(D\),
express a primitive character sum against coefficients by
additive sums at the fractions \(a/D\), \((a,D)=1\), divided
by a Gauss sum of absolute value \(\sqrt D\). Summing squared
values over the primitive characters is at most
\(\varphi(D)/D\) times the sum of squared additive values on
those fractions, by full character orthogonality. Multiplying
by \(D/\varphi(D)\), summing over \(D\) up to the conductor
cap, and applying the additive large sieve proves the
stated composite-bound input.
[Montgomery–Vaughan, Theorem 1](https://personal.science.psu.edu/rcv4/personal/Publications/large_sieve.pdf)

Pass to the \(h\)-fold tensor frame and discard the inactive
slot factors in an absolute Gram row sum. If \(hs_{j,I}\ge1\),
choose Hölder weights of sum one with \(w_i\le h/(2k_i)\).
The displayed slot moments bound that row sum by \(n^{o_X(1)}\).
If \(hs_{j,I}<1\), choose \(w_i=h/(2k_i)\) and put the
remaining weight \(1-hs_{j,I}\) on the constant function
one. There are at most \(n^{\nu j\beta+o_X(1)}\) indices \(g\):
the quadratic case has one character per conductor, and
the nonquadratic case can use up to the conductor size.
In both cases the row sum, hence the Gram operator norm,
is at most \(n^{\nu j\beta(1-hs_{j,I})_++o_X(1)}\).
This proof includes the case \(I\) empty, where the bound
is just the number of indices.

The test vector is the \(h\)-fold tensor of \(|\Omega|\) times
the tuple weight, using the conjugate frame. Its squared
probability-space norm is at most
\((|\Omega|/|\mathcal D|)^h\le n^{h\Delta}\).
Bessel's bound thus bounds the summed \(2h\)-th moment by
\(n^{h\Delta+\nu j\beta(1-hs_{j,I})_++o_X(1)}\).
The factor \(\widetilde\chi(F)\) is a unit phase of modulus
one. Divide by the number \(n^{j\beta-o_X(1)}\) of choices
for the varied primes, average over the frozen choices,
and sum the fixed number of component patterns. This proves
(CCF1).

## Joint residue density in a real window

Let \(W_1,\ldots,W_{2h}\) be independent draws of the same
actual dictionary law, and choose signs
\(\sigma_i\in\{1,-1\}\) in advance. Put
\[
       Z=\prod_{i=1}^{2h}W_i^{\sigma_i},
       \qquad K_H(x)=(1-H|x|)_+\quad(H>0).
\]
The relevant transform of \(Z\) is
\(\mathbb E[\chi(Z)e^{2\pi i\tau\log Z}]\), the product of
one factor \(\Phi_{\boldsymbol q,\chi}(\tau)\) for each plus
sign and its complex conjugate for each minus sign.
Its absolute value is always \(|\Phi_{\boldsymbol q,\chi}(\tau)|^{2h}\).

**CCF2.** Fix \(0<\gamma<\delta_h\). Outside a proportion
\(n^{-\delta_h+\gamma+o_X(1)}\) of the prime tuples, one has
simultaneously for every \(u\in G_{\boldsymbol q}\) and every
real centre \(c\)
\[
\begin{aligned}
 &|G_{\boldsymbol q}|\,
       \mathbb E[{\bf1}_{Z\equiv u\bmod Q_{\boldsymbol q}}
                            K_H(\log Z-c)]\\
 &\quad=\mathbb E[K_H(\log Z-c)]
   +\sum_{\substack{\chi\ne1\\ \chi\text{ not under control}}}
        \overline{\chi(u)}\,
                 \mathbb E[\chi(Z)K_H(\log Z-c)]
   +O(n^{-\gamma}).
\end{aligned}
 \tag{CCF2}
\]
The assertion concerns any prescribed scale \(H>0\), which
may depend on \(X\), or each of a prescribed
\(n^{o_X(1)}\)-size family of scales on a common good set.
It holds for the finitely many sign lists simultaneously.
If the low modes are absent or removed, the window-weighted
residue law has relative error
\(O(n^{-\gamma}/\mathbb E K_H)\) for positive window mass.

To prove this, the Fourier transform of the triangular
kernel is the positive measure
\[
       d\nu_H(\tau)=H^{-1}
         \left(\frac{\sin(\pi\tau/H)}{\pi\tau/H}\right)^2d\tau
\]
of total mass one. Integrate (CCF1) against it and apply
Markov at the prime-tuple variable. At a good tuple
the controlled sum of the absolute Fourier integrals
is at most \(n^{-\gamma}\). The centre contributes the
factor \(e^{-2\pi i\tau c}\) of absolute value one in each
integral, so the good set is independent of the centres.
The triangle inequality and character inversion give
(CCF2). A subpower-sized union over the scales preserves
the exceptional exponent. The same controlled-mode bound
holds if an independent positive random factor is multiplied
onto \(Z\) and \(\log Z\) is replaced by the logarithm of
the corrected product. Its joint character and log transform
has modulus at most one and acts as a multiplier inside the
integrals.

If only the \(S\ne\varnothing\) types are controlled,
all remaining modes are pure quadratic. An independent
product of choices \(1\) or \(\ell_j\), with each choice
uniform, kills these modes exactly if the window stays on
the uncorrected logarithm. A nontrivial pure quadratic mode
has at least one index \(j\) whose marker has character
value \(-1\). For the logarithm of the corrected value
itself, this correction generally changes the weights of its
two choices. The next construction controls that change.

## Pell correction on the corrected logarithm

**CCF3.** Suppose only the nonquadratic-containing types are
used to define \(\delta_h>0\). Fix \(C\ge0\),
\(0<\gamma<\delta_h\), and \(\theta>(C+\gamma)/2\).
One can multiply \(Z\) by an independent positive rational
correction \(P\) so that, up to the same exceptional proportion
as in (CCF2),
\[
 \begin{aligned}
 &|G_{\boldsymbol q}|\,
      \mathbb E[{\bf1}_{ZP\equiv u\bmod Q_{\boldsymbol q}}
                             K_H(\log(ZP)-c)]\\
 &\qquad=\mathbb E[K_H(\log(ZP)-c)]+O(n^{-\gamma}).
 \end{aligned}
 \tag{CCF3}
\]
This is simultaneous over the same prescribed scales
\(H\le n^C\), all \(u,c\), and all sign lists. The correction
has a deterministic set of at most \(2^m\) options with
\(P=1+O(n^{-2\theta})\), of polynomial height at fixed budget.

For each marker take a positive Pell solution
\(a_j^2-\ell_j b_j^2=1\) with \(a_j\asymp_{\ell_j} n^\theta\).
Powers of a fundamental positive Pell unit grow geometrically,
so the first solution with \(a_j\ge n^\theta\) has this bound.
Set \(A_j=a_j^2\), \(B_j=\ell_j b_j^2=A_j-1\). Use
independent uniform choices \(P_j\in\{1,B_j/A_j\}\), and
let \(P=\prod_jP_j\). Remove primes dividing any \(A_jB_j\)
from the test progressions. This removes only \(O_{L,\beta,\theta,m}(1)\)
primes and changes none of the estimates. For every selected
prime, \(B_j/A_j\) has the quadratic signature of \(\ell_j\).
Moreover
\[
       |\log(B_j/A_j)|\ll n^{-2\theta}.
\]
Given any nontrivial pure quadratic mode, choose an active
marker \(j\) for it and condition on all the other corrections
and the original \(Z\). The two terms in the expected
quadratic window transform now have opposite signs.
Their absolute average is at most
\((H/2)|\log(B_j/A_j)|\) by the Lipschitz bound for \(K_H\).
This is \(O(n^{C-2\theta})\). Sum over the at most \(2^m-1\)
such modes. The other modes retain their bound from (CCF2)
by the independent-transform multiplier observation.
This proves (CCF3).

One can instead use integer choices \(A_j,B_j\).
They introduce the deterministic log translation
\(\sum_j\log A_j\) and a deterministic unit dilation at each
good modulus, both harmless for the simultaneous estimate.
Each integer choice has size \(n^{2\theta+o_X(1)}\) and radical
at most \(n^{\theta+o_X(1)}\), with constants allowed to depend
on the marker. The rational choice can have this radical cap
on each of its numerator and denominator, costing up to
\(n^{2\theta+o_X(1)}\) on their product. These are auxiliary
costs and are not absorbed into the original three radical
budgets by the theorem.
[Pell and norm recovery in Erdős #943](ERDOS_943_PELL_RECOVERED_KERNEL_GCD_WEIGHT.md)

## A central window with guaranteed mass

**CCF4.** Take \(h\) plus signs and \(h\) minus signs in
(CCF2) and (CCF3). With or without the small rational Pell
correction, the central window satisfies
\[
                  w_H=\mathbb E K_H(\log(ZP))
                           \gg (1+H\log X)^{-1}.
\]
Here use \(P=1\) if there is no correction. Indeed, the log of an
\(h\)-fold positive product lies in an interval of length
\(O(\log X)\). Divide it into intervals of length \(1/(4H)\).
Two independent products fall in the same interval with
probability at least the reciprocal of the number of intervals.
Their log ratio then has absolute value at most \(1/(4H)\).
In the corrected case \(|\log P|\le1/(4H)\) at sufficiently
large heights by the stipulated choice of \(\theta\). The
triangular weight on this event is at least \(1/2\).

For \(H=n^\kappa\) with \(0\le\kappa<\gamma\), when all
nonprincipal modes are controlled or corrected, this gives
a residue law of relative error
\(O(n^{-\gamma+\kappa+o_X(1)})\) under central window weighting.
The weighted law is supported on positive rational scales
within \(O(n^{-\kappa})\) of one when \(\kappa>0\). It can
supply real control on an auxiliary ratio without needing a
separate hypothesis on the occupancy of a particular window.
Clearing denominators permits use in integer auxiliary
identities where only ratios of positive parameters matter.
Any separate numerator, denominator or support caps needed
for that application still have to be charged.

## Short real scales in recoverable signed histories

**CCF5.** Use the paired singleton and pair blocks in (SS6).
Take a fixed group of good primes as in (CCF4), with local
\(\beta>h_{\rm rad}\), where every original marked radical
\(r_i\) is at most \(n^{h_{\rm rad}}\). Require every
original value in the transformations to be a unit
at every test prime, as holds under the three original
radical caps \(n^{h_{\rm rad}}\) in balanced charts.
Include the fixed odd repair order \(M\) in the
progression conditions and choose the markers off its
support, so \(M\) is coprime to every \(q_j-1\).
Use the known common tree constants and valid split
locus of (SS6). The spectator law may
come from any one actual complete dictionary satisfying the
CCF4 hypotheses. It is independent of the original
histories.

On each connecting bipartite edge of a paired block choose
an independent positive rational \(U_e\) from the
central-window law with \(H=n^\kappa\), \(0<\kappa<\gamma\).
Assign to a leaf \(i\) the product \(S_i\) of its incident
edge units. The two block products agree as positive
rationals, because each edge occurs once on each side.
Also \(S_i=1+O_\ell(n^{-\kappa})\) at every leaf.
The uniform unit law on the edges pushes forward onto
the paired-product subgroup (AT1). A relative density
error \(\rho=O(n^{-\gamma+\kappa+o_X(1)})\) on the
joint CRT unit group stays \(O_\ell(\rho)\) after
pushing forward, by the fixed number of edges. Conditional
on the original rows, the products \(r_i S_i\) therefore
have this law relative to the corresponding CRT coset.

All edge factors and leaf products have polynomial
rational height. Clear their denominators by a common
positive integer \(D_0\equiv1\) modulo the product \(Q\)
of the group primes. This is possible at polynomial
height, since all denominators are units and a multiple
with this congruence can be chosen at most \(Q\) times
their product. At a leaf take
\[
   R_i=\mu_M(|B_i|)(D_0 S_i)^M,\qquad
          |B_i'|=E_i(r_iD_0 S_i)^M
\]
in the exact transformation (BS1). The original defect is
still multiplied by \(R_i\). On dividing the common
scale \(D_0\) out of the powered root coordinate, the
exact \(r_iS_i\) is within a relative \(O_\ell(n^{-\kappa})\)
of \(r_i\). The marked sign and the original three
radical caps have the scope of (SS1) and its signed
orientations. No cap on the auxiliary support is needed.
Choose the split data at each prime by (SS5), and choose
one positive integer transformation parameter \(T\) for
each node by CRT, as in (SS7). The integer decorations
then realize all the finite-field diagrams in this
group. The real window is on the normalized leaf scale,
and imposes no extra constraint on the real node
parameters.

Under the ideal uniform coset the group primes are
conditionally independent given the original histories.
Let \(b_j\) denote the local CT1 or CT2 bound via (SS4)
at \(q_j\) with exact uniform law. For the product of the
signed statistics \(Z_j\) in this group the conditional
bound is
\[
            \left|\mathbf E\prod_{j}Z_j\right|
                          \ll_\ell \prod_j b_j+O_{\ell,m}(\rho).
\]
Indeed the exact uniform expectation factors. For its
difference from the implemented law use the relative
density bound against the uniform coset and
\(\mathbb E_{\rm unif}\prod_j|Z_j|\ll_{\ell,m}1\).
The latter follows from independence across primes
and Cauchy--Schwarz on the uniform marginal diagrams,
with their bounded tree \(L^2\) norms. Use the zero
convention on invalid splits. This compares the joint
law at once and pays its approximation error once.
With several independently sampled groups on the same
original rows, the displayed bounds multiply over the
groups, including their errors. The separate integer
decorations may differ from group to group.

At every prime the normalized tree data still give
the (SS6) list of candidate invariants
\((P_1/P_2,P_3,\ldots,P_b)\) in the original radicals.
This follows from the exact paired-product identities
and is unaffected by the common clearing factor. Combining
the lists across groups of distinct primes by CRT
reconstructs the original rational invariants when the
product of those primes exceeds \(n^{4h_{\rm rad}}\). Assign prime supports as
in (SS6) to get the same \(n^{\sigma+o_X(1)}\) list
price when the original shared radical is at most
\(n^\sigma\). This allows joint real control on actual
spectators with complete original-label recovery.
The signed comparison still requires a lower bound
sensitive to all three original letters.
[Paired-product edge law (AT1)](ABC_SIGNED_AUXILIARY_TREE_COMPARISON.md#the-tree-estimate-and-the-law-of-shared-slots)
· [Exact scale cosets and split reconstruction (SS1–SS7)](ABC_SHARED_SCALE_BILINEAR_TREE_HISTORIES.md#several-auxiliary-primes-in-one-integer-history)
· [Uniform signed cycle bounds (CT1–CT2)](ABC_TWO_SIDED_AUXILIARY_TREE_CANCELLATION.md#hypotheses-and-the-two-sided-projection)

## Comparison with occurrence

The input is one actual complete-value dictionary with cardinality
close enough to its root-box volume. Its internal root coordinates
may be dependent. The convolution and the central-window lower
bound still use independent repeated draws from that dictionary.
A matched three-letter row can have cardinality exponent \(d\)
in a joint box of exponent close to \(3d\). Applying (CCF1) to
that joint law therefore does not give the required independence.

There is also an exponent obstruction to the direct positive
replication using independent spectators. Suppose an occurrence
matching of exponent \(d>0\) uses three marginal dictionaries
each of exponent \(d\). One draw from each marginal gives the
true matched triple with probability \(n^{-2d+o_X(1)}\).
Independently multiply its three letters by three copies of
the same positive spectator law, each a product of \(2h\)
draws from a dictionary of exponent \(d\). The guaranteed
exact-positive histories come from the three spectators being
equal. Their probability is \(n^{-4hd+o_X(1)}\), since each
exact positive product has height-subpower fibre by the
divisor bound. More explicitly, if the spectator product law
is \(\mu\), then
\(\sum_z\mu(z)^3\le\|\mu\|_\infty^2\) gives the upper
exponent and \(\sum_z\mu(z)^3\ge|\operatorname{supp}\mu|^{-2}\)
gives the lower one. The guaranteed exact signal is thus only
\(n^{-(4h+2)d+o_X(1)}\). The independent CRT-uniform modular
additive test has probability of order \(n^{-m\beta+o_X(1)}\).

This positive signal cannot exceed that scale in the frame
range near \(T=d\). In fact \(\delta_h>0\) for the type with
all \(m\) components nonquadratic forces
\(h s_{m,\{1,\ldots,J\}}>1/2\). Since
\(s_{m,\{1,\ldots,J\}}\le T/(4m\beta)\), it follows that
\(hT>2m\beta\). In particular \(d\ge T/2\) already makes
\(4hd>m\beta\). This comparison concerns the subset of guaranteed
exact histories. It does not bound possible additional
identities among the spectator products. Using one common
spectator on the three letters instead preserves exact
occurrence without that probability cost, but the common unit
then cancels from the additive congruence. A signed construction
must retain the original incidence without either loss.

The whole-modulus same-target collision bounds in the external
abc work use an actual two-monomial congruence and exact common
powers after recording complementary factors. Their diagonal
term leaves the \(1/3\) balanced barrier. The inverse-sieve
theorems for additive structure modulo primes require restrictions
on the same set at many prime moduli. Neither input supplies
the joint primitive incidence required here.
[Zeng, §§3–5](https://hexagonmath.org/pdf/2611.00621v1)
· [Croot–Mao–Yip, Theorem 1.11](https://arxiv.org/html/2510.08862v2#S1.SS2)
· [The external collision diagonal floor](ABC_EXTERNAL_MOVING_ROWS_LIMIT.md#the-diagonal-floor-and-an-actual-arithmetic-separation)
· [The arithmetic lower bound needed by the tree comparison](ABC_SIGNED_AUXILIARY_TREE_COMPARISON.md#what-the-positive-input-would-need-to-retain)
