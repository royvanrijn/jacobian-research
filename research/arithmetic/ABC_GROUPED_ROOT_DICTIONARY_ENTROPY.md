# Grouped root dictionary bounds and balance

The root-count exponent can be controlled by joint dictionary sizes on
bounded groups of high-degree coordinates. The cost is the total of the
group cardinality exponents. The numerical root height caps remain the
ones from the original anatomic coordinates. A quarter inequality
forces all three dictionary cost totals to balance at saturation.

## Group dictionary count

Use primitive triples in the signed anatomic equation
\[
 \sum_{T=1}^{3}\epsilon_T t_{T,0}
       \prod_{e} x_{T,e}^{\,e}=0,\qquad \epsilon_T\in\{-1,1\}.
\]
The coefficients \(t_{T,0}\) are fixed in each box and have numerical
exponent \(o(1)\) in the small-budget limit, as in the fixed-degree
cutoff. The roots are pairwise
coprime across all nonzero degree positions of the three letters. Let
\(x_{T,e}\le X^{\alpha_{T,e}+o_X(1)}\), where the raw exponent sum
\(L=\sum_{T,e}\alpha_{T,e}\) is at most \(C_0\lambda\) for fixed
\(C_0>0\). Work in dyadic bins for the complete letter values with maximum
comparable to \(X\). At each fixed budget \(\lambda\), fix the anatomic
truncation and the finite degree lists before taking the height
limit. All bounds below are in this order of limits.

Partition each letter's root positions into groups \(G\), each of
at most \(m\) positions for a fixed \(m\). For each nonsingleton
group require all its degrees to be at least \(a/\lambda\), with
fixed \(a>0\). The group dictionary of root vectors has cardinality
at most \(X^{W_G+o_X(1)}\), where \(W_G\le\sum_{e\in G}\alpha_{T,e}\).
Put \(W_T=\sum_{G\subset T}W_G\) and \(W=\sum_T W_T\). For every \(\eta>0\), at sufficiently small fixed budgets the
count is at most
\[
 X^{(W+\min_T W_T)/4+\eta\lambda+o_X(1)}.
 \tag{GD2}
\]
In particular it is at most \(X^{W/3+\eta\lambda+o_X(1)}\).
The groups and dictionaries may vary. The budget threshold depends
on \(C_0,m,a,\eta\) and the fixed uniform height-bin constants.
The \(o_X(1)\) is taken with the budget and truncation fixed.

Refine all original root ranges to dyadic bins at subpower cost
before proving the bound. Replace each raw cap by the upper
exponent of the refined bin, and replace each group cap by its
minimum with the new sum of raw caps in that group. These changes
can only strengthen the required bound. Thus, inside a refined
bin, each root has logarithmic height
\(\alpha_{T,e}+o_X(1)\) as well as the asserted upper cap.
At fixed budget a height subsequence fixes the limiting raw caps,
group caps and count exponent. The dyadic root loss and arbitrary
inner exponent allowances are taken before the outer budget limit.

### Diagonal moments for a joint dictionary

For purposes of an exponent bound, one may first restrict the
sample by a prime colouring, with one colour for each original
root position. For any row, the probability that every prime in its
roots has the designated position colour is at least
\(c^{-\omega(abc)}=X^{-o_X(1)}\), where \(c\) is the number of
colours at the fixed truncation. Thus some colouring retains the
original exponent. All group dictionaries below can be restricted
to that colouring and hence retain the stated upper caps.

Consider a full modulus \(M\ge X^{w_0+o_X(1)}\) from a locked
collection of coordinates on one letter, for some fixed
\(w_0>0\), with fullness threshold at least \(\kappa/\lambda\) for
fixed \(\kappa>0\). Use the integer fullness parameter
\(k=\lfloor\kappa/\lambda\rfloor\) in (FC2). On a free letter consider a group \(G\)
of size at most \(m\), and restrict its dictionary to unit tuples
modulo \(M\). For a character \(\chi\) of the unit group, put
\[
 S_G(\chi)=\sum_{\boldsymbol x\in\mathcal D_G}
     \chi\left(\prod_{e\in G}x_e^{\,e}\right).
\]
In the \(2h\)-th moment choose the first \(h\) group tuples,
at most \(X^{hW_G+o_X(1)}\) choices, and collapse the second
\(h\) tuples into products on each coordinate. If
\(h\max_{e\in G}\alpha_{T,e}=O(\lambda)\), the collapsed
products have height \(X^{O(\lambda)}\). The colouring gives
pairwise coprimality of the products across positions. Equation
(FC2) then counts the collapsed congruence solutions at
\(X^{o(\lambda)+o_X(1)}\), uniformly in the residue and the
fixed group exponents. A collapsed tuple has divisor-many
preimages at fixed \(h\). Consequently the normalized character
average of \(|S_G|^{2h}\) is at most
\(X^{hW_G+o(\lambda)+o_X(1)}\).
[Bounded-root congruence count](ABC_OMITTED_COORDINATE_GRAPH.md#the-determinant-count-and-amplification)

To combine the groups on the free letter, use Hölder with
\(h_G=\lceil L_T/\alpha_G\rceil\), where
\(\alpha_G=\sum_{e\in G}\alpha_{T,e}\) and \(L_T\) is
the sum of those raw group masses on that letter. Then
\(\sum_G1/h_G\le1\) and the collapsed height cap is at most
\(X^{2L_T+o_X(1)}\). At fixed truncation one can first pass to
limiting raw box caps as \(X\to\infty\), freeze groups with
arbitrarily small total raw caps in the energy copies, and
then take the moment orders fixed for the remaining groups.
All the resulting height-exponent allowances can be absorbed
at this stage. Hölder and the trivial energy bound for a
letter with zero raw caps give residue-energy exponent at
most \(W_T+o(\lambda)\).

If a prescribed collection of groups on letter \(T\) gives
such a full modulus, write \(w_{\rm lock}\) for their summed
cost \(W_G\). Freezing them costs at most \(w_{\rm lock}\).
The free letters are equal up to a fixed unit twist modulo
the modulus. The Cauchy–Schwarz residue bound therefore yields
\[
                    D\le w_{\rm lock}+\tfrac12(W-W_T)+o(\lambda).
 \tag{GD1}
\]
Two complete free values determine the primitive triple.
This uses the actual congruence relation on every lock,
including groups whose variables have different powers.

### Supported letters and the quarter inequality

Suppose a sequence violates (GD2) with a fixed positive normalized
margin. Pass to convergent subsequences of the totals and count
exponents. Use \(\lambda\) as the unit of cost and write \(d=D/\lambda\),
\(w_T=W_T/\lambda\), \(w=\sum_T w_T\) for the limits. Put \(u=w-2d\).
Call a letter supported when \(w_T>u\).

On a supported letter \(T\), any block of groups with a fixed positive
limiting numerical weight and all its degrees at least
\(\gamma/\lambda\), for some fixed \(\gamma>0\), must cost at least
\(\theta_T=(w_T-u)/2>0\). This follows from (GD1). The following
concentration argument applies to any letter with such a lower bound
\(\theta>0\) on the cost of every positive-weight full block. Set aside
all groups of normalized cost at least \(\theta/4\). Their number is
bounded by \(4C_0/\theta\). For any fixed \(0<\gamma<a\), the remaining
groups with at least one degree below \(\gamma/\lambda\) are singletons.
Their total numerical weight is at most \(C_0\gamma+o(1)\). Partition the
other groups greedily into blocks of cost at most \(\theta/2\). Every
closed block except possibly the last has cost greater than
\(\theta/4\), since each next group has cost below \(\theta/4\). Thus
there are boundedly many blocks. If this part has positive weight along
a subsequence, a block of positive weight costs at most \(\theta/2\),
contradicting the lower bound. Sending \(\gamma\) to zero shows that all
numerical weight outside the set-aside groups tends to zero. The same
conclusion holds for a zero-height letter using the empty set of groups.

After passing to subsequences to fix the cardinalities and the needed
weight limits, each such letter consists of boundedly many narrow roots
of height \(X^{O(\lambda)}\) times a wider cofactor of height \(X^{o(1)}\),
apart from its fixed coefficient. In this paragraph and below the
\(o(1)\) in numerical cofactor heights means convergence to zero in
the budget limit after the fixed-budget height limit. A full modulus
on an opposite letter releases one or two such terms by (GW) at
cost \(o(\lambda)\), using the amplification at any fixed order
before sending the budget to zero. The bound applies to the complete
signed monomial of the released roots and cofactors. An exact complete
ratio recovers the pair of primitive terms, and each narrow variable
divides one of them. For amplification, the prime-position colouring
keeps opposite letter supports disjoint across rows, so the same
exact-ratio fibres and ordered factorizations are divisor-bounded
in the amplified tuples. Fixed coefficients are absorbed in the
congruence target. [Full-modulus wide cofactor count](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-bounded-number-of-small-roots)

We now show that if \(w_C<d\), the other two letters both have cost
at most \(u\). Suppose on the contrary that \(A\) is supported.
Freeze \(C\) in full at cost at most \(w_C\). Any positive-weight full
block on \(B\) with fullness of order \(1/\lambda\) gives a congruence
for \(A\) at a fixed unit target. Its exact value together with the
fixed \(C\) determines the triple. The preceding release therefore shows
that every such \(B\)-block has normalized cost at least \(d-w_C>0\).
The concentration argument then applies also to \(B\).

If \(C\) has positive numerical height, its groups with all degrees
at least \(\gamma/\lambda\), for \(\gamma>0\) sufficiently small,
contain a block of positive weight. Indeed the low groups are
singletons of total weight at most \(C_0\gamma+o(1)\), and the fixed
coefficient in the anatomic boxes has numerical exponent \(o(1)\).
Fix this block at cost at most \(w_C\) and release \(A,B\) jointly
using (GW), giving \(d\le w_C\). If \(C\) has zero numerical height,
\(A,B\) are the tall pair. Their complete positive ratio lies in a
strip of relative width \(X^{-\tau}\) about one for any fixed
\(0<\tau<1\). Release both by (RUW1), with their bounded narrow roots
and the ratio of their wider cofactors. Its height is \(X^{o(1)}\)
and an exact signed ratio has divisor-bounded fibres as above
(without amplification). The fixed coefficients enter the
centre of the strip. For each fixed chart parameter \(s>0\) the
wider-height hypothesis holds eventually. Sending \(s\) to zero
after the budget limit gives \(d\le0\). Each case contradicts
\(w_C<d\). [Unequal-height real strip count](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-real-strip-with-unequal-root-heights)

Finally, \(4d>w+\min_Tw_T\) implies \(d>w_C\) for a cheapest
letter \(C\), since \(w\ge 3w_C\). If both other costs were at most
\(u\), we would have \(w-w_C\le2u=2w-4d\), the opposite inequality.
This proves (GD2). At equality \(d=w/3\) it forces
\(w_A=w_B=w_C=w/3\).

## Consequence for joint root entropy

Let the box sample be uniform and suppose its normalized exponent
\(D/L\) tends to \(1/3\) with \(L/\lambda\to1\) and all three height
exponents bounded away from zero. For any allowed group partition,
the sum of the joint marginal entropies of the groups is
\[
            \sum_G H((x_{T,e})_{e\in G})=(L+o(\lambda))\log X .
\]
To prove the lower bound, cap each realized group surprise
\(-\log_X \Pr(\boldsymbol x_G)\) at \(\alpha_G\).
Bin these capped values with total rounding loss at most
\(\xi\lambda\) for any fixed \(\xi>0\). At fixed truncation the
number of joint bin profiles is independent of \(X\), and each
profile supplies dictionaries with caps equal to the bin tops
or the raw caps, whichever is smaller. A bin profile with
total cost below \(L-\delta\lambda\) is smaller than the
full sample by a fixed power when \(\xi\) is small in terms
of \(\delta\). It follows that the typical total capped score
is \(L-o(\lambda)\) in the sequential limits. Taking expectations
and removing the caps proves the lower bound. Each group's
dictionary size is bounded by its raw coordinate ranges,
which gives the matching upper bound.

In particular any matching of pairs of high-degree coordinates,
each pair within one letter, has sum of pair mutual informations
\(o(\lambda\log X)\) in these limits. The other positions may be
singletons. Apply the joint statement with the pair partition
and subtract it from the raw upper bound for the sum of the
individual entropies. The statement allows dependence carried
by larger collections of positions. The same argument can bin with
separate letter totals. In particular (GD2) yields a saving whenever
the cheapest letter total is smaller than the average by a fixed
multiple of \(\lambda\).
