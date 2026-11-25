# Grouped root dictionary bounds and balance

The root-count exponent can be controlled by joint dictionary sizes on
bounded groups or shrinking valuation windows. The cost is the total of
the group cardinality exponents. The numerical root height caps remain the
ones from the original anatomic coordinates. A two-letter inequality
forces all three dictionary cost totals to balance at saturation.

## Group dictionary count

Use primitive triples in the signed anatomic equation
\[
 \sum_{T=1}^{3}\epsilon_T t_{T,0}
       \prod_{e} x_{T,e}^{\,e}=0,\qquad \epsilon_T\in\{-1,1\}.
\]
The coefficients \(t_{T,0}\) are fixed in each box and have numerical
exponent \(o(1)\) in the small-budget limit, as in the fixed-degree
cutoff. The roots are pairwise coprime across all nonzero degree positions of
the three letters. Let
\(x_{T,e}\le X^{\alpha_{T,e}+o_X(1)}\), where the raw exponent sum
\(L=\sum_{T,e}\alpha_{T,e}\) is at most \(C_0\lambda\) for fixed
\(C_0>0\). Work in dyadic bins for the complete letter values with maximum
comparable to \(X\). At each fixed budget \(\lambda\), fix the anatomic
truncation and the finite degree lists before taking the height
limit. All bounds below are in this order of limits.

Partition each letter's root positions into groups \(G\). A group
can have at most \(m\) positions, for fixed \(m\), or its scaled
degrees can lie in at most \(m_0\) windows of diameter
\(\rho_\lambda\to0\), uniformly in the budget limit. The
groups can include degrees of any size, including a shrinking
window at scaled degree zero. Already frozen original positions
can be included as cost-zero singletons. Their raw caps still
enter \(L\). The group dictionary of root vectors has cardinality
at most \(X^{W_G+o_X(1)}\), where \(W_G\le\sum_{e\in G}\alpha_{T,e}\).
Put \(W_T=\sum_{G\subset T}W_G\) and \(W=\sum_T W_T\). For every
\(\eta>0\), at sufficiently small fixed budgets the count is at most
\[
 X^{(W-\max_T W_T)/2+\eta\lambda+o_X(1)} .
 \tag{GD3}
\]
In particular the coefficient \((W+\min_T W_T)/4\) is valid
(GD2), and hence also \(W/3\). The groups and dictionaries may
vary subject to the uniform margins, the fixed dimension bounds,
and the stated rates of convergence of the window and exceptional
coefficient heights. The budget threshold can depend on these
parameters and \(\eta\).
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


For a bounded-window group the moment can be counted without bounding
its number of positions. After collapsing by original position,
multiply together the positions in each window and use a degree
representative in that window. The new roots have heights
\(X^{O(\lambda)+o_X(1)}\), and the product of signed corrections has
height \(X^{O(\rho_\lambda)+o_X(1)}\). Apply (GW) and its
overlap-compatible amplification to these compressed tuples.
The unchanged prime colours and original degrees recover all
original collapsed products from the exact monomial value. The
same holds for coordinatewise products of compressed tuples.
The exact-output multiplication has divisor-many ordered
factorizations. This proves the claimed moment bound also for
the bounded-window groups. The aggregate height hypothesis
\(h\sum_{e\in G}\alpha_{T,e}=O(\lambda)\) suffices for both
types of group and is available from the choice below.

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
such a full modulus from any subset of its original positions,
write \(w_{\rm lock}\) for their summed cost \(W_G\).
The projected dictionaries for the locked positions have
cost at most \(w_{\rm lock}\), since projecting a dictionary
cannot increase its size. One can also fix the entire groups
at this cost. The free letters are equal up to a fixed unit twist
modulo the modulus. The Cauchy–Schwarz residue bound therefore yields
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

For fixed \(\gamma>0\), call the projection of a group onto its
original positions with degrees at least \(\gamma/\lambda\) its
high projection at \(\gamma\). The product of the corresponding
powers is \(\lfloor\gamma/\lambda\rfloor\)-full, irrespective of
the other positions in the group. Across all groups the discarded
positions have total numerical weight at most \(C_0\gamma+o_X(1)\),
because \(\sum_i\alpha_i\le C_0\lambda\). Thus whenever a collection
of groups on a supported letter \(T\) has high projections of
fixed positive limiting numerical weight, it must cost at least
\(\theta_T=(w_T-u)/2>0\), by (GD1).

The following concentration argument applies to any letter with
such a lower bound \(\theta>0\) on the cost of every block whose
high projection supplies a positive-weight full modulus. Set aside
all groups of normalized cost at least \(\theta/4\). Their number
is bounded by \(4C_0/\theta\). For fixed \(\gamma>0\), partition the
other groups greedily into blocks of cost at most \(\theta/2\),
charging the costs of the original groups. Every closed block
except possibly the last has cost greater than \(\theta/4\),
since each next group has cost below \(\theta/4\). Thus there
are boundedly many blocks. If their high projections have
positive limiting weight along a subsequence, one block supplies
a positive-weight full modulus at cost at most \(\theta/2\).
This contradicts the lower bound. Sending \(\gamma\) to zero
bounds the entire numerical weight outside the set-aside groups
by \(o(1)\). The same conclusion holds for a zero-height letter
using the empty set of groups.

Compress every selected window group as in the moment proof.
The sum of its correction heights with the remaining cofactor
height still tends to zero in the budget limit. After passing
to subsequences for the needed weight limits, each such letter
consists of boundedly many narrow roots
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

If \(C\) has positive numerical height, the high projections of
its groups together give a full divisor of positive weight for
\(\gamma>0\) sufficiently small. The other original positions
have numerical weight at most \(C_0\gamma+o_X(1)\), and the fixed
coefficient has numerical exponent \(o(1)\). Fix these projections
at cost at most \(w_C\) and release \(A,B\) jointly using (GW),
giving \(d\le w_C\). If \(C\) has zero numerical height,
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

### The two-letter effective cost

Now suppose a sequence violates (GD3) with positive normalized
margin. The preceding bound gives \(d\le w/3\). A most expensive
letter \(S\) is supported since \(w_S>w-2d=u\). Three supported
letters would give \(w>3u\), contrary to \(d\le w/3\).
If exactly two are supported, the third costs less than
\(w-2u=4d-w\le d\). This contradicts the conclusion of the
quarter-inequality proof for a letter of cost below \(d\).

Suppose \(S\) is the only supported letter. The other two, \(U,V\),
satisfy \(w_U+w_V<2d\). Label them so that \(w_V\le w_U\);
then \(w_V<d\). Fix \(V\) in full. As in the quarter-inequality
proof any positive-weight full block on \(U\) has cost at least
\(b=d-w_V>0\), from (GW) on \(S\). If \(U\) has positive
numerical height, it has a group with a positive-weight high
projection and cost \(\beta_U\ge b-\xi\) for any fixed
\(0<\xi<b\). To see this, pack the groups costing at most
\(b-\xi\) into greedy blocks of that cap. Any two adjacent
blocks have combined cost greater than the cap, so only
boundedly many blocks occur. Their high projections have
vanishing numerical weight, because a block of positive
projected weight would supply a modulus below cost \(b\).
Choose \(\gamma\) so small that \(C_0\gamma\) is below the positive
height of \(U\). A group of cost greater than the cap then has
a positive-weight high projection along a subsequence. Take
\(\xi<2d-w_U-w_V\), giving \(w_U-\beta_U<d\).
For a subsequent release the entire selected group is counted.
In a window group the bounded number of narrow roots and the
vanishing cofactor height include all its original positions.

If both \(U,V\) have positive height, fix the complement of
the selected group on \(U\). Any positive-weight full block
on \(V\) must cost at least \(d-(w_U-\beta_U)>0\), by releasing
that group and \(S\) under the lock. The same greedy argument
gives a group on \(V\) with a positive-weight high projection,
of cost \(\beta_V\). Lock the high projection of each selected
group in turn at the original group cost and fix the complement
of the entire other group. Release that entire group together with
\(S\). Equation (GW) gives
\[
 d\le\beta_U+w_V-\beta_V,\qquad d\le\beta_V+w_U-\beta_U .
\]
Adding contradicts \(w_U+w_V<2d\).

If one of \(U,V\) has height zero, \(S\) and the other letter
are tall. When that other letter has cost below \(d\), fix it
in full and release \(S\) in the real strip using (RUW1).
When its cost is at least \(d\), it is \(U\), so fix only its
complement of cost \(w_U-\beta_U<d\) and release both \(S\) and
the selected group using (RUW1). Both cases contradict \(d\).
For these releases fixed factors enter the centre, and the
two complete tall values are primitive, with divisor-bounded
representation for each exact ratio. This proves (GD3).
The all-singleton count follows upon taking \(m=1\).

### Conditional dictionary inequality

Suppose a prescribed record \(\mathcal R\) fixes any set of
original root positions. Partition the remaining positions
within each letter into allowed groups and take
\(E_T=\sum_{G\subset T}H(\boldsymbol x_G\mid\mathcal R)\), using
natural logarithms for entropy. Then for each letter \(T\),
\[
 (D-o(\lambda))\log X\le H(\mathcal R)
                             +\tfrac12\sum_{U\ne T}E_U .             \tag{GD4}
\]
Include the frozen positions as cost-zero singletons in (GD3).
Bin the record surprise and the conditional group surprises
with total rounding error arbitrarily small compared with
\(\lambda\log X\); the number of profiles is bounded at fixed
truncation. For each profile and each record, the upper
group surprise gives a uniform dictionary cap. The upper
record surprise bounds the number of records. For a given
\(T\), (GD3) and summation over records make every profile
with realized score sum below \((D-\zeta\lambda)\log X\)
power-small at sufficiently small budgets relative to any
fixed \(\zeta>0\). Taking expectations proves (GD4).
The bound retains all three letters' raw cap and height
constraints even in the orientation where one letter's
dictionary exponents do not enter the right-hand side.


## Fixed complete values and prescribed least-repaired targets

Fix a positive integer \(z\) as the value of any specified complete
letter on a primitive shell \(X/2<c\le X\). Suppose
\(\operatorname{rad}(abc)\le X^\lambda\) and the other letters
\(U,V\) have individual radical caps \(R_U,R_V\); the latter may
be taken to be \(X\) if no separate caps are imposed. For every
\(\eta>0\), at all sufficiently small fixed \(\lambda\), uniformly
in \(z,R_U,R_V\), the count is
\[
 \ll_{\lambda,\eta,\varepsilon}
 X^{\eta\lambda+\varepsilon}
 \min\left\{
    (X^\lambda/\operatorname{rad}z)^{1/4},
    R_U^{1/2},R_V^{1/2}
 \right\}.                                                        \tag{GD5}
\]
The case \(\operatorname{rad}z>X^\lambda\) is empty.

To prove this choose the anatomic truncation with radical and
coefficient-fixing losses at most a sufficiently small multiple
of \(\eta\lambda\), then take the height limit at that truncation.
All exact root tuples on the fixed letter have
divisor-many possibilities at fixed truncation, since each root
divides \(z\). Freeze them at inner subpower cost. Its charged cost
is zero and its original raw cap is retained. On the other two
letters charge the raw caps of the whole root dictionary. Their
sum is at most \(\log_X(X^\lambda/\operatorname{rad}z)\) up to
the truncation and dyadic errors, since the three complete values
are pairwise coprime. The individual sums are at most
\(\log_XR_U,\log_XR_V\) up to the same errors. Apply (GD3) and
bound half the smaller outside cost also by one quarter of
their total. The fixed-budget bin and representation losses
go inside the arbitrary \(\varepsilon\). The corresponding bound
holds with \(X^\lambda\) replaced by a smaller given shell radical
budget. At balance it gives fixed complete-value exponent
\((1/6+o(1))\lambda\). A balanced saturating family therefore
has at least \(X^{(1/6-o(1))\lambda}\) distinct values on each
letter and maximum marginal probability at most
\(X^{-(1/6-o(1))\lambda}\) there. Here these \(o(1)\) follow the
same height-first convention.
[Anatomic reduction and truncation](https://arxiv.org/html/2410.12234v2#S2)

For least target repair of order \(M\), prescribe the actual
repaired target \(H=H_M(c)\). Restrict the primitive quotient
by \(\operatorname{rad}(abc)\le c^\lambda\) and
\(\operatorname{rad}(c)\ge c^{\ell\lambda}\) for fixed
\(0<\ell\le1\). With \(t=\min\{1,(M\ell\lambda)^{-1}\}\),
the resulting pointwise count satisfies
\[
 \ll_{\lambda,M,\ell,\eta,\varepsilon}
 H^{\,t((1-\ell)/4+\eta)\lambda+\varepsilon}.                       \tag{GD6}
\]
In fact \(c\mid H\), \(\operatorname{rad}(c)=\operatorname{rad}(H)\),
and \(c\le H^t\) because \(H\ge c\) and
\(H\ge\operatorname{rad}(c)^M\). For each original \(c\) use
(GD5) with \(X=c\) and \(z=c\). Sum over the divisor-many choices
of \(c\) and use \(c\le H^t\). The candidates have one prescribed
least repair per quotient. If \(c\) is also given, the same
argument gives \(\ll c^{\eta\lambda+\varepsilon}
(c^\lambda/\operatorname{rad}(H))^{1/4}\) directly.
[Least target repair](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)

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

In particular any matching of pairs of coordinates, each
pair within one letter, has sum of pair mutual informations
\(o(\lambda\log X)\) in these limits. The other positions may be
singletons. Apply the joint statement with the pair partition
and subtract it from the raw upper bound for the sum of the
individual entropies. The statement allows dependence carried
by larger collections of positions. The same argument can bin with
separate letter totals. In particular (GD2) yields a saving whenever
the cheapest letter total is smaller than the average by a fixed
multiple of \(\lambda\). The stronger (GD3) forces equality of
the effective letter totals whenever their sum attains the
coefficient one third.
