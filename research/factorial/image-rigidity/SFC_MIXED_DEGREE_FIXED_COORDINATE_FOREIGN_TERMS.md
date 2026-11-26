# Mixed-degree factorial windows with foreign terms on fixed coordinates

Let \(b_1,\ldots,b_s\) be distinct nonnegative integer exponent
columns. Write \(d_i=\sum_\nu b_{\nu i}\) and
\(S_i=\{\nu:b_{\nu i}>0\}\). A degree block \(G\) consists of
all columns of the same degree \(d_G\). For each block of at
least two columns, let \(C_G\) be the coordinates where all
columns in \(G\) have the same **positive** exponent, and
write \(b_C\) for this common restriction. Assume that for
each foreign column \(j\notin G\), at least one of the
following conditions holds.
\[
          S_j\subseteq C_G,\qquad
          S_j\setminus S_i\ne\varnothing\quad\text{for every }i\in G .
 \tag{1}
\]
There are no restrictions at singleton degree blocks.

**Theorem.** For all sufficiently large integers \(D\), the
ordinary factorial moments of \(\sum_i c_iX^{D b_i}\) have
the actual support-size window at every positive starting
order, for all coefficients in characteristic zero. There
is also a uniform threshold for ordinary integer \(Db_i\)
at real scale \(D\) on compact normalized shape families
with positive lower bounds on their nonzero entries,
center distances, and distinct degree differences.

The first case of (1) permits a foreign term to survive at
every peak of a repeated-degree block. For example the
theorem applies to \(ZX,ZY,Z^3,Z^4\), with arbitrarily
chosen coefficients after dilation, and to any finite
number of repeated blocks satisfying (1). A foreign
constant is included. At a repeated block the surviving
foreign terms depend only on coordinates with fixed
exponents. The other active coordinates undergo a common
radial change, preserving the angular critical equations.
The resulting critical value is one analytic function of
the original saddle frequency throughout the block.

This extends the
[incomparable-coordinate theorem](SFC_MIXED_DEGREE_INCOMPARABLE_COORDINATES.md),
the [singleton extension](SFC_MIXED_DEGREE_SINGLETON_FOREIGN_BLOCKS.md),
and the [constant theorem](SFC_HOMOGENEOUS_BLOCK_WITH_CONSTANT.md).
The proof applies the constant theorem's finite derivative
transfer to several analytic degree annuli at once.

## The common saddle map of a repeated block

Work over \(\mathbb C\) and discard zero coefficients. The
condition persists under this restriction. The
[exponential-range theorem](SFC_MIXED_DEGREE_TOP_BLOCK_TAIL_WINDOWS.md#all-coefficient-vectors-up-to-an-exponential-starting-order)
reduces to a sequence of possible vanishing windows starting
at \(m>e^{\eta D}\), with \(D\to\infty\) and \(\eta>0\) fixed.
Put \(t=m+\ell\) for the finitely many offsets
\(0\le\ell<s'\), where \(s'\) is the actual support size.
Use the raw scaled coefficient \(q_i(t)\) and the pure peak
\(P_i\) in the
[separated-support saddle proof](SFC_DILATED_SUPPORT_SADDLE_ENVELOPE.md).
Let \(Q=\max_i|q_i(m)|\). Choose eligible labels by
\(|q_i(m)|\ge e^{-\gamma D}Q\), taking \(\gamma>0\)
small after \(\eta\) and all the fixed shape bounds are
chosen. The noneligible boxes and exterior contribution
are \(O(e^{-cmD})\) after the scalar normalization in the
[incomparable-coordinate proof](SFC_MIXED_DEGREE_INCOMPARABLE_COORDINATES.md#two-local-estimates).
A singleton uses the full shifted expansion and bounds of
[the singleton proof](SFC_MIXED_DEGREE_SINGLETON_FOREIGN_BLOCKS.md#the-shifted-expansion).

Fix an eligible repeated-block label \(i\in G\), and write
\(C=C_G\), \(V=C^c\), \(d_V=d_G-\sum_{\nu\in C}b_\nu\).
One has \(d_V>0\), since the columns in the block are
distinct. Factor its rescaled homogeneous block as
\[
 H_G(y)=\sum_{i\in G}\frac{c_i(mD)^{Dd_G}}{Q}y^{Db_i}
                = y_C^{Db_C}\widetilde H_G(y_V).
\]
Let \(\theta_i\) be the block-only saddle with the
inactive coordinates zero. Its critical value
\(w_i=e^{-Dd_G}H_G(\theta_i)\) is nonzero and
\(w_i=(q_i(m)/Q)e^{O(e^{-cD})}\) by peak separation.
In the coordinates \(C\), the active critical equations
give \(\theta_{i,C}=b_C\) exactly. In \(V\),
Euler's identity on the active coordinates gives
\(\sum_\nu\theta_{i,V,\nu}=d_V\).

For a foreign survivor \(j\), meaning \(S_j\subseteq C\),
put
\[
 \kappa_{Gj}=e^{-Dd_G}\frac{c_j(mD)^{Dd_j}}{Q}b_C^{D b_{j,C}},
       \qquad B_j=(\kappa_{Gj}/w_i)(t/m)^{D(d_j-d_G)}.
 \tag{2}
\]
There is no dependence on \(i\) in \(\kappa_{Gj}\).
It equals \((q_j(m)/Q)P_j(b_i)^D\) for any \(i\in G\),
because \(b_{i,C}=b_C\) and \(S_j\subseteq C\). Its modulus
is \(O(e^{-c_0D})\) with \(c_0\) from fixed peak separation.
In particular all \(B_j\) are exponentially small for
eligible \(i\) on a relative complex time disc
\(|t/m-1|<a/D\).

Keep the survivors and temporarily omit the foreign
terms in the second case of (1). The full saddle has
coordinates \(x\) on \(C\) and \(\alpha\theta_{i,V}\) on
\(V\), with \(x,\alpha\) given by the critical equations in
the following function:
\[
 \begin{split}
 \mathcal L_G(x,\alpha;B)
   &=\log\left(T(x)\alpha^{Dd_V}
                       +\sum_j B_j T_j(x)\right)\\
   &\quad-D\left(\sum_{\nu\in C}(x_\nu-b_\nu)
                                      +d_V(\alpha-1)\right),\\
 T(x)&=\prod_{\nu\in C}(x_\nu/b_\nu)^{D b_\nu},\qquad
 T_j(x)=\prod_{\nu\in C}(x_\nu/b_\nu)^{D b_{\nu j}} .
 \end{split}                                                     \tag{3}
\]
Indeed the foreign survivors have no \(V\) dependence,
so the active \(V\) gradient on a rescaled old saddle
has all entries equal. Their common critical equation
is the \(\alpha\) equation. The \(C\) equations are
exactly those for \(x\).

At \(B=0\), the Hessian of \(\mathcal L_G/D\) at
\((b_C,1)\) is diagonal with entries
\(-1/b_\nu\) and \(-d_V\). On a fixed-radius \(1/D\)
polydisc around this point all the displayed powers and
their order-\(h\) derivatives are bounded by \(O_h(D^h)\).
The foreign correction to the gradient of
\(\mathcal L_G/D\) is \(O(\sum |B_j|)\), and to its
derivative \(O(D\sum |B_j|)\). Contraction using the
inverse pure Hessian therefore gives the unique local
solution holomorphic in \(B\) for
\(\sum |B_j|\le a_0/D\) with \(a_0\) small enough. Its
displacement is \(O(\sum |B_j|)\). The critical value
\(g_G(B)\) satisfies
\[
                 |g_G(B)|=O\!\left(\sum_j |B_j|\right).
\]
This follows from the quadratic size \(O(D(\sum |B_j|)^2)\)
of the change of the pure phase and the size
\(O(\sum |B_j|)\) of the foreign logarithmic term. Its
first derivatives are bounded, also directly from (3)
since the critical gradient vanishes. The argument is
uniform on the stated compact shape sets.

Thus the mode exponential at \(t\), before the common
normalization \(e^{Dd_*\ell}\), has base
\[
                 (t/m)^{Dd_G}w_i\exp(g_G(B)).
 \tag{4}
\]
Reinstating any omitted foreign term costs
\(C_L(mD)^{-L}\) times its own saddle-size bound for every
fixed \(L\). In fact each omitted quotient vanishes to
order proportional to \(D\) in an inactive coordinate
of this peak. By Cauchy's formula its modulus on smaller
complex boxes is bounded by \(C_R e^{-cD}\) times the
\(R\)-th power of that coordinate for any fixed \(R\).
On the translated contour the corresponding real linear
endpoint loss absorbs the log-phase perturbation; integrating
the difference of exponentials gains \((tD)^{-R}\) from
the endpoint integral. Taking \(R\) large absorbs the
factor \(t\) from exponentiation. This is precisely the
[mode-relative local comparison](SFC_MIXED_DEGREE_INCOMPARABLE_COORDINATES.md#two-local-estimates),
now at the full surviving saddle (3). None of the omitted
terms changes the value or active gradient at the endpoint
itself.

## Universal analytic functions for the filter

Let \(T_0=\ell/m\), \(\mu_j=D(d_j-d_G)\),
and
\[
 {\cal G}_G(w,T_0)
       =(1+T_0)g_G\!\left(
          \left((\kappa_{Gj}/w)(1+T_0)^{\mu_j}\right)_j\right).
\]
Here the notation means that \(g_G\) is evaluated at the
displayed vector. Put
\[
 R_G=e^{D(d_G-d_*)},\qquad
 U_i=w_i e^{{\cal G}_G(w_i,0)},\qquad
 v_i=R_G w_i e^{\partial_{T_0}{\cal G}_G(w_i,0)} .
 \tag{5}
\]
Thus \(U_i\) is the base at \(m\), while \(v_i\) includes
the first-order phase variation of (4). The map
\(w\mapsto v/R_G\) differs from the identity by
\(O(D\sum_j|\kappa_{Gj}|)\) on
\(e^{-2\gamma D}<|w|<8\), with logarithmic derivative
of its multiplicative correction \(O(D e^{-cD})\).
Contraction near each point in a fixed enlargement
of the eligible normalized annulus gives a holomorphic
inverse \(w=w_G(v/R_G)\), as in the
[constant theorem](SFC_HOMOGENEOUS_BLOCK_WITH_CONSTANT.md#radial-change-of-the-full-saddle).
Define, using this inverse,
\[
 H_G(v,T_0)={\cal G}_G(w,T_0)-{\cal G}_G(w,0)
                           -T_0\partial_{T_0}{\cal G}_G(w,0).
 \tag{6}
\]
This function is \(O(e^{-cD})\) on relative frequency
annuli with a fixed radial enlargement and on
\(|T_0|<a_1/D\), and has a double zero in \(T_0\).
It has the same type of bounds under fixed logarithmic
\(v\) derivatives by relative-disc Cauchy estimates.
For eligible repeated-block labels the shifted mode is
\[
 U_i^m F_{G,\ell}(v_i) A_i(m+\ell),\quad
 F_{G,\ell}(v)=v^\ell
 \exp\left(Dd_G[(m+\ell)\log(1+\ell/m)-\ell]
                      +mH_G(v,\ell/m)\right).
 \tag{7}
\]
At singleton groups set \(H_G=0\), because their phase
variation beyond first order is already in their amplitude.
For all groups, \(F_{G,\ell}(R_G x)/R_G^\ell\)
differs from \(x^\ell\) by \(O(D/m)\) uniformly on
\(a_D\le |x|\le 4\), where \(a_D=\frac1{16} e^{-\gamma D}\).
This uses the double zero in (6) and fixed \(\ell\).

At repeated labels the peeled amplitude is
\(A_i(t)=(tD)^{-\rho_i}{\cal B}_i(t)\),
\(\rho_i=r-\#S_i/2\). Its leading amplitude differs
from the positive pure amplitude by \(O(e^{-cD})\) and the
peeled error is \(O((mD)^{-1/2})\) on a complex time
disc of radius comparable to \(m/D\), by the
[localized sector estimate](SFC_HOMOGENEOUS_BLOCK_WITH_CONSTANT.md#radial-change-of-the-full-saddle)
applied to (3). The small displacements preserve the
negative real quadratic and linear endpoint bounds
uniformly on translated contours and their sides.
Cauchy estimates on the peeled amplitude give
\[
   \frac{m^h\Delta^h A_i(m)}{(mD)^{-\rho_i}}
          =(-1)^h(\rho_i)^{\overline h}K_i^*+o(1)
            \quad\text{for each fixed }h\ge0,\qquad K_i^*>0 ,
 \tag{8}
\]
because \(D^h((mD)^{-1/2}+e^{-cD})\to0\).
For singleton labels the estimate at \(h=0\) and the
weaker \(O((mD)^{-\rho_i}/m)\) bound on each fixed
positive difference suffice. The effective frequencies
lie in separated degree annuli
\(e^{-\gamma D}R_G/4\le |v_i|\le 2R_G\).
A raw maximizing label still has
\(|U_i|\ge\exp(-C e^{-cD})\).


## Analytic interpolation across degree annuli

For at most \(s'-1\) prescribed nodes \(w_j\), each with its
degree group and possible multiplicity, seek a monic filter
vector \((p_0,\ldots,p_n)\), \(n\) the number of nodes, such
that
\[
                      P_G(v)=\sum_{\ell=0}^n
                                      p_\ell F_{G,\ell}(v)
\]
has the prescribed roots in group \(G\). The coefficients
are uniformly bounded and on each expanded group annulus
\[
                   P_G(v)=\prod_{\rm prescribed\ j}(v-w_j)
                                        (1+O(e^{-cD})).        \tag{9}
\]
The error is relative and holomorphic after factoring.

Here is the extension of the
[annular determinant argument](SFC_MIXED_DEGREE_INCOMPARABLE_COORDINATES.md#interpolation-on-separated-annuli).
Normalize each group's nodes and possible test point by
\(R_G\). In a group of size \(n_G\), replace the rows by
Newton divided differences in the normalized coordinate
\(x=v/R_G\). For each exponent \(\ell\) this row retains the
factor \(R_G^\ell\). The entries for \(x^\ell\) in the
divided matrices are bounded polynomials in the nodes;
by the annular Cauchy formula their perturbations due
to (7) are at most \(C D e^{C_s\gamma D}/m\).
Indeed the two contour circles at radii \(a_D,4\)
stay at least \(a_D\) from the normalized nodes. This
formula extends the divided matrices through repeated
nodes and gives their prescribed confluent conditions
as in the
[constant theorem](SFC_HOMOGENEOUS_BLOCK_WITH_CONSTANT.md#interpolation-of-analytic-modes-on-two-annuli).

Partition the exponents among the row groups in the
determinant expansion. The ordered partition assigning
larger exponents to larger \(R_G\) has the largest
radial factor \(\prod_G R_G^{\sum \ell}\). Any other
partition loses at least \(e^{-Dg}\), where \(g\) is
the least distinct normalized degree difference, by
swapping an inverted pair of exponents. For the ordered
partition of consecutive exponents in the denominator,
the monomial divided determinants are
\(\prod_{j\in G}x_j^{L_G}\), with \(L_G\) the least
exponent assigned to the group. They are bounded below
by \(e^{-C_s\gamma D}\). Take \(\gamma\) small in
terms of \(s,\eta,g\), so both the other partitions
and all the analytic corrections are exponentially small
relative to this leading term. Thus the divided
denominator is nonzero and differs from the ordinary
polynomial divided determinant by exponentially small
relative error.

For a Cramer minor omitting one exponent from
\(\{0,\ldots,n\}\), the ordered radial factor is no
larger than the denominator's, since \(R_G\le1\)
and the ordered exponent list only increases.
Its leading divided determinants contain the same
monomials \(\prod x_j^{L_G}\) times bounded polynomials,
as follows by factoring the smallest exponent at
each node. The other partitions and analytic errors
are controlled by the same losses, so the Cramer
coefficients are bounded. Finally augment the determinant
by one evaluation point assigned to its own group.
Apply the same relative comparison after dividing
by the enlarged Vandermonde in that group. The
augmentation contributes precisely its factors
\(x-x_j\) with prescribed same-group nodes. Dividing
the augmented determinant by the denominator yields
(9), including at confluent nodes after canceling
these exact factors. Cauchy's formula on discs of
radius comparable to \(|v|\) gives the derivative
bounds on the holomorphic relative error used below.

## Testing repeated and singleton modes

For a repeated group use the finite
[derivative transfer](SFC_HOMOGENEOUS_BLOCK_WITH_CONSTANT.md#from-shift-powers-to-analytic-derivatives)
with \({\cal X}=(v\partial_v)/m\). Its hypotheses
hold by (6) and (7), since the pure degree factor
in (7) is independent of \(v\). For completeness,
\({\cal X}F_{G,\ell}/F_{G,\ell}=T_0+v\partial_v H_G\).
The ratio for higher derivatives satisfies
\(L_{h+1}=(T_0+v\partial_v H_G)L_h+{\cal X} L_h\),
\(L_0=1\). For each fixed \(J\), the finite Taylor
coefficient matrix of \(L_0,\ldots,L_J\) in \(T_0\)
through degree \(J\) is the identity plus \(o(1)\),
because (6) is exponentially small on a disc of
radius comparable to \(1/D\). Its remainder at
\(|T_0|\le s'/m\) is \(O(m^{-J-1})\), absorbing
fixed powers of \(D\) in that exponential bound.
Consequently testing the exact finite Newton expansion
of an amplitude in \(\ell\), and using (8), gives
for repeated labels
\[
 U_i^m(mD)^{-\rho_i}
 \left[
  \sum_{h=0}^J(\lambda_{ih}+o(1)){\cal X}^hP_{G_i}(v_i)
                         +O(m^{-J-1})\right],                 \tag{10}
\]
where \(\lambda_{ih}=K_i^*(-1)^h(\rho_i)^{\overline h}/h!\)
for \(0\le h\le n\) and vanishes for \(h>n\).
A singleton uses the unmodified Newton-difference formula
\[
         U_i^m\sum_{h=0}^n
                   \frac{v_i^h P_{G_i}^{(h)}(v_i)}{h!}
                                      \Delta^h A_i(m)
\]
and its bounds from the singleton theorem.

Use the same relative padding and maximize as in the
[annular cluster argument](SFC_MIXED_DEGREE_INCOMPARABLE_COORDINATES.md#relative-clusters-and-the-window-contradiction):
\[
 d_{ij}=|v_i-v_j|+|v_i|/m,\quad
 V_i=\prod_{j\ne i}d_{ij},\quad
 W_i=|U_i|^m(mD)^{-\rho_i}V_i,\quad W_p=\max_i W_i .
\]
Pass to subsequences and let \(I\) be the cluster where
\(m(v_i/v_p-1)\) stays bounded, with limits \(\xi_i\).
It lies in a single degree group by annular separation.
Prescribe all outside nodes and \(k\) copies of \(v_p\)
for \(0\le k<|I|\). The normalization of the resulting
test is \(U_p^m(mD)^{-\rho_p}(v_p/m)^k
[P_{G_p}(v)/(v-v_p)^k]_{v=v_p}\), of modulus comparable
to \(W_p(|v_p|/m)^{k-|I|+1}\) by (9).
All approximation errors disappear after normalization:
\(V_i\ge e^{-C_sD}m^{-(s'-1)}\) and \(m>e^{\eta D}\)
allow fixed \(J,L>s'+C_s/\eta\) in (10) and the local
errors. The exponential exterior bound is negligible
relative to the raw-maximizer lower bound on \(W_p\).

For an outside repeated node, a nonzero term of
\({\cal X}^hP_{G_i}(v_i)\) differentiates its own exact
root. Expanding iterated Euler derivatives into ordinary
derivatives, charge at least \(|v_i|/m\) to this removed
factor, and use the padded distances for other factors.
Cauchy bounds from (9) give a cost \(O(|v_i|^{-1})\)
per derivative on the holomorphic relative factor.
The resulting bound for (10) relative to the
normalizer is at most a constant times
\[
 \frac{W_i}{W_p}\,
 \frac{|v_i|/m}{d_{ip}}
                  \left(\frac{|v_p|/m}{d_{ip}}\right)^{|I|-1-k},
\]
which tends to zero by the relative separation.
An outside singleton has one prescribed root in its
group and needs only the one \(1/m\) loss, so the
[singleton filter estimate](SFC_MIXED_DEGREE_SINGLETON_FOREIGN_BLOCKS.md#filtering-singleton-modes)
gives the same conclusion.

If \(p\) is a singleton alone, its \(k=0\) test has a
surviving positive prefactor and all its positive
differences vanish in the normalizer. If the cluster
lies in a repeated group, \(U(v)/(v/R_G)\) from (5)
has logarithmic derivative \(O(D e^{-cD})\). Hence
\((U_i/U_p)^m\to e^{\xi_i}\) for \(i\in I\). The
padding and maximality imply \(\rho_i\ge\rho_p\)
within \(I\), with strict inequality discarding the
term. The holomorphic factor of \(P_{G_p}\) after
dividing by \((v-v_p)^k\) has ratio tending to one
from \(v_p\) to \(v_i\); its derivatives on the
relative scale \(v_p/m\) vanish by the outside
distances in (9). Thus in (10) only derivatives
of the \(k\)-fold factor survive, giving the
limit equations
\[
    \sum_{\substack{i\in I\\\rho_i=\rho_p}}
       K_i^* e^{\xi_i}
          \sum_{h=0}^k\binom{k}{h}\xi_i^{k-h}
                      (-1)^h(\rho_p)^{\overline h}=0
                \qquad(0\le k<|I|).
\]
These are the consecutive-degree monic Vandermonde
equations. Coincident \(\xi_i\) combine with weights
\(e^{\xi_i}\) times positive sums, so they are
impossible. This contradicts the vanishing window
and proves the theorem. Characteristic zero follows
by embedding the finitely generated coefficient
field into \(\mathbb C\).

The window convention is that of
[Edo and van den Essen, Conjecture 2.8](https://arxiv.org/abs/1304.3956).
Condition (1) leaves foreign terms that survive on
varying coordinates of repeated blocks outside the
scope of the common saddle map.


