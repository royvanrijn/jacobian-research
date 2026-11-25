# Upper quantiles after fixed root locks

Use the anatomic endpoint conventions and the three-factor uniqueness observation in [An explicit positive-degree mass in terms of numerical height](ABC_QUANTITATIVE_DEGREE_WINDOWS.md#an-explicit-positive-degree-mass-in-terms-of-numerical-height). The limiting root masses of the letters are \(l_W\), their numerical heights are \(h_W\), and the normalized count exponent is \(d\). Let
\[
 H_*=\min_{P\ne Q}(h_P+h_Q).
\]
The product of fixed factors on the three letters determines the primitive triple if its numerical weight exceeds \(H_*\). We use this fixed-lock count throughout.

## The residual quantile

On each letter \(W\) prescribe a subset of the original root coordinates to lock. Call its limiting normalized root mass \(\mu_W\) and numerical weight \(w_W\), and write \(\mu_0=\sum_W\mu_W\), \(w_0=\sum_W w_W\). The chosen subsets are taken in the box before counting their actual root values. Let \(\sigma_{W,L}\) be the prelimit measure of the remaining coordinates on \(W\), and let \(\sigma_W\) be its weak limit after passing to a subsequence. Thus
\[
 m_W:=\sigma_W([0,\infty))=l_W-\mu_W,\qquad
 J_W:=\lim_L\int r\,d\sigma_{W,L}(r)=h_W-w_W .
\]
The residual escaped first moment is
\[
 e_W^{\rm res}=J_W-\int r\,d\sigma_W(r)\ge0.
\]
The residual moment includes just the numerical weight still available after the prescribed locks.

Choose any collection \(I\) of the two or three available letters (the formula also allows a single letter), and pool their remaining measures:
\[
 \sigma_L=\sum_{W\in I}\sigma_{W,L},\qquad
 \sigma=\sum_{W\in I}\sigma_W,\qquad
 M=\sum_{W\in I}m_W,\qquad e=\sum_{W\in I}e_W^{\rm res}.
\]
For \(0\le t\le M\), let \(F_\sigma(t)\) denote the largest numerical weight of mass \(t\) in \(\sigma\), allowing a fraction of a cutoff atom. This is the top integral of \(r\) under mass constraint \(t\). Its prelimit counterpart \(F_L(t)\) allows a fraction of the last original coordinate in the pooled list. For every fixed \(0<t<M\),
\[
 \lim_{L\to0} F_L(t)=F_\sigma(t)+e .
 \tag{QL1}
\]

To see this, use the upper-quantile identity
\[
 F_\sigma(t)=\inf_{y\ge0}
      \left(ty+\int(r-y)_+\,d\sigma(r)\right).
\]
Take all mass strictly above a quantile threshold \(y\) and enough of the mass at \(y\) to reach \(t\) for equality. The same identity holds for the prelimit \(\sigma_L\). If \(t\) is fixed positive, the infima can be restricted to a fixed compact \(y\)-range, since evaluating at \(y=0\) gives the bounded first moment, and the objective is at least \(ty\). Rewrite the integral as the full first moment minus \(\int\min(r,y)\,d\sigma_L\). The latter integrals converge uniformly on that compact range by weak convergence and their uniform Lipschitz bound in \(y\). The full first moment tends to \(\int r\,d\sigma+e\). This proves (QL1), including when the threshold has atoms. It also shows that \(F_\sigma\) is continuous up to the endpoints for a finite first moment. The term \(e\) enters for every positive \(t\), because the escaped weight has vanishing mass available in the top quantile.

## One-coordinate rounding and the resulting inequalities

Put
\[
 g_L=\max_{\substack{W\in I\\ i\ {\rm remaining}}}\beta_{W,i,L},
 \qquad g=\limsup_L g_L .
\]
Choosing whole coordinates from a prelimit quantile costs at most \(t+g_L\): rank all remaining root coordinates by their scaled degree and round up the last one. The numerical weight is at least \(F_L(t)\). One can instead use the actual overshoot in place of \(g_L\) at the particular \(t\). Ties may be ordered arbitrarily. In particular,
\[
 g\le \max_{W\in I}\sup_{r\ge0}\sigma_W(\{r\}).
 \tag{QL2}
\]
Indeed a sequence of individual coordinates with positive normalized mass has bounded scaled degrees by the first-moment bound. On a subsequence its scaled degree converges and its limiting mass contributes to an atom of the residual measure on that letter. Thus in (QL2) it is the largest atom on an individual letter, including an atom at zero, rather than the sum of atoms at the same scale across letters.

If \(0<t<M\) and \(\mu_0+t+g<d\), the three-factor count gives
\[
 w_0+F_\sigma(t)+e\le H_* .
 \tag{QL3}
\]
Otherwise (QL1) supplies a numerical weight exceeding \(H_*\) on a subsequence. Rounding up and adding the original locks would leave a total normalized freeze cost strictly below \(d\), a contradiction. By continuity in \(t\), (QL3) holds at \(t=\min(M,d-\mu_0-g)\) if this number is positive. Without any round-up condition, if \(\mu_0<d\), locking only a growing scaled tail of the residual pool costs \(o(1)\) and collects numerical weight tending to \(e\). Hence one also has \(w_0+e\le H_*\) in that case.

### Three letters, without initial locks

In a balanced \(1/3\) endpoint choose all three letters and no initial locks. If \(g\) bounds the largest atom on each letter and \(g<1/3\), then
\[
 F_{\nu_A+\nu_B+\nu_C}(1/3-g)+e_A+e_B+e_C
       \le 1+h_{\rm short},
 \tag{QL4}
\]
since the two largest numerical heights are one. The same formula with two letters pooled instead of three is valid, with just the corresponding escaped moments and largest atom. For a nonatomic profile the mass in this restriction is exactly \(1/3\).

### Two letters, after a tail lock on the third

Let \(T\) be the third letter and fix a positive continuity cutoff \(\theta\) for \(\nu_T\). Lock all \(T\) coordinates with scaled degree at least \(\theta\). Write
\[
 \mu=\nu_T([\theta,\infty)),\qquad
 R_T(\theta)=\int_{[0,\theta)}r\,d\nu_T(r).
\]
The full prelimit numerical weight of the lock tends to \(h_T-R_T(\theta)\), which includes \(e_T\). Pool the unaltered measures of the other two letters \(U,V\). If \(u=\min(h_U,h_V)\) and \(g\) bounds the largest atom of either \(\nu_U\) or \(\nu_V\), then
\[
 F_{\nu_U+\nu_V}(1/3-\mu-g)+e_U+e_V
       \le u+R_T(\theta)
 \tag{QL5}
\]
whenever the argument of \(F\) is positive. This uses \(H_*-h_T=u\) when there are two tall letters. No assumption about which letter is short is required.

## Two explicit diffuse sectors

Suppose the two tall measures are both uniform of mass \(1/3\) on \([2,4]\). Their pool has \(g=0\), no escaped weight, and
\[
 F(t)=4t-\frac32t^2,\qquad 0\le t\le2/3 .
\]
By (QL3) with no initial locks and only these two letters, \(F(1/3)=7/6\le1+h_{\rm short}\). Thus the numerical short height is at least \(1/6\), regardless of its degree profile.

There is a uniform finite-box version with an explicit saving. If the short numerical height is at most \(1/10\) and the two tall tails with scaled degrees at least \(61/20\) have combined normalized root mass at most \(8/25\) and numerical weight at least \(111/100\), the frozen factors determine the triple. Indeed \(111/100>1+1/10\), so the short–tall determinant vanishes after freezing. The height-box count has exponent at most \((8/25+o(1))L\), with the usual arbitrarily small truncation exponent losses. At the uniform tall model these two quantities are respectively \(19/60\) and \(893/800\), strictly inside the sector. One may use a fixed upper cutoff above the uniform support if a bounded-scale neighborhood is preferred. This yields a saving on an open neighborhood of the uniform tall pair.

At all numerical heights one, the same uniform pair as the opposite letters in (QL5) gives, whenever \(R_T(\theta)\le1/100\),
\[
 \mu\ge \frac{\sqrt{10-6R_T(\theta)}-3}{3}>\frac1{20}.
\]
For example, take a formal completely nonatomic \(T\) measure of mass \(22/75\) uniformly on \([0,1/50]\), and its remaining mass \(1/25\) uniformly on a short interval centered at \(3739/150\). Its first moment is one. For a cutoff \(\theta=10\), it has \(R_T=11/3750\) and \(\mu=1/25\), so the formal profile with the two uniform opposite measures is excluded even though all three letters are nonatomic. For a fixed-box realization near this formal profile, locking the two opposite tails from scaled degree \(63/20\) and the \(T\) tail from \(10\) has mass \(97/300\) and numerical weight \(2431/2400+3739/3750>2\) at its center. Thus the coefficient \(13/40\) is valid throughout a suitable neighborhood of that box sector.

## A positive scaled-degree cutoff with an opposite-atom allowance

At a balanced limiting abc endpoint, \(d=l_A=l_B=l_C=1/3\). For a letter \(T\), write \(h=h_T>0\) for its numerical height, scaled so the tallest height is one. Let \(U,V\) be the opposite letters, of heights \(u,1\) where \(0<u\le1\). Let \(g\) bound the limsup of the largest individual normalized prelimit root mass on either \(U\) or \(V\). One permitted choice is
\[
 g=\max_{W=U,V}\sup_{r\ge0}\nu_W(\{r\}).
\]
Then for every \(0<\theta\le h/100\),
\[
 \boxed{\quad \nu_T([\theta,\infty))+g>\frac{h}{178}-\frac{\theta}{9}.\quad}
 \tag{PC1}
\]
In particular
\[
 \nu_T([h/100,\infty))+g>\frac{361}{80100}h>\frac h{225}.
\]
The estimate allows a weak limiting atom at \(\theta\) and includes escaped first moments.

Use the fixed-factor count and exact quantitative determinant parameters from [An explicit positive-degree mass in terms of numerical height](ABC_QUANTITATIVE_DEGREE_WINDOWS.md#an-explicit-positive-degree-mass-in-terms-of-numerical-height). For the quantile rounding convention see [Upper quantiles after fixed root locks](research/arithmetic/#one-coordinate-rounding-and-the-resulting-inequalities).

### Prelimit tails, round-up, and median dispersion

Put \(\mu^*=\nu_T([\theta,\infty))\) and \(s=\mu^*+g\), and suppose \(s\le h/178-\theta/9\). The prelimit root mass of the \(T\) coordinates with scaled degree at least \(\theta\) has limsup at most \(\mu^*\) by weak convergence on a closed set. Their complementary numerical weight is at most \(\theta(1/3+o(1))\), since all the omitted coordinates have scaled degree less than \(\theta\). Thus this lock has limiting cost at most \(\mu^*\) and numerical weight at least \(h-\theta/3\). These one-sided bounds remain valid when roots approach an atom at \(\theta\) from either side. Similarly, a prelimit opposite root coordinate whose mass is bounded below has bounded scaled degree and contributes a weak limiting atom along a subsequence. Its limsup mass is therefore bounded by the largest atom of its own letter, justifying the displayed choice of \(g\).

Let \(F_L(t)\) be the prelimit top quantile integral on \(U,V\), allowing a fraction of the last coordinate, and put \(t_0=1/3-s\). For any fixed \(t<t_0\) with \(t>0\), freeze those top coordinates, rounding the last one up. This costs at most \(t+g+o(1)\). If \(\limsup F_L(t)>u+\theta/3\), combine with the \(T\) tail lock above. The combined weight exceeds \(h+u\), so the determinant for the pair \(T,U\) vanishes by the three-factor count, whereas the total root mass is strictly below \(1/3\). This contradicts the endpoint. Therefore
\[
 \limsup F_L(t_0)\le u+\theta/3.
\]
Here \(t\uparrow t_0\) is allowed since top-quantile slopes on a positive mass range are uniformly bounded by the prelimit first moment divided by the left mass.

Write \(M_L\to2/3\) for the pooled free-letter mass and \(J_L\to1+u\) for its full prelimit first moment. By the decreasing top-quantile average,
\[
 \limsup F_L(M_L/2)\le\frac{u+\theta/3}{1-3s}.
\]
Choose a pooled prelimit median and pass to a subsequence on which it tends to \(a\). It is bounded, and it is bounded away from zero because the average degree of the bottom half is at least \(3(1-3h/175)\) in the limit, as follows from the bound below. Group every opposite root at common unweighted base exponent \(q=\lfloor a/L\rfloor\). Its narrow ratio \(z\) has normalized height at most \(1/3\). The absolute median dispersion is \(2F_L(M_L/2)-J_L\), and the imbalance between the two signed correction weights is \(1-u+o(1)\), since the opposite root masses both tend to \(1/3\). Thus the wider factor has exponent
\[
 E\le \frac{u+\theta/3}{1-3s}-u
    \le\frac{3s+\theta/3}{1-3s}
    \le\frac{3h}{175}.
\]
The same bound implies the asserted bottom-half moment bound because \(F_L(M_L/2)\le u+3h/175+o(1)\). This argument used the full prelimit first moments, including any escape.

### The full modulus and exact margins

On \(T\) now freeze the coordinates with scaled degrees at least \(\kappa=(3/2)h\), and take \(k=\lfloor\kappa/L\rfloor\). Their root mass \(\mu\) is at most \(\mu^*\). The weight below \(\kappa\) is at most \(\theta/3+\kappa\mu^*\): below \(\theta\) use the direct bound above, and between \(\theta\) and \(\kappa\) use the total prelimit tail mass. Hence the locked weight \(w\) satisfies
\[
 \begin{split}
 w&\ge h-\theta/3-\kappa s\\
  &\ge h\left(1-\frac{3h}{356}
             -\frac{\theta}{h}\left(\frac13-\frac h6\right)\right)
   \ge \frac{52861}{53400}h .
 \end{split}
\]
For the last inequality use \(\theta/h\le1/100\), \(h\le1\), and the positivity of \(3/356-1/600\).

The complete shallow record capacity is at least \(w/(2\kappa)>161/500\). Choose a selective shallow record with normalized weight \(b\ge161/500\), adding one copy of an original root at a time until this threshold. Its mass is at most \(161/500+\mu\). The prime-power depth condition holds because the allowed number of copies of each root of degree \(e\) is \(\lfloor e/k\rfloor\).

Use the arbitrary-rank, first-coordinate reconstruction version of the one-root count with
\[
 \alpha=\frac{543}{1000},\qquad M=6,\qquad r=\frac{37}{100}h .
\]
Then
\[
 B\le\frac{19}{84}h,\qquad
 \frac{b\alpha(M+1)r}{2}
   \ge\frac{22642557}{10^8}h>\frac{19}{84}h,
 \qquad
 \alpha(M+1)r=\frac{140637}{10^5}h<\kappa.
\]
The median gives \(r<a\). For terminal reconstruction,
\[
 \frac{2B}{1-\alpha}\le\frac{9500}{9597}h
       <\frac{52861}{53400}h\le w .
\]
The last comparison is exact, since \(52861\cdot9597-9500\cdot53400=7017>0\). Thus the retained modulus determines \(z\), and the original full modulus then determines \(v\). The chart cost is the actual record mass \(b\). Restoring the lock gives
\[
 d\le \mu+b\le\frac{161}{500}+\frac{2}{178}<\frac13 ,
\]
which proves (PC1).

## Packing selected atoms and the nonatomic part

One can restrict the eligible coordinates in a quantile before
rounding. Ignoring a coordinate leaves its value free, so it has
no freeze cost. This permits a separate choice of the atomic
factors to lock. Write \(\nu_{W,\mathrm{na}}\) for the nonatomic
part on letter \(W\). For a selected set of atoms across the
letters, let \(p\) be their total root mass and \(w\) their
total numerical weight. If \(p<d\), put
\[
 \sigma=\sum_W\nu_{W,\mathrm{na}},\qquad
 t=\min\{\sigma([0,\infty)),d-p\}.
\]
Then, whenever \(t>0\),
\[
 w+F_\sigma(t)+\sum_W e_W\le H_* .
 \tag{QL6}
\]
In particular, at the minimum nonatomic endpoint this uses
\(d=1/3,\ H_*=2\) and no escaped first moments. The quantile
is taken only over the nonatomic parts, including when the
selected atoms have very different masses on the three letters.

For a finite selected set, freeze its shrinking atomic windows.
Keep a finite list containing the selected atoms and omit the
unselected windows in that list from the eligible pool. For each
fixed list, choose the shrinking windows slowly enough to
recover their full atomic masses. Discarding these windows from
the pool leaves the escaped first moments there, since their
degrees stay bounded for each fixed list. Enlarging the list
through all atoms makes the largest remaining atom tend to zero.
The remaining weak measures exceed their nonatomic parts by
positive atom measures whose first moments tend to zero.
Adding a measure to a quantile pool can increase the top
weight at a fixed mass by at most the added first moment.
Thus the top weight tends to \(F_\sigma\) at every mass up to
the nonatomic mass. Apply (QL1) to this pool, and round as
in (QL2)–(QL3). Any positive mass strictly below \(t\) has
strict freeze margin for a sufficiently large list, proving
(QL6) by continuity. For an infinite selected set apply the
finite statement along increasing finite subsets and use
convergence of the selected weights. If the limiting eligible mass
vanishes, a slowly growing cutoff on the prelimit eligible tail
collects its escaped weight at root cost \(o(1)\). This gives
\(w+\sum_W e_W\le H_*\) when \(t=0\).

### Tracking individual coordinates at a limiting atom

A weak atom may comprise several large coordinates or many
small ones. To use the actual rounding granularity, rank the
original coordinates on each letter by decreasing normalized
root mass before taking the small-budget limit. Pass to a
diagonal subsequence on which the mass at each rank converges.
Whenever its limit \(\beta_{W,j}\) is positive, the first-moment
bound keeps that coordinate's scaled degree bounded, so require
that degree to tend to \(r_{W,j}\) as well. The measures
\[
 \lambda_W=\sum_{j:\beta_{W,j}>0}\beta_{W,j}\delta_{r_{W,j}},
 \qquad \omega_W=\nu_W-\lambda_W
\]
are nonnegative, since every finite partial sum of \(\lambda_W\)
is dominated by \(\nu_W\) on nonnegative continuous tests.

Select any subcollection of the positive limiting coordinates,
with total mass \(p<d\) and weight \(w\). Then (QL6) holds
with \(\sigma=\sum_W\omega_W\). To see this, fix a finite
selected subcollection and take \(K\) large enough to contain its
ranks. Omit from the eligible pool all other coordinates among the
first \(K\) positive limiting ranks on each letter. Do not remove any rank whose limiting mass is
zero. Its scaled degree could tend to infinity and carry escaped
weight. The removed positive ranks have bounded scaled degrees
for fixed \(K\), so the eligible escaped moment is unchanged.
As \(K\) grows, the largest remaining coordinate mass tends
to zero, while the difference of the eligible weak measure
from \(\sum_W\omega_W\) has first moment tending to zero.
The same quantile and continuity argument proves the claim,
and finite approximations handle an infinite selected
subcollection. In this form, mass contributed to a weak atom
by coordinates whose individual masses vanish remains available
for fractional quantile selection.

### A diffuse example that requires skipping an intervening atom

Here is a formal \(\Gamma=1/2\) profile with all numerical heights one. Put \(c=1/6\).
- \(U\): an atom of mass \(c\) at \(5\), mass \(1/50\) uniformly on \([599/100,601/100]\), and mass \(11/75\) uniformly on \([0,7/11]\).
- \(V\): an atom of mass \(c\) at \(0\), and mass \(c\) uniformly on \([4,8]\).
- \(T\): an atom of mass \(c\) at \(59/10\), and mass \(c\) uniformly on \([0,1/5]\).

Each total mass is \(1/3\) and each first moment is one. For \(U\), the atomic, high diffuse and low diffuse weights are \(5/6,3/25,7/150\), respectively. For \(V\) only its diffuse part carries numerical weight. For \(T\) the atomic and diffuse weights are \(59/60,1/60\).

Lock the \(T\) atom, the high diffuse \(U\) band, and the diffuse \(V\) tail above \(91/20\). Their combined mass is
\[
 \frac16+\frac1{50}+\frac{23}{160}=\frac{793}{2400}<\frac13,
\]
and their numerical weight is
\[
 \frac{59}{60}+\frac3{25}+\frac{17319}{19200}
 =2+\frac{103}{19200}>2 .
\]
This is an explicit use of (QL6), with the \(U\) atom at \(5>91/20\) excluded from the quantile pool. An open neighborhood of these fixed box windows therefore has normalized count exponent at most \(331/1000\), by the strict fixed-factor uniqueness margin. 
For comparison, take the original representative of each listed atom to be a single macroscopic coordinate, so the prelimit round-up cap in (QL4) and (QL5) is \(g=c\). The three-letter unstructured quantile (QL4) then tests mass \(c\) in a pool supported in \([0,8]\), with weight at most \(8/6<2\). The third-letter tail variants (QL5) also hold.
- For a lock on \(T\) up to its atom scale, or on \(U\) up to its atom scale, the locked mass is at least \(c\), so the residual quantile argument is nonpositive. Above the \(T\) atom there is no tail and the residual numerical weight is one. Above \(5\) on \(U\) the residual numerical weight is at least \(22/25\), whereas the opposite top integral is at most \(8/6\). Thus the positive-budget inequalities on these two letters hold at every permissible cutoff.
- For a cutoff on \(V\) at most \(4\), its entire diffuse part is locked, so again the quantile argument is nonpositive. For a cutoff in \((4,8]\), the available opposite quantile mass \(t\) is precisely the mass of the diffuse \(V\)-part below the cutoff. Its residual weight is at least \(4t\). The opposite pool is supported in \([0,601/100]\), hence its top integral is at most \((601/100)t\le1+4t\) for \(t\le c\). Above \(8\) the residual weight is one and the bound is immediate.

Even using the actual overshoot after locking the \(T\) atom, ordinary greedy selection from the whole opposite pool encounters the \(U\) atom at \(5\). Stopping before it collects only the \(V\) tail above \(5\) and the narrow high diffuse \(U\) band, whose combined weight is \(39/48+3/25<1\). Rounding up the atom exceeds the root budget, since positive mass was already used before reaching it and this atom alone has the full remaining budget \(c\). Selecting the eligible diffuse factors fills the available mass below the atom's scale.
