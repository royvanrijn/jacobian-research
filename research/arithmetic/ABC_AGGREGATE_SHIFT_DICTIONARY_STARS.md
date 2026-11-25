# Common-row CRT bounds with a varying shift dictionary

Use boxes with \(n=X^L\), taking the fixed-budget height limit before the small-budget limit. Let \(U,T,K\) be primitive complete integer letters, with \(U,T>0\), fixed signs, and
\[
 U=Fxy,\qquad K=Gk,\qquad T=\epsilon Fxy+\eta Gk .
\]
Take \(k>0\) and absorb its sign into \(\eta\). The positive full factors \(F,G\) are fixed, possibly unit. The original residual roots \(U/F\) have normalized root cap \(a\). The row block \(x\) uses cap \(p\) and numerical height at most \(S_xp+o(1)\), and the remaining column factor \(y\) uses cap \(a-p\) and numerical height \(H_y+o(1)\). Use the actual height cap of \(y\), which may be smaller than that of the original residual \(U/F\). The available \(k\)-dictionary has size at most \(n^{b+o(1)}\). It need not have a prescribed factorization. Fix a choice of original root representation per primitive triple. Suppose the fibre has normalized count exponent \(r=a+b-\delta\).

On \(T\) choose a compact positive degree band of lower endpoint \(s\) with enough negligible-size root coordinates to select a fixed block of root cap \(\beta\) and numerical weight at least \(s\beta-o(1)\). This applies if the band has nonatomic mass greater than \(\beta\). One row may include a disjoint tail block of root cap \(\zeta\) and numerical weight at least \(e-o(1)\). Take \(e=\zeta=0\) if no tail is used. Let \(m\ge2\) be fixed, \(q=m\beta+\zeta\), and assume \(p>\delta\).

The graph on row \(x\) and column \((y,k)\) is simple after prescribing the representation; \(x,y,k\) determine the triple. Its two part caps are \(p\) and \(C=a-p+b\). Discard columns of degree below \(n^{p-\delta-\varepsilon}\) with \(0<\varepsilon<p-\delta\), losing at most \(n^{r-\varepsilon+o(1)}\) edges. On the surviving columns the distinct ordered \(m\)-stars number \((1-o(1))\sum d_j^m\), since every degree tends to infinity. By the power mean this is at least \(n^{mr-(m-1)C+o(1)}\). Divide by \(n^{mp+o(1)}\), the row-tuple cap, to get fixed distinct \(x_1,\ldots,x_m\) with at least \(n^{C-m\delta+o(1)}\) common columns.

On row \(i\), record the selected target block (also the tail on row one if used), whose corresponding complete powers form \(M_i\mid T_i\). The joint roots cost \(q\). For every shared prime and each full multiplicity,
\[
 \gcd(M_i,M_j)\mid x_i-x_j ,
\]
because \(T_i-T_j=\epsilon Fy(x_i-x_j)\) and \(Fy\) is a unit at every prime dividing the gcd by primitivity. This gives a deterministic lower numerical weight
\[
 W=e+s m\beta-\binom m2 S_xp
\]
for the lcm \(M\) of the recorded locks, up to \(o(1)\). Indeed the loss of sum of prime exponents upon taking their maximum is at most the sum of all pairwise minima.

There are three ways to count the columns for each record.

- If the \(k\)-dictionary lies in an interval of numerical width \(X^{h+o(1)}\), and \(W>h\), each fixed \(y\) allows at most one \(k\). Congruences for any two possible \(k\)'s give \(M_i\mid G(k-k')\) in every row, and \(G\) is a unit at those primes. Thus there are at most \(n^{a-p+o(1)}\) columns, which contradicts the star exponent if \(q<b-m\delta\). More generally if the dictionary is covered by \(n^{\theta+o(1)}\) such intervals, replace this last condition by \(q+\theta<b-m\delta\). One can likewise use any uniform bound \(n^{\theta+o(1)}\) on its residue occupancy at moduli of weight at least \(W-o(1)\).
- If \(W>H_y\), each fixed \(k\) allows at most one \(y\), since \(Fx_i\) is a unit modulo \(M_i\). There are at most \(n^{b+o(1)}\) columns per record, so obtain a contradiction if \(q<a-p-m\delta\). When both single-coordinate conditions hold the smaller per-record cap applies.
- If \(k\) has numerical height at most \(h_k\) and \(W>H_y+h_k\), for any two columns the determinant \(k_1y_2-k_2y_1\) is divisible by every \(M_i\). For example multiply the two target congruences by the opposite \(y\) values and subtract; the remaining coefficient \(G\) is a unit. Thus the determinant vanishes, and \(\gcd(k,y)=1\) on each primitive edge makes the coprime positive pairs equal. There is one column per record, which contradicts the star exponent if \(q<C-m\delta\).

The last criterion uses the actual absolute height of \(k\), whereas the first uses only the width of a containing interval or residue occupancy. All estimates are uniform over the fixed records and shifts, and apply with arbitrary small fixed margins before the inner limit.

## Removing the overlap penalty for certain row blocks

Suppose the \(x\)-block includes an original indexed root \(z\) whose normalized cap is \(p_1>\delta\). For each fixed column, fixed other row \(x_j\), and fixed \(\xi>0\), offending rows with \(\log_X\gcd(M_i,M_j)\ge\xi\) have a full gcd divisor \(d\mid M_j\) of numerical weight at least \(\xi\), and \(x_i\equiv x_j\bmod d\). Fix the other roots of \(x_i\), at cost \(p-p_1\), and enumerate divisor-many \(d\) at fixed truncation. The equation on \(z\) is \(z^\ell\equiv A\bmod d\) with \(A\) a unit and \(\ell\) the original index. On each prime power the number of unit roots is at most \(2\ell\) (cyclic odd-primary unit group and at most two cyclic factors at 2 with the smaller factor order two). Since \(\ell\) is fixed before the inner limit, all unit classes number \(X^{o_X(1)}\) uniformly. For small \(L\), the original \(z\)-range is shorter than \(d\), so each class gives at most one \(z\). Thus offending rows number at most \(n^{p-p_1+o(1)}\) for a fixed other row and column. This is negligible relative to the regularized degree \(n^{p-\delta-\varepsilon}\), with \(\varepsilon<p_1-\delta\). Fixing the other \(m-1\) row slots in a star and summing over pairs shows that good stars give the same lower bound. Each pairwise overlap now has numerical weight below arbitrary \(\xi\), so in the three criteria one may use \(W=e+s m\beta\) with a strict margin. This includes an exact-power row block using a single root.

There is a diffuse variant. Suppose the selected \(T\) band is of bounded scaled degrees \([s,S]\) and its individual prime powers all have numerical weights tending to zero, as holds for whole coordinates chosen from a nonatomic band. At most one tail is used. For any pairwise gcd, one side comes from the bounded band, so each prime contribution to the gcd is \(o(1)\), even for a pair meeting the tail. Require the optional tail to lie above the band. All prime powers in the gcd then have scaled orders at least \(s-o(1)\). Choose \(k_0 L\to\kappa=s/2\), and fix \(0<\varepsilon_1<\kappa/16\). An offending gcd of numerical weight at least \(\varepsilon_1p\) yields a \(k_0\)-full sub-divisor \(d\mid M_j\) of weight \(\varepsilon_1p+o(1)\), by adding its full prime powers until that threshold. Again there are divisor-many \(d\).

Suppose the row block of cap \(p\) lies in a compact positive band with width \(D\). Group it as \(x=z^\ell v\), with \(z\le n^{p+o(1)}\), \(v\le X^{Dp+o(1)}\), and \(\ell\) at the lower edge. Use (UW) from ABC_QUANTITATIVE_DEGREE_WINDOWS.md with dimension one, narrow cap constant \(C_0=\max(1,2\kappa/\varepsilon_1)\), and full modulus \(d\). The wider exponent relative to \(d\) is \(\le (D+o(1))/\varepsilon_1\), so \(D<c_1\varepsilon_1/(2C_0)\) suffices, where \(c_1\) is the dimension-one constant from (UW). Take one narrow/wider tuple per exact \(x\) and hence exact-fibre bound one. The number of offending \(x_i\)'s is at most \(n^{2\varepsilon_1p/\kappa+o(1)}\), which is at most \(n^{p/4+o(1)}\). If \(\delta<3p/4\), this is negligible per regularized column when the regularization loss is small enough. Thus the three criteria hold with \(W=e+s m\beta-\binom m2\varepsilon_1p\). The band width is chosen after \(\varepsilon_1\), independently of \(p\). The auxiliary divisor \(d\) is used only to count bad stars and is not recorded on the retained stars.

## Numerical separation from the two elementary fixed-parameter methods

Take \(a=3/20\), \(b=1/10\), height \(3/5\) on the varied residual \(U\) letter with degrees uniformly on \([39/10,41/10]\), and arbitrary shift dictionary of height \(h=1/4\). Give \(T\) mass \(1/200\) uniformly on \([18/5,37/10]\) and mass \(2327/9400\) uniformly on \([23/10,12/5]\). The target first moment is \(3/5\), consistent with \(T=U\pm K\) and the shorter shifts. The target mass need not equal \(a\) or \(b\) for this comparison.

In any fixed-shift reconstruction from the high target band alone, using all affordable cap \(a\) across arbitrary many rows gives nominal maximum weight at most \((37/10)a<3/5\), so even before overlap loss the corresponding (FS1) condition fails when the row block is small. Even with the reduced column height and any choice of \(p<a\), the residual \(U\) coordinates all have scaled degrees at least \(39/10\), so their height is at least \((39/10)(a-p)\). Any viable nominal fixed-row packing costs at most \(a-p\) and weighs at most \((37/10)(a-p)\). Conversely, if \(U\) is fixed, no single target block with cap \(b\) can reach \(h\), since its weight is at most \( (1/200)(37/10)+(b-1/200)(12/5)=493/2000<1/4\).

In the aggregate first-coordinate criterion take \(m=20,\ \beta=1/250,\ p=1/10^6,\ \delta=1/10^9,\ \zeta=e=0\). Then \(q=.08<b-m\delta\), \(p>\delta\), and the deterministic weight is at least
\[
     (18/5)(.08)-190(41/10)10^{-6}=.287221 >.25.
\]
The aggregate saving is therefore strict below \(a+b-\delta=1/4-10^{-9}\), with no arithmetic restriction on the shift dictionary beyond size and height. This formal sector compares the indicated CRT profiles and does not assert that such primitive fibres exist at the unsaved count.

## Scope at the balanced endpoint

The first reconstruction requires density near the product of the
row-root capacity and the shift dictionary capacity. In an
unsaved affordable pair record, if the residual letter cap is
\(a=1/3-f\) and the other recorded block costs \(g\), then the
unsaved fibre exponent is \(r=a-g\). For a shift dictionary cap
\(b\) this sets \(\delta=b+g\). The interval criterion
\(q<b-m\delta\) is then vacuous. The other reconstructions have
different thresholds and must be checked individually. In
particular, the estimate applies to genuinely near-product
multiplicative splits, and an unrestricted use would need an
additional arithmetic argument connecting such splits to the
surviving affordable records.
