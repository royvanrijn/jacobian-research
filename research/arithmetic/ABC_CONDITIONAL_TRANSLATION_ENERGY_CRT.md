# Quantitative additive-energy exclusion in fixed unit classes

This note gives a quantitative scalar-energy test in the actual residue cells, and its application at the smaller conditional dictionary scale. Put \(n=X^L\) and always take the height limit with \(L\) fixed. All fixed margins below precede the small-budget threshold.

## TE1–TE2. A uniform one-dictionary lemma

Fix positive constants \(\nu,H_0,r,v,s,\beta_0\) and \(w\ge0\). Let \(\mathcal S\) be a set of \(N\) positive integers \(z\le X^{H_0}\) with \(\operatorname{rad}(z)\le n^r\), all in one **unit** residue class \(a\bmod D\), with \(X^\nu\le D\le X^{H_0}\). Assume their common containing real interval has length at most \(D X^w\). For each \(z\in\mathcal S\), require that the full prime powers \(p^{v_p(z)}\) with \(Lv_p(z)\ge s\) and \(\log_n p\le \omega\) carry total normalized radical mass at least \(\beta_0\). Thus their squarefree support has logarithm at least \(\beta_0\log n\). One may require the powers to lie in a specified compact valuation band instead. Suppose
\[
                           s v>w,\qquad \beta_0>r-2v.       \tag{TE1}
\]
Then there are constants \(\omega_0,e,L_0>0\) in terms of these fixed bounds and margins such that, for each fixed \(0<L<L_0\), all sufficiently large \(X\), and \(N\ge n^v\), uniformly in the residue class, modulus, interval and dictionary, if \(0<\omega\le\omega_0\) then
\[
                            E_+(\mathcal S)<N^3n^{-e}.      \tag{TE2}
\]
The same statement for \(N\ge n^{v-\tau}\) holds for sufficiently small fixed \(\tau>0\), decreasing \(v\) using the strict margins. The bound may count a subdictionary whose elements arise on genuine rows with any additional eligibility conditions.

### Proof

Choose \(0<b<\beta_0\) with \(v>(r-b)/2\), and choose \(t\) with \(w/s<t<v\). Fix an integer \(m\ge2\) with \(t/m<\beta_0-b\). Fix \(\omega_0\) small enough that \(t+m\omega_0<v\). Choose \(e>0\) so small that \(e<v/2\) and
\[
            v-me>\max\{(r-b)/2,\ t+m\omega_0\}.
\]
All inequalities have positive fixed margins. Fix sufficiently small positive normalized counting allowances inside those margins.

Suppose \(E_+(\mathcal S)\ge N^3 n^{-e}\). Write \(r_h=\#\{z\in\mathcal S:z+h\in\mathcal S\}\). Omit the zero shift, whose energy contribution is \(N^2\). Shifts with \(r_h<N n^{-e}/4\) contribute at most \(N^3 n^{-e}/4\), since \(\sum_h r_h=N^2\). Thus the surviving nonzero shifts \(\mathcal H\) have
\[
             |\mathcal H|\gg Nn^{-e},\qquad
             r_h\ge Nn^{-e}/4\quad(h\in\mathcal H).
\]
The first bound follows because each shift contributes at most \(N^2\) to the energy. Use the bipartite graph joining \(z\in\mathcal S\) to \(h\in\mathcal H\) if \(z+h\in\mathcal S\). The mean \(z\)-degree is at least \(|\mathcal H|n^{-e}/4\), which tends to infinity since \(e<v/2\). The power mean, after deleting degrees below half the mean, counts \(\gg_m N(|\mathcal H|n^{-e})^m\) ordered \(m\)-stars with distinct neighbors. Some fixed ordered list of distinct \(h_1,\ldots,h_m\) therefore has \(\gg_m Nn^{-me}\) common bases \(z\in\mathcal S\). All shifts are multiples of \(D\), and their magnitudes are at most \(X^{H_0}\).

Take the prescribed prime pool \(P\) to support \(D\prod_{i<j}|h_i-h_j|\). It has uniform polynomial height for the fixed list size. For each position \(i\), the number of common bases with at least \(b\) normalized radical mass of \(z+h_i\) on this pool has bound
\[
                         O(n^{(r-b)/2+\eta} X^\varepsilon)
\]
for arbitrary fixed normalized allowance \(\eta>0\) at sufficiently small fixed \(L\), with \(\varepsilon>0\) subsequently arbitrarily small for that \(L\). Indeed fix the exact full supported part \(Q\) of \(z+h_i\) at height-subpower cost. It is a unit modulo \(D\). The quotient \((z+h_i)/Q\) has outside radical at most \(n^{r-b}\), height at most \(X^{H_0}\) and belongs to one fixed unit class modulo \(D\). Apply the absolute residue bound (RB2) with a list of size one, uniform in \(Q\), and then use injectivity of \(z\mapsto z+h_i\). The power-height specified prime pool count and the residue estimate both remain uniform when \(P,D,a\) depend on the previously chosen shifts. The finitely many deletions leave \(\gg_m N n^{-me}\) bases.

On each surviving base, for every \(i\) the marked band on \(z+h_i\) has at least \(\beta_0-b\) radical mass outside \(P\). Take full powers off \(P\) by a greedy selection of normalized prime cost between \(t/m\) and \(t/m+\omega_0\). Denote their product by \(F_i\). The products are pairwise coprime, because a shared prime would divide \(h_i-h_j\), and all are coprime to \(D\). They satisfy \(\prod_i F_i\ge X^{st}>X^w\). Their joint squarefree support has at most \(n^{t+m\omega_0}\) choices. On each support, label the primes by their positions \(i\) and enumerate their full powers at height-subpower cost, uniformly since each \(F_i\le X^{H_0}\). Thus there are at most \(n^{t+m\omega_0}X^{o_X(1)}\) packet records. Each record prescribes \(z=a\bmod D\) and \(z=-h_i\bmod F_i\) for every \(i\). The common modulus is \(D\prod_iF_i>D X^w\), so at most one \(z\) in the original interval is possible. This contradicts the surviving count by the choice of \(e\). All constants or fixed-support losses are absorbed by the original strict margins, with no exchange of the limits in \(X\) and \(L\).

[Absolute residue count RB2](ABC_BALANCED_DICTIONARY_RESIDUE_CELLS.md#regular-balanced-cohorts-and-absolute-residue-pricing)
· [Prescribed support enumeration](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

## TE3–TE4. All three original caps and the conditional scale

Suppose a common three-register cell on genuine primitive original rows \(A+B=C\) has fixed pairwise coprime full locks \(D_A,D_B,D_C\) and fixed signed register values on those moduli, in the sense of RB6–RB7 and RT6. Assume that the row count is at least \(n^{2v-\tau}\) and each of the three complete dictionaries is at most \(n^{v+\tau}\). Assume an original role \(T\) has numerical upper height \(X^{h_T}\). Let \(\gamma_W\) be fixed lower bounds for \(\log_X D_W\), and suppose that the product of the two locks not belonging to \(T\) has a fixed positive numerical floor. Apply the scalar lemma to the dictionary \(z=T/D_T\). The two opposite signed registers and the additive relation prescribe a complete residue for \(T\) modulo \(D_{T'}D_{T''}\). This residue is a unit by primitivity. The fixed coefficient \(D_T\) is a unit there, so \(z\) is in one unit class with \(D=D_{T'}D_{T''}\). Its containing interval has length at most \(DX^{w}\) with
\[
                  w=\max(0,h_T-\gamma_A-\gamma_B-\gamma_C).
\]
Any improved numerical interval length for \(T\) can replace the height \(h_T\). A negative actual leftover exponent gives at most one dictionary element directly. Let \(r_T\) bound the **remaining radical** of \(T/D_T\), and let \(\beta_0\) be as in the scalar lemma on the actual full powers of that quotient. The other two original radical caps and all eligibility conditions are retained. In particular, all input rows and the local density assertion still come from the balanced three-cap graph.

Indeed, every signed pair of complete labels determines its edge. The size \(N_T\) of each original dictionary lies between \(n^{v-2\tau}\) and \(n^{v+\tau}\). Cauchy--Schwarz on the marked-role output of the actual additive rows gives mixed energy at least \(n^{3v-3\tau}\) on the other two dictionaries. The mixed energy is bounded above by the geometric mean of their individual energies, each at most \(n^{3v+3\tau}\). Switching roles supplies for whichever original role \(T\) is chosen
\[
             E_+(\mathcal D_T)\ge n^{3v-9\tau}
                              \ge N_T^3 n^{-12\tau}.
\]
Passing from its values to \(T/D_T\) preserves exact additive energy. Consequently the local cell is impossible at small enough tolerances and budget whenever
\[
    \boxed{\quad s v>\max(0,h_T-\gamma_A-\gamma_B-\gamma_C),
                         \qquad \beta_0>r_T-2v .\quad}       \tag{TE3}
\]
There is a fixed positive separation from the full local exponent when both inequalities and the positive opposite-lock depth have fixed margins. To make the order of choices explicit, reserve half the margins to replace \(v\) in TE1 by \(v-\xi\) with fixed small \(\xi>0\), get \(e\) from the scalar lemma, and then choose \(2\tau<\xi\), \(12\tau<e\) and \(L\) sufficiently small. Because RT6 has near-square density, a suitable non-atomic band in the original **post-own-lock** profile may be used whenever compactness supplies the uniform lower marked mass and the small individual prime condition in that band. This step concerns actual prime valuations. A raw coordinate cap without uniformly occupied marked mass is insufficient.

For a single RB4 cell on a lock \(D_W\), the free dictionaries have local exponent \(u=1/3-\alpha\), the edge exponent is \(2u\), and the marked dictionary may have exponent as large as \(\alpha\). The mixed energy may be used directly, avoiding the loss of passing to individual energy by Cauchy--Schwarz. Write \(\rho=\alpha-u\). At tolerance \(\tau\), the true sums give mixed energy at least \(n^{3u-\rho-3\tau}\). Each dictionary has size \(N_T\) between \(n^{u-2\tau}\) and \(n^{u+\tau}\). Assume \(u>\rho\) with a fixed margin, otherwise this lower bound need not beat the diagonal. Drop the zero shift. Drop also shifts with either representation count below \(n^{u-\rho-5\tau}/4\), using \(\sum_h r_T(h)=N_T^2\) for the other dictionary. The remaining shifts have cardinality at least a constant times \(n^{u-\rho-5\tau}\), since the product of two representation counts is at most \(n^{2u+2\tau}\). On either dictionary their mean degree in the translate graph is at least a constant times \(|\mathcal H|n^{-\rho-6\tau}\). This tends to infinity if \(u>2\rho\) with a fixed margin. It suffices to ensure this here: the required criterion below implies \(u>2\rho\) once \(m\ge2\), because its right side is nonnegative. The normalized \(m\)-star size on the chosen role is hence at least \(u-m\rho-(6m+2)\tau\). The same prescribed-pool deletion and CRT argument apply. One directly excludes such a cell if, for the packet parameters in the scalar proof,
\[
 \boxed{\quad
 u-m(\alpha-u)>\max\{(r-b)/2,\ t+m\omega_0\},\qquad st>w,
 \qquad t/m<\beta_0-b .\quad}       \tag{TE4}
\]
Here \(r\) is the counted free letter's radical cap, \(0<b<\beta_0\), and \(m\ge2\) is fixed. This uses the third-role radical cap through the balanced cell and genuine target dictionary. No switch to three registers is needed. For comparison, passing through individual energy alone would give \(E_+(\mathcal D_U)\ge N_U^3n^{-4(\alpha-1/6)-12\tau}\) and would lose twice as much of the \(\alpha-u\) margin in the star size.

## Scope and comparison with compact packing and translation cubes

At the exact minimum scale \(v=1/6\), with \(r=1/3\)
and unit-width parameter exponent \(w=1-\gamma\), the
condition is \(s>6(1-\gamma)\) with positive high diffuse
mass. That cutoff is already in the filler-only branch
of MPK1, including its all-role refinement MPK2.
Those results also treat mixed atomic and diffuse
packing. TC1 and TC5 give simultaneous original-row
translations and fresh support prices in fixed fibres.
[Compact packet packing MPK1–MPK2](ABC_CROSS_ROLE_COMPACT_PACKING.md)

TE2 is a scalar energy inequality for any actual fixed
unit residue dictionary. The pruning of the pool of
shift differences supplies the exact conditional slack
\(r_T-2v\) in TE3. That application uses the cap after
stripping the role's own full divisor and the numerical
depths of all three locks. The profile inequalities and
positive conditional dictionary exponent still have to
be supplied by the actual cell. At occurrence
\(\alpha=1/3\), a residual exponent zero provides no
power margin for this argument.
