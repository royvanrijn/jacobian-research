# No graded full Long jet at degree 25, with arbitrary parity

Over every characteristic-zero field, the complete degree-25 graded
Long family fails the commutator equations through \(h^5\). This
includes the odd directions, both quadratic-field branches, every
admitted first-image deformation, and every second and third
homogeneous freedom.

No nonlinear source or target normalization is presumed. Only the
uniform linear source \(XQ\) and target \(RS\) shears are used.
The proof retains four entire field strata, including the zero locus
of the parameter inverted on one stratum.

## Statement and complete coefficient spaces

Use source order \((X,Q,P,Z)\), brackets
\(\{P,X\}=\{Z,Q\}=1\), standard Moyal ordering, and weights
\[
 \operatorname{wt}(X,Q,P,Z,h)=(-1,1,4,2,3),\qquad
 \operatorname{wt}(R,T,D,S)=(-1,2,4,1).
\]
Let \(f=(R,T,D,S)\) be the full canonical Long polynomial symbols from
[the exact transcription](../dixmier-audit/verify_long_pbw.py).

**Theorem.** There is no scalar full Weyl jet
\[
 F=f+hF_1+h^2F_2+h^3F_3+h^4F_4+O(h^5),\qquad
                         \deg F_{m,i}\le25,                 \tag{1}
\]
with coefficient weights \(w_i-3m\), satisfying all six canonical
commutators through \(h^5\). The coefficients may have either parity.

Hence an arbitrary-parity full Long lift in this grading, with a
common spatial cap, requires cap at least 26. This gives no cap-26
conclusion, normalization into the grading, nongraded obstruction
or DC1/DC2 resolution.

The ENTIRE coefficient space at index \(m\) consists of
\[
 X^\alpha Q^\beta P^\gamma Z^\delta,\qquad
 \alpha=\beta+4\gamma+2\delta+3m-w_i\ge0,\quad
                 \alpha+\beta+\gamma+\delta\le25.            \tag{2}
\]
Here exponent names are unrelated to the scalar parameters below.
The linear map \(A_m\) sends \(U\) to
\(\{U_i,f_j\}+\{f_i,U_j\}\) on all six pairs.

The commutator coefficient of order \(h^N\) is
\[
 \sum_{\substack{a+b+\ell=N\\\ell\ge1\ {\rm odd}}}
             \frac{\Pi^\ell(F_{a,i},F_{b,j})}{2^{\ell-1}\ell!}.
                                                               \tag{3}
\]
The unknown is \(A_{N-1}F_{N-1}\). Write its negative known right
side as \(B_N\). At \(h^5\) this includes
\(\Pi^3(F_{1,i},F_{1,j})/24\); the even-jet formula is insufficient.

## The nine-dimensional first kernel

Put
\[
\begin{split}
 t&=XQ,\quad B=3X^2P+2(1-3XQ)Z,\quad A_0=Q+3XQ^2-XB/3,\\
 Y&=RA_0-3t,\\
 H_4&=Y^4+4tY^2-\tfrac{16}{3}tY+10t^2,\\
 H_6&=R^4D+
 \tfrac1{25}(9Y^5+45tY^3+45t^2Y-60tY^2+110tY-720t^2),\\
 M&=X^2Z,\\
 H_5&=9Y^5+45tY^3-60tY^2+135t^2Y-10tY-555t^2,\\
 J&=(RS)^2 .
\end{split}                                                    \tag{4}
\]
These are polynomials in commuting Weyl-symbol variables. The first space has
component counts \((102,141,171,127)\), totaling 541. A retained
532-square minor is \(930421\) modulo \(1000003\). Nine exact
independent kernel vectors are the Hamiltonian variations of
\[
                    t,Y,Y^2,R^2T,H_4,H_6,M,H_5,J.             \tag{5}
\]
Their retained nine-square coefficient minor has determinant
\(-109117440\) over \(\mathbb Q\). Thus the rational rank is exactly
532 and (5) is complete.

The new cap-25 direction is \(J\), with variation
\[
                         (0,2R^2S,-2RS^2,0),
\]
of component degrees \((-,17,25,-)\). It is retained. Its finite
abstract target flow does not provide a uniform spatial-cap
normalization for all corrected tuples.

The identity \(RS=-6t-3Y^2-8Y\) permits replacing \(Y^2\) by \(RS\).
The two exact shears have the formulas
\[
\begin{gathered}
 P\mapsto P-\lambda hQ,\quad Z\mapsto Z-\lambda hX,
                          \quad X,Q\text{ fixed},\\
 T\mapsto T+\mu hR,\quad D\mapsto D-\mu hS,
                          \quad R,S\text{ fixed}.
\end{gathered}
\]
They are linear in the source generators or in the four images, so
both preserve any common cap and the grading. After these two
justified first gauges,
the COMPLETE first Hamiltonian is
\[
                  G=aY+bR^2T+cH_4+dH_6+eM+fH_5+\eta J,
 \qquad F_1=\{G,f\}.                                        \tag{6}
\]
The parameter \(\eta\) is distinct from the formal variable \(h\).

## Exact h3 field classification and all second coefficients

The entire \(F_2\) space has counts \((71,102,127,91)\), totaling
391. Its rank is 387: a 387-square minor is \(88338\) modulo
\(1000003\), and the four exact independent kernels are
\[
                  K_{R^3},K_{R^4S},K_{X^3},K_{R^4a_0},
 \qquad K_H=\{H,f\},\quad R^4a_0=R^3(Y+2).                   \tag{7}
\]
Here \(a_0=A_0+X^{-1}\), while \(R^4a_0=R^3(Y+2)\) is globally polynomial. The \(X^3\) kernel has a nonzero
momentum-linear \(S\)-component, while the last has nonzero positional
\(S\)-component. Together with the nonzero \(T\)-component of
\(K_{R^4S}\) and the \(D\)-only \(K_{R^3}\), this proves independence.

Two retained rational row functionals, supported on 369 and 71 rows,
annihilate all 391 columns and evaluate the complete \(h^3\) right
side to
\[
                  d^2-250df+35625f^2,\qquad e\eta.             \tag{8}
\]
The full row set has 2333 entries.

The witness supplies an admitted polynomial second particular
\(\Psi\) with the following explicit reduction convention:
\[
               d^2\longmapsto250df-35625f^2,\qquad
               e\eta\longmapsto0.                            \tag{9}
\]
The first rewrite strictly lowers the \(d\)-exponent; the second
sets every divisible monomial to zero. Each replacement subtracts
a multiple of the corresponding polynomial in (8).
The replay reconstructs \(B_3\), performs these declared reductions,
and verifies
\[
                         A_2\Psi=\operatorname{red}(B_3)
\]
coefficientwise. Its constant tuple is exactly the retained
59-coefficient \(g_0\) in
[the earlier witness](graded_long_cap23_witness.json).
No unretained algebraic solve is used.

Thus (8) is necessary and sufficient through \(h^3\), and the
ENTIRE second family on that locus is
\[
 F_2=\Psi(a,b,c,d,e,f,\eta)
          +uK_{R^3}+vK_{R^4S}+\gamma K_{X^3}+\delta K_{R^4a_0}.
                                                               \tag{10}
\]
All four additional freedoms remain. The conic in (8) retains the
two branches \(d=(125\pm100\sqrt{-2})f\) after constant extension.
No reality condition or separate real/imaginary parameter equations
are imposed.

## Exact h4 consequences and all field strata

Allow every \(F_3\) coefficient in (2): its counts are
\((47,71,91,62)\), totaling 271. Its rank is 270, certified by a
270-square minor equal to \(167431\) modulo \(1000003\), and by
the exact nonzero kernel
\[
                         K_{R^6}=(0,0,-6R^5,0).              \tag{11}
\]

First reconstruct the ENTIRE \(B_4\) using (6) and (10), with the
valid reductions (9). Two rational functionals, supported on 217
and 49 specified rows, annihilate all 271 columns and give exactly
\(f^3\) and \(e^3\). Over a field this forces
\[
                              f=e=d=0,                       \tag{12}
\]
where \(d=0\) then follows from the conic in (8). Both complex
branches are covered by this individual cube argument.

After precisely (12), define
\[
 \ell=a-\frac{544}{27}c,\qquad
 \mathcal K=976c^2+5184c\eta+27\delta.
\]
Three further retained rational functionals, each with support 209,
annihilate every \(A_3\) column and give
\[
                \ell\mathcal K,\qquad
                (c+36\eta)\mathcal K,\qquad
                864\gamma\eta-\ell^3.                        \tag{13}
\]
These necessary identities are checked on the full reduced
right-side coefficient set, with 1844 support rows.

All field points of (13) are covered by the following four strata.
The parameters \(b,u,v,w\) remain arbitrary throughout.

| Stratum | Parameter conditions |
|---|---|
| \(A_\times\) | \(\eta\ne0,\ a=544c/27+\ell,\ \gamma=\ell^3/(864\eta),\ \delta=-(976c^2+5184c\eta)/27\) |
| \(A_0\) | \(\eta=0,\ a=544c/27,\ \delta=-976c^2/27,\ \gamma\) arbitrary |
| \(B_\gamma\) | \(a=-2176\eta/3,\ c=-36\eta,\ \gamma=0,\ \delta\) arbitrary |
| \(B_0\) | \(a=c=\eta=0,\ \gamma,\delta\) arbitrary |

To see coverage, if \(\mathcal K=0\), split into \(\eta\ne0\)
and \(\eta=0\) and use the third equation of (13). If
\(\mathcal K\ne0\), the first two equations give
\(\ell=0,c=-36\eta\). Then \(\gamma\eta=0\), giving the last
two strata. Their overlaps and all denominator-zero points are
included. No nonreduced ideal equality is asserted.

For EACH stratum the witness retains every coefficient of actual
\(F_1,F_2,F_3\) frames. Their \(F_1,F_2\) are verified to be the
stated parameter substitutions in (6) and (10). Their entire
\(F_3\) families are verified coefficientwise as solutions of
\(A_3F_3=B_4\), including an arbitrary \(wK_{R^6}\).
Rank \(270/271\) makes these families exhaustive. Their respective
\(F_3\) coefficient counts, including the kernel, are
\(664,259,231,67\).

On \(A_\times\) these are Laurent polynomials in \(\eta\),
with only the explicitly allowed \(\eta^{-1}\) denominators.
The other three strata cover \(\eta=0\) without inversion.
These substitutions parameterize compatibility points; they are
not operator automorphisms or gauge normalizations.

## Every entire stratum fails at h5

The full \(F_4\) space has counts \((29,47,62,40)\), totaling 178.
In particular all 29 first-image coefficients remain unknown.

For each entire stratum, use all the retained \(F_1,F_2,F_3\) in
the order-\(h^5\) coefficient formula (3).
The four right sides have respectively \(31,12,12,9\)
nonzero parameter monomials. Every row set has 1446 entries.

The witness supplies four exact rational functionals, with supports
\(142,138,139,140\), one for each stratum. They satisfy
\[
                   \lambda A_4=0
          \quad\hbox{on ALL 178 coefficient columns},
 \qquad \lambda B_5=1                                     \tag{14}
\]
as full parameter identities. On \(A_\times\), (14) is an identity
in the Laurent parameter ring with \(\eta\ne0\); every Laurent
coefficient is checked. Every nonconstant coefficient evaluates to
zero on all four strata. No free parameter is discarded or sampled.

Thus each entire field stratum would require \(0=1\).
Equations (5)--(14) cover every admitted jet in (1) after only the
two uniform first gauges, proving the theorem.

## Retained evidence and scope

[The witness](graded_long_cap25_full_jet_witness.json) stores three
rank minors, the nine-kernel independence minor, the full second
family with its explicit reductions, the staged necessity
functionals, and all four complete third families and final
functionals. The [checker](check_graded_long_cap25_full_jet.py)
uses Python Fraction polynomial arithmetic, bare modular
determinants, and rational row evaluations. It verifies the
parameter pullbacks and every coefficient identity above, including
the \(\Pi^3(F_1,F_1)\) contribution and rejected corruptions of the
four final functionals.

A separate normal-PBW multiplication control verifies all six
commutators through \(h^4\) at the new \(A_\times\) point
\[
 (b,c,\ell,\eta,u,v,w)=(0,0,1,1,0,0,0),
 \qquad a=1,\quad\gamma=1/864,\quad\delta=0.
\]
It converts the retained Weyl symbols to normal order and multiplies
them using the PBW rule, separately from the Moyal right-side engine.
This single positive control checks conventions and the admitted
family; the universal obstruction uses the exact parameter identities.

The [receipt](graded_long_cap25_full_jet_controls.json) pins
repository-relative inputs and the canonical Long producer, which
supplies polynomial input only. The theorem is a field-valued
classification in the specified grading. It proves no nonreduced
scheme equality, no bounded normalization into the grading,
no cap-26 obstruction, no full nongraded nonlifting theorem, and
no DC1/DC2 resolution.
