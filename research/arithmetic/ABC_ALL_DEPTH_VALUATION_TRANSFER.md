# All-depth valuation certificates and prime-power fibres

Hanxin Zhang's Theorem 1.1 gives the claimed uniform exponent
\(34/59\) for \(0<\lambda\le1\). Corollary 1.2 combines
it with the two earlier bounds to give
\[
                   \min\{2\lambda/3,(23\lambda+3)/40,34/59\}.
\]
The constant is the smallest of the three terms for
\(\lambda>1183/1357\). The estimate therefore addresses
a different range from the small-budget coefficient \(1/3\).
[Zhang, Theorem 1.1 and Corollary 1.2](https://doi.org/10.5281/zenodo.22151481)

## Analytic inputs and exact certificate

The paper retains each prime's original valuation, writing
\(T_m\) for the radical mass at valuations above \(m\).
The claimed all-depth inequality is \(59D-12T_m\le34\),
where \(D\) is the logarithmic count of a box in units of
\(\log X\). The numerical height moments give
\(T_m\le3/(m+1)\). The arithmetic rows supplement the
BBLT estimates with multilayer primitive slopes, the
Heath–Brown square-full estimate and a fourth-moment
estimate on an exposed fifth-power layer using Salberger's
diagonal-surface bound. The proof transports the weak
duality multipliers from depth eight and verifies later
layer columns by a period-210 calculation. Appendix A's
critical family is for the relaxation before the
exposed surface row, as the appendix specifies.
[Primary paper, §§2–8 and Appendix A](https://doi.org/10.5281/zenodo.22151481)

The published standard-library checker was replayed at
commit 4bf4eeb9dfa7b40ef5dca9f24fb005a74e468c02
with command \(\texttt{python3 -B -I -S verify.py}\).
It returned exit code zero, complete coverage, and the
claimed inequality \(59D-12T_m\le34\) for every \(m\ge8\)
on its encoded anatomical program. The checked inputs had
SHA-256 values

- verify.py — 01254809f2d43fd4a15a2af80924d3ad78b091654e3e6c2b23ef28ed72edc94f
- proof_data.txt — 617071475b8a23fa6bd6171c39e082581c618337d73f1139b920b834ff4a6fd1

This is a replay of the finite proof checker, with its
tail extension and regressions. The analytic transfer from
integer triples to every row has a separate audit. In
particular a changed selector family requires renewed
proof for that family.
[Public pinned checker](https://github.com/hanxinzhang/abc-universal-bound-34-59/tree/4bf4eeb9dfa7b40ef5dca9f24fb005a74e468c02)

Hongmin Li's revised radical-budget calculation reports a
smaller exponent than Zhang's pointwise envelope for
\(0.70\le\lambda<682839/724933\), with equality at the
upper endpoint of that comparison. The November versions
must be distinguished from the original archived claims
and from the aiXiv withdrawal of the first posting.
These are preprint statements about the stated ranges.
[Li, radical budget, version 4](https://rstation.ai/problems/math/papers/17/v4)
· [Li, squares and twisted energy, version 4](https://rstation.ai/problems/math/papers/28/v4)

## The divisibility row available at every depth

The prime-power fibre theorem gives an additional arithmetic
row on the original valuation coordinates. For a fixed
prime \(\ell\), \(d=\ell^j\) with
\((\ell-1)(d-2)/2\ge2\), and any permutation of the three
letters, set \(V_a(d)=\sum_{d\mid e}a_e\), and define
\(V_b(d),V_c(\ell)\) in the same way. With total radical
mass \(L=a+b+c\), it gives
\[
                        D\le L-V_a(d)-V_b(d)-V_c(\ell).
 \tag{AT1}
\]
For a finite visible subset of the indicated layers one
can use its mass in place of the full sum. It also follows
by applying the primitive fibre theorem with just those
visible factors put into the powers, and counting the
coefficient radicals by the complement. Thus tail
truncation with nonnegative mass preserves the row.
[Three-cap defect inequality](ABC_PRIME_POWER_FERMAT_DESCENT.md#a-three-cap-coefficient-defect-inequality)

This confronts a degeneration of fixed exposed valuation
indexes. Take a model whose three radical allocations are
concentrated on the same large prime degree \(j\) and whose
total radical mass is \(\lambda\) with \(j\lambda/3\le1\)
near equality. Any fixed finite list of low-index selectors
misses this mass. In the exact single-layer model, taking
\(d=\ell=j\) in (AT1) forces \(D\le0\). The fibre estimate also controls
nearby profiles with variable nontrivial coefficient
defects, by the radical mass of those defects. The proof supplies
the row for each fixed \(j\), and its fixed-coefficient
height estimate is uniform for variable large common orders
by enumeration of the two roots. In the small-budget
height limit the degree and residue selectors can be fixed
at each \(\lambda\) before taking \(X\) to infinity.

The same row family does not force a uniform saving on
arbitrary diffuse profiles. Equal mass spread over the
distinct prime valuations in \([J,2J]\) on every letter
gives total recovered mass \(o(L)\) for every choice in
(AT1), even when the prime selector moves with \(J\).
These mass data describe only the relaxation. The remaining
arithmetic step must couple distinct valuation classes
across the three letters.

For any fixed finite set of the new divisibility rows,
their valuation columns have periodic indicators, with
period dividing the least common multiple of the orders
used and \(210\). They are compatible with the
periodic–affine tail form of the all-depth certificate.
A new set of dual multipliers would still have to be
checked. At a fixed small budget \(\lambda\), on balanced
boxes with each letter's radical cap at most
\((1/3+\tau)\lambda\) and \(0<\tau<1/3\), only
\(d=O_\tau(1/\lambda)\) can recover more than
\(2\lambda/3\). Indeed the third term's mass is at most
\((1/3+\tau)\lambda\), and the first two satisfy
\(d(V_a(d)+V_b(d))\le2\) by their height moments.
The useful selectors can thus be fixed before taking
the height limit at this budget.


## Three unequal orders at one prime

For \(\ell\)-power orders \((m_a,m_b,m_c)\) with
\(1/m_a+1/m_b+1/m_c<1\), there is an all-depth row
\[
                         D\le L-\sum_i V_i(m_i).
\]
It follows from a count depending only on the joint defect
support. The cover uses an \(\ell\)-group, so its
specialization fields have a bounded cyclic-\(\ell\) tower.
For example the orders may be \((2,4,8)\), \((4,4,4)\)
or \((3,3,9)\) and any coordinatewise multiples. The
period argument above still applies, with the three orders
included in the least common multiple.
[Three-order prime-primary support bound](ABC_MISMATCHED_GRID_NUMBER_FIELDS.md#a-support-only-bound-from-a-single-prime-in-the-three-orders)
