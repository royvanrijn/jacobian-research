# Repeated-root quartic binary GVC jet search

## 1. Status and outcome

This note records an **experiment**, not a theorem and not a counterexample.
It attacks the genuinely nonhomogeneous part of the binary degree-five GVC
frontier left open in
[`BINARY_DEGREE_FIVE_GVC_FRONTIER.md`](BINARY_DEGREE_FIVE_GVC_FRONTIER.md).
The subsequent exact
[quadruple-root theorem](BINARY_QUARTIC_QUADRUPLE_ROOT_GVC.md),
[triple-plus-simple theorem](BINARY_QUARTIC_TRIPLE_SIMPLE_ROOT_GVC.md),
and [double-root theorem](BINARY_QUARTIC_DOUBLE_ROOT_GVC.md)
close all four repeated-root quartic orbits with arbitrary higher jets.
The data below are retained as discovery history, not as a dependency of
those proofs.

The first architecture capable of moving a mixed defect to depth proportional
to the moment is
\[
 \Lambda=\Lambda _4+\Lambda _5,\qquad P=P_5+P_4.             \tag{1.1}
\]
Bounded modular searches on the repeated-root quartic orbits
\((4),(3+1),(2+2)\) found no counterexample.  On the \((4)\) orbit, a
conditioned fifth-moment search left two points, but both lie on an explicit
support-separated face and therefore satisfy the GVC conclusion by a direct
degree estimate.  Adding the complete defect-two data
\[
 \Lambda=\Lambda _4+\Lambda _5+\Lambda _6,\qquad
 P=P_5+P_4+P_3                                             \tag{1.2}
\]
left one delayed point through moment three in 500,000 conditioned samples
over \(\mathbf F_{29}\); its fourth pure moment is nonzero.

The experiment itself does not exclude characteristic-zero points outside
the samples, later failure of a bounded survivor, the \((2+1+1)\) orbit,
deeper jets, or additional lower pieces.  The later exact theorems do.

## 2. Why (1.1) is the minimum migrating ansatz

Let \(j\) be the number of \(\Lambda _5\) selections in a term of
\(\Lambda^m\), and let \(\ell\) be the number of \(P_4\) selections in a term
of \(P^m\).  Relative to the maximum output degree \(m\), that term has defect
\[
 k=j+\ell.                                                  \tag{2.1}
\]
Thus one fixed output layer receives contributions from all pairs
\((j,\ell)\) with \(j+\ell=k\).  For \(k\) proportional to \(m\), these
contributions can in principle cancel while escaping every fixed-layer
vanishing theorem.

With homogeneous \(P_5\), the \(\Lambda_5\)-selection number determines the
output degree by itself, so this first cancellation is absent.  With
homogeneous \(\Lambda_4\), the split-symbol theorem already proves GVC for an
arbitrary nonhomogeneous \(P\).  Both off-leading pieces in (1.1) are
therefore essential to the first genuinely new face.

At defect two, the new primitive selections are \(\Lambda_6\) and \(P_3\);
they interact with two defect-one selections.  This gives (1.2).

## 3. Normalized repeated-root slices

The search uses the following monomial representatives and differential-unit
quotients:
\[
\begin{array}{c|c|c|c}
\text{orbit}&\Lambda_4&P_5&\text{retained }\Lambda_5
\\ \hline
(4)&X^4&y^2C_3(x,y)&
 \langle Y^5,XY^4,X^2Y^3,X^3Y^2\rangle\\
(3+1)&X^3Y&x^5\text{ or }y^3C_2(x,y)&
 \langle Y^5,XY^4,X^2Y^3,X^5\rangle\\
(2+2)&X^2Y^2&x^4L(x,y)\text{ or }y^4L(x,y)&
 \langle Y^5,XY^4,X^4Y,X^5\rangle .
\end{array}                                                  \tag{3.1}
\]
The omitted quintic terms are precisely the products of \(\Lambda_4\) with
a linear differential unit.  For the defect-two \((4)\) search, the retained
sextic jet is
\[
 \Lambda_6\in
 \langle Y^6,XY^5,X^2Y^4,X^3Y^3\rangle.                    \tag{3.2}
\]

All arithmetic is exact in the declared finite field.  The primes satisfy
\[
 p>5m_{\max},                                                \tag{3.3}
\]
so none of the falling-factorial differentiation coefficients in the pure
moment window degenerates merely because an exponent reaches the
characteristic.

## 4. Search results

The unconditioned defect-one search first solves moment one for the
coefficient of \(P_4\) dual to \(\Lambda_4\).  At \(p=23\), 100,000
normalized samples on each orbit give:
\[
\begin{array}{c|r|r|r|r}
\text{orbit}&\text{samples}&M_2=0&M_2=M_3=0&M_2=M_3=M_4=0\\ \hline
(4)&100000&35&1&1\\
(3+1)&100000&84&0&0\\
(2+2)&100000&17&2&2.
\end{array}                                                  \tag{4.1}
\]
These are bounded survivor counts, not point counts of an algebraic scheme.

The \((4)\) orbit was then conditioned on moment one and the two degree-one
coefficients of moment two.  Among 500,000 normalized samples at \(p=29\),
\[
\begin{array}{c|r}
\text{stage}&\text{survivors}\\ \hline
\text{solvable defect-one affine system}&1805\\
M_2=0&55\\
M_2=M_3=0&2\\
M_2=M_3=M_4=0&2\\
M_2=\cdots=M_5=0&2.
\end{array}                                                  \tag{4.2}
\]
Both final points have
\[
 \max_x\operatorname{Supp}(P_5)=1,\qquad
 \max_x\operatorname{Supp}(P_4)=3,\qquad
 \min_x\operatorname{Supp}(\Lambda_5)=2.                    \tag{4.3}
\]
They are GVC-safe for an all-order reason.  Let \(h\) be the number of
\(P_5\) selections and \(j\) the number of \(\Lambda_5\) selections in a
putative contribution to \(\Lambda^m(QP^m)\), where
\(\deg Q=q\) and \(\deg_xQ=q_x\).  Total degree requires
\[
 j\le h+q.                                                   \tag{4.4}
\]
The operator \(x\)-order and input \(x\)-degree satisfy
\[
\begin{aligned}
 \operatorname{ord}_x&\ge4(m-j)+2j=4m-2j,\\
 \deg_x(\text{input})&\le h+3(m-h)+q_x=3m-2h+q_x.
\end{aligned}                                                \tag{4.5}
\]
Their difference is at least \(m-2q-q_x\), hence is positive for
\(m>3q\).  Therefore every mixed contraction vanishes beyond that bound.

For the defect-two search, the early affine equations were solved rather
than sampled blindly.  Among 500,000 \((4)\)-orbit samples at \(p=29\),
\[
\begin{array}{c|r}
\text{stage}&\text{survivors}\\ \hline
\text{solvable through defect one}&1722\\
\text{solvable through }M_2&1108\\
M_2=M_3=0&1\\
M_2=M_3=M_4=0&0.
\end{array}                                                  \tag{4.6}
\]
The unique delayed point is, modulo \(29\),
\[
\begin{aligned}
P_5={}&y^5+xy^4,\\
P_4={}&24y^4+20xy^3+25x^2y^2+18x^3y+4x^4,\\
P_3={}&25y^3+12xy^2+2x^2y+14x^3,\\
\Lambda_5={}&26Y^5+11XY^4+24X^2Y^3+19X^3Y^2,\\
\Lambda_6={}&26Y^6+25XY^5+23X^2Y^4+27X^3Y^3.
\end{aligned}                                                \tag{4.7}
\]
It satisfies the first three pure equations but
\[
 \Lambda^4(P^4)=20+17y\ne0\pmod {29}.                       \tag{4.8}
\]
This is a delayed failure, not a lift candidate.

## 5. Reproduction and next gate

Run the pinned search with:

```bash
.venv/bin/python \
  scripts/search_binary_repeated_quartic_gvc_jets_mod_p.py
```

It takes roughly two minutes on the development machine.  A small regression
run that prints its result instead of replacing the pinned artifact is:

```bash
.venv/bin/python \
  scripts/search_binary_repeated_quartic_gvc_jets_mod_p.py --quick
```

The generated record is
[`binary_repeated_quartic_gvc_jet_search.json`](../artifacts/generated-results/binary_repeated_quartic_gvc_jet_search.json),
with whole-file SHA-256
`299297fcf936c9041cb2da7ab4f7d124271504e09af02066f92d0cfc4e930f14`.

The later exact calculations carry out all three proposed next steps:
they compute the \((3+1)\) defect-one radical, close the \((2+2)\)
threshold face, and treat the previously omitted \((2+1+1)\) orbit.  No
binary GVC counterexample occurs through polynomial degree five.
