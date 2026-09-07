# Rank-ten Gale-matroid survivors

## 1. Outcome and scope

The cyclic-complement lemma stated in `OP-QHNW10` is false.  It is false
literally because loops are allowed, and it remains false after adding
looplessness or even simplicity.  Thus the rank-ten frontier must use the Hessian trace
identities on the surviving matroids; the proposed purely matroidal route to
Waring rank eleven cannot be used as stated.

The first exact checker is
[`verify_quartic_hn_rank10_matroid_survivors.py`](../scripts/verify_quartic_hn_rank10_matroid_survivors.py).
It also freezes a complete **loopless, nonsimple** catalogue slice relative
to `matroid-database==0.3`: 37 coloured isomorphism types survive the
combinatorial constraints, 35 of them are rationally representable, and none
of those 35 supports any rank-six Gram matrix satisfying even the first two
Hessian trace identities.  This is not a catalogue of simple rank-four
matroids on ten elements.

The complete Gale-loop slice contains 115 abstract coloured types, of which
111 are characteristic-zero realizable.  The complementary-splitting
argument closes all 66 types with at least two loops, and a matroidal
self-square support obstruction closes 16 of the 45 one-loop types.  Exact
realization-open saturations close 25 more, one disconnected active matroid
splits into low-dimensional summands, and a loop--triple Witt obstruction
closes the final three.  Thus **all 111 characteristic-zero Gale-loop types
are impossible**.

The final simple census contains five abstract survivors, all rationally
representable.  The checker
[`verify_quartic_hn_rank10_simple_survivors.py`](../scripts/verify_quartic_hn_rank10_simple_survivors.py)
closes all five by a universal six-point rank-two-flat obstruction using
only the first two traces.  Consequently the complete rank-ten Gram branch
is empty and every essential six-variable quartic HN counterexample has
Waring rank at least eleven.

## 2. The literal loop survivor

Let

\[
 M=U_{4,6}\oplus U_{0,4}.
 \tag{2.1}
\]

This matroid is representable over every field.  Its cocircuits are exactly
the three-element cocircuits of \(U_{4,6}\); the four loops occur in no
cocircuit.  Every nonloop parallel class is a singleton.  A basis uses four
of the six nonloops, so its six-element complement consists of four loops
and two nonloops.  Both nonloops are coloops of that restriction.  Hence no
basis has cyclic complement.

Therefore `OP-QHNW10`, read literally, already has a characteristic-zero
survivor.

It does not survive the Gram equations.  The four Gale-loop coordinates are
free coordinates on \(\ker K\).  The first two traces make the Gram block on
the corresponding four Waring covectors zero coefficientwise.  Those
covectors are independent by Gale duality, so they would span a
four-dimensional totally isotropic space, again exceeding Witt index three.

## 3. A loopless rational survivor

Looplessness does not repair the statement.  Consider the rank-four
seven-point simplification represented by

\[
 K_{\mathrm{simp}}(a,b)=
 \begin{pmatrix}
 1&0&1&0&0&0&1\\
 0&1&1&0&0&0&a\\
 0&0&0&1&0&1&1\\
 0&0&0&0&1&1&b
 \end{pmatrix},
 \qquad ab(a-1)(b-1)\ne0.
 \tag{3.1}
\]

Give the first three projective points multiplicity two and the last four
multiplicity one.  Equivalently, duplicate columns 1, 2, and 3 of (3.1).
The resulting \(4\times10\) Gale matrix has three parallel pairs and no
loops or triple class.

The normalized realization ideal has no equations.  Its realization scheme
over \(\mathbb Q\) is the open set

\[
 \boxed{
 \operatorname{Spec}
 \mathbb Q[a,b,(ab(a-1)(b-1))^{-1}].
 }
 \tag{3.2}
\]

Thus \((a,b)=(2,3)\) is an exact rational realization.  Direct rank-oracle
enumeration verifies:

1. every cocircuit has size at least three;
2. the parallel multiplicities are \((2,2,2,1,1,1,1)\); and
3. the complement of every basis has a coloop.

In the Matsumoto--Moriyama--Imai--Bremner reverse-lex catalogue used by the
checker, the simplification is the seven-element rank-four entry with local
index 19 and rankline

```text
00***0******0000******************0
```

This is a characteristic-zero, loopless counterexample to the proposed
cyclic-complement lemma.

## 4. Realization-ideal census below ten simple points

Filtering the complete rank-four catalogues on five through nine elements,
then assigning total multiplicity ten subject to the double/triple bound,
gives 37 loopless nonsimple coloured isomorphism types.  They lie over 17
simple matroids:

\[
 5\text{ on 7 points},\qquad
 7\text{ on 8 points},\qquad
 5\text{ on 9 points}.
 \tag{4.1}
\]

For every underlying matroid the checker chooses a basis, writes a
fundamental matrix \([I_4\mid D]\), normalizes a spanning forest in the
nonzero bipartite support of \(D\), and constructs the realization data

\[
 I_M=\bigl(\det K_J:J\text{ a nonbasis}\bigr),
 \qquad
 \prod_{B\text{ a basis}}\det K_B\ne0.
 \tag{4.2}
\]

Two eight-point simplifications, local indices 269 and 845, have unit
normalized ideals over \(\mathbb Q\).  The other 15 underlying matroids have
exact rational points on (4.2), giving exactly 35 characteristic-zero
coloured survivor types.  The catalogue completeness here is relative to
the external `matroid-database==0.3` files.  The checked-in script freezes
and independently replays the 17 returned ranklines, their realization
ideals, and the 37 colour orbits; it does not vendor or regenerate the
external catalogue.  The source wheel used for the extraction had SHA-256

```text
85d0304575784ceb4797014fb5a761443d2a0bbf01aafc38ef293d2c64a1b5ce
```

Here is the complete coloured census.  A digit string records the
multiplicity of each ordered point of the simplification.  Types in the same
row have the same realization scheme (4.2).

| catalogue key | coloured multiplicity representatives | realization over \(\overline{\mathbb Q}\) |
| --- | --- | --- |
| (7, 15) | `1231111`, `2221111` | rational point |
| (7, 19) | `1111231`, `1112221` | rational point |
| (7, 26) | `1111123`, `2111113`, `2111122`, `3111112` | rational point |
| (7, 32) | `1111123`, `1112113`, `1112122`, `1113112` | rational point |
| (7, 54) | `1111123`, `1111222` | rational point |
| (8, 269) | `11111113` | unit ideal |
| (8, 581) | `11131111`, `11221111` | rational point |
| (8, 586) | `11131111`, `11221111`, `21121111`, `31111111` | rational point |
| (8, 587) | `11131111`, `11221111` | rational point |
| (8, 588) | `11131111`, `11221111`, `21121111`, `31111111` | rational point |
| (8, 589) | `11111113`, `11111122` | rational point |
| (8, 845) | `11111113` | unit ideal |
| (9, 188841) | `111121111` | rational point |
| (9, 188846) | `111121111`, `211111111` | rational point |
| (9, 188847) | `111121111` | rational point |
| (9, 188848) | `111121111`, `211111111` | rational point |
| (9, 188849) | `111121111` | rational point |

For a machine-readable construction of every scheme, run

```bash
.venv/bin/python scripts/verify_quartic_hn_rank10_matroid_survivors.py --details
```

Each `QHNW10_REALIZATION_RECORD` is exact JSON containing the rankline, the
normalized matrix, generators of \(I_M\), every nonconstant basis minor to
invert, and either a rational realization matrix or a certified unit-ideal
flag.  This makes the realization construction inspectable without
expanding the often large product in (4.2) in the note.

## 5. The six-point rank-two-flat trace obstruction

Every one of the 35 characteristic-zero types contains a rank-two Gale flat
\(R\) carrying exactly six elements, counted with parallel multiplicity,
and

\[
 \operatorname{rk}_K(E\setminus R)\ge2.
 \tag{5.1}
\]

The multiplicity profiles on \(R\) are

\[
 (3,2,1),\ (2,2,2),\ (3,1,1,1),\
 (2,2,1,1),\ (2,1,1,1,1).
 \tag{5.2}
\]

Let \(L_R=\ker K_R\), a four-plane in \(k^6\).  The first Hessian trace
restricted to \(L_R\) is

\[
 \sum_{i\in R}g_{ii}l_i^2=0.
 \tag{5.3}
\]

For every profile in (5.2), the coordinate-square restriction map

\[
 k^6\longrightarrow\operatorname{Sym}^2(L_R^*),
 \qquad
 (D_i)\longmapsto\sum_iD_il_i^2
 \tag{5.4}
\]

is injective.  In the same projective normalizations used below, the checker
certifies six-by-six minors

\[
 4,\quad4,\quad4\lambda,\quad4\lambda,\quad4\lambda.
 \tag{5.5}
\]

For the five profiles in their displayed order, the projective-direction
open set makes these determinants nonzero.  Consequently
\(g_{ii}=0\) for every \(i\in R\), without assuming that the full Gram
diagonal is isotropic.

The second trace restricted to \(L_R\) now becomes

\[
 \sum_{i<j\in R}2g_{ij}^{2}l_i^2l_j^2=0.
 \tag{5.6}
\]

The square-pair restriction map

\[
 k^{15}\longrightarrow \operatorname{Sym}^4(L_R^*),
 \qquad
 (A_{ij})\longmapsto\sum_{i<j}A_{ij}l_i^2l_j^2
 \tag{5.7}
\]

is injective for every profile in (5.2).  For three projective directions,
the checker certifies determinants \(-256\) and \(-512\) for profiles
\((3,2,1)\) and \((2,2,2)\).  For four directions normalized to
\(0,\infty,1,\lambda\), it certifies

\[
 -256\lambda^4,\qquad -512\lambda^2
 \tag{5.8}
\]

for \((3,1,1,1)\) and \((2,2,1,1)\).  These are nonzero on the required
open set.  With at least five directions, the triangle-circuit argument of
Section 5.1 of
[`QUARTIC_HN_WARING_RIGIDITY.md`](QUARTIC_HN_WARING_RIGIDITY.md)
annihilates every edge.

Consequently (5.3) and (5.6) force \(G_R=0\).  Gale duality gives

\[
 \operatorname{rk}_V(R)
 =|R|-4+\operatorname{rk}_K(E\setminus R)
 \ge4.
 \tag{5.9}
\]

The Waring vectors indexed by \(R\) would therefore span a totally
isotropic subspace of dimension at least four in a nondegenerate
six-dimensional quadratic space.  Its Witt index is three, a contradiction.

Thus

\[
 \boxed{
 \text{none of the 35 loopless nonsimple characteristic-zero survivors
 supports a rank-six Gram matrix satisfying the first two HN traces.}
 }
 \tag{5.10}
\]

For the explicit survivor of Section 3, the six-element flat is the union
of the three parallel pairs and (5.9) has rank five.

## 6. The Gale-loop census

The same external catalogues give a complete loop census because a
ten-element Gale matroid with at least one loop has at most nine nonzero
elements.  After quotienting parallel-multiplicity colourings by the
automorphism group of the simple nonzero-column simplification, the counts
are

| Gale loops | abstract coloured types | characteristic-zero types |
| ---: | ---: | ---: |
| 1 | 47 | 45 |
| 2 | 49 | 48 |
| 3 | 18 | 17 |
| 4 | 1 | 1 |
| **total** | **115** | **111** |

The 115 types lie over 66 simplifications.  Three normalized realization
ideals are the unit ideal over \(\mathbb Q\), removing four coloured types;
the other 63 simplifications have exact rational realization matrices.

The external extraction is performed by
[`enumerate_quartic_hn_rank10_loop_survivors.py`](../scripts/enumerate_quartic_hn_rank10_loop_survivors.py)
and frozen in
[`quartic_hn_rank10_loop_survivors.json`](../artifacts/generated-results/quartic_hn_rank10_loop_survivors.json).
The catalogue wheel, catalogue-file bundle, and generated artifact have
respective SHA-256 hashes

```text
85d0304575784ceb4797014fb5a761443d2a0bbf01aafc38ef293d2c64a1b5ce
7c30ea1ef0a8ccf7548e7fbbdb64b636c889fe2185dae59f9a46c6571c6bb006
a7f5c04a05cdd6d32c723ae4ecfdfbb8c2018cf2072f00254f7c3a50c15e71b9
```

### 6.1 Two or more loops are impossible

Let \(Z\) be the set of \(z\) Gale loops and let \(O\) be its complement.
Gale duality gives

\[
 \operatorname{rk}_V(Z)=z,\qquad
 \operatorname{rk}_V(O)=6-z.
 \tag{6.1}
\]

The two Waring spans are therefore complementary.  After a linear source
change the quartic splits into essential summands in \(z\) and \(6-z\)
variables, so its Hessian determinant is, up to a nonzero scalar, the
product of the two summand Hessian determinants.  For \(2\le z\le4\), both
summands use at most four variables.  The characteristic-zero
low-dimensional Hesse theorem makes the Hessian determinant of each
essential summand nonzero.  Their product cannot vanish, contradicting HN.

Thus all 66 characteristic-zero types with at least two loops are closed.

### 6.2 The one-loop self-square obstruction

Let \(f\) be the Waring vector indexed by the unique Gale loop, let \(K_A\)
be the active \(4\times9\) Gale matrix, and put \(L_A=\ker K_A\).  If

\[
 a_j=\langle f,v_j\rangle ,
\]

then \(GK^{\mathsf T}=0\) gives \(a\in L_A\).  The free coordinate in the
first trace gives \(\langle f,f\rangle=0\), while its squared coefficient in
the second trace is

\[
 \sum_{j=1}^9a_j^2l_j^2=0
 \qquad(l\in L_A).
 \tag{6.2}
\]

Nondegeneracy forces \(a\ne0\).  Hence a one-loop HN candidate requires a
nonzero **self-square vector**: a vector \(a\in L_A\) whose coordinatewise
square is a diagonal quadratic relation on \(L_A\).

There is a realization-independent support sieve.  Put
\(T=\operatorname{supp}(a)\).  The full-support dependence \(K_Aa=0\) makes \(T\) cyclic in
the active Gale matroid.  The coordinate projection of \(L_A\) to \(T\) is
totally isotropic for the nondegenerate diagonal form in (6.2), so

\[
 \operatorname{rk}_V(T)
 =|T|-4+\operatorname{rk}_{K_A}(E_A\setminus T)
 \le\left\lfloor\frac{|T|}{2}\right\rfloor.
 \tag{6.3}
\]

Sixteen of the 45 characteristic-zero one-loop types have no subset \(T\)
satisfying both conditions and are therefore closed for every realization.
Together with Section 6.1, this immediately closes 82 of the 111
characteristic-zero loop types.

For all 45 one-loop types, the checker
[`verify_quartic_hn_rank10_loop_survivors.py`](../scripts/verify_quartic_hn_rank10_loop_survivors.py)
constructs an exact rational realization.  At ten points the nine
coordinate squares are independent; at the other 35 points all five
projective charts of the self-square scheme have unit Gröbner basis.  Since
the self-square incidence is projective over the realization base, its image
is closed.  Thus every type has a nonempty Zariski-open neighbourhood on
which the first two traces are impossible.

### 6.3 Exact closure of the special one-loop loci

After the support sieve, 29 types could still have lived on special closed
realization loci.  The same SymPy checker saturates the realization and
self-square equations by the product of the required basis minors in all
five projective kernel charts.  Twenty-one of the 29 saturated schemes have
unit ideal over \(\mathbb Q\).  Thus Sections 6.1--6.3 already close 103 of
the 111 characteristic-zero loop types.

The remaining eight types have three independent closures.

First, the active Gale matroid of catalogue type `(7,70)` with colour
`1111113` is disconnected.  Its components have Waring dimensions three
and two; the Gale loop supplies a one-variable component.  The quartic
therefore splits into essential summands in dimensions \(3,2,1\).
Hessian-determinant factorization and the characteristic-zero Hesse theorem
in at most four variables exclude it.

Second, the following four exact realization-plus-self-square saturations
have unit ideal over \(\mathbb Q\):

| catalogue key | colour |
| --- | --- |
| `(7,3)` | `1111113` |
| `(8,581)` | `11121111` |
| `(9,188841)` | `111111111` |
| `(9,188846)` | `111111111` |

The first, second, and fourth rows are checked in all five projective kernel
charts.  The third has a unique support passing (6.3), so its forced-support
chart suffices.  The checker
[`verify_quartic_hn_rank10_loop_closure.py`](../scripts/verify_quartic_hn_rank10_loop_closure.py)
reconstructs the ideals from the frozen ranklines, replaces the realization
open set by a squarefree basis-minor product, and asks Singular for the exact
characteristic-zero standard bases.

It remains to close catalogue types `(7,0)`, `(7,1)`, and `(7,2)`, with
colours `1111113`, `1111113`, and `3111111`.  Each has exactly one support
passing (6.3), and that support is disjoint from its active triple parallel
class \(C\).  Let \(f\) be the Waring vector of the Gale loop.  Then (6.2)
gives

\[
 \langle f,f\rangle=0,\qquad
 \langle f,v_c\rangle=0\quad(c\in C).
 \tag{6.4}
\]

The three Waring vectors indexed by \(C\) are independent.  On their
three-dimensional slice, write their linear forms as \(y_1,y_2,y_3\) and
put \(s=y_1+y_2+y_3\).  Gale parallelism makes every outside linear form a
multiple of \(s\) on this slice.  The first trace is therefore a linear
combination of

\[
 s^2,\ y_1^2,\ y_2^2,\ y_3^2,
 \tag{6.5}
\]

which are independent.  Hence the three Gram diagonals on \(C\) vanish.
The second trace is a linear combination of the seven independent quartics

\[
 s^4,\quad s^2y_i^2\ (1\le i\le3),\quad
 y_i^2y_j^2\ (1\le i<j\le3).
 \tag{6.6}
\]

The coefficient of \(y_i^2y_j^2\) after the preceding diagonal vanishing is
twice \(\langle v_i,v_j\rangle^2\), so every mutual Gram entry on \(C\)
vanishes.  Thus \(C\) spans a totally isotropic three-plane.  Equation (6.4)
adjoins the independent isotropic vector \(f\), producing a totally
isotropic four-plane in a nondegenerate six-space.  This contradicts Witt
index three.

Consequently all 45 characteristic-zero one-loop types, and hence all 111
characteristic-zero Gale-loop types, are closed.  No loop type reaches the
higher Gram traces or the collision equations.

## 7. Complete simple ten-point census

The remaining branch is finite and nonempty at the matroid level.  The exact
enumerator
[`enumerate_quartic_hn_rank10_simple_survivors.py`](../scripts/enumerate_quartic_hn_rank10_simple_survivors.py)
starts from all 190,214 rank-four nine-element catalogue entries, of which
185,981 are simple.  Every simple rank-four ten-element matroid has a
non-coloop deletion which is a simple rank-four nine-element matroid, so this
loses no ten-point type.

A rank-preserving one-element extension is encoded by a modular cut in the
flat lattice of its deletion.  In rank four, simplicity fixes the cut on
rank-zero and rank-one flats, leaving Boolean variables only on lines and
hyperplanes.  The modular-pair axioms, the cocircuit bound, and the
cyclic-complement condition reject all but 13 deletion types before the
residual Boolean solve:

| exact scan outcome | deletion types |
| --- | ---: |
| coloop deletion | 68 |
| forced line--plane modular conflict | 185,661 |
| forced selected hyperplane of size at least seven | 239 |
| forced minimal cut already survives | 8 |
| optional cut still required | 5 |

Z3 then enumerates every modular cut on the 13 residual types.  Five have no
model; the other eight give 23 catalogue-extension models.  Exact abstract
isomorphism testing reduces them to precisely five ten-point matroids.  Their
basis counts

\[
 114,\quad110,\quad108,\quad105,\quad90
 \tag{7.1}
\]

are already pairwise distinct, so the five displayed representatives are
nonisomorphic.  The generated census is
[`quartic_hn_rank10_simple_survivors.json`](../artifacts/generated-results/quartic_hn_rank10_simple_survivors.json).
Finite-field representability is not used anywhere in this enumeration.
The pinned nine-element catalogue file has SHA-256
`c664af9c4a4edf406bbe6b973db88d6d0d411d9a2052801e685bb2701fd94819`;
the generated census has SHA-256
`9304dff4987e85e9e26694a62ea38fd5549441689cb5a297b3659f5616aac725`.
The enumerator and verifier have respective SHA-256 hashes
`ec393864a7054a5cf0d6fe486a558ccb2cd568871cf43633791fc8324b5e8f3b`
and
`4776d87d288fe8561664c476cef95730dd6f3817426a221b4764350158f3d63e`.

Thus the cyclic-complement lemma remains false even after imposing
simplicity.  The gain comes from the trace equations, not from repairing the
matroid statement.

## 8. Realization schemes of the five survivors

For each survivor, the verifier chooses a basis, writes a normalized
fundamental matrix \([I_4\mid D]\), and forms the exact realization scheme

\[
 \operatorname{Spec}
 \mathbb Q[x_0,\ldots,x_{s-1}]/I_M
 \setminus V\!\left(\prod_{B\in\mathcal B(M)}\det K_B\right).
 \tag{8.1}
\]

The complete schemes have the following small presentations.  The last
column gives one exact rational point in the displayed variable order.

| survivor | bases | variables | \(I_M\) | nonconstant open minors | rational point |
| ---: | ---: | ---: | --- | ---: | --- |
| 0 | 114 | 7 | \((x_3-x_2)\) | 73 | \((-2,-1,2,2,3,-2,-1)\) |
| 1 | 110 | 5 | \((0)\) | 41 | \((-2,-2,-1,2,-1)\) |
| 2 | 108 | 5 | \((0)\) | 50 | \((-2,-1,2,3,-2)\) |
| 3 | 105 | 4 | \((0)\) | 38 | \((-2,-2,-1,2)\) |
| 4 | 90 | 4 | \((0)\) | 40 | \((-2,-1,2,-2)\) |

Consequently all five abstract survivors are representable over
\(\mathbb Q\), not merely over a finite field or over
\(\overline{\mathbb Q}\).  Run

```bash
.venv/bin/python scripts/verify_quartic_hn_rank10_simple_survivors.py --details
```

to print every normalized matrix, ideal generator, open-minor count, rational
matrix, and trace determinant as exact JSON.

## 9. Universal six-point-flat trace obstruction

Every one of the five survivors has the same distinguished rank-two flat

\[
 F=\{0,1,2,3,4,9\}.
 \tag{9.1}
\]

The six columns are simple, hence are six distinct points of a projective
line.  After row operations, column scalings, and a relabelling they have the
form

\[
 K_F=
 \begin{pmatrix}
 1&0&1&1&1&1\\
 0&1&1&a&b&c
 \end{pmatrix},
 \tag{9.2}
\]

where \(0,1,a,b,c\) are pairwise distinct.  A basis of
\(\ker K_F\) gives the six value forms

\[
 -u_1-u_2-u_3-u_4,\quad
 -u_1-au_2-bu_3-cu_4,\quad
 u_1,\quad u_2,\quad u_3,\quad u_4.
 \tag{9.3}
\]

The verifier constructs their coefficient matrices directly.  A selected
\(6\times6\) minor of the six squares is

\[
 -4c(a-b),
 \tag{9.4}
\]

and a selected \(15\times15\) minor of the fifteen off-diagonal square
products is

\[
 -128bc^5(a-b)^3(b-c).
 \tag{9.5}
\]

Both are nonzero on the six-distinct-point realization open set.  The
checked matrices for the five normalized schemes may use different selected
coordinates, but every irreducible factor of their determinants is verified
to occur among the required nonzero basis-minor factors.

Restrict the first Hessian trace to the four-dimensional Waring-value
subspace supported on \(F\).  Independence in (9.4) forces

\[
 g_{ii}=0\qquad(i\in F).
 \tag{9.6}
\]

After (9.6), the restricted second trace is

\[
 2\sum_{i<j,\ i,j\in F}g_{ij}^2l_i^2l_j^2=0.
 \tag{9.7}
\]

Independence in (9.5) gives \(g_{ij}=0\) for all \(i,j\in F\).  Hence the
Waring span indexed by \(F\) is totally isotropic.  Gale duality gives

\[
 \operatorname{rk}_V(F)
 =|F|-\operatorname{rk}K+\operatorname{rk}_K(E\setminus F)
 =2+\operatorname{rk}_K(E\setminus F).
 \tag{9.8}
\]

The complement rank is three on the first four survivors and two on the
fifth, so this totally isotropic span has dimension five or four.  Both
contradict Witt index three for a nondegenerate symmetric form on a
six-space.

Therefore none of the five characteristic-zero simple survivors supports a
rank-six Gram matrix satisfying even the first two HN traces.  This argument
also covers the special realization loci on which the unrestricted
ten-coordinate first trace could admit a non-isotropic diagonal: the six
diagonal entries on \(F\) still vanish by (9.4).

## 10. Rank-ten closure

Sections 2--6 close all Gale-loop and loopless nonsimple types.  Sections
7--9 enumerate the complete simple branch, construct all five realization
schemes, and close every one by the six-point-flat trace obstruction.  Thus
there is no essential six-variable quartic HN polynomial of Waring rank ten.
Combined with the rank-nine invertibility theorem, this proves

\[
 \boxed{\text{every essential six-variable quartic HN counterexample has
 Waring rank at least eleven.}}
 \tag{10.1}
\]

No rank-ten trace component survives, so the higher traces and the
nonradial collision equations are not reached.  The calculation does not
construct a counterexample and has not received external mathematical
review.
