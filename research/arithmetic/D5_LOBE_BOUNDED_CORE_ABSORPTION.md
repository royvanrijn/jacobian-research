# A root-preserving thirteen-sheet core absorbs D5 lobes

A canonical wreath-leaf tuple with one appended D5 lobe admits a
uniform conversion back to the canonical wreath-leaf orbit. Two
secondary-edge sweeps reduce every degree to the same thirteen-sheet
core. The proved degree-thirteen polynomial connectivity supplies
one fixed controller, and a final elementary slide restores the
whole star while retaining the common root. Iteration gives a quadratic
total upper bound for the unbuffered many-lobe bouquet.

This settles a particular absorption family. It does not give arbitrary-tree
connectivity or a linear procedure for an entire bouquet.

## The precise domain and bound

Write
\[
 \begin{split}
 W_r&=(A,D_{r-1},\ldots,D_1),\\
 A&=(x\,a_1)(01),\qquad
 D_j=(0\,a_j)(1\,b_j),\\
 a_j&=2j,\quad b_j=2j+1,\quad x=2r .
 \end{split}                                                   \tag{1}
\]
This is the canonical tuple of the
[wreath-leaf normal-form theorem](IMPRIMITIVE_LEAF_DELETION_NORMAL_FORM.md).
Let \(r\geq4\), add four new labels
\(\alpha,\beta,\gamma,\delta\), and append
\[
                F=(0\,\alpha)(\beta\,\gamma),\qquad
                B=(\alpha\,\beta)(\gamma\,\delta).             \tag{2}
\]
The D5 lobe in (2) is attached specifically at the canonical sheet 0.

**Theorem 1.** There is one fixed thirteen-sheet braid word \(w_*\),
independent of \(r\), such that the tuple \((W_r,F,B)\) can be taken
to \(W_{r+2}\), up to one simultaneous sheet relabeling, using
\[
                         |w_*|+14r-38                       \tag{3}
\]
ordinary signed Hurwitz moves. One may take
\[
                         |w_*|\leq446265984.                 \tag{4}
\]
The final relabeling fixes sheet0. The endpoint has a wreath deletion
at the new sheet \(\alpha\), which need not remain a leaf at intermediate stages.

The fixed controller exists by the retained exact finite
degree-thirteen connectivity theorem. Its literal word has not
been recovered or replayed here. Equations (3)-(4) are a certified
existence bound, rather than a claim of a short practical controller.

Throughout,
\[
 Q_i(a,b)=(b,b^{-1}ab),\qquad
 Q_i^{-1}(a,b)=(aba^{-1},a),
\]
and sequences act from left to right.

## 1. The two sweeps

For seven distinct letters, direct substitution gives
\[
 Q^{3}\bigl((u\,v)(p\,q),(v\,w)(s\,t)\bigr)
 =Q^{-3}\bigl((u\,v)(p\,q),(v\,w)(s\,t)\bigr)
 =\bigl((u\,v)(s\,t),(v\,w)(p\,q)\bigr).                    \tag{5}
\]
This is a secondary-edge exchange preserving both primary edges.
Checking (5) on one labeled seven-element set proves it for every
choice of seven distinct letters by simultaneous conjugation.

First move the secondary edge \((\beta\,\gamma)\) left through
\(D_1,\ldots,D_{r-1}\), using in that order
\[
                         Q_r^{-3},Q_{r-1}^{-3},\ldots,Q_2^{-3}.
\]
Formula (5) gives
\[
 A,\quad F_1=(0\,a_{r-1})(\beta\,\gamma),\quad
 E_{r-1},\ldots,E_1,\quad B,
\]
where
\[
 E_j=(0\,a_{j-1})(1\,b_j),\qquad a_0=\alpha.
\]
Apply \(Q_{r+1}^{3}\) to \(E_1,B\). Their new values are
\[
 F_2=(0\,\alpha)(\gamma\,\delta),\qquad
 B'=(\alpha\,\beta)(1\,b_1).
\]
Move \(F_2\)'s secondary edge left through \(E_2,\ldots,E_{r-1}\),
using \(Q_r^{-3},\ldots,Q_3^{-3}\). The resulting tuple is
\[
 A,F_1,F'_2,E'_{r-1},\ldots,E'_2,B',                         \tag{6}
\]
where
\[
 F'_2=(0\,a_{r-2})(\gamma\,\delta),\qquad
 E'_j=(0\,a_{j-2})(1\,b_j)\ (j\geq3),\qquad
 E'_2=(0\,\alpha)(1\,b_2).
\]
These two sweeps and their intervening exchange cost exactly
\(6(r-1)\) ordinary moves. The distinctness hypotheses of (5)
hold at every step.

## 2. Separating a fixed core from the outside swaps

The last three factors in (6) are \(E'_3,E'_2,B'\).
Move this three-factor block left across \(E'_{r-1},\ldots,E'_4\)
by positive moves. Each crossing costs three ordinary moves.
For \(K=E'_3E'_2B'\), the block identity is
\[
 (\eta,E'_3,E'_2,B')\longmapsto
                  (E'_3,E'_2,B',K^{-1}\eta K).
\]
The core factors remain unchanged, and
\[
                         K^{-1}(0)=a_1,\qquad K^{-1}(1)=b_3.
\]
All other labels of these outside factors are fixed by \(K\).
Thus after \(3(r-4)\) moves the tuple is
\[
       (\mathcal C_r,\ H_{r-1},\ldots,H_4),\qquad
       H_j=(a_1\,a_{j-2})(b_3\,b_j),                         \tag{7}
\]
where the six-factor core is
\[
 \mathcal C_r=\bigl(
 (01)(x\,a_1),\
 (0\,a_{r-1})(\beta\,\gamma),\
 (0\,a_{r-2})(\gamma\,\delta),\
 (0\,a_1)(1\,b_3),\
 (0\,\alpha)(1\,b_2),\
 (\alpha\,\beta)(1\,b_1)\bigr).                              \tag{8}
\]
It has exactly thirteen distinct sheets. All outside swaps in (7)
are anchored at the same pair \(\{a_1,b_3\}\), and their other
two sheets are outside the core.

Relabel the thirteen core sheets by
\[
 \begin{array}{c|ccccccccccccc}
 \text{core label}&0&1&2&3&4&5&6&7&8&9&10&11&12\\ \hline
 \text{original sheet}&0&1&a_1&b_1&a_{r-2}&b_2&a_{r-1}&b_3&
 x&\alpha&\beta&\gamma&\delta .
 \end{array}                                                \tag{9}
\]
Then (8) becomes the fixed tuple
\[
 \mathcal C=\bigl(
 (01)(28),\
 (06)(10\,11),\
 (04)(11\,12),\
 (02)(17),\
 (09)(15),\
 (9\,10)(13)\bigr),                                         \tag{10}
\]
with product
\[
                C=(0\,9\,11\,12\,10\,8\,2\,4\,6\,1\,3\,5\,7).
 \tag{11}
\]
It generates \(A_{13}\). Its product makes it transitive.
In odd degree, a nontrivial block would have size at least three,
and a double transposition cannot move two whole such blocks;
thus this transitive group is primitive. The square of the product
of the first two factors of (10) is \((0\,1\,6)\).
Jordan's 3-cycle consequence, in the retained
[Jones formulation, Theorem1.1](sources/jones-primitive-cycle-1209.5169v1.pdf),
therefore gives \(A_{13}\).

Omitting sheet9 from (11) gives the antipodal pairs
\[
             \{0,4\},\{6,11\},\{2,7\},\{1,12\},\{3,10\},\{5,8\}.
 \tag{12}
\]
In particular the anchoring pair \(\{2,7\}\) of every outside
swap is one block of the required lower wreath action.
This is an explicit check of the interface conversion, rather
than an application of abstract attachment transitivity.

## 3. A controller which retains the common root

Let \(W_6\) denote (1) on thirteen sheets, let \(P_6\) be its
product, and define the unique cyclic conjugator \(q\) by
\[
                        q(P_6^i(12))=C^i(9),\quad0\leq i<13.
\]
Its array on labels0 through12 is
\[
                    q=(0,4,11,6,12,1,10,3,8,5,2,7,9).
\]
Set \(\mathcal Z=qW_6q^{-1}\). Then
\[
 \begin{split}
 \mathcal Z=\bigl(&
 (9\,11)(04),\
 (02)(47),\
 (45)(08),\\
 &(0\,10)(34),\
 (0\,12)(14),\
 (46)(0\,11)\bigr).
 \end{split}                                                \tag{13}
\]
It has product \(C\) and a wreath deletion at sheet9, which is
\(\alpha\) in (9). Its root block is \(R=\{0,4\}\),
its attachment block is \(U=\{6,11\}\), and the outside
anchoring pair is the third block \(V=\{2,7\}\).
Its first star swap is
\[
                              J=(02)(47),                  \tag{14}
\]
swapping \(R,V\).

The choice of the leaf is important. For the canonical tuple
\(W_6\), \(P_6(0)=12\). The target leaf was chosen to be
\(9=C(0)\). Hence
\[
                       q(0)=C^{-1}(q(12))=0.               \tag{15}
\]
Thus the final canonical relabeling can retain the common root
of all still unprocessed D5 lobes.

The [complete degree-thirteen polynomial census](DEGREE_FOURTEEN_HURWITZ_ORBIT.md)
proves a single absolute \(D^6\), full-cycle orbit with
\(446265625\) vertices. Both (10) and (13) lie in it.
Consequently they are connected by an absolute path of length at
most \(446265624\).

With product \(C\) fixed, the only sheet normalizations are
powers of \(C\), since a full cycle has cyclic centralizer.
Common sheet conjugation commutes with every Hurwitz generator.
The normalizations along an absolute path therefore accumulate
into a single final conjugation by \(C^j\).
A central full twist on six factors conjugates the tuple
by \(C^{-1}\) and has length \(6\cdot5=30\).
At most twelve twists correct the final power. Thus one fixed
literal word satisfies
\[
                 w_*(\mathcal C)=\mathcal Z,\qquad
                 |w_*|\le446265624+12\cdot30=446265984.       \tag{16}
\]
The same word, embedded by (9), works in every tuple (7).
It fixes every outside factor because it uses only the first
six factor positions. No literal word for this fixed controller
is recovered here; its existence and bound follow from the
complete finite theorem, not from selected seed loops.

## 4. The final root and center adjustment

After \(w_*\), the outside factors have the form
\[
                              H_j=(2\,u_j)(7\,v_j).
\]
Each commutes with the four core swaps after \(J\) in (13):
their supports are disjoint. Move all the outside factors left
across these four swaps, at a cost \(4(r-4)\).
Then move \(J\) right past the outside factors by inverse moves:
\[
                   (J,H_j)\longmapsto(JH_jJ,J)
                       =\bigl((0\,u_j)(4\,v_j),J\bigr).      \tag{17}
\]
This costs another \(r-4\) moves.

Every factor except the first now swaps the root block
\(R=\{0,4\}\) with a distinct two-sheet block.
The first factor is \((9\,11)(04)\), and the final factor
remains \((46)(0\,11)\), swapping the root with its attachment
block \(U=\{6,11\}\).

There is an explicit final relabeling to \(W_{r+2}\).
Give sheets0,4 labels0,1 and sheet9 label \(2r+4\).
Read the remaining factors from right to left and label their
other sheet-pairs \(2j,2j+1\), orienting each pair by its
images of root labels0,1. The last factor gives sheet11 label2,
so the first factor is the required
\((2r+4\,\,2)(01)\). This relabeling fixes the original sheet0.

The total cost outside \(w_*\) is
\[
                 6(r-1)+3(r-4)+4(r-4)+(r-4)=14r-38,
\]
which proves Theorem1. \(\square\)

## Corollary 2. A quadratic upper bound for the whole bouquet

Let
\[
 \mathcal T_t=(A_0,B_0,\ldots,A_{t-1},B_{t-1}),\qquad t\geq2,
\]
where
\[
 A_i=(0\,a_i)(b_i\,c_i),\quad B_i=(a_i\,b_i)(c_i\,d_i),\qquad
 (a_i,b_i,c_i,d_i)=(4i+1,4i+2,4i+3,4i+4).
\]
This is the unbuffered family of the
[lobe distance theorem](UNBOUNDED_HURWITZ_LEAF_NORMALIZATION_DISTANCE.md).
It can be taken to \(W_{2t}\), up to one common sheet relabeling
fixing0, in at most
\[
 6656+(t-2)\,446265984+14t^2-52t+48                         \tag{18}
\]
ordinary signed Hurwitz moves. In particular all these bouquets
lie in the canonical wreath-leaf orbit, uniformly in \(t\).

**Proof.** Start with the first two lobes, on nine sheets.
They have full alternating group and full-cycle product \(P\).
For this degree-nine premise, use explicitly the retained
degree-nine row of the
[finite-census receipt](degree_fourteen_hurwitz_verification.json):
the complete character bound gives59049 literal tuples of fixed
9-cycle product, and the retained degree-nine BFS visits6561
absolute tuples. The cyclic centralizer has order9 and acts freely:
a nonidentity rotation has all its sheet orbits of length3 or9,
so the four-sheet support of a fixed DT cannot be invariant.
Thus the independent absolute total is59049/9=6561, matching
the reached orbit. The retained cyclic degree-nine seed has full
A9, as certified by its explicit3cycle and the odd-degree block
argument. Consequently this is complete single-orbit connectivity,
not merely the presence of a6561-element orbit in a larger space.
No fresh census is performed here.

Choose the conjugate of \(W_4\) with product \(P\) and target
leaf \(P(0)\). Since the canonical product satisfies
\(P_4(0)=8\), its unique cyclic conjugator has \(q(0)=0\).
As in (16), an absolute path of length at most6560 and at most
eight full twists of length12 give a fixed literal base
controller of length at most
\[
                              6560+8\cdot12=6656.
\]
Its endpoint relabeling fixes0.

All remaining lobe letters are outside these nine sheets.
Consequently this relabeling leaves the remaining lobes in their
original form, still attached at0. Apply Theorem1 to each successive
lobe. At every stage its final relabeling fixes0 and can be extended
to fix all future lobe letters, so the induction hypothesis is
retained literally.

The old factor counts at these stages are \(r=2j\),
\(2\leq j\leq t-1\). Therefore
\[
             \sum_{j=2}^{t-1}(14(2j)-38)=14t^2-52t+48.
\]
This gives (18). Intermediate common relabelings commute with
Hurwitz moves and can all be pushed to a single final common
relabeling, without changing the ordinary word length. \(\square\)

## Exact evidence and limits

The [checker](verify_d5_lobe_core_absorption.py) verifies (5),
all prefix words and outside formulas, the thirteen-sheet product,
the target conjugator fixing0, every antipodal pair, the
controller endpoints, the commuting outside moves, the literal
final inverse slide, and the explicit canonical relabeling.
One-lobe controls use every \(4\leq r\leq12\), and \(r=25,50\).
Whole-bouquet endpoint controls use \(t=2,3,4,5,8,12\).

The [receipt](d5_lobe_core_absorption_verification.json) distinguishes
these actual word replays from substitutions at the two imported
finite-controller endpoints. No new large census, replay of a
recovered fixed controller word, or optimal length claim is made.
The finite source proves that those words exist. The written
identities and root-preserving induction prove uniformity.

The one-lobe domain is canonical \(W_r\) with a lobe attached at
canonical sheet0. The whole-bouquet statement is for the displayed
unbuffered family. No arbitrary old-root absorption, general
arbitrary-tree connectivity, all-degree orbit classification,
linear total bouquet procedure, or conversion between distinct
extension interfaces by assumption is asserted.
