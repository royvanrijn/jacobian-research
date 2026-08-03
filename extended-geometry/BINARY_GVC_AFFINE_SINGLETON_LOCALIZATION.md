# Binary GVC Hall shell: affine singleton localization

## 1. Scope and outcome

At orders `p+r`, the two-digit Hall transform has an empty or singleton
high quotient on each selection side.  This note proves the exact
state-level promotion theorem for those singleton fibres.

The result has three parts.

1. Every state sharing an affine fibre with `p e_h+a` differs from it in
   one fixed integer kernel.  Its difference conformally decomposes into
   the fixed Graver basis of that kernel.  If the differences are unbounded
   as `p` varies, one fixed primitive move occurs with unbounded
   multiplicity along a subsequence.
2. If the marked column is an exposed vertex, the whole fibre stabilizes
   after the signed recentering `n=p e_h+b`.  The finite set of corrections
   `b`, their carry states, and their factorial units are independent of
   `p`.
3. Globally in \(p\), all singleton fibres form a finite module over one
   fixed affine return semigroup.  A Hilbert basis gives a finite generator
   alphabet; after any fixed additive packet projection, the Graver basis of
   the projected generator matrix gives all primitive factorization moves.

Thus an exposed singleton mark cannot wander through genuinely new
prime-dependent packets.  Every remaining unbounded singleton ambiguity is
positive-density repetition of a fixed primitive affine relation.  This is
the precise justification for reducing **state relations** to Graver moves.

It is not yet the full Hall-shell inheritance theorem.  A vanishing moment
is a linear sum of specialized state monomials; it does not canonically pair
its summands into binomial semigroup relations.  The remaining step is to
promote one of the primitive state relations below, or a complete fibre
containing it, to a pure-zero characteristic-zero packet.  Section 6 records
that distinction explicitly.

## 2. Homogenized affine fibres

Let

\[
 A=(A_1,\ldots,A_s)\in\mathbb Z^{d\times s}
\]

be an integer configuration whose first row is \((1,\ldots,1)\).  The first row
records selection order; the other rows may simultaneously record radial
coordinates, translation-source coordinates, side markings which have been
lifted to integers, or any other additive shell data.

Fix a mark \(h\) and a bounded low state \(a\in\mathbb N^s\), with
\(|a|=r\).  Put

\[
 x_p=p e_h+a
\tag{2.1}
\]

and consider its nonnegative affine fibre

\[
 \mathcal F_p=\{y\in\mathbb N^s:Ay=Ax_p\}.
\tag{2.2}
\]

Let

\[
 L=\ker_{\mathbb Z}A.
\]

For \(u,v\in\mathbb Z^s\), write \(u\sqsubseteq v\) when they lie in the
same orthant and \(|u_i|\leq |v_i|\) for every \(i\).  The Graver basis
\(\mathcal G(A)\) is the finite set of nonzero \(\sqsubseteq\)-minimal
elements of \(L\).

## 3. Singleton-fibre Graver localization

> **Theorem 3.1 (affine singleton Graver localization).**
> For every prime \(p\) and every \(y\in\mathcal F_p\), put \(z=y-x_p\).
> Then:
>
> 1. \(z\in L\), and \(z\) has a conformal decomposition
>    \[
>      z=g_1+\cdots+g_t,
>      \qquad g_j\in\mathcal G(A),\quad g_j\sqsubseteq z;
>      \tag{3.1}
>    \]
> 2. the negative support of every \(g_j\) is contained in
>    \[
>      \{h\}\mathbin\cup\operatorname {supp}(a);
>      \tag{3.2}
>    \]
> 3. for any sequence \(p_n\) and \(y_n\in\mathcal F_{p_n}\), either the
>    differences \(y_n-x_{p_n}\) range over a finite set, or, after passage
>    to a subsequence, one fixed \(g\in\mathcal G(A)\) occurs in (3.1) with
>    multiplicity tending to infinity.

### Proof

The fibre equality gives \(Az=0\).  If nonzero \(z\in L\) is not a Graver
element, there is a nonzero \(g\in L\) with \(g\sqsubsetneq z\).  Then
\(z-g\) lies in \(L\), is conformal to \(z\), and has smaller
\(\ell_1\)-norm.  Repeating
terminates and proves (3.1).

If \(i\) is outside \(\{h\}\cup\operatorname{supp}(a)\), then
\((x_p)_i=0\).  Nonnegativity of \(y\) gives \(z_i\geq0\).  A conformal
summand cannot be negative in that
coordinate, proving (3.2).

The Graver basis is finite.  If the differences are unbounded, the number
of summands in (3.1) is unbounded because every Graver element has bounded
\(\ell_1\)-norm.  Pigeonhole then gives a fixed element whose multiplicity is
unbounded along a subsequence.  If the differences are bounded, their
integer values form a finite set.  This proves the dichotomy.  \(\square\)

The theorem is insensitive to the prime.  The apparent change of low digits
with `p` is therefore not a change of the primitive relation alphabet.  It
is either a bounded signed correction or an increasing number of copies of
one of finitely many fixed moves.

### 3.1 A global Hilbert-module normal form

The subsequence dichotomy can be strengthened at the level of state sets.
Put

\[
 B=(-Ae_h\mid A)\in\mathbb Z^{d\times(s+1)}
\]

and define

\[
 \mathscr M_a=\{(p,y)\in\mathbb N^{1+s}:B(p,y)=Aa\},
 \qquad
 \mathscr H_h=\{(N,z)\in\mathbb N^{1+s}:B(N,z)=0\}.
\tag{3.3}
\]

The slice of \(\mathscr M_a\) at first coordinate \(p\) is exactly
\(\mathcal F_p\).

> **Theorem 3.2 (finite Hilbert-module normal form).**
> There are finitely many base solutions
> \(m_1,\ldots,m_c\in\mathscr M_a\) and a finite Hilbert basis
> \(h_1,\ldots,h_g\) of \(\mathscr H_h\) such that
> \[
>  \mathscr M_a=\bigcup_{i=1}^c(m_i+\mathscr H_h),
>  \qquad
>  (p,y)=m_i+\sum_{j=1}^g k_jh_j
>  \quad(k_j\in\mathbb N)
> \tag{3.4}
> \]
> for every singleton-fibre state.  Fix any integer additive packet map
> \(\pi:\mathbb Z^{1+s}\to\mathbb Z^r\), and let \(\Pi\) have columns
> \(\pi(h_j)\).  Within one translate, every equality between two Hilbert
> factorizations with the same \(\pi\)-data lies in the toric ideal of
> \(\Pi\), and every such equality is a conformal sum of its fixed Graver
> basis.

### Proof

Order \(\mathbb N^{1+s}\) componentwise.  Dickson's lemma makes the set of
componentwise-minimal members of \(\mathscr M_a\) finite; call them
\(m_1,\ldots,m_c\).  If \(x\in\mathscr M_a\) is not minimal, choose
\(x'<x\) in \(\mathscr M_a\).  Then \(x-x'\in\mathscr H_h\).  Repeating
strictly decreases the \(\ell_1\)-norm and ends at one of the \(m_i\), which
proves the first equality in (3.4).

The set \(\mathscr H_h\) is the lattice-point semigroup of the rational
polyhedral cone

\[
 \{x\in\mathbb R_{\geq0}^{1+s}:Bx=0\}.
\]

Gordan's lemma gives its finite Hilbert basis and hence the second equality
in (3.4).  Finally, two coefficient vectors \(k,k'\in\mathbb N^g\) have the
same projected packet data exactly when \(\Pi(k-k')=0\).  The standard
conformal decomposition argument from Theorem 3.1 reduces that kernel vector
to the Graver basis of \(\Pi\).  \(\square\)

The pure scaling return \((1,e_h)\) always belongs to \(\mathscr H_h\); the
other Hilbert generators are finitely many primitive affine returns.
Fixed additive marks, fixed congruence classes, and fixed carry cases can be
included by adjoining integer rows and nonnegative slack variables before
applying the theorem.

For the sharp levels \(0,1,2\), with \(h=1\) and \(a=0\), the Hilbert basis
is

\[
 h_0=(1;0,1,0),\qquad h_1=(2;1,0,1).
\tag{3.5}
\]

Indeed every state in the fibre has the unique form

\[
 (p;t,p-2t,t)=(p-2t)h_0+t h_1.
\tag{3.6}
\]

After projecting only to the order coordinate, the primitive toric move is
\(2h_0\leftrightarrow h_1\).  Taking \(t=(p-1)/2\) gives exactly the
positive-density centered family in (5.1).  Thus the Hilbert-module and
Graver descriptions meet on the first sharp obstruction.

Theorem 3.2 is a **set-theoretic** normal form.  Factorial weights are not
multiplicative in a Hilbert factorization, different translates may overlap
through additional finite module syzygies, and a specialized linear sum need
not vanish on one translate.  Thus (3.4) justifies computing Hilbert and
Graver bases once per support, but it does not prove the weighted inheritance
lemma in Section 6.

## 4. Exposed marks have stable signed digits

Delete the homogenizing coordinate from the columns and write the remaining
vectors as \(v_i\).  Suppose an integral linear form \(\ell\) exposes \(v_h\)
uniquely:

\[
 \ell(v_h)-\ell(v_i)\geq\delta>0
 \qquad(i\ne h).
\tag{4.1}
\]

> **Theorem 4.1 (exposed-singleton stabilization).**
> Under (4.1), every \(y\in\mathcal F_p\) satisfies
> \[
>   \sum_{i\ne h}y_i
>   \leq
>   \frac{D_a}{\delta},
>   \qquad
>   D_a=\sum_i a_i\bigl(\ell(v_h)-\ell(v_i)\bigr).
>   \tag{4.2}
> \]
> Consequently there is a finite set
> \[
>  \mathcal B=\{b\in\mathbb Z^s:
>     Ab=Aa,\ |b|=r,\ b_i\geq0\ (i\ne h)\}
>  \tag{4.3}
> \]
> such that, for all sufficiently large \(p\),
> \[
>   \mathcal F_p=\{p e_h+b:b\in\mathcal B\}.
>   \tag{4.4}
> \]
> In particular the marked high quotient, the signed low corrections, and
> every additive marking of those corrections are independent of \(p\).

### Proof

Apply \(\ell\) to the non-homogenizing rows of \(Ay=Ax_p\) and use
\(|y|=p+r\):

\[
 \begin{aligned}
 \sum_i y_i\bigl(\ell(v_h)-\ell(v_i)\bigr)
 &=(p+r)\ell(v_h)-
   \left(p\ell(v_h)+\sum_i a_i\ell(v_i)\right)\\
 &=D_a.
 \end{aligned}
\tag{4.5}
\]

Every off-mark summand on the left is at least \(\delta y_i\), proving
(4.2).  Put \(b=y-p e_h\).  Its off-mark entries are nonnegative and bounded,
while \(b_h=r-\sum_{i\ne h}b_i\) is bounded as well.  Thus only the finite
set (4.3) occurs.  Conversely every member of (4.3) gives a nonnegative
\(p e_h+b\) once \(p\) is sufficiently large, proving (4.4).  \(\square\)

The same proof with a nontrivial exposed face bounds only the mass outside
that face.  Positive-density motion can then occur inside the face, which is
exactly where another face descent or a primitive return packet is needed.

### 4.1 Exact carry and unit

Write \(b_h=s\), put \(b_i\geq0\) off the mark, and assume \(p\) is larger
than all bounded entries.  Since \(|b|=r\), the state \(p e_h+b\) has
multinomial

\[
 M_p(b)=\frac{(p+r)!}{(p+s)!\prod_{i\ne h}b_i!}.
\tag{4.6}
\]

If \(s\geq0\), it has no carry and

\[
 M_p(b)\equiv
 \frac{r!}{s!\prod_{i\ne h}b_i!}
 \pmod p.
\tag{4.7}
\]

If \(s=-t<0\), it has exactly one carry and Wilson's theorem gives

\[
 \frac{M_p(b)}p\equiv
 (-1)^{t-1}
 \frac{r!(t-1)!}{\prod_{i\ne h}b_i!}
 \pmod p.
\tag{4.8}
\]

Thus crossing the standard digit boundary changes only a fixed carry bit
and a fixed signed factorial unit.  It does not change the underlying
singleton mark or create a new prime-dependent packet.

### 4.2 A long radial-carry interval is one Hasse row

There is a second binary simplification which applies before any individual
carry positions are separated.  Fix \(0\leq R<p-1\).  If the two low radial
digits add to \(p+R\), write them as

\[
 t,\quad p+R-t,\qquad R<t<p.
\]

Wilson's theorem gives the pointwise identity

\[
 \begin{aligned}
 t!(p+R-t)!
 &\equiv
 (-1)^{t-R}\frac{t!}{(t-R-1)!}\\
 &=(-1)^{t-R}t^{\underline{R+1}}
 \pmod p.
 \end{aligned}
\tag{4.9}
\]

Consequently, for

\[
 F(T)=\sum_{t=R+1}^{p-1}c_tT^t,
\]

the complete carry interval satisfies

\[
 \boxed{
 \sum_{t=R+1}^{p-1}c_t\,t!(p+R-t)!
 \equiv-F^{(R+1)}(-1)\pmod p.}
\tag{4.10}
\]

Indeed
\(F^{(R+1)}(-1)=\sum_t c_t t^{\underline{R+1}}
(-1)^{t-R-1}\).  If the two residues instead add to \(R\), there are only
the fixed \(R+1\) positions \(0,\ldots,R\), with kernel \(t!(R-t)!\).

Thus an unbounded binary radial-carry interval is not an unbounded family of
independent observations.  After the high radial quotient and bounded
residue are fixed, it is one derivative (equivalently one Hasse row, up to
the invertible factor \((R+1)!\)).  Any remaining prime dependence lies in
the coefficient array \((c_t)\), i.e. in the Cartier/selection fibre, not in
the radial factorial kernel.

## 5. Sharpness at an interior mark

Take the one-dimensional configuration with levels `0,1,2` and the interior
mark `h=1`.  For every odd prime `p`, set

\[
 x_p=(0,p,0),\qquad
 y_p=\left(\frac{p-1}{2},1,\frac{p-1}{2}\right).
\tag{5.1}
\]

Both states have mass `p` and total level `p`, but the first has a singleton
high digit and the second has an empty high digit.  Their difference is

\[
 y_p-x_p=\frac{p-1}{2}(1,-2,1),
\tag{5.2}
\]

an unbounded repetition of the centered-triple Graver move.  Hence the
exposed-vertex hypothesis in Theorem 4.1 is necessary.

There is also an exact two-side warning.  The packets `(x_p,y_p)` and
`(y_p,x_p)` have the same radial factorial and the same product of their two
multinomial weights at every odd prime.  Their carry vectors are `(0,1)` and
`(1,0)`.  Scalar data after side exchange cannot select one of them.  This is
not a new obstruction: (5.2) shows that the whole family is made from the
already known centered-triple atom.  It does show why state-level Graver
localization must not be silently promoted to a vanishing packet identity.

### 5.1 The complete primitive three-level ghost

The sharp family has an exact extension to every primitive triple of
collinear levels.  Let \(u,v\geq1\) be coprime, put \(s=u+v\), and use the
levels \(0,v,s\), with the middle level marked.  The nonnegative states of
mass \(p\) and total level \(pv\) are exactly

\[
 n_t=(ut,p-st,vt),
 \qquad 0\leq t\leq\left\lfloor\frac p s\right\rfloor.
\tag{5.3}
\]

Give the three channels coefficients \(A,B,C\), and set

\[
 X=\frac{A^uC^v}{B^s}.
\]

For a prime \(p>s\), delete the pure singleton \(t=0\), divide the remaining
coefficient by \(pB^p\), and reduce modulo \(p\).  The complete affine ghost
is

\[
 G_{p;u,v}(X)=
 \sum_{t=1}^{\lfloor p/s\rfloor}
 \frac{(p-1)!}{(ut)!(vt)!(p-st)!}X^t
 \in\mathbb F_p[X].
\tag{5.4}
\]

Since \(st<p\), Wilson's theorem gives the termwise form

\[
 \frac{(p-1)!}{(ut)!(vt)!(p-st)!}
 \equiv
 \frac{(-1)^{st-1}}{st}\binom{st}{ut}
 \pmod p.
\tag{5.5}
\]

This is a truncated logarithmic constant term.  Normalize
\(a=A/B\), \(c=C/B\), so that \(X=a^uc^v\), and put

\[
 f(z)=1+a z^{-v}+cz^u.
\]

Then, as a formal series in \(X\),

\[
 \operatorname{CT}_z\log f(z)
 =\sum_{t\geq1}
   \frac{(-1)^{st-1}}{st}\binom{st}{ut}X^t,
\tag{5.6}
\]

and equivalently

\[
 G_{p;u,v}(X)
 \equiv\frac{\operatorname{CT}_z f(z)^p-1}{p}\pmod p.
\tag{5.7}
\]

For \((u,v)=(1,1)\), the non-support root \(X=1\) is the familiar centered
triple.  The exact bounded census in Section 7 tests every coprime
\(1\leq u\leq v\leq6\), every rational number of numerator at most \(40\)
and denominator at most \(20\), the primes through \(43\), and every
root-of-unity order through \(80\).  It also excludes all 2,139 primitive
irreducible minimal-polynomial candidates of degrees two and three and
coefficient height at most four.  It finds only \(X=1\) for the centered
triple, besides the support root \(X=0\).  This is bounded evidence, not an
all-prime classification of algebraic roots.  Formula (5.7) also explains
why such a classification is a genuinely arithmetic problem rather than a
formal consequence of Graver localization.

> **Theorem 5.2 (all-width cyclotomic separation).**
> Let \(u,v\geq1\) be coprime and let \(s=u+v\).  Suppose a root of unity
> \(\zeta_n\) is a zero of \(G_{p;u,v}\) after reduction at every prime of
> \(\mathbb Q(\zeta_n)\) above every rational prime
> \(p>\max\{s,3\}\) with \(p\nmid n\).  Then
> \[
>  (u,v,\zeta_n)=(1,1,1).
> \tag{5.8}
> \]
>
> Equivalently, among everywhere-good Galois-stable torsion candidates, every
> primitive three-level affine ghost is character-separated except for the
> already-safe centered root.

### Proof

List the primes larger than \(\max\{s,3\}\) as
\(q_1<q_2<\cdots\), and let \(q_j\) be the first one not dividing \(n\).
Because \(q_j\nmid n\), the reduction of \(\Phi_n\) modulo \(q_j\) is
separable.  Vanishing at every prime above \(q_j\) therefore implies

\[
 \Phi_n(X)\bmod q_j\quad\text{divides}\quad G_{q_j;u,v}(X).
\tag{5.9}
\]

If \(j=1\) and \(s\geq3\), Bertrand's postulate gives \(s<q_1<2s\).
Thus (5.4) has only its \(t=1\) term and is a nonzero multiple of \(X\),
which has no cyclotomic factor.  If \(s=2\), coprimality forces
\((u,v)=(1,1)\), \(q_1=5\), and

\[
 G_{5;1,1}(X)=X(X-1).
\tag{5.10}
\]

Since \(5\nmid n\), the only cyclotomic polynomial dividing (5.10) is
\(\Phi_1=X-1\).  This gives the centered root in (5.8).

It remains to exclude \(j\geq2\).  Then every \(q_i\) with \(i<j\) divides
\(n\), and hence

\[
 \varphi(n)\geq\prod_{i<j}(q_i-1).
\tag{5.11}
\]

For \(j=2\), Bertrand gives \(q_2<2q_1<4s\) when \(s\geq3\), while
\(q_1-1\geq4\); for \(s=2\), the pair is \((q_1,q_2)=(5,7)\).
In both cases

\[
 q_1-1>\left\lfloor\frac{q_2}{s}\right\rfloor.
\]

The inequality propagates: Bertrand gives \(q_{i+1}<2q_i\), while
multiplying the left side by \(q_i-1>2\) more than compensates for this
factor of two.  Formula (5.5) shows that the last coefficient in (5.4) is
nonzero, since \(s\lfloor q_j/s\rfloor<q_j\).  Therefore

\[
 \varphi(n)>\left\lfloor\frac{q_j}{s}\right\rfloor
 =\deg G_{q_j;u,v}.
\tag{5.12}
\]

This contradicts the divisibility (5.9), and proves the theorem.
\(\square\)

The Galois-stability clause is the one supplied by an algebraic
characteristic-zero packet: applying all embeddings and reducing at every
good prime above \(p\) forces the full separable minimal polynomial to divide
the ghost.  The displayed theorem assumes that no further rational primes are
discarded.  It does not by itself show that an arbitrary finite exceptional
set in a spread-out Hall packet is harmless, classify a general non-torsion
algebraic cross-ratio, or inherit the complete three-level ghost from the
linear Hall shell.

## 6. What this closes, and the remaining lemma

Theorems 3.1--4.1 close five pieces of the \(p+r\) Hall-shell problem.

1. **No new primitive move can appear indefinitely with the prime.**  The
   affine fibre uses one fixed Graver basis, and every unbounded branch repeats
   one member of it.
2. **The full state family has one finite semigroup alphabet.**  Theorem 3.2
   writes every affine singleton state in a finite module over a fixed
   Hilbert-basis return semigroup; its nonfree factorizations have a fixed
   Graver basis after any fixed additive packet projection.
3. **Exposed singleton marks are promoted.**  After signed recentering, their
   complete correction alphabet, carries, units, and markings are fixed.
4. **The genuinely moving case has positive density.**  It lies in a
   nontrivial convex face and contains an unbounded fixed primitive move.  In
   the one-dimensional level model the first such move is exactly the
   centered triple; (5.4) gives the complete ghost for every primitive
   three-level face, and Theorem 5.2 removes every new everywhere-good
   torsion survivor.
5. **Individual radial carry positions compress to one bounded-order row.**
   Formula (4.10) replaces the full interval by a derivative at \(-1\); the
   unresolved object is the selection-fibre coefficient polynomial on which
   that row acts.

Together with the all-span consecutive-residue theorem, this removes
prime-dependent *state combinatorics* as an independent mystery.  It was
tempting to ask for the following coefficient-blind linear-to-toric
implication:

> **Insufficient module-only inheritance statement.**  Write the complete
> affine state set in the form (3.4).  If its factorially weighted linear shell
> vanishes, then either a whole module translate/profile vanishes, or a
> conformal Graver block for one fixed additive projection, with a common high
> quotient, inherits a pure-zero identity (up to the already-safe
> beta/centered decompositions).

This statement is false under the displayed module hypotheses alone.
[Proposition 2.1 of the translation-tangent note](BINARY_GVC_TRANSLATION_TANGENT_RIGIDITY.md)
takes two free module translates, multiplies both by the common factorial
\((N!)^2\), and uses the isoperiodic pair
\(z+z^{-1},z^2+z^{-2}\).  Their complete pure shell vanishes at every scale,
neither translate vanishes, both internal Graver bases are empty, and a fixed
Laurent multiplier gives a nonzero odd tail.  This is not a rank-one Cartesian
Hall packet or a GVC counterexample; it proves that a module idempotent cannot
follow from (3.4) and factorial weighting alone.

The corrected remaining target must use the Cartesian Taylor down-set and its
adjacent-channel derivative identities.  The same translation-tangent note
proves that, for a translated monomial at primitive slope, the complete
linearized diagonal-period kernel is exactly the flat torus direction, and it
deduces exact flatness of every \(q^a\)-order character collision once the
underlying prime \(q\) is sufficiently large.  What remains is to expose a
nonflat exceptional-small-prime, mixed-prime, or two-dimensional curvature
row over one common high quotient, or to construct a nonterminal rank-one
Cartesian Hall packet where that fails.  Prime-power
signatures and Graver bases solve the packet after this Hall-specific
inheritance; they do not construct it from a specialized linear shell.

### 6.1 The multiplier can be reduced to shifted powers

There is no need to carry an arbitrary multiplier \(Q\) through the final
Hall argument.  De Bondt's
[*A few remarks on the Generalized Vanishing
Conjecture*](https://arxiv.org/abs/1206.2836), Theorem 1 and Corollary 2,
gives the following equivalence for a fixed \((\Lambda,P)\):

\[
\begin{aligned}
&\Lambda^m(QP^m)=0
   &&\text{for every }Q\text{ and all }m\gg_Q 0,\\
\Longleftrightarrow\quad
&\Lambda^m(P^{m+d})=0
   &&\text{for every fixed }d\geq1\text{ and all }m\gg_d 0.
\end{aligned}
\tag{6.1}
\]

For completeness, the input is de Bondt's commutator lemma

\[
 \Lambda^{M-D}\widetilde f=0,\quad \deg Q\leq D
 \quad\Longrightarrow\quad
 \Lambda^M(Q\widetilde f)=0.
\tag{6.2}
\]

Its proof filters the commutator \([\Lambda,Q]\widetilde f\) by the degree of
the polynomial multiplying a constant-coefficient derivative of
\(\widetilde f\), then inducts on \((M,D)\); commuting one copy of
\(\Lambda\) lowers that multiplier degree.  If the shifted-power statement
in (6.1) holds, take \(D\geq\deg Q\).  With \(n=m+D\), it gives
\(\Lambda^{n-D}(P^n)=0\) for all large \(n\); (6.2), with
\(\widetilde f=P^n\), gives \(\Lambda^n(QP^n)=0\).  The converse is obtained
by taking \(Q=P^d\).

This reduction removes arbitrary affine multiplier weights from the theorem
target.  In a promoted finite trace the extra factor is the complete power
\(P^d\), so it shifts the common exponent rather than inserting an unrelated
character weight.  It does **not** prove promotion of the pure shell: the
premise \(\Lambda^m(P^m)=0\) must still yield a fixed pure-zero packet before
(6.1) can be used.

Accordingly, the narrowest remaining binary statement is:

> Pure Hall-shell vanishing, together with the full Cartesian translation
> tower, either gives a flat scale-compatible twist or exposes a nonflat
> curvature block with one common high quotient; after that exposure, prove
> the shifted-power rows of that fixed packet.

The second clause is covered by the existing finite-trace and Laurent
separator machinery.  The primitive one-direction tangent and
large-prime-power torsion cases of the first clause are now proved in the
translation-tangent note.  Exceptional small primes, mixed-prime torsion,
two-dimensional Taylor down-sets, and their factorial-compatible Hall
exposure remain open.

## 7. Exact replay

Run

```bash
.venv/bin/python scripts/research_binary_gvc_frobenius_carry.py \
  --radial-limit 3 --order-limit 2 --bridge-limit 29 --residue-limit 2
.venv/bin/python scripts/research_binary_gvc_ghost_shell.py
```

The carry check enumerates endpoint fibres at two successive
large primes, verifies stabilization and (4.7)--(4.8), replays (5.1)--(5.2)
through every odd prime at most `29`, and verifies (4.9)--(4.10).  The ghost
check verifies (5.3)--(5.7), searches the stated rational window, and tests
cyclotomic divisibility through order `80` and primitive irreducibles through
degree three and height four.  These computations are
regressions for the exact identities above; their finite survivor census is
not a proof beyond its stated bounds.
