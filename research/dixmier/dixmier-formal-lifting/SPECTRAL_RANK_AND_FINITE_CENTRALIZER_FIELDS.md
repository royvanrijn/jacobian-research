# Spectral rank and finite actual centralizer fields

**QSF1.** The actual joint centralizer of an algebraically independent
commuting m-tuple in a scalar rank-m Weyl algebra has finite fraction-field
degree over its position field. No finite generation of the centralizer
and no individual conjugates are needed. This is an application of
Braverman--Etingof--Gaitsgory's quantum-integrable-system theorem, rather
than a new centralizer-finiteness theorem.

In the Long deformation setting of QTB1, actual individual conjugates
therefore leave only a finite non-Galois extension of degree at least
three, or the polynomial centralizer. A nonzero spectral rank one or two
forces the polynomial alternative and finite full completion; a finite set
of positive ranks with greatest common divisor one or two also suffices.
No spectral rank or individual conjugates with these properties are supplied for a
new Long pair.

## The actual scalar field theorem

Let L be a characteristic-zero field, let A=A_m(L), and choose commuting
algebraically independent elements F_1,...,F_m in A. Put
\[
 P=L[F_1,\ldots,F_m],\quad Q=\operatorname{Frac}P,\quad
 B=\operatorname{Cent}_A(F_1,\ldots,F_m).
 \tag{1}
\]
Then B is commutative and
\[
 S^{-1}B=H=\operatorname{Frac}B,\qquad
 S=P\setminus\{0\},\qquad [H:Q]<\infty.
 \tag{2}
\]
The localization in (2) is central inside B. It is not a claim that
A is module-finite over P, that B is module-finite over P, or that
B is finitely generated as an L-algebra.

The precise external input is
[Braverman--Etingof--Gaitsgory, arXiv:alg-geom/9607012v4](https://arxiv.org/pdf/alg-geom/9607012v4),
Section 1.2(3), Section 2.2 (Lemma, Theorem and Corollary), with the
holonomicity argument in Appendix A.1--A.3 and the commutativity proof
in Sections 6.1--6.3. Printed pages 4, 7, 17--19 and 20--21 contain
these statements and arguments. The source assumes
algebraically closed constants; descent to L is justified below.
The cached PDF and exact source pin are recorded in
[the source audit](spectral_centralizer_source_audit.json).

Here is the application with its module structure made explicit.
Give A a commuting right P action by the F_i, and form
\[
 M=A\otimes_P Q
 \quad\hbox{as a left }A_Q=A\otimes_L Q\hbox{-module}.       \tag{3}
\]
The external Q action is right multiplication by the rational functions
of F; it is not their original left action on A. The module is cyclic:
the external denominators make every a\otimes q an A_Q multiple of 1.
It is nonzero, since A is torsion-free as a right P module.

For completeness, the dimension argument behind the cited holonomicity
does not require any independence of leading symbols. Regard A as a
module over A\otimes_L P by left A and right P multiplication.
The filtered-dimension inequalities in the cited Appendix give
\[
 2m=\dim_A A\ \ge\ \dim_{A\otimes P}A\
       \ge\ \dim_{A_Q}M+m.
 \tag{4}
\]
Bernstein's inequality supplies the opposite bound
\(\dim_{A_Q}M\ge m\), so M is holonomic. Finite-dimensionality of
endomorphisms of a holonomic algebraic D-module is the remaining
external D-module input used in Section 2.2.

An endomorphism of M is determined by its value on 1. Embed M in the
Weyl division ring by a\otimes s^{-1}\mapsto a\,s(F)^{-1}.
The equations imposed on that value are
\[
 (F_i-u_i)v=0,\qquad u_i\hbox{ acts by right }F_i,
\]
which say exactly that its numerator commutes with every F_i.
Consequently
\[
 \operatorname{End}_{A_Q}(M)\simeq(S^{-1}B)^{\mathrm{op}}.
 \tag{5}
\]
The cited centralizer theorem makes B commutative. A finite-dimensional
commutative Q-domain is a field, giving (2). In the formal Long setting
commutativity also follows directly from QSC1's classical centralizer
calculation, so the extra general commutativity theorem is unnecessary
there.

The cited commutativity argument can also be seen directly from the
finite-dimensional algebra in (5). Each b in B is algebraic over the
central subalgebra P. Taking leading terms of such a relation makes
gr(B) algebraic over gr(P); these are domains with the induced Weyl
filtration. If B were noncommutative, take its smallest nonzero
commutator degree drop to obtain a nonzero Poisson bracket on gr(B).
It annihilates the central gr(P). Differentiating a minimal algebraic
equation in characteristic zero then forces it to annihilate every
element of gr(B), a contradiction. This is the argument of the cited
Sections 6.1--6.3, with its algebraicity and domain hypotheses verified.

To descend the cited statements, extend constants from L to an algebraic
closure \(\overline L\). The centralizer is the intersection of the
kernels of the finitely many L-linear maps \(\operatorname{ad}F_i\).
A field extension preserves those kernels and their finite intersection:
\[
 B\otimes_L\overline L
 =\operatorname{Cent}_{A_m(\overline L)}(F_1,\ldots,F_m).
\]
The corresponding localization is
\((S^{-1}B)\otimes_Q\overline L(u_1,\ldots,u_m)\).
Its finite dimension over the larger field implies the same finite
dimension over Q by faithful flatness; commutativity descends as well.
Algebraic independence survives this extension. Thus every hypothesis
of the algebraically closed theorem has been checked.

## The field degree divides every nonzero spectral rank

Fix a Weyl polarization \(A=L[x]\langle p\rangle\) and put
\[
 V=Q(x_1,\ldots,x_m)\otimes_{Q[x_1,\ldots,x_m]}M,\qquad
 r=\dim_{Q(x)}V.
 \tag{6}
\]
Holonomicity makes r finite, although it can be zero. It is the generic
rank in the original source positions, not a Weyl-division-ring index.

If r>0 and d=[H:Q], then
\[
                              d\mid r.                    \tag{7}
\]
Indeed the endomorphism action (5) commutes with the left source
coordinates and hence acts on V. Extending this action by Q(x) gives
a unital action of
\[
 H\otimes_Q Q(x)=H(x).
\]
The tensor product is a field of degree d over Q(x), since the latter
is purely transcendental over Q. A nonzero vector space with this
field action has Q(x)-dimension a positive multiple of d. This proves
(7), without assuming that B is finitely generated or constructing a
maximal spectral variety. Compare the source's Section 2.3, Corollary 2
and its proof.

The condition r>0 is essential. For m=1 and F=x, the module (3) has
relation (x-u)1=0. Source-coordinate localization inverts the nonzero
polynomial x-u and gives V=0, although H=Q. Fourier transforming this
same example to F=p gives the rank-one equation (p-u)1=0 and r=1.
Thus rank zero says nothing about d, and rank depends on polarization.
Every exact finite Weyl automorphism preserves the field degree d, but
can change r.

For the control F=p^k in rank one, H=L(p), Q=L(p^k), and d=k.
The monic equation (p^k-u)1=0 gives the spectral rank r=k, showing
that nontrivial field degree and nonzero rank can agree in every
positive degree. These are all-degree algebraic controls, not a
finite numerical scan.

## Low spectral rank is a finite-completion certificate

Now retain the deformation hypotheses of
[QTB1](TAYLOR_BRANCHES_AND_CENTRALIZER_RIGIDITY.md):
the actual positions belong to the finite-PBW h-Weyl algebra over
K[[h]], have independent classical reductions f_i with full polynomial
Poisson centralizer K[f], and admit actual generic individual conjugates.
Set L=K((h)). QET1 removes their parameter poles, and QTB1 excludes
every nontrivial finite Galois subfield of H/Q.

Together with (2), this gives the exact alternatives
\[
 H=Q,\quad\hbox{or}\quad 3\le[H:Q]<\infty
                         \hbox{ and }H/Q\hbox{ is non-Galois}. \tag{8}
\]
The previously permitted infinite algebraic centralizer field cannot
occur for these actual scalar Weyl positions. This strengthens only
the field conclusion; it does not make the centralizer module-finite
or finitely generated.

Suppose there is a polarization, obtainable by an exact finite Weyl
coordinate change if needed, in which 1<=r<=2. Equation (7) gives d<=2.
The quadratic possibility can be excluded directly, without the
characteristic-zero affine-space covering theorem or branch purity.

If H=Q(theta), theta^2=a in Q, QTB1's Taylor-branch argument gives,
at each irreducible base divisor pi, an extension valuation with value
group Z and v(pi)=1. Hence
\[
                     2v(\theta)=v_\pi(a)
\]
makes every irreducible-divisor valuation of a even. Factor numerator
and denominator in the UFD L[t_1,...,t_m]; this gives a=c*b^2 for
c in L^times and b in Q^times. But theta/b is algebraic over L,
so all spectral derivations kill it. QTB1's joint-constant calculation
puts it in L, making a a square in Q. Thus d=2 is impossible, and H=Q.

QTB1's differential-overring localization theorem and the scalar-unit
condition now give B=P on the generic fiber. The valuation descent
and explicit polynomial curvature homotopy in QFC1 give a full finite
integral completion directly; that application does not use QFC1's
finite-etale-cover argument. If the prescribed classical individual
conjugates already form a flat tuple, the corrections may preserve
that tuple, using QET1's classical-shift adjustment. Therefore this
low-rank completion criterion needs QTB1's Taylor/flatness inputs,
but not its full Galois-factor covering theorem.

There is a stronger certificate using several polarizations. Let
sigma_1,...,sigma_l be exact finite Weyl automorphisms of A over L,
and compute positive spectral ranks r_j for the transformed tuples
sigma_j(F). Their centralizer field extensions are all isomorphic
to H/Q, so
\[
             d\mid\gcd(r_1,\ldots,r_l).                   \tag{9}
\]
If this greatest common divisor is one or two, the same argument
forces H=Q and finite full completion. For example, ranks (3,4)
give gcd one, while ranks (4,6) give gcd two; neither certificate
requires any individual rank to be at most two. Gcd one even avoids
the quadratic exclusion. The automorphisms used for these rank
measurements need not preserve the classical graph or be h-regular:
the actual field isomorphism over L transports d, and QTB1's
constant and divisor restrictions are applied to the original pair.

Thus a finite collection of independently certified positive ranks can
be more effective than minimizing a single rank. This is still a
conditional search criterion. No such rank collection is computed for
a new Long pair here.

More generally, a nonpolynomial survivor under these hypotheses must
have a finite non-Galois field degree d>=3 dividing every certified
positive spectral rank in every such polarization. QTB1's divisor
inertia and normal-stabilizer restrictions continue to apply to this
single finite extension.

This certificate needs both actual individual conjugates and one or
several rigorously certified positive ranks. Rank information alone does not
supply the individual conjugates. For example F=p^2 has r=d=2 but
no polynomial individual conjugate: if [p^2,T]=1 then expanding
T=\sum_j a_j(p)x^j and comparing the largest x power forces j=1
and 2p*a_1(p)=1, impossible in L[p]. Its quadratic Galois field
also fails QTB1's necessary condition. There is no contradiction.

## Remaining Long construction gate and evidence scope

For actual bounded commuting Long positions, field degree is therefore
automatically finite. A computation of their full spectral module can
now supply a finite field-degree bound or a positive-completion
certificate from positive ranks with gcd one or two. No particular
commuting Long deformation with such a rank collection, actual individual mates, or finite full
completion has been constructed here. The theorem gives no uniform
spatial-degree bound and does not identify the centralizer field with
the completed division degree three of QDI1.

The cached primary source was retrieved and its stated sections and
Appendix proof inspected. The identification (5), constant-field
descent, divisibility argument and low-rank application are written
arguments; the examples are direct algebraic controls. No numerical
experiment or checker is offered as proof of those universal statements.
Finite generation, a new Long radical certificate, a bounded Long lift,
DC1/DC2, independent theorem replay, formal verification, external
review and novelty priority remain unasserted.
