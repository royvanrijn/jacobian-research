# A quadratic residue model closes the duplicated-three-cycle family

For every integer \(b\ge3\), put \(m=b+2\), \(N=2m\).
The minimal full-cycle finite passport
\[
               (3^2 1^{2b-2}),\qquad (b^2 1^4),\qquad T      \tag{1}
\]
has exactly two full Hurwitz components. One has actual
monodromy \(S_N\); the other has \(S_2\wr A_m\) for odd \(b\)
and \(S_2\wr S_m\) for even \(b\). All block systems, deck
groups and separate exact orbit sizes are given in W3.

The proof supplies an explicit quadratic plane model for
the primitive component before the final finite label
quotient. A simple zero of its discriminant proves geometric
irreducibility for all \(b\). This closes the
\(\min(a,b)=3\) part of the
[first-deficit duplicated-cycle problem](POLYNOMIAL_FIRST_DEFICIT_COMPONENT_REDUCTION.md).
The \(\min(a,b)=2\) part has its separate
[rational model](POLYNOMIAL_DT_DUPLICATED_CYCLE_COMPONENTS.md).

## Primary comparison and conventions

Wajnryb, [Theorem 1, retained p2](sources/wajnryb-two-special-fibers-1995.pdf),
has a connectivity theorem with at most two nonsimple
critical values and fixed monodromy. Its definition and
precise theorem were inspected in the preceding residue
work. Here infinity and both heavy finite values are
nonsimple, so that theorem does not apply. Bogatyrëv--Gendron,
[published passport and Corollary 2, pp1485--1486](sources/bogatyrev-gendron-pell-abel-2025.pdf),
treat two finite involution factors. Neither heavy factor
in (1) is an involution. These are hypothesis comparisons,
not imported component statements. No novelty claim is made.

The inputs used are U1--U3 of the first-deficit note:
the complete imprimitive sector, the actual groups, and the
geometric-factor component criterion. Their nonemptiness
and count inputs are the corrected primary
Goulden--Jackson cactus formula already audited there.
Full braids allow every distinct factor-type placement.
Absolute classes quotient by powers of the fixed product.
Polynomial covers and geometric components are over
\(\mathbb C\); the equations below are defined over \(\mathbb Q\).

## Theorem W1. A reduced-degree residue equation for arbitrary \(a,b\)

Let \(a,b\ge2\), \(n=a-1\), \(m=a+b-1\). For the passport
\((a^2 1^{2b-2}),(b^2 1^{2a-2}),T\), normalize the first
critical pair to \(0,1\). For the other pair use
\[
 S=x+y,\quad V=xy,\quad U=S-1,\quad W=V-U,\quad
 t=W/V,\quad Q=z^a(z-1)^a(z^2-Sz+V)^b .
 \]
On the primitive locus \(V,W,U\ne0\), so \(t\ne0,1\).
Define the polynomial
\[
 H_{a,b}(S,V)=
 \sum_{\ell=0}^{n}\ \sum_{j=0}^{\lfloor\ell/2\rfloor}
 (-1)^j
 \binom{2a-\ell-2}{a-\ell-1}
 \binom{b+\ell-j-1}{\ell-j}\binom{\ell-j}{j}
 S^{\ell-2j}V^{a-\ell+j-1}.                                \tag{2}
\]
Then the top residue condition is
\[
 \Delta=(-1)^a\left(
       \frac{H_{a,b}(S,V)}{V^m}
              -\frac{H_{a,b}(2-S,W)}{W^m}\right)=0.           \tag{3}
\]
After forgetting the order of \(x,y\), the primitive equation is
\[
 F_{a,b}(V,t)=
 \frac{t^mH_{a,b}(1+V(1-t),V)
              -H_{a,b}(1-V(1-t),tV)}{t-1}=0.               \tag{4}
\]
It is a polynomial of degree exactly \(a-1\) in \(V\),
primitive over \(\mathbb C[t]\). There are no vertical
components omitted by removal of \(t-1\).
Exchanging the normalized first critical pair induces
\[
                            (t,V)\longmapsto(1/t,tV).        \tag{5}
\]
For \(a=b\), exchange of the two heavy branch colors is
also forgotten after quotienting both within-pair orderings.
Thus the remaining primitive component problem can be
studied using factors of a polynomial whose degree in \(V\)
is \(\min(a,b)-1\), with these finite label quotients.

**Proof.** The residue at zero is the coefficient of \(z^{a-1}\)
in
\[
 (-1)^a(1-z)^{-a}
                    (V-Sz+z^2)^{-b}.
 \]
Expanding the two factors and multiplying by \(V^m\)
gives exactly (2). Reflection gives the second residue
in (3). Substituting \(S=1+V(1-t)\), \(W=tV\) gives
the numerator in (4), denoted \(R(V,t)\).
It vanishes at \(t=1\) for every \(V\).

Put
\[
                  C=\binom{b+n-1}{n},\qquad
                  C'=\binom{b+n-2}{n-1}.
 \]
The first two coefficients of \(R\), as a polynomial in \(V\), are
\[
 \begin{split}
 R_0(t)&=C(t^m-1),\\
 R_1(t)&=Cn(1-t)(t^m+1)+2C'(t^m-t).
 \end{split}                                                \tag{6}
\]
At a root \(\zeta^m=1\), \(\zeta\ne1\), the second expression is
\[
                         R_1(\zeta)=2mC'(1-\zeta)\ne0.
 \]
The root \(1\) of \(R_0\) is simple. Hence the content
of \(R\) in \(\mathbb C[t]\) is exactly \(t-1\), up to a
nonzero constant. This proves primitivity of (4) and
excludes any other vertical component.

The expansion has degree at most \(n\) after substitution.
At \(t=0\) it gives
\[
                      F_{a,b}(V,0)=C(1-V)^n,
 \]
so its \(V\)-degree is exactly \(n\).
All changes of variables are birational on the indicated
configuration open, and the CRT exclusions of U3 specify
the good part of the curve. Reflection interchanges \(V,W\),
giving (5). If the profiles agree, the within-pair group
must first be quotiented before the color exchange acts;
no unsupported action on a partially labeled curve is
assumed. Exchanging \(a,b\) allows the smaller degree
in the last assertion. \(\square\)

## Theorem W2. Uniform geometric irreducibility when \(a=3\)

For every \(b\ge3\), the primitive polynomial \(F_{3,b}\)
is geometrically irreducible. Its admissible open is
nonempty, and gives one primitive full component.

An explicit quadratic model is as follows. Put
\(c=b(b+1)/2\), \(M=b+2\), and define
\[
 \begin{split}
 r_0(t)&=c(t^M-1),\\
 r_1(t)&=b\left(((b+3)-(b+1)t)t^M
                          -((b+3)t-(b+1))\right),\\
 r_2(t)&=\left(6+3b(1-t)+c(1-t)^2\right)t^M\\
       &\hspace{9mm}
                -\left(6t^2-3bt(1-t)+c(1-t)^2\right).
 \end{split}                                                \tag{7}
\]
Then
\[
             F_{3,b}(V,t)=
                  \frac{r_2(t)V^2+r_1(t)V+r_0(t)}{t-1}.
                                                               \tag{8}
\]
If \(\mathcal D_b(t)\) is its quadratic discriminant, then
\[
             \mathcal D_b(0)=0,\qquad
             \mathcal D_b'(0)=2b^2(b+1)\ne0.                 \tag{9}
\]
In particular the quadratic model has no splitting over
\(\mathbb C(t)\).

**Proof.** Formula (2) specializes to
\[
                  H_{3,b}=6V^2+3bSV+cS^2-bV.
 \]
Substitution in (4) gives (7)--(8). Before division by
\(t-1\), the coefficient values and first derivatives at
zero are
\[
 \begin{array}{c|ccc}
       &r_0&r_1&r_2\\ \hline
 t=0&-c&2c&-c\\
 \partial_t|_0&0&-b(b+3)&b(b+4).
 \end{array}                                                \tag{10}
\]
Thus the discriminant \(r_1^2-4r_0r_2\) has value zero
and derivative \(4cb=2b^2(b+1)\). Division by \((t-1)^2\)
does not change those two assertions at zero, proving (9).
A rational-function square has even valuation at every
point; the simple zero at \(t=0\) therefore proves that
the discriminant is not a square even in \(\mathbb C(t)\).
The quadratic is irreducible there, and W1's content
calculation proves geometric irreducibility of the plane
polynomial.

The point \(t=0\) is a collision boundary, not a point
claimed to belong to the Hurwitz space. It is used only
to prove the global nonsquare assertion.
For nonemptiness of the good open, U2 counts the entire
imprimitive component and gives the positive remaining
literal total
\[
                    \omega m^2\bigl(5(2b-1)-1\bigr)>0,
 \]
where \(m=b+2\) and \(\omega=3(2-\boldsymbol1_{b=3})\).
Every remaining tuple is primitive by U1 and gives a point
on (8). Hence the good locus is nonempty. A nonempty
open of this irreducible curve is irreducible. Forgetting
the remaining finite labels preserves connectedness,
so there is exactly one primitive full component.
This argument uses the independent total count, not a
finite experiment, for uniform existence. \(\square\)

## Theorem W3. Complete groups, blocks, deck and separate sizes

For every \(b\ge3\), passport (1) has exactly two absolute
and two literal full-braid orbits. Put
\[
 m=b+2,\qquad \omega=3(2-\boldsymbol1_{b=3}),\qquad
 K=\begin{cases}A_m,&b\text{ odd},\\S_m,&b\text{ even}.\end{cases}
 \]
Their invariants are
\[
 \begin{array}{c|c|c}
 &\text{primitive}&\text{imprimitive}\\ \hline
 \text{actual group}&S_{2m}&S_2\wr K\\
 \text{all block-count divisors}&\{1,2m\}&\{1,m,2m\}\\
 \lvert\mathrm{Deck}\rvert&1&2\\
 \text{literal orbit}&2\omega m^2(5b-3)&\omega m^2\\
 \text{absolute orbit}&\omega m(5b-3)&\omega m
 \end{array}.                                               \tag{11}
\]
Here \(\omega\) counts all distinct placements of the three
factor types; the equal-profile value \(b=3\) is included.

**Proof.** U2 gives exactly one imprimitive component,
with the displayed group, all blocks, deck and size.
W2 gives exactly one primitive component; U1 identifies
its group as \(S_{2m}\) and its deck as trivial.
The primitive aggregate in U2 is
\[
            \omega m^2\bigl(5(2b-1)-1\bigr)
                              =2\omega m^2(5b-3),
 \]
which is now its separate component size. Dividing
literal counts by \(2m/\lvert\mathrm{Deck}\rvert\)
gives the absolute counts. \(\square\)

## Exact support and remaining scope

Run from the canonical checkout:

    python3 research/arithmetic/verify_polynomial_duplicated_three_cycle.py

The [checker](verify_polynomial_duplicated_three_cycle.py)
compares the general residue expansion with exact local
residues, verifies the quadratic coefficient identities
and the discriminant's first jet, and takes complete
bounded colored-braid absolute closures. Both actual
groups, every block divisor and deck order are checked
against independent hook totals. The quotient alternating
group is certified by the supports of literal three-cycle
conjugates, rather than by enumerating a factorial-sized
group. The [receipt](polynomial_duplicated_three_cycle_verification.json)
and [source ledger](polynomial_duplicated_three_cycle_source_ledger.json)
record exact domains, inputs and assurance limits.

W1 gives a lower-degree algebraic route for all duplicated
cycle lengths. The cases \(\min(a,b)\le3\) now have complete
component classifications. Geometric factors and their
finite label action remain to be determined for
\(\min(a,b)\ge4\), and for the other unresolved first-deficit
passports, including higher-\(p\) two-factor profiles.
Deeper deficits and zero-simple passports remain in the
overall goal outside already resolved subfamilies.
No genus or rationality assertion is made for the
quadratic model or its final finite quotient. Independent
review, replay and formal verification remain pending.
