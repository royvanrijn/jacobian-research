# Two rational components for a double transposition and a duplicated cycle

For every integer \(b\ge2\), put \(m=b+1\), \(N=2m\).
Consider minimal full-cycle factorizations with finite passport
\[
             (2^2 1^{2b-2}),\qquad(b^2 1^2),\qquad T.        \tag{1}
\]
There are exactly two full Hurwitz components, one with actual
group \(S_N\), the other with actual group \(S_2\wr S_m\).
Their affine-quotient moduli curves are rational over
\(\mathbb Q\), and explicit rational coefficient charts and
uniform admissible rational points are given below.

This evaluates the
[first-deficit determinant criterion](POLYNOMIAL_FIRST_DEFICIT_COMPONENT_REDUCTION.md)
for one entire infinite higher-index family. It does not
evaluate the determinant for arbitrary repeated cycle lengths
\(a,b\ge3\), or for the remaining first-deficit passports.

## Primary comparison and conventions

Bogatyrëv--Gendron, [published Corollary 2 and its passport,
pp1485--1486](sources/bogatyrev-gendron-pell-abel-2025.pdf),
already classify the \(b=2\) case through their two-involution
family. The actual statement and passport were inspected.
Their two distinguished finite factors are products of
transpositions, so that result does not directly apply to
\(b>2\). No novelty claim is made.

Wajnryb, [Theorem 1, retained preprint p2](sources/wajnryb-two-special-fibers-1995.pdf),
classifies covers with at most two nonsimple critical values,
fixed branching data and fixed monodromy. His actual definition
and theorem were inspected. In (1), infinity and both heavy
finite values are nonsimple, so there are three such values.
That theorem is comparison, not a component input here.

The proof uses the residue criterion and monodromy/count
results U1--U3 of the first-deficit note. Full braids allow
all distinct factor-type placements. Absolute classes quotient
by simultaneous powers of the fixed product \(P_N\).
Critical-point labels used in coefficient charts are forgotten
before counting these components.

## Theorem V1. The primitive residue curve has a rational quotient

Normalize the two critical points over the double-transposition
value to \(0,1\). Write the critical points over the other
heavy value as \(x,y\), and put
\[
 S=x+y,\qquad V=xy,\qquad U=S-1,\qquad W=V-U,\qquad
 Q=z^2(z-1)^2(z^2-Sz+V)^b.
 \]
The configuration requires \(V,W\ne0\), \(x\ne y\), and
\(x,y\notin\{0,1\}\). The top residue condition is
\[
 \Delta=\operatorname{res}_0\frac{dz}{Q}
              +\operatorname{res}_1\frac{dz}{Q}
       =\frac{2V+bS}{V^{b+1}}
             -\frac{2V-(b+2)U+b}{W^{b+1}}=0.                \tag{2}
\]
The imprimitive locus is exactly \(U=0\).
On its complement put \(t=W/V\), so that \(t\ne0,1\).
Define
\[
 \begin{split}
 B_b(t)&=\sum_{j=0}^{b}t^j,\\
 A_b(t)&=\sum_{j=0}^{b-1}(j+1)(b-j)t^j .
 \end{split}
 \]
After quotienting by interchange of \(x,y\), every primitive
configuration satisfying (2) lies on the single rational graph
\[
 V(t)=\frac{bB_b(t)}{(t-1)^2A_b(t)},\qquad
 U(t)=-\frac{bB_b(t)}{(t-1)A_b(t)},\qquad S(t)=1+U(t).        \tag{3}
\]
The exact admissible locus is a nonempty open subset of
this graph, as proved in V2.

**Proof.** Differentiate at the two double poles. At zero,
the residue is \((2V+bS)/V^{b+1}\). Reflecting \(z\mapsto1-z\)
gives the second term in (2). U1--U2 show that every
imprimitive cover has a quadratic deck involution. It swaps
\(0,1\), hence is \(z\mapsto1-z\); it also swaps \(x,y\).
Thus \(S=1\), or \(U=0\). Conversely at \(U=0\) the two
critical pairs have that common involution, so their unique
CRT polynomial modulo constants is invariant under it.
Every admissible such polynomial is imprimitive.

Away from \(U=0\), substitution \(W=tV,\ U=V(1-t)\) in (2)
gives an equation linear in \(V\):
\[
 D_b(t)V=b(1-t^{b+1}),\qquad
 D_b(t)=(b+2-bt)t^{b+1}+b-(b+2)t.                           \tag{4}
\]
The elementary polynomial identities
\[
 D_b=-(t-1)^3A_b,\qquad
 1-t^{b+1}=-(t-1)B_b
 \]
give (3). There is no omitted vertical component at \(t\ne1\):
if \(t^{b+1}=1\), then
\[
                         D_b=2(b+1)(1-t).
 \]
Consequently the two coefficients in the linear equation
have gcd exactly \(t-1\), even over \(\mathbb C\).
After removal of that factor, (4) is geometrically irreducible
in the \((t,V)\)-plane, by the degree-one equation and Gauss's
lemma. The change of variables is birational on \(V\ne0,t\ne1\).
It proves the rational quotient assertion.

This argument does not claim that the equation in the
ordered pair \(x,y\) is irreducible: any possible splitting
of that quadratic cover is merged when this pair is
forgotten. The quotient itself is the irreducible graph (3).
\(\square\)

## Lemma V2. Exact coefficient charts and uniform nonemptiness

For either \(S,V\) on (3) or \(S=1\), put
\[
 P=z(z-1),\quad B=z^2-Sz+V,\quad
 E=P^2\bigl((P^2)^{-1}\bmod B^b\bigr).
 \]
Assume \(P,B\) are coprime. Then \(E\) is the unique polynomial
of degree less than \(2b+4\) with
\[
                         E=0\bmod P^2,\qquad E=1\bmod B^b.
 \]
Condition (2) makes its coefficient of degree \(N+1\) zero.
On the open set \(L=[z^N]E\ne0\), set
\[
 f=E/L,\qquad v=1/L,\qquad
 A=\frac{f'}{NP B^{b-1}}=z-\alpha,\qquad w=f(\alpha).
                                                               \tag{5}
\]
The necessary and sufficient remaining good conditions are
\[
 \operatorname{disc}(B)\ne0,\qquad
 \gcd(A,PB)=1,\qquad w\ne0,\quad w\ne v.                    \tag{6}
\]
Together with coprimality and \(L\ne0\), they give the exact
passport (1), and every polynomial in the two charts arises
this way after the stated affine normalizations.

Both good loci are nonempty for every \(b\ge2\), with
explicit rational points:
\[
 \begin{array}{c|c|c}
 &S&V\\ \hline
 \text{primitive}&1-v_b&v_b\\
 \text{imprimitive}&1&2
 \end{array},
 \qquad
 v_b=\frac{b(2^{b+1}-1)}{(b-2)2^{b+1}+b+4}.
                                                               \tag{7}
\]
The primitive point is \(t=2\) in (3).

**Proof.** CRT gives \(E\), and the top coefficient vanishes
exactly when (2) holds, by the residue argument in U3.
The congruences imply \(PB^{b-1}\mid f'\). Its degree is
\(2b=N-2\), so the quotient \(A\) is monic linear.
Conditions (6) give the exact local orders, the sole
remaining simple critical point, and distinct critical
values. All derivative zeros have been accounted for.
This also proves the converse.

For the rational points, \(A_b(2)=(b-2)2^{b+1}+b+4\) and
\(B_b(2)=2^{b+1}-1\). Writing the positive denominator
of \(v_b\) as \(d_b\), one has
\[
 v_b-1=\frac{2^{b+2}-2b-4}{d_b}>0,\qquad
 3-v_b=\frac{(2b-6)2^{b+1}+4b+12}{d_b}>0.
 \]
For the second inequality the numerator is \(4\) when
\(b=2\), and is positive term by term for \(b\ge3\).
Thus \(1<v_b<3\) for every \(b\ge2\). Hence for the primitive
point
\[
 \operatorname{disc}(B)=(1-v_b)^2-4v_b
                         =v_b^2-6v_b+1<0.
 \]
At the imprimitive point it equals \(-7\). In both cases
\(B\) is strictly positive on the real line, and (2) holds.

We verify all exclusions uniformly. If \(\deg E<N\), its
derivative must be \(cPB^{b-1}\) for a real constant \(c\).
Its integral from \(0\) to \(1\) is nonzero unless \(c=0\),
since \(P<0\) and \(B>0\) on that interval. This contradicts
\(E(1)=E(0)=0\); the alternative \(c=0\) contradicts the
two different CRT values. Thus \(L\ne0\).
After making \(f\) monic, the integral of
\(NPB^{b-1}(z-\alpha)\) from \(0\) to \(1\) is zero.
The positive weight \(-PB^{b-1}\) shows that \(0<\alpha<1\).
All prescribed local orders are therefore exact.
The derivative is positive on \((0,\alpha)\), so \(w>0\).
Also \(v\ne0\). If \(v=w\), degree and multiplicities force
\[
                       f-v=(z-\alpha)^2 B^b .
 \]
Its right side is nonnegative for real \(z\), contradicting
\(f(0)=0<w\). Thus \(v\ne w\).
This proves every open condition in (6) for both points.
The primitive point has \(U=-v_b\ne0\); the other has \(U=0\).
\(\square\)

## Theorem V3. Complete components, actual groups and exact sizes

Passport (1) has exactly two absolute and two literal
full-braid orbits for every \(b\ge2\). Put
\[
                       W=3(2-\boldsymbol1_{b=2}).
 \]
Their invariants are
\[
 \begin{array}{c|c|c}
 &\text{primitive}&\text{imprimitive}\\ \hline
 \text{actual group}&S_{2m}&S_2\wr S_m\\
 \text{all block-count divisors}&\{1,2m\}&\{1,m,2m\}\\
 \lvert\mathrm{Deck}\rvert&1&2\\
 \text{literal orbit}&2Wm^2(3b-2)&Wm^2\\
 \text{absolute orbit}&Wm(3b-2)&Wm
 \end{array}.                                               \tag{8}
\]
Here \(W\) is the number of distinct placements of the three
factor types; these are full-orbit counts.

**Proof.** U1--U2 give one imprimitive component with the
displayed group and size, and identify every other tuple
as \(S_N\). V1 and V2 show that, after quotienting the
unordered \(b\)-critical pair, the entire primitive locus
is a nonempty irreducible rational graph with finitely
many algebraic exclusions. Forgetting the remaining
finite labels and taking connected affine quotients
preserve its connectedness. Hence it supplies exactly
one primitive component. The two actual groups prevent
the two components from merging.

For \(a=2\), U2's combined primitive literal size is
\[
 Wm^2\bigl(3(2b-1)-1\bigr)=2Wm^2(3b-2).
 \]
Since it is one component, this is its separate size.
Dividing each literal size by \(N/\lvert\mathrm{Deck}\rvert\)
gives (8). Full-twist rotations identify the absolute
and literal orbit counts as in K2. \(\square\)

## Theorem V4. Rational moduli coordinates and exceptional label action

Both components, modulo source and target affine changes,
are rational curves over \(\mathbb Q\).
For \(b>2\), their function fields have the explicit coordinates
\[
           q=t+t^{-1}\quad\text{(primitive)},\qquad
           V\quad\text{(imprimitive)}.                       \tag{9}
\]
For \(b=2\), when the two heavy branch profiles agree, one
also quotients by exchange of these branches. The coordinates
can be chosen as
\[
 \left(q+2+\frac1{q+2}\right)
       \quad\text{(primitive)},\qquad
              \frac{4V^2}{4V-1}
       \quad\text{(imprimitive)}.                           \tag{10}
\]
In all cases the exact good exclusions are those of V2
and their finite quotient images.

**Proof.** The \(x,y\) interchange has already been removed
by using \(S,V\). Exchanging the normalized critical points
\(0,1\) reflects the source \(z\mapsto1-z\), sending
\((S,V)\mapsto(2-S,W)\). On (3), palindromicity of
\(A_b,B_b\) gives \(t\mapsto1/t\). On \(U=0\) it fixes
\(V\). For \(b>2\) the branch profiles differ, so there
is no further color exchange. These observations give
the exact invariant fields (9).

For \(b=2\), exchange the two heavy branches and normalize
the new critical pair to \(0,1\). On the imprimitive curve,
\[
                  V\longmapsto\frac{V}{4V-1},
 \]
whose invariant \(V+V/(4V-1)\) is the second coordinate
in (10). On the primitive curve put
\(L_0=S-S^2+2V\). The corresponding action on \(q\) is
\[
 q\longmapsto\frac{L_0^2}{VW}-2
                   =\frac1{q+2}-2.
 \]
Indeed the new two critical points are
\(-x/(y-x),(1-x)/(y-x)\), and direct symmetric evaluation
gives the first expression; substituting (3) at \(b=2\)
gives the second. Thus \(q+2\) is inverted. Taking its
sum with its inverse proves the primitive coordinate
in (10). This also checks the extra finite action in
the equal-profile case. \(\square\)

## Exact support and scope

Run from the canonical checkout:

    python3 research/arithmetic/verify_polynomial_dt_duplicated_cycle.py

The [checker](verify_polynomial_dt_duplicated_cycle.py)
verifies polynomial identities, CRT congruences, every
critical-value exclusion at rational chart points, and
the source and color actions. It also takes complete
colored-braid absolute closures in a specified bounded
range, checking both actual groups, all blocks and deck
orders against independent hook-character totals.
The [receipt](polynomial_dt_duplicated_cycle_verification.json)
and [source ledger](polynomial_dt_duplicated_cycle_source_ledger.json)
give exact domains and source pins.
The uniform claims use V1--V4, not finite extrapolation.

The case \(\min(a,b)=3\) is also closed by the
[quadratic residue theorem](POLYNOMIAL_DUPLICATED_THREE_CYCLE_COMPONENTS.md).
The entire duplicated-cycle family is now closed by the
[four-core path normalization](POLYNOMIAL_DUPLICATED_CYCLE_PATH_NORMALIZATION.md),
which gives exactly two components for all lengths.
The residual determinant problem includes
the other unresolved first-deficit passports. Deeper deficits and
the unrestricted component classification remain open
outside already resolved families. Independent replay,
mathematical review and formal verification remain pending.
