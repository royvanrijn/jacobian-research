# Elliptic factors beyond the two doubled Mordell projections

The doubled genus-five curve has two further types of elliptic
factor. One comes from the cubic quotient on the parameter line
and gives a third rank-positive quadratic twist over \(\mathbb Q\).
The other appears through a pure cubic constant extension.
Its rational order-three direction pulls back into the span of
the two existing norm classes.

## The remaining Jacobian factors

**Theorem DJ1.** Let \(b\in\mathbb Q^\times\) be noncubic,
\(m=b^2\), and \(C_m\) be the curve of
[DB1](ABC_CUBIC_TANGENT_DOUBLED_BRANCH_FAMILY.md#classification-of-the-doubled-branch-model).
Put \(L=\mathbb Q(\rho)\) with \(\rho^3=m\) and let \(M\) be
its normal closure. Let
\[
 H_3:U^2=X^3+X^2-X,\qquad
 H_4:V^2=X^3-\tfrac34X^2+6X .
\]
There is a rational isogeny decomposition
\[
 \operatorname{Jac} C_m\ \sim_{\mathbb Q}\
 E_b\times E_{b^2}\times H_3\times
 \ker\!\left(\operatorname{Res}_{L/\mathbb Q}H_4
                    \xrightarrow{\operatorname{Tr}}H_4\right).
 \tag{1}
\]
The quadratic twist of \(H_3\) by the common field polynomial
has generic ranks \(1,1,1\) over \(\mathbb Q(s),L(s),M(s)\),
respectively. For \(H_4\) those ranks are \(0,1,2\),
respectively.

The map to \(H_3\) is
\[
 w=ms^3,\qquad X=(w+1)/(w-1),\qquad U=Y/(w-1)^2.
\]
It is the quotient of \(C_m\) by \(s\mapsto\zeta_3s\) over
the algebraic closure. For the other factor, set \(t=\rho s\),
so \(w=t^3\). The involution
\[
 (t,Y)\longmapsto(-1/t,-Y/t^6)
\]
has invariants \(u=t-1/t\) and \(v=Y/t^3\). Its quotient is
\[
 G:\quad v^2=(u^3+3u)(u^3+3u+4).
\]
Indeed \(u^3+3u=t^3-t^{-3}\) and the displayed equation
is \(Y^2/t^6\) from DB1. There is a degree-three map
from \(G\) to \(H_4\) given by
\[
 D=(u+1)(u^2+3),\qquad
 X=\frac{6u}{D},\qquad V=\frac{3(2u^2+3u+3)v}{D^2}.
 \tag{2}
\]
Substitution verifies the elliptic equation. Let \(g_0,g_1,g_2\)
be the three resulting maps on \(C_m\) to \(H_4\), obtained by
using \(t,\zeta_3t,\zeta_3^2t\) in (2). Each has degree six
and the Galois group of \(M/\mathbb Q\) permutes them.

The two \(j\)-invariants of \(H_3,H_4\) are \(16384/5\) and
\(804357/500\); in particular both curves are non-CM.
They are also geometrically nonisogenous.
If non-CM elliptic curves over \(\mathbb Q\) were geometrically
isogenous, the Galois action on their rank-one integral Hom
space would give a quadratic character, so one would be
rationally isogenous to a quadratic twist of the other.
Both have good reduction at \(13\) and their traces there
are \(2,-4\). The character would have to be unramified at a
prime of common good reduction, with the two traces agreeing
up to sign, which fails at \(13\).

The trace map \(g_0+g_1+g_2\) on the curve is invariant under
\(s\mapsto\zeta_3s\) and factors through \(H_3\) over the
algebraic closure. It induces zero on Jacobians since
\(\operatorname{Hom}(H_3,H_4)=0\).
Thus \(g_0^*+g_1^*+g_2^*=0\), and the nonzero pullbacks span
the two-dimensional standard representation of \(S_3\):
the permutation representation modulo its invariant line
is irreducible over \(\mathbb Q\).
Over the algebraic closure the images of two independent
pullbacks are distinct elliptic subvarieties because
\(\operatorname{End}(H_4)=\mathbb Z\).
The sum of the pullbacks, restricted to triples with sum zero,
is defined over \(\mathbb Q\) and gives a dimension-two image.
This supplies the last factor of (1).

Each of \(E_b,E_{b^2}\) occurs by DB1, and they are rationally
nonisogenous by GT5 since \(b\) is noncubic.
They are CM and neither is geometrically isogenous to the
two non-CM factors. With the distinction between \(H_3,H_4\)
already proved, the dimensions sum to five. This proves (1).

For a constant elliptic target, rational maps from \(C_m\)
are controlled by homomorphisms from its Jacobian, up to constant
translations. There are rational points over \(s=0\), so this
identification applies over each constant field in the statement.
The hyperelliptic involution acts as \(-1\) on the Jacobian.
Consequently the rank of the quadratic twist is the rank of
the corresponding Hom group, because the only anti-invariant
constant points are two-torsion. For \(H_3\), (1) has one copy
even over the algebraic closure. For \(H_4\), take invariants
of the standard two-dimensional representation under \(S_3\),
a transposition, and the trivial group to obtain \(0,1,2\).

## The coefficient-extension norm relation

**Proposition DJ2.** The rational order-three point of \(H_4\)
gives a marked norm class on each of the three projections (2).
Over \(L(C_m)\) the class from \(g_0\) is \([\beta/\alpha]\)
in the notation of DB2. Over \(M(C_m)\) this holds for all
three projections. Thus these additional maps give no new
marked cubic norm direction.

The point \((1,5/2)\) on \(H_4\) has order three and
tangent \(V=\tfrac32X+1\). Its Kummer function is
\(\xi=V+\tfrac32X+1\), with norm \(-(X-1)^3\).
The rational cyclic subgroup is unique: the division polynomial
is \(3(X-1)(X^3+12X+12)\) with no rational root in the cubic.
Define
\[
 r_1=1+\frac{Y}{(w-1)^2},\qquad
 r_2=1+\frac{Y}{(w+1)^2}.
\]
Directly in the doubled curve's function field after adjoining
\(t\) with \(t^3=w\), one has
\[
 \xi\frac{r_1}{r_2}
 =
 \left[
 -\frac{(t+1)(t^2-3t-1)}{2(t^2+t+1)}
 +
 \frac{Y}{2(t-1)(t^2+t-1)(t^2+t+1)}
 \right]^3 .
 \tag{3}
\]
Since \(\alpha=br_1,\ \beta=b^2r_2\) and
\(b=(\rho^2/b)^3\) is a cube in \(L\), this gives the
stated relation. The same identity holds after substituting
\(\zeta_3^it\). The Jacobi factors in (1) therefore give more
elliptic points over the extension without adding independent
marked order-three functions.
The [Sage checker](verify_mordell_doubled_jacobian_complement.py)
replays (2), (3), and the target curve invariants.

## Three simultaneous marked quadratic twists

**Corollary DJ3.** At least
\[
 X^{1/6}\exp\!\left(-O\!\left(\frac{\log X}{\log\log X}\right)\right)
\]
distinct imaginary quadratic fields of discriminant of absolute
value at most \(X\) give positive rank on all three quadratic
twists of \(E_6,E_{36},H_3\). The primes \(2,3,5\) split,
and the conductors of the three twists are
\(972d^2,243d^2,20d^2\). Their finite \(3\)-Selmer dimensions
are even, even, and odd, respectively.
The cubes of the split prime classes over \(2,3\) in DB5
continue to hold.

Use \(b=6\) and the interval of
[DM2](MORDELL_CUBIC_DOUBLED_SIMULTANEOUS_TWISTS.md#counting-distinct-simultaneous-twists)
with the additional congruence \(s\in5\mathbb Z_5\).
This retains \(\gg H^2\) rational parameters of height at most
\(H\), by counting primitive pairs with numerator divisible
by \(5\) and denominator a unit at \(5\).
It makes \(\Delta(s)\equiv1\bmod5\), so \(5\) splits.
The first two twists have the asserted properties by DM1.
The third anti-invariant point on \(H_3\) has nontrivial
\(3\)-isogeny evaluation by DB6 and the valuation row at \(2\)
of DB2. In the twist coordinates \(x'=d(X-1),\ y'=d\sqrt d\,U\)
the equation is \(y'^2=x'^3+d(2x'+d)^2\).
The [Cohen–Pazuki Kummer function](sources/cohen-pazuki-2009-three-isogeny.pdf)
is \(y'-\sqrt d(2x'+d)=(\sqrt d)^3(U-(2X-1))\).
Its class is the inverse of the class of \(\gamma\), by the cube norm.
On \(H_3\) there is no anti-invariant three-torsion
over a quadratic field: the division polynomial of DB6 has
only its rational marked root over such fields because the
other factor is irreducible cubic. Thus the point has infinite
order, as in DM1.

The conductor of \(H_3\) is \(20\). The extra twist is locally
isomorphic to it at \(2,5\), and each prime dividing \(d\)
contributes exponent two, since it is a ramified quadratic
twist of good reduction at a prime at least seven.
For its rational \(3\)-isogeny, the local ratios of the base
curve at \(2,3,5\) are \(1/3,1,1\), respectively, using the
minimal invariant differentials and Tamagawa numbers. At primes
dividing \(d\), both source and target twists have type \(I_0^*\)
and the Tamagawa numbers agree: the original odd isogeny
identifies their 2-torsion and the number of roots of the
Tate polynomial is the number of nonzero 2-torsion points
in the residue field. The differential ratio is a unit
there. At infinity the order-three kernel is nonreal on
the negative twist and the isogeny induces an isomorphism
on the real component groups. Thus the global Selmer ratio
is \(1/3\). The composed-isogeny sequence and the alternating
Cassels–Tate pairing give odd finite \(3\)-Selmer dimension
with no rational three-torsion, by the parity argument in
[Bhargava–Elkies–Shnidman, Proposition 49](sources/bhargava-elkies-shnidman-mordell-isogeny-selmer.pdf).
The count follows by the uniform per-field multiplicity
from DM2 in the restricted congruence family.
