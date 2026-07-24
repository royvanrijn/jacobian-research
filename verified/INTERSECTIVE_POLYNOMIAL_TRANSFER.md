# Intersective-polynomial transfer and a rank-one fixed Keller pencil

This note advances two arithmetic Keller-fiber problems:

1. it proves an unconditional quadratic-stabilized transfer theorem for
   intersective polynomials over number fields;
2. it proves that one fixed quintic Keller pencil contains infinitely many
   rational quadratic-cubic factorizations with matching discriminant
   squareclass, and gives two new Hasse-failing fibers in that pencil.

The remaining gap for infinitely many Hasse failures is isolated to local
solubility at the varying ramified primes.

## 1. Quadratic-stabilized transfer

Let \(K\) be a number field. Let \(f\in K[X]\) be squarefree, have no
\(K\)-rational root, and have a root in every completion of \(K\), including
every real completion.

### Theorem

There are \(\alpha\in K^\times,\beta\in K\) and an irreducible quadratic
\(R\in K[W]\) such that

\[
 P(W)=f(\beta+\alpha W)R(W)
\]

is squarefree and tangent-admissible:

\[
\begin{aligned}
P(1)-P(0)&=P'(0),\\
P'(1)&\ne P'(0),\\
P''(1)&\ne2(P'(1)-P'(0)).
\end{aligned} \tag{1.1}
\]

Consequently \(P\) is the inverse polynomial of a complete regular fiber of
an explicit polynomial Keller map of \(\mathbb A^3_K\). That fiber has a
point over every completion of \(K\) but no \(K\)-point. Its degree is

\[
 \deg f+2.
\]

Scheme-theoretically, it is the disjoint union of the original finite etale
scheme, after its affine change of generator, and one quadratic etale
scheme.

### Proof

Put \(g(W)=f(\beta+\alpha W)\), and abbreviate

\[
 u=g(0),\quad v=g(1),\quad d=g'(0).
\]

Write

\[
 R(W)=AW^2+BW+C.
\]

The tangent equation for \(P=gR\) is the single linear equation

\[
 vA+(v-u)B+(v-u-d)C=0. \tag{1.2}
\]

Because \(f\) has no \(K\)-root, \(u,v\ne0\). Thus (1.2) is an affine
two-parameter linear system of quadratics.

Solving (1.2) for \(A\), the discriminant \(B^2-4AC\) is a binary quadratic
form in \(B,C\). It is a square in \(K[B,C]\) exactly when

\[
 vu-vd-u^2=0. \tag{1.3}
\]

The left side of (1.3), viewed as a polynomial in \(\alpha,\beta\), is not
identically zero. Choose \(\beta\) with \(f'(\beta)\ne0\). For fixed
\(\beta\), its leading term as a polynomial in \(\alpha\) comes from

\[
 f(\beta+\alpha)\bigl(f(\beta)-\alpha f'(\beta)\bigr)
 -f(\beta)^2
\]

and is nonzero. Hence a generic affine change makes the quadratic
discriminant nonsquare in \(K(B,C)\).

The two weighted nondegeneracy conditions restrict to two further linear
forms on (1.2). They are not identically zero after a generic affine change.
Indeed, on the direction \(C=0\), where

\[
 A=-\frac{v-u}{v}B,
\]

they are, up to nonzero scalar factors,

\[
 u(v+g'(1))-v^2 \tag{1.4}
\]

and

\[
 u(2g'(1)+g''(1))-2vg'(1). \tag{1.5}
\]

For fixed \(\beta\), the highest \(\alpha\)-degree terms in (1.4) and (1.5)
are respectively \(-v^2\) and \(-2vg'(1)\), so neither expression is
identically zero. We may therefore choose \(\alpha,\beta\) avoiding
(1.3)-(1.5), \(u=v\), and the finitely many other endpoint degeneracies.

The resultant \(\operatorname{Res}(g,R)\) is not identically zero on (1.2).
The direction \(C=0\) has roots \(0\) and \(v/(v-u)\), both in \(K\);
neither is a root of \(g\), since \(g\) has no \(K\)-root. Thus coprimality,
\(A\ne0\), and the two nondegeneracy conditions define a nonempty open
subset of the parameter plane.

Finally, \(Y^2-(B^2-4AC)\) is irreducible over \(K(B,C)\). Hilbert
irreducibility supplies a \(K\)-point of the parameter plane at which
\(R\) is irreducible, while simultaneously avoiding the preceding proper
closed subsets. This proves (1.1), squarefreeness, and irreducibility of the
auxiliary factor.

Every local root of \(f\) transfers through the affine change and remains a
root of \(P\). Neither factor has a \(K\)-root. The exact weighted transfer
theorem in
[ARITHMETIC_KELLER_FIBER_ENGINEERING.md](../ARITHMETIC_KELLER_FIBER_ENGINEERING.md)
then supplies the asserted complete Keller fiber.

## 2. The fixed quintic pencil

Use the degree-five seed

\[
H(W)=-531441W^5+1003833W^4-758889W^3+286497W^2
\]

from
[MINIMAL_HASSE_PRINCIPLE_KELLER_FIBER.md](MINIMAL_HASSE_PRINCIPLE_KELLER_FIBER.md).
Its fixed \(C=1\) pencil is

\[
 E_{s,t}(W)=H(W)-sW+t.
\]

For

\[
 q(W)=W^2+aW+b,
\]

division of \(H\) by \(q\) has a linear remainder. Removing that remainder
selects a unique \((s,t)\) for which \(q\mid E_{s,t}\). Requiring the
quadratic field of \(q\) to equal the discriminant field of the cubic
quotient gives a surface. On the slice \(a=-7/9\), it is the quartic

\[
\begin{aligned}
Y^2={}&3(57395628b^3-26749197b^2\\
     &\qquad+4181544b-219512)(49-324b). \tag{2.1}
\end{aligned}
\]

## 3. Rank-one theorem

The Jacobian of (2.1) has the Weierstrass model

\[
\begin{aligned}
y^2={}&x^3-7996592727x^2\\
&+21315161806412546304x\\
&-18938742282153164303832264192. \tag{3.1}
\end{aligned}
\]

Exact two-descent and saturation in PARI/GP give:

\[
\boxed{\operatorname{rank}E(\mathbb Q)=1,\qquad
       E(\mathbb Q)_{\mathrm{tors}}=0.}
\]

The conductor is

\[
 3202146=2\cdot3^3\cdot19\cdot3121.
\]

It follows that (2.1) has infinitely many rational points. Since the map to
the \(b\)-line has degree two, infinitely many distinct \(b\)-values occur.
The fixed Keller pencil therefore contains infinitely many rational
quadratic-cubic factorizations with matching discriminant squareclass.

This does not by itself prove that infinitely many specializations have both
factors irreducible: rational reducibility is a thin arithmetic condition,
not a finite Zariski-closed locus on the elliptic curve. Even after
irreducibility, local solubility at every ramified prime remains separate:
a local \(S_3\) decomposition group can act without a fixed point on both
orbits.

## 4. Two further Hasse-failing fibers of the same map

Besides the original point \(b=37/243\), (2.1) contains

\[
\left(\frac{3139}{19764},\frac{137805}{3721}\right),
\qquad
\left(\frac{31687}{204363},\frac{23067311}{707281}\right).
\]

### First new fiber

The quadratic and cubic factors are

\[
\begin{aligned}
q_2={}&19764W^2-15372W+3139,\\
c_2={}&-533628W^3+592920W^2-216099W+25429.
\end{aligned}
\]

The corresponding pencil parameters are

\[
s_2=\frac{3207688047}{59536},\qquad
t_2=\frac{239464893}{59536}.
\]

Thus the Keller target is

\[
\left(
\frac{79821631}{6857475552},
\frac{3207688047}{59536},
1
\right).
\]

Both factors are irreducible. At every unramified prime, the common
\(S_3\) group is covered by the cubic point stabilizers and the quadratic
\(A_3\) stabilizer. The ramified primes are

\[
2,3,5,61,9187,7765337,
\]

and exact \(p\)-adic factorization gives a root in at least one factor at
each of them.

### Second new fiber

The factors are

\[
\begin{aligned}
q_3={}&204363W^2-158949W+31687,\\
c_3={}&-1839267W^3+2043630W^2-751770W+89959.
\end{aligned}
\]

Here

\[
s_3=\frac{38120229081}{707281},\qquad
t_3=\frac{2850530833}{707281},
\]

and the Keller target is

\[
\left(
\frac{2850530833}{244398120426},
\frac{38120229081}{707281},
1
\right).
\]

The ramified primes are

\[
2,3,19,29,389,751,3121,458701.
\]

Again both factors are irreducible and one has a root over \(\mathbb Q_p\)
at every listed prime. Thus both displayed targets give complete regular
Hasse-failing fibers of the same fixed Keller map.

## 5. What remains for an infinite Hasse theorem

Rank one supplies infinitely many rational factorizations. To promote
infinitely many of them to Hasse failures, one must first prove simultaneous
irreducibility of the quadratic and cubic factors for infinitely many
Mordell-Weil multiples, and then exclude local decomposition group \(S_3\)
at every newly ramified prime.

Three plausible routes remain:

1. combine an elliptic Hilbert-irreducibility argument with the two factor
   covers to obtain infinitely many irreducible specializations;
2. prove that the pencil identities force a local root at every prime on
   that irreducible subset, or identify a Mordell-Weil progression on which
   the possible local \(S_3\) primes are excluded;
3. apply a squarefree or power-free sieve to the cubic discriminant values,
   reducing all ramification to transposition type.

The first route is now testable by eliminating \(b,Y\) together with the
local irreducible-Eisenstein conditions. The second can be explored using
the formal group of (3.1). The third is likely conditional without further
factorization of the discriminant divisor.

## Verification

Run

```bash
.venv/bin/python scripts/verify_intersective_transfer_identity.py
.venv/bin/python scripts/verify_minimal_hasse_keller_fiber.py
"$(brew --prefix pari)/bin/gp" -q scripts/verify_fixed_hasse_resolvent_rank.gp
```

The first Python checker verifies the universal tangent, discriminant, and
nondegeneracy identities. The second verifies the quartic derivation. The
PARI/GP certificate verifies its Jacobian, conductor, torsion, exact rank,
the three rational points, irreducibility of the four new factors, and all
ramified-prime local roots.
