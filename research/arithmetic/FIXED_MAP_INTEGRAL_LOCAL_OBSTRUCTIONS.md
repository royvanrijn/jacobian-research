# Local obstructions for integral targets of the fixed quintic Keller map

These are necessary conditions for Hasse failures at integral targets of the fixed map `Phi`, together with exact local statements at three. They do not assert that an integral Hasse-failing target exists, or that the displayed conditions are sufficient for global Hasse failure. The no-integral result in [the moving translated-pure-cubic branch, Theorem 4.1](FIXED_MAP_MOVING_QUADRATIC_FAMILY.md) remains a separate restricted theorem; its proof has been reviewed independently of the arguments below.

Write a target as `(u,b,c)`, put `pi=-u`, and use the scalar inverse polynomial

\[
 E(S)=\pi^5S^5-\frac32\pi^4S^4+\frac32\pi S^3
       -\left(\frac54\pi+\frac9{32}b\right)S^2
       +\frac9{16}S-\frac9{32}c.                        \tag{1}
\]

For `pi!=0` its derivative-unit root scheme is the entire fibre. By [the ambient arithmetic criterion FMA1–FMA3](FIXED_MAP_AMBIENT_LOCAL_CLASSIFICATION.md), a Hasse failure has squarefree degree-five inverse polynomial, irreducible factorization `2+3`, and matching quadratic resolvents. In particular its discriminant is a nonzero rational square. Multiplying a quintic by a nonzero scalar, or scaling its variable, multiplies its discriminant by a square, so this restriction applies equally to `32E` and the monic polynomial in `W=pi S`,

\[
 P(W)=W^5-\frac32W^4+\frac3{2\pi^2}W^3
       +\beta W^2+\frac9{16\pi}W-\frac9{32}c,
 \qquad \beta=-\frac5{4\pi}-\frac{9b}{32\pi^2}.       \tag{2}
\]

The `pi=0` plane and the `c=0` plane have rational points by the existing reconstruction: at `c=0`, `S=0` is a simple root because `E'(0)=9/16`. They cannot contain a Hasse failure.


The imported local results are J. S. Milne, *Algebraic Number Theory*, version 3.00, 11 February 2008: Theorem 7.32 (strict Newton/Hensel lifting), printed p.114, and Theorem 7.33 with its proof (unique coprime monic factor lifting), printed pp.114–115. These exact statements were read in the [author-hosted primary text](https://www.jmilne.org/math/CourseNotes/ANT300.pdf), retained as [milne-ant300.pdf](sources/milne-ant300.pdf). Proposition 7.44, printed pp.118–119, describes roots and factors from Newton polygons; its coefficient order is reversed from our exponent-index convention, so its slopes have the opposite sign. Lemma ILO1 supplies its own short supporting-line argument for the particular root-existence implication used at three.
## 1. Lemma ILO1: an adjacent Newton edge gives a local root

Let `f(S)=sum a_i S^i` be a nonzero polynomial over `Q_p`. Suppose the lower convex hull of the finite points `(i,v_p(a_i))` has a maximal edge whose endpoints have consecutive integer indices. Then `f` has a root in `Q_p`.

**Proof.** Write its endpoints as `(j,v_j)` and `(j+1,v_(j+1))`, and put `s=v_(j+1)-v_j`, an integer. Substitute `S=p^(-s)T` and divide by `p^(v_j-js)`. The resulting polynomial has integral coefficients, and its reduction modulo `p` consists exactly of two nonzero terms `A T^j+B T^(j+1)`. Maximality of the edge means that no other coefficient lies on its supporting line. The nonzero residue root `T=-A/B` is simple, since the derivative at that root is `B T^j!=0`. Ordinary simple-root Hensel lifting produces the claimed `Q_p`-root. ∎

The argument uses a supporting-line rescaling and Hensel's lemma, rather than assuming that every horizontal segment of a diagram has this consequence. Edges of horizontal length greater than one require additional analysis.

## 2. Theorem ILO2: restrictions when pi is a three-adic unit

Let `pi in Z_3^*` and `b,c in Z_3`, and suppose `Disc(P)` is a nonzero square in `Q_3`. Then

\[
 P\text{ has a root in }\mathbb Q_3
 \quad\Longleftrightarrow\quad
 \pi\equiv7\text{ or }8\pmod9.                         \tag{3}
\]

The square-discriminant assumption also forces

\[
 c\equiv3\pmod9\qquad\text{or}\qquad c\equiv9\pi\pmod{27}.
                                                        \tag{4}
\]

In the first branch `c=3 mod9`, the discriminant is always a nonzero square in `Q_3`, independently of `b`; its valuation is six. In the second branch (4) is only necessary, and higher digits can still obstruct a square discriminant.

Consequently an integral Hasse-failing target with `3` not dividing `u` must satisfy

\[
 u\equiv1\text{ or }2\pmod9,\qquad
 c\equiv3\pmod9\ \text{or}\ c\equiv-9u\pmod{27}.       \tag{5}
\]

**Proof.** Polynomial (2) is monic integral over `Z_3`. Modulo three it equals `W^2(W+pi)^3`, with coprime factors. Coprime Hensel lifting gives unique monic factors

\[
 P=QH,\quad Q=W^2+q_1W+q_0,\quad
 H=W^3+h_2W^2+h_1W+h_0,
\]

of degrees two and three with these reductions. Modulo nine (2) is exactly divisible by `W^2`, and the remaining cubic still has unit constant term. Uniqueness of coprime lifting therefore gives

\[
 Q\equiv W^2\pmod9,\qquad
 H\equiv W^3-\frac32W^2+\frac3{2\pi^2}W-\frac5{4\pi}
           \pmod9.                                    \tag{6}
\]

Shift `H` by `W=Z-pi`, writing the result as `Z^3+A Z^2+B Z+C`. Then `A,C in 3Z_3`, and

\[
 B\equiv3\left(\pi^2+\pi+\frac1{2\pi^2}\right)\pmod9,
 \quad v_3(B)=1,\quad B/3\equiv\pi\pmod3,              \tag{7}
\]

\[
 C\equiv-\pi^3-\frac32\pi^2-\frac{11}{4\pi}\pmod9.
                                                        \tag{8}
\]

For the six unit classes `pi=1,2,4,5,7,8 mod9`, the right side of (8) is respectively `6,6,3,3,0,0 mod9`.

The cubic discriminant formula now has one uniquely least-valuation term, `-4B^3`; all other terms have valuation at least four. Thus

\[
 v_3(\operatorname{Disc}H)=3,\qquad
 \operatorname{Disc}(H)/27\equiv-\pi\pmod3.             \tag{9}
\]

In particular this discriminant is nonsquare. The factors `Q,H` are coprime modulo three, so their resultant is a unit. The product discriminant identity and the assumed square discriminant force `Disc(Q)` to be nonsquare as well. Thus `Q` supplies no `Q_3`-root.

If `C` has valuation one, the shifted cubic is Eisenstein and has no `Q_3`-root. One can also see the root obstruction directly: for a proposed root of nonnegative integer valuation, either its cubic term (valuation zero) or its constant term (valuation one) is uniquely least; at negative valuation the cubic term is uniquely least. If `C in 9Z_3`, substitute `Z=3T` and divide by nine. This gives

\[
 3T^3+A T^2+(B/3)T+C/9\in\mathbb Z_3[T],
\]

whose reduction is a linear polynomial with nonzero slope. Hensel gives a `Q_3`-root. Combining this with (8) proves (3).

It remains to prove (4). Exact coefficient comparison gives

\[
 q_0=-\frac{9c}{32h_0},\qquad
 q_1h_0+q_0h_1=\frac9{16\pi},
 \quad h_0\equiv\pi^{-1}\pmod3,\quad h_1\in3\mathbb Z_3.
\]

Since `q_0 in9Z_3`, the second summand has valuation at least three. Hence

\[
 v_3(q_1)=2,\qquad q_1/9\equiv1\pmod3,
 \qquad v_3(q_0)=2+v_3(c)\quad(c\ne0).                 \tag{10}
\]

If `c` is a unit, `Disc(Q)=q_1^2-4q_0` has valuation two; together with (9), this makes `Disc(P)` have odd valuation five. If `v_3(c)>=3`, or `c=0`, it has valuation four, giving odd total valuation seven. Both cases are excluded.

If `v_3(c)=1`, write `c=3c_1`. Then

\[
 \operatorname{Disc}(Q)/27\equiv-\pi c_1\pmod3.
\]

Multiplying by (9) shows that `Disc(P)` has valuation six and normalized unit congruent to `c_1 mod3`, since the squared unit resultant is one modulo three. A three-adic unit is a square exactly when its residue modulo three is one, by simple-root Hensel lifting. This is precisely `c=3 mod9`, proving the stated sufficient square-discriminant assertion.

Finally if `v_3(c)=2`, write `c=9c_2`. Equation (10) gives

\[
 \operatorname{Disc}(Q)/81\equiv1-\pi c_2\pmod3.
\]

Unless `pi c_2=1 mod3`, this has valuation exactly four and again gives odd full valuation seven. Therefore `c=9pi mod27` is necessary in this last case. This proves (4), and `pi=-u` gives (5). ∎

## 3. Theorem ILO3: targets with pi divisible by three are locally soluble at three

If `pi in3Z_3` is nonzero and `b,c in Z_3`, the entire fibre has a `Q_3`-point. In particular divisibility of the first integral target coordinate by three is not itself an obstruction to local solubility.

**Proof.** Use

\[
 F(S)=32E(S)=32\pi^5S^5-48\pi^4S^4+48\pi S^3
               -(40\pi+9b)S^2+18S-9c.                \tag{11}
\]

If `c=0`, its simple root zero proves the assertion. Otherwise put `e=v_3(pi)>=1`, `t=v_3(c)>=0`, and `k=v_3(40pi+9b)`, allowing `k=infinity` if that coefficient vanishes. The coefficient valuations at indices zero through five are

\[
 (2+t,\ 2,\ k,\ 1+e,\ 1+4e,\ 5e),\qquad
 k\ge\min(e,2).                                        \tag{12}
\]

The segment joining indices three and five has slope `2e-1/2` and is a maximal lower supporting edge: index four lies strictly above its line, as do indices zero, one and two. At index two that line has value `3/2-e<k`; at indices one and zero its values are `2-3e` and `5/2-5e`, also strictly below the actual valuations.

The lower hull between indices zero and three cannot consist of a single segment. If `e+2t>1`, index one lies strictly below the chord, whose value there is `(5+2t+e)/3>2`. The only remaining case is `e=1,t=0`; here `k=1` exactly, while the chord has value two at index two. Thus some interior index, one or two, is a genuine vertex. A hull spanning three horizontal units with an interior integer vertex has at least one maximal edge of horizontal length one. Lemma ILO1 supplies a simple `Q_3`-root of `F`. Simplicity also follows from the simple lifted residue root in that lemma, so it survives the derivative localization and gives a point of the entire fibre. ∎

## 3a. Corollary ILO3a: the quadratic resolvent must ramify at three

Suppose an integral target `(u,b,c)` of the fixed map has a Hasse-failing fibre, and write its rank-five algebra as `K x L_3`, with quadratic field `K` and cubic field `L_3`. Put `pi=-u`. Then

\[
 K\otimes\mathbb Q_3\simeq
 \begin{cases}
 \mathbb Q_3(\sqrt{-3}),&\pi\equiv7\pmod9,\\
 \mathbb Q_3(\sqrt3),&\pi\equiv8\pmod9\text{ or }3\mid\pi.
 \end{cases}                                          \tag{3a.1}
\]

These cases exhaust the possibilities by ILO2 and `pi!=0`. The entire local fibre algebra is `Q_3 x (K tensor Q_3) x (K tensor Q_3)`. In particular the global quadratic field `K` must ramify at three. If `K=Q(sqrt(-3))`, necessarily `pi=7 mod9`, or equivalently `u=2 mod9`; the case `3|pi` is excluded for this resolvent regardless of how the cubic field is generated.

**Proof.** By FMA2, the splitting field of the global fibre algebra has group `S_3`, acting on its quadratic and cubic components. Its decomposition group `D_3` is proper because the fibre has a local point. If a local irreducible quadratic orbit occurs, its degree forces `2` to divide `|D_3|`. A proper subgroup of `S_3` with even order is exactly `C_2`. In the combined action on two and three points, such a subgroup has orbit sizes `2,2,1`. Both size-two orbits have the same local quadratic splitting field. Thus it suffices to identify one irreducible quadratic local factor.

When `pi` is a unit, ILO2 gives `pi=7 or8 mod9`. Its Hensel cubic factor `H` has a local linear factor and discriminant squareclass `-3pi`, by (9). Its remaining quadratic factor is irreducible with this squareclass, so it gives the required orbit. A unit congruent to one modulo three is a square in `Q_3` by simple-root Hensel lifting. For `pi=7 mod9`, `pi` is such a square, giving field `Q_3(sqrt(-3))`; for `pi=8 mod9`, `-pi` is a square, giving `Q_3(sqrt(3))`.

Now let `e=v_3(pi)>=1`. Choose `theta` with

\[
 \theta^2=-\frac3{2\pi^2},\qquad v_3(\theta)=\frac12-e,
 \qquad F=\mathbb Q_3(\theta).
\]

This is a ramified quadratic extension. Its squareclass is `-6`, equal to the squareclass of three over `Q_3` because the unit `-2=1 mod3` is a local square. Normalize (2) by

\[
 \frac{P(\theta T)}{\theta^5}
 =T^5-\frac3{2\theta}T^4-T^3
    +\frac\beta{\theta^3}T^2
    +\frac9{16\pi\theta^4}T-\frac{9c}{32\theta^5}.
                                                        \tag{3a.2}
\]

The four coefficients besides the displayed `T^5-T^3` have strictly positive valuation. Indeed their valuations are, respectively,

\[
 e+\frac12,\quad
 \ge\min(2e-3/2,e+1/2)>0,\quad
 3e,\quad\ge5e-1/2>0,
\]

using `v_3(beta)>=min(-e,2-2e)` and `b,c in Z_3`. The reduction over the residue field of `F` is therefore `T^3(T^2-1)`. Both unit roots `1,-1` are simple in characteristic three and lift in `F`. A resulting root `W=theta T` has nonintegral valuation `1/2-e`, so it is not in `Q_3`; its field of definition is exactly `F`, of degree two. This gives the required irreducible quadratic orbit. The preceding `D_3=C_2` argument identifies both quadratic local factors with `F=Q_3(sqrt(3))` and proves (3a.1).

Both fields displayed in (3a.1) are ramified quadratic extensions of `Q_3`, proving the final assertions. ∎

## 4. Two-adic restrictions

Use the coefficient-index convention `(i,v_2(a_i))` for `F(S)=sum a_iS^i`.
The slopes in this convention are the negatives of the root valuations.
For `pi!=0`, put `e=v_2(pi)`, and write `k=v_2(40pi+9b)`, allowing
`k=infinity` if that coefficient vanishes. For integral `pi,b,c`, the
coefficient valuations of (11) at indices one through five are

\[
 (1,k,4+e,4+4e,5+5e),\qquad v_2(a_0)=v_2(c)\ge0,              \tag{13}
\]

with the constant coefficient omitted when `c=0`.

### Theorem ILO4: the complete odd-b dichotomy at two

Let `pi,b,c` be integers, `pi!=0`, and `b` odd.

- If `e=0,1,3`, the target cannot be Hasse-failing.
- If `e=2` or `e>=4`, its entire fibre always has a `Q_2`-point.

The first assertion uses the necessary irreducible `2+3` factorization and
matching quadratic-resolvent condition. It does not claim that all targets
in its excluded cases have no `Q_2`-point. The second assertion is unconditional
and allows arbitrary `c`, including repeated inverse polynomials.

**Proof.** Here `k=0`. For `0<=e<=3`, the rightmost lower Newton edge joins

\[
 (2,0)\quad\hbox{to}\quad(5,5+5e).                              \tag{14}
\]

The intermediate point at index three is above its line because
`4+e>(5+5e)/3` in this range; index four is also strictly above. The points
at zero and one are strictly above the supporting line on the left,
independently of `c`. If `e=0,1,3`, the slope has reduced denominator three.
The Newton-polygon theorem gives a degree-three factor over `Q_2`, all
of whose roots have valuation `-(5+5e)/3`. Every root field therefore has
ramification index divisible by three. The factor is irreducible, its
field has degree three, and this extension is totally ramified.

We give the needed local resolvent argument explicitly. In any totally
ramified cubic extension `T/Q_2`, the residue field is `F_2` and cubing is
bijective on the unit group by simple-root Hensel lifting. If `varpi` is a
uniformizer, `2/varpi^3` is a unit, so multiplying `varpi` by a suitable unit
gives an element whose cube is two. Hence `T=Q_2(cuberoot(2))`. Its quadratic
resolvent is `Q_2(sqrt(-3))`, a nonsplit unramified quadratic extension:
`X^2+X+1` is irreducible modulo two, and the cubic discriminant has
squareclass `-3`.

If the rational inverse quintic had the irreducible quadratic-times-cubic
factorization required by a Hasse failure, the local irreducible cubic
factor just found would have to be its cubic factor. A quadratic factor
cannot contain an irreducible local factor of degree three. Matching the
quadratic resolvent would then make the quadratic factor nonsplit over
`Q_2` as well. Neither factor would have a `Q_2`-root, contradicting local
solubility. This proves the exclusions.

For `e=2`, set `S=2^(-5)Z`. The polynomial

\[
 2^{10}F(2^{-5}Z)\in\mathbb Z_2[Z]
 \quad\hbox{reduces to}\quad Z^2(Z^3+1).                        \tag{15}
\]

Its residue root `Z=1` is simple. For `e>=4`, set `r=e+4`; then

\[
 2^{2r}F(2^{-r}Z)\in\mathbb Z_2[Z]
 \quad\hbox{reduces to}\quad Z^2(Z+1).                          \tag{16}
\]

Again `Z=1` is simple. Hensel gives a simple `Q_2`-root of `F` in both
cases. The nonzero derivative survives the derivative localization in the
entire-fibre formula, giving the asserted actual local point. \(\square\)

### Proposition ILO5: further even-b restrictions and automatic local points

For integral `pi,b,c` with `pi!=0`, the following additional statements hold.
The valuation `v_2(b)=infinity` is permitted where `4|b` is specified.

| Conditions | Consequence |
|---|---|
| `(v_2(pi),v_2(b))=(0,1),(1,2),(2,1)` | no Hasse failure |
| `v_2(pi)=0`, `v_2(b)=2`, `c` odd | inverse quintic irreducible over `Q_2`; no Hasse failure |
| `(v_2(pi),v_2(b))=(1,1)` | entire fibre has a `Q_2`-point |
| `pi` odd, `4|b`, `2|c` | entire fibre has a `Q_2`-point |

These cases do not exhaust all even-`b` targets.

**Proof.** In the first three cases, `k=v_2(b)` and the rightmost edge
joins `(2,k)` to `(5,5+5e)`. Its slope is respectively `4/3,8/3,14/3`.
Every other coefficient is strictly above its supporting line. The same
irreducible totally ramified cubic and resolvent argument used in ILO4
excludes Hasse failure, for every `c`.

If `pi` is odd and `v_2(b)=2`, the integral polynomial `F(Z/2)` reduces,
when `c` is odd, to

\[
 Z^5+Z^4+Z^2+Z+1\in\mathbb F_2[Z].                             \tag{17}
\]

It has no linear factor, and division by the only irreducible quadratic
`Z^2+Z+1` has remainder one. A reducible quintic without a linear factor
would have a quadratic factor, so (17) is irreducible. Its leading
coefficient is a unit. Thus `F`, and hence `E`, is irreducible over `Q_2`,
which precludes the rational `2+3` factorization required for Hasse failure.

If `(e,v_2(b))=(1,1)`, the integral polynomial `2^5F(2^(-3)Z)` reduces
to `Z^2(Z^3+1)` and has simple root one. Finally, if `pi` is odd, `4|b`
and `2|c`, the integral polynomial `F(Z/2)` reduces to

\[
 Z^5+Z^4+\epsilon Z^2+Z,\qquad
 \epsilon=(b/4)\bmod2.                                        \tag{18}
\]

Its residue root zero is simple. Simple-root Hensel lifting and derivative
localization prove both automatic local-point assertions. \(\square\)

The Newton-polygon and rational factor statements used here are
[J. S. Milne, *Algebraic Number Theory*, v3.00, Proposition 7.44,
printed p. 118, PDF p. 126](https://www.jmilne.org/math/CourseNotes/ANT300.pdf).
Milne orders monic coefficients in descending degree; reversing the indices
gives the slope convention above. The root-field ramification deduction
uses the extended value group in Remark 7.43 and `n=ef` in Corollary 7.42.
The simple-root lifts, including cubing on local units, use Proposition 7.31,
printed p. 113, PDF p. 121. The exact retained source is
[sources/milne-ant300.pdf](sources/milne-ant300.pdf); these applied statements
were read directly. The tame cubic's resolvent identification is proved above,
rather than imported as a claim about all ramified decomposition groups.

## 5. Verification and remaining scope

Run the checker from the repository root:

```sh
.venv/bin/python research/arithmetic/verify_fixed_map_integral_local_obstructions.py --output research/arithmetic/fixed_map_integral_local_obstructions_checks.json
```

The checker verifies exact coefficient and discriminant identities, all unit classes modulo nine used above, the finite normalized-unit residue computations, the proper `S_3` orbit argument for the ramified resolvent, and the stated two-adic residue and coefficient-valuation supporting-line controls. Such residue enumerations certify the explicitly reduced finite steps. They do not replace the unrestricted Hensel-factor, valuation, or arithmetic proofs. No bounded target search is used as a proof of nonexistence.

Theorem ILO2 is a necessary-condition result for Hasse failures, with an exact local statement under a square-discriminant hypothesis. Theorem ILO3 gives automatic local solubility at one place; it does not imply global solubility. General moving quadratic factors and general cubic generators still leave the integral Hasse-existence question unresolved. The results have written proofs and scoped exact controls; no external review, formalization or unrestricted counting theorem is claimed.

**An explicit local escape case.** The integral target `(u,b,c)=(1,4,12)`,
with `pi=-1`, passes both low-prime local tests. Its inverse polynomial scaled
by 32 is

\[
 F(S)=-32S^5-48S^4-48S^3+4S^2+18S-108,
 \qquad\operatorname{Disc}(F)=2^{24}3^6\cdot50260753.
\]

At two, `v_2(F(2))=3>2v_2(F'(2))=2`; at three,
`v_3(F(5))=4>2v_3(F'(5))=2`. Strict Hensel lifting therefore supplies
simple local roots at both places, which survive derivative localization.
The discriminant is also a square in both `Q_2` and `Q_3`: its two-adic
normalized unit is one modulo eight, and its three-adic normalized unit
is one modulo three. It is not a rational square, because 50260753 is not
an integer square. In fact the fibre has no `Q_11`-point: its leading
coefficient is an eleven-adic unit, but `F` has no root modulo eleven.
This example is not Hasse-failing. It demonstrates that the local tests at
two and three leave genuine integral targets. The global matching-resolvent
condition and the other primes still require separate checks.
