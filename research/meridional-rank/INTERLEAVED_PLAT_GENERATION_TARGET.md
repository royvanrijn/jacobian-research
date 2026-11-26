# An interleaved four-bridge link for a three-meridian test

The explicit 8-plat below has bridge number four. Its two components are an unknot and a trefoil, so their separate knot bridge numbers sum to three. The larger bridge bound comes from the distance of the bridge sphere of the whole link. A generating tuple with one meridian of the unknot and two of the trefoil would supply the desired separation. Two triples of cap meridians reached the limit of a bounded coset enumeration, and generation remains unverified.

## The plat and its components

Number braid positions from left to right. Attach top and bottom caps to `(1,2)`, `(3,4)`, `(5,6)`, and `(7,8)`. Read the following 33 rows from top to bottom. Odd-numbered rows have boxes with braid generators `σ_2, σ_4, σ_6`. Even-numbered rows have boxes with generators `σ_1, σ_3, σ_5, σ_7`. Every box has exponent 4 except those in this table.

| Row | Changed exponents |
| --- | --- |
| 1 and 33 | `σ_2^3`, `σ_6^3` |
| 2, 4, 6, and 8 | `σ_3^3` |
| 7 | `σ_4^(-3)` |
| 10 and 12 | `σ_5^3` |
| 11 | `σ_4^3` |

Label the four top caps `A,B,A,B`, in order. The first row changes the lane sequence from `AABBAABB` to `ABABABAB`. Rows 2 and 4 temporarily interchange lanes 3 and 4. Rows 6 and 8 repeat that interchange. Rows 10 and 12 temporarily interchange lanes 5 and 6. The last row restores `AABBAABB`. Between these episodes, all crossings involve one strand of each label. The intermediate rows 3, 7 and 11 are the only rows with same-label crossings. Deleting the other component and renumbering the four surviving lanes gives the subbraids

```text
A: σ_2^4 σ_2^(-3) σ_3^4 = σ_2 σ_3^4
B: σ_1^4 σ_1^4 σ_2^3    = σ_1^8 σ_2^3.
```

The end twists `σ_3^4` for `A` and `σ_1^8` for `B` lie next to the respective closing caps and can be undone by rotating those caps. The 4-plats of `σ_2` and `σ_2^3` are the unknot and a trefoil, up to the choice of mirror convention. Each is connected. The [exact strand checker](verify_interleaved_plat.py) verifies these braid restrictions and the connected components of the full plat by tracking the strands and cap pairings. A separate internal Spherogram 2.4.1 run using the braid word and consecutive-cap closure also produced two components, which simplified to the empty unknot diagram and a three-crossing trefoil diagram.

## Independent bridge number

The exhibited bridge sphere has four bridges. The two component bridge numbers give a lower bound of three, because a generic height function must have at least one maximum on the unknot and at least two on the trefoil. All 33 rows have twists of absolute value at least three, and the length exceeds `4·4·(4−2)=32`. The high-distance minimality theorem therefore gives bridge number four. This uses the link version, via [Moriah and Purcell, *Diagram uniqueness for highly twisted plats*, Corollary 2.8](https://users.monash.edu/~jpurcell/papers/DiagramUnique082017.pdf), combining the distance formula and the bridge comparison theorem for links. [Johnson and Moriah, Theorem 1.2](https://msp.org/agt/2016/16-6/agt-v16-n6-p08-p.pdf), [Johnson and Tomova, Theorem 4.4](https://msp.org/agt/2011/11-4/agt-v11-n4-p05-p.pdf).

The literal diagram has even-column gcds `1,1,1`, so the sufficient nonunit-gcd criterion does not apply to this presentation. That observation gives no conclusion about arbitrary meridional tuples or other diagrams. [Blair, Kjuchukova and Pfaff, Proposition 1.4 and Section 2.2](https://doi.org/10.1017/S030500412610214X).

## Constructive generation test

The four top cap meridians `x_1,x_2,x_3,x_4` generate the presented link group. For a three-element tuple, use one meridian on `A` and two on `B` to be compatible with the induced tuple after deletion of `A`. Deletion sends the `A` meridian to the identity, and the trefoil group is noncyclic, as witnessed by its quotient onto `S_3` from two meridians sent to adjacent transpositions.

In the disk-meridian presentation start with the ordered labels `x_1,x_1^(-1),...,x_4,x_4^(-1)`. Each cap connects one upward orientation around its puncture to the other, so the two consistently oriented disk meridians are inverse elements. At a positive crossing, incoming labels `(a,b)` become `(a b a^(-1), a)`, and at a negative crossing they become `(b,b^(-1) a b)`. The script introduces one new generator at each underpass and the corresponding relation. At each bottom cap the product of adjacent labels is set to the identity. This retains the actual meridians for all generators and gives a presentation with 452 generators and 452 relators. A constructive certificate may express all four `x_i` in the proposed three meridians. Each equality can be checked by a product of conjugates of the retained relators and free reduction.

For the two allowed triples of top caps, a Felsch Todd–Coxeter enumeration in GAP 4.15.1 reached its limit of 20,000 defined cosets. The four-cap control closed at index one. An enumeration closing at index one for a proposed triple would certify generation. Reaching the limit leaves the index undetermined. The trial tests two tuples of actual meridians of this 8-plat.

```sh
python3 research/meridional-rank/verify_interleaved_plat.py
python3 research/meridional-rank/verify_interleaved_plat.py --gap-input | gap -q
```

## Why the annular-twist component bound does not give the target

An annular twist of the fixed six-strand, thirteen-row Purcell–Zupan parent, with both annular boundary components retained, illustrates an allocation problem. The six cap meridians of its parallel components and one meridian of the unknot generate the parent group. The two parallel components each have bridge number three. At 37 annular twists, the catching disk has Euler characteristic `−1` and gives `b(K_37) ≥ ceil(37/72+1/2)=2` for the twisted component. Hyperbolicity of the fixed exterior rules out the essential-annulus alternative, and the nontrivial boundary knots rule out putting the annulus on a sphere. The resulting link has at least eight bridges by adding its component bounds. Its original allocation of one generator on the twisted knot cannot yield seven generating meridians. Filling the other components along the changed meridians would give a cyclic group for the nontrivial twisted knot. [Purcell and Zupan, Lemmas 3.2 and 3.4 and Corollary 4.3](https://arxiv.org/pdf/1512.03869v2), [Baker, Gordon and Luecke, Definition 1.2 and Theorem 1.3](https://msp.org/agt/2016/16-1/agt-v16-n1-p01-p.pdf).

More generally, if a link has a generating tuple of total length smaller than the sum of its component knot bridge numbers, one component gets fewer meridians than its own bridge number. Filling all other components sends the tuple onto that component's knot group. Any separation certified solely by this sum would already exhibit a single-component counterexample. The interleaved plat gives a whole-link bridge bound beyond that sum and leaves the constructive upper bound as the unresolved step.
