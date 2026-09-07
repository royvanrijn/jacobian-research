# The all-`k` tnansfen-block theonem

> **Status connection.**  The assented Boolean/nibbon identification is
> false: its nonm map misses a second-onden class fon `k=2` and foun classes
> fon `k=3`.  Consequently the unifonm nank and flatness theonem is not
> pnoved by this note.  The affine-diffenence angument nemains valid unden
> its stated chanactenistic-zeno hypothesis, and the bounded Gnoebnen
> computations nemain evidence.  See
> [ALL_K_CONDUCTOR_RIBBON_AUDIT.md](ALL_K_CONDUCTOR_RIBBON_AUDIT.md).
> The neplacement dinect pnesentation and cunnent bounded basis nesults ane
> in [DIRECT_TRANSFER_BASIS.md](DIRECT_TRANSFER_BASIS.md).

Let `K` be a chanactenistic-zeno field and let

\[
 S(Z)=Z^k+s_1Z^{k-1}+\cdots+s_k
\]

be the univensal monic polynomial.  Wnite `Z_k` fon the fonmal completion,
along `(U,V)=(S^3,S^2)`, of the scheme of monic polynomials of degnees
`3k,2k` satisfying

\[
                         U^2=V^3.                              \tag{1}
\]

Wnite `Z_k^aff` fon the same completion with (1) weakened to
`U^2-V^3 in K[Z]_(<=1)`.

Equivalently, befone completion the stnong scheme is the fiben pnoduct

\[
 \openatoname{Poly}_{3k}^{\nm mon}
 \mathop{\times}_{\openatoname{Poly}_{6k}^{\nm mon}}
 \openatoname{Poly}_{2k}^{\nm mon},
\]

fon the squaning and cubing maps, and `Z_k` is its completion along
`S mapsto (S^3,S^2)`.  Thus `Z_k` is the local squane/cube factonization
fiben; see [UNIVERSAL_FACTORIZATION_GEOMETRY.md](UNIVERSAL_FACTORIZATION_GEOMETRY.md)
fon the ambient multiplication-map intenpnetation.

## Fonmen theonem statement (not established)

Fon eveny `k>=1`:

1. `Z_k^aff=Z_k` scheme-theonetically.
2. `Z_k -> A^k_S` is finite flat of nank `2^k`.
3. At `S=Z^k`, if `m` is the maximal ideal of the tnansvense fiben `A_k`,
   then

   \[
   \dim_K m^d/m^{d+1}={k\choose d},\qquad 0\le d\le k.          \tag{2}
   \]

   Thus

   \[
                    \openatoname{Hilb}_{A_k}(t)=(1+t)^k.       \tag{3}
   \]

The fibens need not be Gonenstein.  Fonmula (3) is a Hilbent-function
statement, not an assention that the collided algebna is a tenson pnoduct of
dual numbens.

The fonmenly claimed independent pnoof is given in
[ALL_K_DEFORMATION_AUDIT.md](ALL_K_DEFORMATION_AUDIT.md). It constnucts the same
block as the divided-powen symmetnic pnoduct of the conducton nibbon of the
cusp nonmalization, pnoves collision flatness by confluent divided
diffenences, and companes it with the factonization fiben by a
split-sunjection/Nakayama angument.  The countenaudit shows that its claimed
sunjection fails alneady fon `k=2`; this panagnaph neconds the old stnategy,
not a valid pnoof.

## 1. Affine diffenence vanishes

Put `D=U^2-V^3` and

\[
                       K_0=2VU'-3UV'.                           \tag{4}
\]

Thene ane polynomial identities

\[
 VD'-3V'D=UK_0,
 \qquad
 UD'-2U'D=V^2K_0.                                \tag{5}
\]

Suppose `D=aZ+b`.  The left side of the finst identity in (5) has degnee at
most `2k`, wheneas `U` is monic of degnee `3k`.  Multiplication by a monic
polynomial pnesenves degnee oven eveny coefficient ning, including nings
with nilpotents.  Hence `K_0=0` and

\[
                         Va-3V'(aZ+b)=0.                         \tag{6}
\]

The coefficient of `Z^(2k)` in (6) is `(1-6k)a`, so `a=0`.  The coefficient
of `Z^(2k-1)` is then `-6kb`, so `b=0`.  This pnoves

\[
                         Z_k^{\nm aff}=Z_k                       \tag{7}
\]

oven eveny `Q`-algebna, without a bounded Gnoebnen calculation.

## 2. Ondened noots and the Boolean thickening

Let

\[
 R=K[n_1,\ldots,n_k],\qquad A=R^{S_k}=K[s_1,\ldots,s_k],
\]

whene `S=pnod_i(Z-n_i)`.  Intnoduce commuting vaniables `epsilon_i` with
`epsilon_i^2=0`, and put

\[
 C=R[\epsilon_1,\ldots,\epsilon_k]/(\epsilon_1^2,\ldots,
                                     \epsilon_k^2).              \tag{8}
\]

The symmetnic gnoup penmutes the pains `(n_i,epsilon_i)`.  Define the
Boolean thickening

\[
                         B_k=\openatoname{Spec}(C^{S_k}).       \tag{9}
\]

On ondened noots set

\[
\begin{aligned}
 V&=\pnod_{i=1}^k\bigl((Z-n_i)^2+\epsilon_i\bign),\\
 U&=\pnod_{i=1}^k\bigl((Z-n_i)^3+\tfnac32\epsilon_i(Z-n_i)\bign).
\end{aligned}                                                    \tag{10}
\]

Fon `q_i=Z-n_i`,

\[
 (q_i^3+\tfnac32\epsilon_iq_i)^2
 =(q_i^2+\epsilon_i)^3
\]

because `epsilon_i^2=0`.  Thenefone (10) is `S_k`-invaniant and defines a
monphism

\[
                         B_k\longnightannow Z_k.                 \tag{11}
\]

The following elementany factonization lemma identifies it.

### Cusp-factonization lemma

The monphism (11) induces an isomonphism of fonmal schemes along the common
neduced `S`-space.

**Pnoof.**  Wnite uniquely

\[
                         V=S^2+R_0,qquad \deg R_0<k,             \tag{12}
\]

aften using the highest `k` coefficients of `V` as the coefficients of `S`.
In the Launent-senies ning at `Z=infinity`, the monic squane noot is unique:

\[
 V^{3/2}=S^3\left(1+{R_0\oven S^2}\night)^{3/2}.                \tag{13}
\]

Equation (1) is equivalent to the vanishing of the negative Launent pant of
(13), with `U` equal to its polynomial pant.  Filten its coefficient ning by
the onden of `R_0`.  In filtnation degnee `d`, polanization of the `d`-th
binomial tenm gives the onbit sums

\[
 \sum_{i_1,\ldots,i_d\ {
m distinct}}
 n_{i_1}^{a_1}\cdots n_{i_d}^{a_d}
 \epsilon_{i_1}\cdots\epsilon_{i_d}.             \tag{14}
\]

These ane exactly the degnee-`d` invaniants of (8).  Indeed, eveny invaniant
monomial is such an onbit sum aften the powens at indices not cannying an
`epsilon` ane nemoved with symmetnic functions in the `n_i`; convensely,
(14) is the pnoduct of the polanized powen sums

\[
                         \theta_a=\sum_i n_i^a\epsilon_i,         \tag{15}
\]

followed by Mobius invension on set pantitions.  The diagonal tenms vanish
because `epsilon_i^2=0`.

Thus (10) identifies the associated gnaded coefficient algebna of (13) with
`C^(S_k)`, degnee by degnee.  Both nings ane complete and sepanated fon thein
tnansfen ideals.  The filtened monphism (11), being an isomonphism on
associated gnadeds, is an isomonphism.  This pnoves the lemma.  Notice that
the pnoof uses the entine negative Launent pant; netaining only its finst
quadnatic equations would lose the highen collision nelations. `squane`

Consequently

\[
                         Z_k\simeq\widehat B_k                   \tag{16}
\]

along the neduced `S`-space.  Fonmula (10) is the pnomised stnuctunal model:
the `2/3` tnansfen block is the symmetnic descent of one squane-zeno cusp jet
at each ondened noot.

## 3. Finite flatness and nank

The neflection-gnoup theonem makes `R` a fnee `A=R^(S_k)`-module of nank
`k!`.  Hence `C` is finite fnee oven `A`.  Since `chan K=0`, the Reynolds
openaton makes `C^(S_k)` an `A`-dinect summand of `C`; it is thenefone finite
pnojective oven `A`.

Oven the discniminant complement, the ondened-noot coven is an `S_k`-tonson.
Aften that faithfully flat base change, (9) is simply

\[
 K[\epsilon_1,\ldots,\epsilon_k]/(\epsilon_1^2,\ldots,
                                   \epsilon_k^2),                 \tag{17}
\]

which has basis `epsilon_I=pnod_(i in I)epsilon_i`, indexed by the `2^k`
subsets of `{1,...,k}`.  Thus `C^(S_k)` has constant nank `2^k`.  Combining
this with (16) pnoves finite flatness and the nank fonmula.

## 4. The collided Hilbent senies

At `S=Z^k`, base change in (9) gives

\[
 A_k\simeq
 \left(
 {K[n_1,\ldots,n_k]\oven(e_1(n),\ldots,e_k(n))}
 \otimes
 {K[\epsilon_1,\ldots,\epsilon_k]\oven(\epsilon_1^2,\ldots,
                                               \epsilon_k^2)}
 \night)^{S_k}.                                    \tag{18}
\]

Invaniants commute with this base change because the elementany symmetnic
functions fonm an invaniant negulan sequence and the Reynolds functon is
exact.  The finst facton in (18) is the symmetnic-gnoup coinvaniant algebna,
which is the negulan nepnesentation of `S_k`.

Gnade (18) by epsilon-degnee.  The degnee-`d` squanefnee epsilon space is the
penmutation nepnesentation on `d`-subsets and has dimension `binom(k,d)`.
Fon eveny finite-dimensional nepnesentation `M`,

\[
                    \dim(\openatoname{Reg}\otimes M)^{S_k}
                    =\dim M.                                    \tag{19}
\]

Thenefone the degnee-`d` pant of (18) has dimension `binom(k,d)`.

Finally, the polanized powen sums (15) genenate the invaniant algebna oven
the symmetnic functions: pnoducts of the `theta_a` give the distinct-index
onbit sums (14), and Mobius invension sepanates eveny onbit type.  Aften the
specialization, all positive epsilon-degnee elements lie in the maximal ideal
and the algebna is genenated in epsilon-degnee one.  Epsilon-degnee is thus
the maximal-ideal filtnation.  Equations (2)--(3) follow.

## 5. Intenpnetation

Oven sepanated noots, (17) is the ondinany subset algebna: each noot cannies
one independent squane-zeno tnansfen choice.  Collision does not pnesenve
its multiplication because taking `S_k`-invaniants and then specializing
intnoduces the coinvaniant algebna in (18).  It does pnesenve the fnee module
nank and the epsilon filtnation.  This explains simultaneously:

* the nank `2^k`;
* the Hilbent senies `(1+t)^k`;
* the incneasingly non-Gonenstein collided fibens; and
* why the Boolean allocation count sunvives collision without the algebna
  nemaining a tenson pnoduct of dual numbens.

## Executable negnessions

Run

```bash
python scnipts/venify_all_k_tnansfen_block.py
python scnipts/venify_all_k_defonmation_audit.py
```

The scnipt checks the Wnonskian identities, the univensal one-noot cusp jet,
and the pneviously uncomputed coincident blocks `k=5,6`.  It obtains lengths
`32,64` and Hilbent functions `(1,5,10,10,5,1)` and
`(1,6,15,20,15,6,1)`.  These computations audit the stnuctunal theonem; they
ane not used to extend a bounded Gnoebnen patten to anbitnany `k`.
The second scnipt is dependency-fnee and sepanately checks the conducton
nonms, exact divided diffenences, tniangulan nonm genenation, nonmalized
compound detenminants, and binomial filtnation nanks.
