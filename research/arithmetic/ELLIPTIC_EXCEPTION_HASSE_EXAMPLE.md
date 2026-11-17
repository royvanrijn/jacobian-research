# A Hasse failure with a different quadratic resolvent on the elliptic exception

The elliptic exceptional fibre at `(pi,ell)=(128/25,9/10)` contains a genuine Hasse-failing target of the same fixed map `Phi`. Its quadratic component is `Q(sqrt(-81275190))`, distinct from `Q(sqrt(-3))`. This is an explicit existence result, not a claim that infinitely many points on the elliptic curve satisfy all local conditions.

## 1. Theorem EEH1: an explicit full Hasse-failing fibre

The rational target

\[
 y=\left(-\frac{128}{25},
 \frac{95820387033679170183314}{398297113921800046875},
 \frac{4611224705434499552560201273232737}
      {962454086534234237899900108800000}\right)          \tag{1}
\]

of the fixed degree-five determinant-one map `Phi` has a full finite étale fibre of rank five. That fibre has points over every completion of `Q` and no rational point. Its quadratic component is

\[
 K=\mathbb Q(\sqrt{d}),\qquad d=-81275190=-30\cdot2709173.
\]

**Proof.** Use `W=pi S` and the factorization chart of [the resolvent fibration](FIXED_MAP_RESOLVENT_FIBRATION.md), with

\[
 \pi=128/25,\quad\ell=9/10,\quad
 v=39976591281/36871680400.
\]

\[
q(W)=W^2+\ell W+v,
\]

and the monic cubic

\[
 h(W)=W^3-\frac{12}5W^2+
       \frac{85557356255259}{75513201459200}W-
       \frac{346044364289717257088531}{278429863018443603968000}.
                                                        \tag{2}
\]

Exact coefficient multiplication gives the inverse polynomial at (1):

\[
 P(W)=q(W)h(W),\quad
 P(W)=W^5-\frac32W^4+\frac3{2\pi^2}W^3+
       \left(-\frac5{4\pi}-\frac{9b}{32\pi^2}\right)W^2
       +\frac9{16\pi}W-\frac9{32}c.
                                                        \tag{3}
\]

The two factor discriminants have the common nonsquare class `d`:

\[
 \operatorname{Disc}(q)=d\left(\frac2{9601}\right)^2,
\]

\[
 \operatorname{Disc}(h)=d
 \left(\frac{684863008168809711969}
            {890975561659019532697600}\right)^2.       \tag{4}
\]

For an integral model, put

\[
 D_0=278429863018443603968000,\qquad
 H(T)=D_0^3h(T/D_0).
\]

Its exact integer coefficients are

\[
 \begin{aligned}
 H(T)={}&T^3-668231671244264649523200T^2\\
       &+87834695638335309390968046192139663650324480000T\\
       &-26826462523882131227875243325944959137572304355400352703388319744000000.
 \end{aligned}                                        \tag{5}
\]

Reduction of `H` modulo seventeen has no root. A cubic over a field is irreducible if and only if it has no root, so `H` and `h` are irreducible over `Q`. The negative discriminant in (4) makes the quadratic irreducible as well. The factors are coprime because their irreducible degrees differ, and both are separable in characteristic zero. The entire-fibre reconstruction therefore gives precisely

\[
 \Phi^{-1}(y)\simeq\operatorname{Spec}
        \bigl(K\times\mathbb Q[W]/h(W)\bigr),         \tag{6}
\]

with rank five, equal to the map's generic degree. Neither field factor is `Q`, so there is no rational point.

The exact discriminant factorization of the integral cubic is

\[
 \operatorname{Disc}(H)=
 -2^{65}3^7 5^{15}7^2 13^2 353^2 1321^2 2591^2
       9601^{16}\cdot2709173\cdot230703199^2.           \tag{7}
\]

Every prime on the right is prime, as verified by trial division in the accompanying certificate checker. It remains to exhibit local points at these eleven primes. All other finite primes are unramified in the cubic splitting field. By (4) the quadratic field is its unique quadratic resolvent, so the `S_3`-covering argument in [FMA2/FMA3](FIXED_MAP_AMBIENT_LOCAL_CLASSIFICATION.md) supplies a local point at every unramified prime.

At the following primes the quadratic field supplies a local point. Each table entry satisfies `z^2=d mod p` and `2z!=0 mod p`, so simple-root Hensel lifting applies.

| `p` | `z` |
|---|---:|
| 7 | 2 |
| 13 | 2 |
| 1321 | 412 |
| 9601 | 2888 |
| 230703199 | 89353083 |

At each remaining prime, the cubic supplies the point. The integers below have the displayed **exact** valuations for the polynomial (5).

| `p` | `x` | `v_p(H(x))` | `v_p(H'(x))` |
|---|---:|---:|---:|
| 2 | 2129920 | 47 | 21 |
| 3 | 2 | 5 | 2 |
| 5 | 3775 | 12 | 5 |
| 353 | 13 | 1 | 0 |
| 2591 | 2376 | 1 | 0 |
| 2709173 | 2124817 | 1 | 0 |

Every row satisfies `v_p(H(x))>2v_p(H'(x))`. The strict Newton/Hensel criterion gives a `Q_p`-root of `H`, hence of `h` by `W=T/D_0`. The exact imported local lifting statement is Milne, *Algebraic Number Theory*, version 3.00, Theorem 7.32, printed p.114, in [the author-hosted primary text](https://www.jmilne.org/math/CourseNotes/ANT300.pdf), with [retained copy](sources/milne-ant300.pdf). The cubic has a real root because its degree is odd. These certificates prove local solubility at every place and complete the proof. ∎

## 2. Location on the elliptic exception and limits

The base `(pi,ell)=(128/25,9/10)` is the noncentered elliptic degeneration of Proposition 3.7 in the resolvent-fibration note. Equations (4) give a nonzero rational square root of `Disc(q)Disc(h)` and therefore a rational point of its matching-resolvent cover. Dividing by the repeated linear factor gives its elliptic normalization; the checker verifies the corresponding exact point on

\[
 T^2=X^3-10758382X+6600759744,
 \qquad X=6400v/3-3984.
\]

No assertion about a Mordell–Weil basis, full rank calculation, or the density of Hasse points is needed for this example. The quadratic field differs from `Q(sqrt(-3))` because `d/(-3)=27091730` is not a rational square; its two-adic valuation is odd. The first coordinate in (1) is neither `-1` nor `2`, so the target lies outside the two previously counted Hasse-failing lines. Its nonintegral first coordinate means it does not settle the integral-target Hasse-existence problem.

## 3. Exact verification

The data are [elliptic_exception_hasse_certificate.json](elliptic_exception_hasse_certificate.json). Run

```sh
python3 research/arithmetic/verify_elliptic_exception_hasse_example.py --output research/arithmetic/elliptic_exception_hasse_checks.json
```

The verifier uses only Python's standard library: integer arithmetic, rational fractions, polynomial coefficients, trial division, and modular evaluations. It derives the factorization and integral model from the supplied chart parameters; checks both discriminant squareclasses and the entire discriminant factorization; proves cubic irreducibility by exhaustive reduction modulo seventeen; and checks every local witness and complete coverage of the discriminant primes. No PARI/GP computation or unverified change of primitive element is used in verification.

These checks certify this explicit example under the retained whole-fibre reconstruction and classical local/Galois criteria. They do not formally verify those general inputs and are not external mathematical review. No claim about infinitely many Hasse-failing points on this elliptic fibre follows from the example.
