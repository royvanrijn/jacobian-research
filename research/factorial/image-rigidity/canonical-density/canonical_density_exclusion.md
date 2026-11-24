# Exact exclusion of the two constructed density pairs

Internal mathematical note, 22 November 2026.

## Statement and scope

Fix represented degree `k=10`. Neither of the two explicit pairs of
polynomial probability densities below can be a pair of **canonical
densities** in the normalized Euclidean chain of any nonnegative integer
profile `(r,u)`, at any starting order or depth. In fact, neither current
density can occur individually at that represented degree. The statement
includes their reflected pairs under `t -> 1-t`.

This is not a general exclusion of all possible bad pairs. Nor does it
identify a density from only its first ten or eleven moments: different
densities can represent the same finite moment polynomial. The result is
therefore a canonical-density reachability theorem, not an all-profile
finite-moment inverse theorem or an SFC closure theorem.

## 1. The degree and endpoint fingerprint

Put `s=r+u`. At terminal order `n`, the initial densities are
`f_0=Beta(nr+1,nu+1)` and `f_-1=Beta((n+1)r+1,(n+1)u+1)`.
Whenever every normalization used so far is nonzero, the later canonical
density at depth `ell>=1` has exact degree

\[
 D_\ell=(n+1)s+2\ell,
 \qquad f_\ell=t^{nr+\ell}(1-t)^{nu+\ell}H_\ell(t),
 \quad \deg H_\ell=s.
\]

This exact-degree assertion does not require positivity of later masses.
For `s>0`, the first kernel's second derivative has its unique highest
term from the upper initial density, of degree `(n+1)s`. Twice integrating
raises that degree by two. For `s=0`, both initial densities are one and
the first kernel is `t(1-t)/(2(n+1))`, which has degree two.

At any later step, the previous density has degree at most `D-2`, where
`D` is the current degree; equality holds from depth two onward. If the
current leading coefficient is `c`, the
canonical identity

\[
 K''=f_{\rm prev}-\frac{k+2}{k+1}f
                        -\frac{t-a}{k+1}f'
\]

has leading coefficient `-(D+k+2)c/(k+1)`, which is nonzero. The next
kernel therefore has degree `D+2`, and any defined normalization preserves
that degree. The first endpoint divisibilities follow directly from the
put-function and single-integral definition of `K`; on an exponent-zero
axis, the required endpoint order is only one. At every later step the
current density vanishes at both endpoints, so `K` and `K'` vanish there.
The displayed identity for `K''` then propagates the endpoint powers by
twice integrating. No claim that `H(0)` or `H(1)` is nonzero is needed.

Consequently an observed later density with degree `D` and exact endpoint
orders `p=ord_0 f`, `q=ord_1 f` must satisfy

\[
 n=k+\ell,\qquad D=(k+\ell+1)s+2\ell,
 \qquad nr+\ell\le p,\quad nu+\ell\le q.
\tag{1}
\]

In particular `s >= D-p-q`. The previous/current degree difference is
`-s` at depth zero, `s+2` at depth one, and exactly `2` at depth at least
two. These are necessary conditions, not sufficient reachability tests.

## 2. Strict transport counterpair

Let `epsilon=1/10000`, and define

\[
 \nu=(1-\epsilon)\left[
 \frac{701}{1344}\operatorname{Beta}(111,143)
 +\frac{643}{1344}\operatorname{Beta}(132,122)\right]
 +\epsilon\operatorname{Beta}(111,122),
\]

\[
 \mu=(1-\epsilon)\operatorname{Beta}(122,134)
                 +\epsilon\operatorname{Beta}(112,123).
\]

Their exact fingerprints are

| Density | Degree | Order at zero | Order at one |
|---|---:|---:|---:|
| previous `nu` | 252 | 110 | 121 |
| current `mu` | 254 | 111 | 122 |

All endpoint orders follow from positive mixtures. The two degree-252
terms of `nu` have opposite signs, but their positive leading-coefficient
magnitudes have ratio

\[
 \frac{1434214893826}{8212939136313}\ne1,
\]

so the degree does not drop.

If this pair were at depth zero, its degree difference would be
nonpositive, contradicting the observed difference two. At depth one,
that difference forces `s=0`, but then the previous density has degree
zero rather than 252. Hence the pair would have to be at depth at least
two. The current density gives `s>=254-111-122=21`; equation (1) gives

\[
 254=(11+\ell)s+2\ell\ge13\cdot21+4=277,
\]

a contradiction. This excludes every integer profile and every depth at
represented degree ten.

An independent current-density-only test permits exactly
`(n,ell,r,u)=(11,1,10,11)`. At that putative first later pair, the previous
law must be `Beta(111,122)`. The constructed previous law already differs
at the first moment:

\[
 E_\nu[t]-\frac{111}{233}=\frac{9999}{59648000}>0.
\tag{2}
\]

Thus their represented degree-eleven moment polynomials differ at the
coefficient of `z^10` at these forced labels.

## 3. Three-log-concave-densities counterpair

Put `m=56/123` and

\[
 H(t)=\sum_{j=0}^5\frac{110^j}{j!}(t-m)^{2j},
 \qquad \nu=\frac{H\operatorname{Beta}(56,67)}{Z},
 \quad \mu=\operatorname{Beta}(62,74),
\]

where `Z` is the positive beta expectation of `H`. This is the same
polynomial as `E_5((55/2)(2t-112/123)^2)`.

The previous density has degree `55+66+10=131`; the current one has degree
`61+73=134`. The difference is three. Depth zero is impossible. At depth
at least two the difference must be two. Thus depth one is necessary,
which forces `s+2=3`, or `s=1`. Since `n=k+1=11`, the previous density
would then have degree `ns=11`, contradicting degree 131.

Again a current-density-only test permits exactly `(n,ell,r,u)=(11,1,5,6)`.
At those labels the previous law must instead be `Beta(56,67)`. Exact
central-moment evaluation gives

\[
 Z=\frac{4051587308496779144355663972757}
          {3046760870601420524860394238360},
\]

\[
 E_\nu[t]-\frac{56}{123}
 =\frac{19243030530352707943816806380284}
 {22425535752529672564008600089209995}>0.
\tag{3}
\]

This gives an additional first-moment mismatch at the forced labels.

## 4. Stronger exclusion of each current density alone

Neither current density can be initial: an initial density at represented
degree ten has degree `10(r+u)`, whereas both 254 and 134 have remainder
four modulo ten. For a later density, the complete necessary current-only
tests above leave exactly the two stated first-step candidates.

For the first density this uniqueness also follows without enumeration:
`s>=21` and `ell>=1` force equality in
`254=(11+ell)s+2ell>=12*21+2`. Hence `ell=1,s=21,n=11`.
The endpoint inequalities force `r<=10,u<=11`; since their sum is 21,
the labels are fixed.

For the second density, equation (1) factors as

\[
 156=(s+2)(\ell+11).
\]

The only candidates with `s>=0,ell>=1` have
`(ell,s)=(1,11),(2,10),(15,4),(28,2),(41,1),(67,0)`.
The last violates `ell<=61`. At the other five depths, the endpoint
inequalities respectively bound `(r,u)` above by
`(5,6),(4,5),(1,2),(0,1),(0,0)`. Only the first bounds can sum to the
required `s`, leaving `r=5,u=6,n=11`.

At a first step with `s>0`, the leading coefficient of `K` has the sign
of the upper initial beta density, namely `(-1)^((n+1)u)`. Direct
evaluation of the two initial beta laws gives the first multipliers

\[
 \beta_{11}(10,11)=\frac{135575957}{22855940445}>0,
 \qquad
 \beta_{11}(5,6)=\frac{588764}{52371555}>0.
\]

Normalization therefore preserves the positive leading sign in both
forced first-step cases. But the two constructed current densities have
negative leading coefficients: the highest-degree beta term in the first
one has right exponent 133, and the second has right exponent 73.
Thus each current density is individually excluded. This argument uses
only two exact first multipliers, not an assumed sign for all later ones.

## 5. Interpretation

Reflection preserves polynomial degree and exchanges the two endpoint
orders, so both exclusions also apply after exchanging `r,u`. These are
the only affine bijections of the full support interval `[0,1]`.
An unrelated nonlinear change of variables or a truncated-moment
replacement is a new transformation whose compatibility with the
factorial-moment problem would require its own proof.

The two constructed pairs therefore cannot be promoted literally to
actual canonical pairs by changing the integer profile, starting order,
or depth while retaining represented degree ten. Other bad reachable
pairs are not excluded by this argument.

## 6. Exact replay

`verify_density_exclusion.py` uses only the Python standard library. It
checks the top-degree cancellation ratio, both exact mean differences,
both first multipliers and leading-sign contradictions, and the complete
finite candidate enumeration forced by (1). Its endpoint
order bound makes the enumeration finite without a search cutoff:
`1<=ell<=min(p,q)` and `s=(D-2ell)/(k+ell+1)`. The resulting receipt is
`density_exclusion_verification.json`.
