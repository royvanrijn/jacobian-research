# GENUSREVIEW1205: independent review of actual ramification nefness and bigness

Verdict: **CONFIRMED_AT_NAMED_FOUNDATIONS**
Review UTC: 2026-11-07T12:13:18.489274+00:00
Task: `actual-ramification-nef-big-review-sol56-20261107T1205Z`
Configured reviewer: gpt-5.6-sol / ultra. Actual hosted identity: **UNVERIFIED**.
Lifecycle effect: independent mathematical review only; no novelty finding, promotion, successor, or JC2 conclusion.

## Scope, custody, and reading

The pinned startup completed at 2026-11-07T12:03:24.008597+00:00 and wrote the immutable mode-0444 ACK. The copied reader SHA-256 was exactly `47e3aea5fd494d0966aad35312aa3eea689a000b05289bd03646b9e5d577fbe4`. Preflight verified all 17 manifested inputs, physical nonsymlink ancestry, hashes, LF counts, modes, and both original report seals.

Whole required reading completed through verified packet-reader windows: 1836/3500 non-STATE LF and 98/140 STATE LF. This includes MANIFEST 269; nine governors 1342; ADMISSION 45; METHOD 5; PROPOSAL 83; and GENUS_INTERFACE 92 LF. The two historical files physically displayed 175/240 LF. This is four more than the prompt's stated 171 (=79+92), because the verified manifest and reader report PROPOSAL as 83 LF; I charged all four final seal-metadata lines. Help output is charged conservatively as 30/100 LF; repeats and report previews were zero. No command-format correction was needed, and no reader output was truncated.

I accept `MINIMAL-RESOLUTION-RAMIFICATION-GENUS-1`, namely `p_a(R)>=1`, only at the supplied MANUAL/PROVISIONAL named-foundation scope. I did not re-audit that proof. The review used no web/source/history/repository search, CAS, secondary model, control family, degree farm, shared Git/docs, settings/auth/adapters, or protected-tree access.

## Exact statement checked

Let `f:Y -> P^2` be generically finite of degree `N>1`, with `Y` smooth irreducible projective and containing the whole `U=A^2` as an open set. Assume `f` is etale on `U`, `f(U)` misses a line `L`, and `f` contracts no smooth rational `(-1)`-curve. Put

`H=f^*L`, `K=K_Y`, `R=div(det(df))`, and `Q=2K+3H=K+R`.

The frozen proposal claims on this same model that `Q` and the full determinant divisor `R` are nef and big; their common null curves are exactly the `f`-contracted smooth rational `(-2)`-curves; `p_a(R)>=2`; and `h^0(Y,Q)>=2`.

These hypotheses apply to `g=h o q` for the finite projective normalization `h:X->P^2` of a hypothetical Keller map and its minimal resolution `q:Y->X`: `q` is an isomorphism over the smooth source `U`, `h` is finite, and a `g`-contracted curve is `q`-exceptional. Minimality therefore gives the required absence of contracted rational `(-1)`-curves. No smoothness, Q-Gorenstein property, or rational singularities of `X`, and no finiteness of `g`, are used.

## 1. Basic divisor identities

Because `Y` contains dense `A^2`, it is rational. Hence `chi(O_Y)=1`, `Pic^0(Y)=0`, and `K^2+rho(Y)=10`. The pullback `H` is Cartier and nef, with `H^2=N>0`. The determinant of the differential is a nonzero section of `O_Y(K+3H)`, so its full effective Cartier divisor satisfies

`R ~ K+3H`, `Q=K+R=2K+3H`, and `2R=Q+3H`.

This constructs `R` directly on smooth `Y`; it does not pull back a canonical divisor from singular `X`.

## 2. Nefness of Q and R

If `Q` were not nef, the surface cone theorem supplies a `K`-negative extremal ray on which `Q` is negative: on the `K`-nonnegative cone, `Q=2K+3H` is already nonnegative.

For a birational ray, its curve `E` is a smooth rational `(-1)`-curve. If `H.E=0`, then `f` contracts `E`, contrary to the hypothesis. Otherwise `H.E` is a positive integer and

`Q.E=-2+3H.E>=1`,

again impossible.

For contraction to a point, the smooth rational Picard-rank-one surface is `P^2`. Writing `H=m[ell]`, negativity of `Q=(-6+3m)[ell]` forces `m=1`, hence `N=H^2=1`, a contradiction.

For an extremal contraction `pi:Y->B` to a curve, rationality gives `B=P^1`. A general fibre `F` has `F^2=0` and `K.F=-2`. Hodge index excludes `H.F=0`; thus `H.F` is a positive integer, and `Q.F<0` forces `H.F=1`.

Every component of every fibre lies in the one contracted ray and is numerically proportional to `F`, hence has square zero. A connected reducible fibre would have a component meeting another, contradicting its zero intersection with the whole fibre. A multiple fibre `mC` is excluded by `1=H.F=mH.C`. Flatness then gives integral fibres of arithmetic genus zero, hence smooth `P^1` fibres. Thus `Y=F_n`.

With standard section `S`, write `H=S+bF`. Then

`R ~ K+3H = S+(3b-n-2)F`.

Adjunction gives `p_a(R)=0` for this divisor class, independently of the representative's reducibility or multiplicities. This contradicts the accepted `p_a(R)>=1` interface. All extremal cases are therefore excluded, so `Q` is nef.

Since `2R=Q+3H`, `R` is nef. Moreover

`R^2=(Q^2+6Q.H+9N)/4>0`,

so `R` is big. This uses only nefness of `Q,H` and `N>0`.

## 3. The numerically-zero endpoint

Assume `Q` is numerically zero. Then `K=-(3/2)H` numerically and `K^2=9N/4`. Integrality makes `4|N`; while `K^2=10-rho<=9` makes `N<=4`. Thus `N=4`, `K^2=9`, `rho=1`, `Y=P^2`, and `H=2[ell]`.

Divisor-class localization for `U=A^2` is

`O(Y)^* -> O(U)^* -> direct_sum_boundary Z[C_i] -> Pic(Y) -> Pic(U) -> 0`.

Both unit groups are `C^*` and `Pic(U)=0`, so the boundary-prime group maps isomorphically to `Pic(P^2)=Z`. There is exactly one boundary prime `C`, and its effective primitive class makes it a line. Since both `f^*L` and the determinant divisor are supported outside `U`, their classes force `H=2C` and `R=3C` as effective divisors.

Now `H` is ample, so `f` has no positive-dimensional fibre and is finite. Base-changing away from `L` gives a connected finite degree-four map

`P^2\C = A^2 -> P^2\L = A^2`.

Its determinant never vanishes there because `R` is supported on `C`, so it is etale. Complex `A^2` has no nontrivial connected finite etale cover. This contradiction excludes `Q` numerically zero. Notice that finiteness and ampleness are derived only inside this endpoint; they were not assumed for the resolved map.

## 4. Excluding Q^2=0

Suppose the nonzero nef divisor `Q` has `Q^2=0`. Hodge index gives `Q.H>0`, and

`K.Q=-(3/2)H.Q<0`.

For all large `n`, `(K-nQ).H<0`; since `H` is nef, `K-nQ` cannot be effective. Serre duality gives `H^2(O_Y(nQ))=0`. Riemann--Roch on rational `Y` then yields

`h^0(nQ) >= chi(O_Y(nQ)) = 1-(n/2)K.Q`,

which is unbounded.

Choose `n` with `h^0(nQ)>=2` and write `|nQ|=Z+|M|`, with `Z` fixed and `|M|` having no fixed component. Nefness and `Q.(M+Z)=nQ^2=0` imply `Q.M=Q.Z=0`. Two general members of `|M|` share no component, so `M^2>=0`. Hodge index on `Q^perp` gives `M^2<=0`, with equality only when `M` is numerically a positive multiple of `Q`. Hence `M^2=0`.

Two general members consequently do not meet, so `|M|` has no base points. Its morphism has one-dimensional image. After Stein factorization, a smooth connected general fibre `F` has `F^2=0` and is numerically a positive multiple of `Q`. Thus `K.F<0`; adjunction forces `F=P^1` and `K.F=-2`. But `Q.F=0` would give

`0=2K.F+3H.F=-4+3H.F`,

impossible because `H.F` is an integer. Therefore `Q^2>0`; the nef divisor `Q` is big. This establishes the movable-pencil step without importing semiampleness.

## 5. Common null curves

Let `C` be irreducible with `Q.C=0`. Since `Q^2>0`, Hodge index gives `C^2<0`. The identity `2K.C+3H.C=0` gives `K.C<=0`, while adjunction gives

`K.C=2p_a(C)-2-C^2>=-1`.

If `H.C>0`, integrality in `2K.C=-3H.C` in fact makes `H.C` even and `K.C<=-3`, a contradiction. Thus `H.C=K.C=0`. Adjunction then gives `p_a(C)=0` and `C^2=-2`, so `C` is a smooth rational curve contracted by `f`. Conversely, an `f`-contracted smooth rational `(-2)`-curve has `H.C=K.C=0`, hence `Q.C=0`.

Because `2R=Q+3H` and all three divisors are nef, `R.C=0` is equivalent to `Q.C=H.C=0`. Thus `Q` and `R` have exactly the same stated null curves.

## 6. Arithmetic genus and sections

The nef and big divisors `Q,R` satisfy `Q.R>0`. Cartier adjunction gives

`Q.R=(K+R).R=2p_a(R)-2`,

so `p_a(R)>=2`.

For the possibly nonreduced or disconnected effective Cartier divisor `R`, adjunction supplies

`0 -> O_Y(K) -> O_Y(K+R) -> omega_R -> 0`.

Rationality gives `H^0(O_Y(K))=H^1(O_Y(K))=0`, hence `H^0(Y,Q)=H^0(R,omega_R)`. Duality on the Cartier curve yields

`h^0(omega_R)=p_a(R)+h^0(O_R)-1>=2`.

No reducedness, SNC, or connectedness of `R` is used.

## 7. Numerical and blowup checks

For a smooth general member of `|H|` of genus `gamma`, adjunction and `K^2+rho=10` give exactly

`p_a(R)=1+K^2+(9/2)(H^2+K.H)=2+9gamma-rho(Y)`,

and direct expansion gives

`Q^2=16+24gamma-4rho(Y)-3N>0`.

For a boundary blowup `beta:Y'->Y` with exceptional curve `E`, determinant and canonical formulas give

`H'=beta^*H`, `R'=beta^*R+E`, and `Q'=beta^*Q+2E`.

Consequently `R'.E=-1`, `Q'.E=-2`, and adjunction gives `p_a(R')=p_a(R)-1`. Thus extra boundary blowups destroy nefness and change the genus exactly as claimed. They cannot replace the same relatively minimal model in this implication.

## Limitations and disposition

Every requested attack survives: the cone argument covers all extremal contraction types; the numerical-zero endpoint derives rather than assumes finiteness; the `Q^2=0` argument genuinely produces a basepoint-free pencil; null-curve adjunction uses the proven bigness; and the blowup formulas explain the strict minimal-model scope. The numerical identities also check exactly.

The conclusion remains relative to the accepted genus interface and the named standard foundations listed in the frozen inputs. It proves a necessary positivity constraint only. It does not prove properness or surjectivity of a Keller map, identify `U` with a whole finite inverse image, establish novelty, or resolve JC2.

**Final verdict: CONFIRMED_AT_NAMED_FOUNDATIONS.**
<!-- BODY-END -->

## Seal

- Body definition: every byte through the unique standalone `<!-- BODY-END -->` line,
  including its terminating newline; this seal is outside the body.
- Body bytes: `9887`.
- Body SHA-256:
  `036861c502752e326eea8852cb30dd46ec3f13c42c15845b938e5e6b3b279859`.
- Frozen basis: `27c4ab29317fa30501ddcecb9ba25b0ae6bca09c`.
