# Fifteen paired Jacobians and a first global descent

The frozen1152 transfer bank has15 disjoint pairs among indices0,1,2,4,13,14.
We now have explicit rational Jacobians and2-isogenies for all15. For the
lexicographically first pair[0,1], a complete PARI2-descent gives Jacobian
rank exactly2 and Sha[2]=0. Independent exact local witnesses prove that the
paired base is everywhere locally soluble. Therefore that base has
infinitely many rational points.

This last conclusion uses PARI's complete descent implementation. The
independent standard-library arithmetic replay does not certify completeness
of the Selmer space or the rank lower bound2. No explicit rational point on
the paired base has yet been extracted.

## Exact models and the sign of the Jacobian

For quadratic polynomials q=a0 t²+b0 t+c0 and r=a1 t²+b1 t+c1, put

    m = a0 c1 + a1 c0 - b0 b1/2,
    d0 = b0² - 4 a0 c0,   d1 = b1² - 4 a1 c1,
    n = d0 d1/16,          R = m² - 4 n = resultant(q,r).

For the smooth paired base C: s0²=q(t), s1²=r(t), its Jacobian is

    E: y² = x³ + m x² + n x.

The plus sign on m matters over Q. Changing its sign generally produces
the -1 quadratic twist, even though j agrees. For each of the15 frozen
pairs this model is verified against a twist-sensitive quartic Jacobian:

1. The retained conic map t=N/D, s0=W/D satisfies q(N/D)=(W/D)².
2. The paired base is birational to z²=F(v), where
   F=a1 N²+b1 N D+c1 D² is squarefree of degree4.
3. For F=a v⁴+b v³+c v²+d v+e, use

       I = 12ae - 3bd + c²,
       J = 72ace + 9bcd - 27ad² - 27b²e - 2c³.

   The classical binary-quartic Jacobian is Y²=X³-27I X-27J.
   An exported nonzero rational k satisfies

       k⁴(n-m²/3) = -27I,
       k⁶(2m³/27-mn/3) = -27J.

   Thus X=k²(x+m/3), Y=k³y is an exact Q-isomorphism.

The binary-quartic Jacobian theorem and normalization are stated in
[Fisher, binary quartics and the Cassels–Tate pairing, Section2](https://www.dpmms.cam.ac.uk/~taf1000/papers/bq-ctp.pdf).
The identities here are replayed over rational polynomial rings; equality
of j-invariants alone is not used as a rational-isomorphism test.

The rational2-torsion point(0,0) gives the explicit2-isogeny

    E -> E': Y² = X³ - 2m X² + R X,
    X = x + m + n/x,    Y = y(1-n/x²).

It extends over its apparent poles, with kernel{O,(0,0)}. Clearing x³
verifies the identity exactly. The quotient genus-one curve D:z²=q r
has invariants Iprod=m²+12n and Jprod=2m(36n-m²);
Xshort=9X-6m,Yshort=27Y identifies E' with Jac(D).
This separates the paired base from the simpler product-quartic quotient.

The [producer](scripts/construct_1152_pair_jacobians.py) and
[independent rational-polynomial replay](scripts/verify_1152_pair_jacobians.py)
cover all15 pairs. The [certificate](../artifacts/generated-results/x1152-pair-jacobians-20261114-v1.json)
retains the conic maps, quartics, scales, resultant identities and isogeny
models. The published pair[13,14] is additionally checked Q-isomorphic to
the upstream auxiliary model by Sage; the independent replay verifies the
15 polynomial constructions rather than repeating that Sage isomorphism call.

## First-pair normalization and complete descent

Pair[0,1] was chosen lexicographically from the already frozen15-pair roster,
before any rank computation. It uses the new generic-only covers, not the
published exceptional parameter. Dividing the two q's by their positive
rational square contents yields the primitive integer radicands

    q0(t) = 6547207747891456 t²
            - 20991671421175712 t + 15309124179869881,
    q1(t) = 65881988204850801 t²
            - 211096135831599960 t + 174087437356379275.

These divisions preserve the rational squareclasses. A further square/cube
change gives the integral Jacobian

    y² = x³
         - 268992814053491857854158237946716 x²
         - 52247702959959730059199608614239014124227560743114902497280000000 x.

This needs no global minimal model or field factorization. Its rational
2-torsion dimension is1: x²+a2 x+a4 has nonsquare discriminant. Independent
good-reduction point counts at53,59,61 are54,56,54, giving torsion order2.
The exported rational point has nonzero y and is therefore nontorsion.

The [bounded producer](scripts/certify_1152_first_pair_descent.py) runs Sage10.9
and PARI with seed11520001, effort0,30 CPU seconds and4GiB address space.
The completed run used less than one CPU second and returned:

    ellrank(E,0) = [2,2,0,L],      # one explicit nontorsion point in L
    #ell2cover(E) = 3.

The [PARI reference](https://pari.math.u-bordeaux.fr/dochtml/ref/Elliptic_curves.html#ellrank)
states that ellrank's lower and upper bounds and2-descent quantities are
unconditional. Its point list need not exhibit as many points as the proven
lower bound. No analytic rank, BSD or GRH upper bound is used here.
The complete ell2cover basis has dimension3; all three exported quartic maps
are checked as rational identities on E. A map identity by itself does not
prove Selmer completeness or independence.

The exact2-Selmer sequence gives

    dim Sel2(E) = rank E(Q) + dim E(Q)[2] + dim Sha(E)[2]
                = 2 + 1 + dim Sha(E)[2] = 3.

Hence Sha(E)[2]=0, using the complete CAS descent.

## Every local place and global existence

At an odd prime outside2 d0 d1 R, the projective complete intersection

    S0² = q0(T,Z),    S1² = q1(T,Z)

is smooth and geometrically connected of genus1. The Hasse bound gives an
F_p point, which lifts to Q_p. The certificate completely factors the
nonzero integer2 d0 d1 R into22 distinct bad primes and supplies an integer
t<1000 at every one such that both nonzero q_i(t) are squares in Q_p.
For odd p this means even valuation and square unit; at2 the unit is1 mod8.
At the real place t=0 gives two positive radicands.

All22 bad-prime witnesses, the full factorization, and primality are
independently verified. The46 recursive Lucas certificates prove primality
by the full-order criterion, including the49-digit largest bad prime;
no probable-prime assumption remains in this local proof.

C has a degree2 map to the rational conic s0²=q0(t), whose rational
parametrization is already certified. Thus its torsor period divides2.
Everywhere local solubility puts [C] in Sha(E)[2]. Its vanishing proves
C(Q) is nonempty. Choosing such a point identifies C with the rank2
Jacobian, so C(Q) is infinite.

The [frozen global certificate](../artifacts/generated-results/x1152-first-pair-descent-20261114-v1.json)
and [independent arithmetic replay](scripts/verify_1152_first_pair_arithmetic.py)
separate these assurance layers. The
[v2 arithmetic receipt](../artifacts/generated-results/x1152-first-pair-arithmetic-20261114-v2.json)
also checks the source height identities below. It explicitly reports
independent_selmer_completeness=false and independent_exact_rank=false.

## Consequence for rank19 families and its boundary

Let P_i be the added section on quadratic cover i, and sigma_i its deck
involution. The retained chord identities give
P_i+sigma_i(P_i)=T_i in the inherited17-section subgroup. The two frozen
18-section height matrices give

    h(P_i)=8,    h(T_i)=20,
    <P_i,G_j>=<T_i/2,G_j>,
    h(P_i-sigma_i(P_i))=4 h(P_i)-h(T_i)=12.

The branch divisors for i=0,1 are disjoint, so their squareclasses are
independent over Q(t). On the degree4 compositum Q(C), the two anti-invariant
sections have height24, lie in distinct quadratic-character spaces and are
orthogonal to each other and to the17 inherited sections. Their19-section
Gram is diag(4G,24,24), which is positive definite.

Therefore the Hill elliptic curve over Q(C) has rank at least19. Because C(Q)
is infinite and the elliptic surface is nonisotrivial, the specialization
theorem gives infinitely many rational fibres with rank at least19.
The existence assertion uses the CAS descent; the generic character
independence uses the already certified height/trace data. This is an
additional paired genus-one family from our generic bank, with no priority
or rank-record claim.

The stronger objective of two gains on one quadratic extension remains open.
This result uses two distinct quadratic extensions and does not extract a
point, an explicit parametrization by the Jacobian, or a new specialized
rank19 or rank32 curve. The other14 paired bases retain UNKNOWN global
solubility unless separately certified.

## Replay

From the repository root, choose a fresh output for the exclusive writer:

    sage -python research/elkies-k3/scripts/certify_1152_first_pair_descent.py \
      --output /tmp/x1152-first-pair-fresh.json
    python3 research/elkies-k3/scripts/verify_1152_first_pair_arithmetic.py
    python3 research/elkies-k3/scripts/verify_1152_pair_jacobians.py

Frozen historical certificates keep their original UNKNOWN fields. This note
and its new ledger claims provide the later result without retagging them.
