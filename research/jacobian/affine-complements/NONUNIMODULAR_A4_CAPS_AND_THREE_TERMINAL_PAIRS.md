# Nonunimodular A4 caps and three terminal pairs

15 November 2026. PF2RULED89 closes every remaining A4 incoming
trajectory whose first simple-S4 arrival is the FINAL cable stage.
This holds at arbitrary companion depth, without prefix eta, in
both source pencils (the common pencil was already closed by PF87).
For complementary sources with later simple tails, an unrestricted
cap/last-gcd reduction leaves only the ordered source terminal pairs

    (winding,free-tail)=(2,4),(2,8),(4,10).              (1)

These are necessary cases, not compatible finite maps. The first
pair requires that the A4 companion stayed solid until the arrival;
the other two require an earlier nonsolid cyclic3 cap. The full
quartic/higher-degree goal remains OPEN.

## 1. Two source lemmas valid at every characteristic depth

Retain the actual smooth exact TWO-line complementary whole charts
of PF77/79, including the unramified D chart with sigma>=1, pole p,
support j>=sigma+1 and gcd(sigma,support)=1. Its full contact formula is

    n_D=S_D+2p-sigma+1,
    S_D=sum_(sigma conjugates other than1)(2p-j_detected).

Since j_detected<=p, the GENERAL lower bound is

    n_D >= (sigma+1)*p-sigma+1.                        (2)

Also n_D is EVEN. For a characteristic gcd drop d->e=gcd(d,j),
(d-e)*j equals d-e modulo2: an odd difference requires j odd.
Summing over the drops gives S_D=sigma-1 modulo2, proving the claim.
No conductor theorem or bounded extraction is needed.

Now suppose the source has at least two nodes and its last residual
gcd is a>=2. Truncate the support strictly ABOVE its final
gcd-dropping index j_*=sigma-1+kappa, and divide sigma,p and all
retained indices by a.
The resulting primitive chart has valid indices above sigma'=sigma/a,
gcd1, pole p'=p/a and diagonal n_D'. Lower noncharacteristic terms
can be retained; they do not affect the identity. Counting the
conjugates before and at the last drop proves EXACTLY

    n_D=a^2*(n_D'-1)+a-(a-1)*kappa.                    (3)

Before the last drop their total contact is a^2*S_D'; the a-1
remaining conjugates have contact2p-j_*. Substitution gives(3).
For ANY retained two-line source parameters n_R,q, put
delta=n_R*n_D-q^2*p^2 and delta'=n_R*n_D'-q^2*p'^2. Then

    delta=a^2*delta'-n_R*(a-1)*(a+kappa).               (4)

The stripped chart is an actual primitive LOCAL chart. We do NOT
assume that it separately extends to an affine source or finite map.
Its diagonal is even by the preceding parity argument.

## 2. Keep the A4 cap coefficient2 through every solid level

By PF87 all A4 common sources, later A4 births and initial non-trefoil
A4 complementary origins are closed. The remaining origin is T(2,3).
After its true fixed-sheet meridian fill the active companion is an
ACTUAL solid torus. Its closed index3 meridian fill has H1=Z/2.
Thus in the genuine basis m=M^3,L its primitive cap is

    eta=A*m+2L, A odd.                                (5)

No particular value or recurrence for A is used. PF87 shows that
any cyclic3 self loop with3|k adds the unmarked summand(Z/r)^2.
The final R mark is zero and a two-line source with zero R has cyclic
H1, so EVERY trajectory containing such a loop is excluded. Hence
all remaining cyclic3 loops have r,k both units modulo3 and add no H1.

Their unfilled active piece is an annulus with one cone r and true
primitive fibre H=k*m+3rL. Capping by(5) adds order

    J=|3Ar-2k|>0, gcd(r,J)=gcd(r,2).                  (6)

J is nonzero because3 does not divide k. If J=1, r is odd and
the result is solid again; its meridian fill still has H1=Z/2, so
it has form(5) with another odd integer. This permits ALL solid
repetitions without computing A_new. Otherwise disk(r,J) is a
nonsolid terminal piece. The only potentially nonnegative-Euler
pair(2,2) is impossible: r=2, A and k odd make J a multiple of4.
All genuine edges after that first nonsolid outcome persist by
PF83's actual fibre determinants and the retained full recognition.

## 3. Every A4 arrival with no later simple tail is closed

The true7->3 arrival has r_final=2u, k_final=3K, K odd,
gcd(2u,3K)=1. Its drilled main piece is annulus(u,2u), omitting1.
If a nonsolid companion already exists, u>=2 gives nonconsecutive
two-cone root(u,2u), impossible for the source corner fusion; with
later tails it becomes a forbidden internal two-cone annulus.
Only u=1 can remain.

If that nonsolid arrival is FINAL, its root has cone2 and normal
R-meridian/fibre distance1, H1=Z/2 and R=0. The possible source
roots are rho=3 (n_R=7,q=2) or the rho=2 collar (n_R=3,q=1).
For the former, R=0 forces n_D=2B and the determinant2 equation is
7B-2p^2=1. It would require p^2=3 modulo7, impossible.
For the collar root, cone2 and distance1 imply p=2g,sigma=g.
By(2), delta>=2g^2+3g+3>2. Both cases fail.

If the companion is STILL SOLID at the arrival, capping its root
adds J_f=2|uA-K|. For u>=2 it is positive and even. The whole
disk(u,2u,J_f) cannot be(q,q+1,2): u=2 gives all-even orders;
u>2 leaves the nonconsecutive pair(u,2u) if J_f=2. Later tails
retain its three terminal cones. Thus again u=1.

Then J_f=2|A-K| is0 or a positive multiple of4. At0 the fibre cap
gives solid-torus connected-sum L(2,1), with group Z*C2. The torsion
factor survives all later boundary attachments. The source exterior
is aspherical with torsion-free group, so these trajectories fail
at ANY later length, without using boundary-incompressible recognition
on that reducible target. For J_f>=4 the whole exterior is disk(2,J_f).
A FINAL arrival must match the whole collar/single-D model. The
only positive-sigma assignment is p=J_f,kappa=2,sigma=J_f-1.
Its actual R-meridian/fibre distance is2J_f-3>=5, whereas the target
distance is1. Thus every FINAL-stage A4 arrival is closed.

## 4. Later simple tails force a small marked source

Every remaining case has u=1, a nonsolid negative-Euler arrival
exterior and at least one later simple-state tail. Even-k simple
tails fail their equal-two-cone root/internal budget. All passing
simple tails therefore have odd winding and odd k. PF87 preserves

    H1(end)=Z/2, [R]=0, final R-fibre distance2.

PF83 forces the rho=2 collar source, n_R=3, delta=2. Its actual
first D fibre is -b*h_x+Q*h_f with g=gcd(sigma,p), Q=p/g,
b=sigma/g and gcd(Q,b)=1. Matching the normal R meridian
h_x-3h_f gives |Q-3b|=2. The branch Q=3b-2 contradicts(2):

    delta >= g^2*Q*(3b-Q)+3g*(Q-b)+3 >2.

Hence Q=3b+2. Coprimality makes b odd; zero R makes p even,
so g is even. In particular the FINAL simple-tail winding is

    Q=3b+2=5 mod6.                                    (7)

This is an all-depth marked-root constraint, not just a cone-order test.

## 5. Still-solid arrival: only terminal(2,4)

Here the arrival exterior is disk(2,J_f), J_f=0 mod4, and the later
tails leave it terminal. Its cap coefficient at order J_f is
+/-(1+J_f/2) modulo J_f. Indeed for the arrival fibre(K,2) choose
a primitive section(c,d) with2c-Kd=1, so d is odd. The old cap
(A,2) has coefficient2c-Ad=1+(K-A)d; since K-A=+/-J_f/2,
this is1+J_f/2 modulo J_f. The other cone has order2.

Whether J_f is assigned to the source free-tail or winding leg,
its source coefficient must be +/-1: the free-tail coefficient is
+/-1, and winding coefficient is inverse(2-1)=1 if the free order2
is used instead. Thus J_f=4. The assignment(winding4,free2) fails
PF83's full last-gcd congruence modulo16: delta2 would equal
3*(4-3*2)=10 modulo16. Only(winding2,free4) remains.

In fact ONE later simple tail is impossible. If its winding is Q,
the source D path has exactly two nodes, with windings(Q,2).
Then p=2Q,sigma=2b,j_*=sigma+3=2b+3. The internal fibre
intersection is D=p-j_*=4b+1. To see this without extrapolating
from a cluster, remove preceding noncharacteristic integral terms
in the weighted local coordinates and write
u=w^p*(1-c*w^(p-j_*)+...), x=w^sigma.
At the first weighted node, a Bezout monomial has order2 and
u^b/x^Q-1 has order p-j_*. Their determinant-one coordinate change
and the Euclidean resolution transport the next fibre with transverse
coefficient p-j_*; this is exactly the fibre intersection.
The first-node order-Q coefficient is +/-inverse(b).

The actual target tail cone has coefficient +/-inverse(2k), and
its edge intersection is |Q*r_0*k_0-k| with r_0=2,k_0=3K.
Matching the normalized coefficient and edge therefore requires
2D=+/-b modulo Q. With D=4b+1 and Q=3b+2 this forces Q to divide
b-2 or4, impossible for odd b>=1 and Q>=5. Thus this case needs
AT LEAST TWO later simple tails. Extra tails are not declared closed.

## 6. Earlier nonsolid cyclic3 cap: only(2,8) or(4,10)

Let its terminal orders be(r,J) from(6). Put(a,kappa) for the
source winding/free assignment, allowing BOTH assignments and BOTH
coefficient signs. There is an internal arrival cone2 and a later
tail, so the full last-gcd formula applies; a divides p and sigma.

The cone of order r has coefficient +/-inverse(3k). The cap cone
coefficient beta_J satisfies3r*beta_J=2 modulo J: choose a section
(c,d) with3r*c-k*d=1 and compute beta_J=2c-Ad. These are integral
identities, without dividing by2 when r,J are even.

The source coefficients are +/-inverse(kappa-1) at winding a and
+/-1 at free kappa. Combining these with delta2=3*kappa mod a
gives, for EITHER assignment,

    a divides8, kappa divides3a-2 OR3a+2.               (8)

Explicitly, for a=r, multiplying3k=+/-(kappa-1) by2 and using
2k=+/-kappa mod a forces2kappa=0 or4 modulo a; with3kappa=2 this
gives a|4 or a|8. The free coefficient gives kappa|3a+/-2.
For a=J, the cap coefficient and source winding inverse give
a|(kappa+2) or a|(5kappa-2); with3kappa=2 this again gives a|8
or a|4. The free coefficient on r and(6) give kappa|3a+/-2.
All independent signs have been allowed as weaker necessary conditions.

Now a>=2, gcd(a,kappa)=2 and

    2=3*(a-(a-1)*kappa) mod a^2.

The COMPLETE finite divisor lists a=2,4,8 in(8) give exactly(1).
There is no bounded parameter or characteristic-depth search here.
For(2,4), (4) gives delta'=5. But the internal arrival cone2 is an
EARLIER D winding, so p'=p/a is even; n_D' is even by Section1.
Then delta'=3n_D'-p'^2 is even, contradiction. Therefore only
(2,8) and(4,10) remain in this earlier-reset case; both strip to
delta'=8. No claim that all determinant8 charts are impossible is made.

If there is exactly ONE later simple tail and NO cyclic3 loop between
the nonsolid cap and arrival, the source winding path is(Q,2,a).
After stripping a, the two-node chart has p'=2Q,sigma'=2b and
delta'=8. Its final free order is2, by direct substitution of
Q=3b+2 in the contact formula. The original last index would then
lie AFTER its predecessor only if kappa<=a; both8>2 and10>4 violate
that necessary index order. Such a trajectory also closes. Other
reset histories are left open.

## 7. Evidence and explicit limits

The [checker](verify_nonunimodular_a4.py),
[receipt](nonunimodular_a4_receipt.json) and
[manifest](sources/nonunimodular-a4/manifest.json) retain integral cap
and marking controls, exact gcd-contact/strip identities, complete
divisor lists, and true full target CW controls. Larger diagnostic
charts use exact gcd contact arithmetic, not symbolic series expansion.
The completed small cluster and the stopped larger exploratory
calculation are retained under research/.scratch/wright89; no unfinished control
is counted as evidence. Retained primary topology audits are imported.

Countercontrols matter: sigma8,p28,support{28,26,15} has delta2 and
windings(7,2,2), free8, but root distance1. Sigma28,p140,
support{140,48,46,35} has delta2, root distance2, windings(5,7,2,2),
free8. These are primitive local chart/lattice data, not affine
realizations or full marked-topology matches. They rule out excluding
all delta2 charts, or all delta2 charts with an internal cone2, by
those invariants alone. The new ordered/cap constraints remain essential.

Remaining A4 complementary histories: solid until arrival then
terminal(wind2,free4) and at least TWO simple tails; or an earlier
nonsolid cyclic3 cap with terminal(wind2,free8)/(wind4,free10),
at least one simple tail, and the restrictions above. All prior
PF81/83 seed exclusions remain. Incoming2(S3) remains separately.
No full quartic/higher-degree/degenerate exclusion, finite Keller map,
formal verification, external review or novelty/priority is asserted.
The full goal remains OPEN.
