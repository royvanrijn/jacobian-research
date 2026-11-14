# A generic quartic discriminant gap closes the elliptic pole(6,9) profile

Current handoff: [PF2RULED21](UNIBRANCH_POLAR33_OBSTRUCTION.md) closes both
higher-polar33 profiles retained below as PF2RULED20 necessary data. The local
quartic gap and exact two-cusp theorem here remain valid. Current noncontact
partitions are72,54; every retained pole69 common nonic has support<=2.

14 November 2026. PF2RULED20 proves a complete local obstruction for
pencils A-c*z^6 with A(0,w)=w^4*unit. Generic polar order17 is impossible.
If local finiteness forces it to be at least17, it is exactly18, and the
generic fibre has two normalized cusp branches with simple Q poles.

This excludes the elliptic(2,2,2)/contact3 family of pole(6,9), as well as
the noncontact partitions(5,2,2) and(3,2,2,2). The surviving noncontact
partitions are now

 (7,2), (5,4), (4,3,2).

The sixteen marked two-root leading candidates of PF2RULED19 remain
unresolved. Two higher polar profiles, both requiring a degree-eight
boundary branch with polar intersection33, are the next local gates.
No whole height-nine pair or actual compatible completion is proved.

## A formal normal form independent of the pencil value

Let A belong to C[[z,w]], with A(0,w)=w^4*unit. A formal coordinate W,
with z unchanged and invertible derivative in w, puts it in the form

 A=W^4+p(z)*W^2+q(z)*W+r(z), p,q,r in z*C[[z]].          (1)

The coordinate change is independent of c, so the pencil remains

 F_c=W^4+p*W^2+q*W+r-c*z^6.                             (2)

Here is a direct all-order construction. First absorb the fourth root
of the boundary unit into W at z=0. Inductively suppose the normal form
is obtained modulo z^n. Write the next residual coefficient as

 R(W)=R0+R1*W+R2*W^2+4*W^3*H(W).

A change W_old=W_new-z^n*H(W_new) removes its last part modulo z^(n+1),
while the first three terms adjust r,q,p. The existing nonleading
coefficients are divisible by z, so their coordinate-change errors
begin in order n+1. Quadratic errors begin in order2n>=n+1. The
z-adically convergent induction constructs an invertible formal change
with no division by z or c. This includes every higher w,z coefficient.

The distinction from multiplying a Weierstrass polynomial by a unit
matters: no unit depending on w is divided into c*z^6. Keeping c only
in the constant term makes the discriminant coefficient argument valid.
A formal transverse change preserves the boundary point, z and the
unit character of the volume away from E.

## The complete discriminant calculation

Work with c as the generic pencil parameter over C((z)), and put
Delta(c)=Disc_W(F_c). The polar intersection is ord_z Delta. Indeed
the completed curve algebra is free of rank4 over C[[z]], and the norm
of multiplication by F_W is the discriminant. Its determinant order
is the sum of the normalized branch orders of F_W, equivalently the
intersection I(F_c,F_W). The invertible transverse coordinate change
preserves the original polar intersection.

The four coefficients of Delta as a polynomial in c are

 [c^0] =256r^3-128p^2r^2+144pq^2r-27q^4+16p^4r-4p^3q^2,
 [c^1] =-16z^6*(48r^2-16p^2r+9pq^2+p^4),
 [c^2] =128z^12*(6r-p^2),
 [c^3] =-256z^18.                                      (3)

In particular the generic order is at most18. Assume it is at least17.
Then every coefficient in(3) has order at least17. Write ord(0)=infinity.
The c^2 coefficient gives

 r=p^2/6+H, ord(H)>=5.

The bracket in the c^1 coefficient becomes exactly

 -p^4/3+9pq^2+48H^2,                                   (4)

and must have order at least11.

If ord(p)=1, either ord(q)=1 and pq^2 uniquely leads(4) in order3,
or ord(q)>=2 and p^4 uniquely leads it in order4. Both are impossible.

If ord(p)=2, the same comparison forces ord(q)=3. Write the leading
coefficients as p=P2*z^2+..., q=Q3*z^3+..., with P2!=0. Equation(4)
requires Q3^2=P2^3/27. The order12 coefficient of [c^0]Delta is

 P2^6+(Q3^2-P2^3/27)*(19P2^3-27Q3^2)=P2^6!=0.           (5)

This exact nonzero certificate contradicts order>=17.

Thus ord(p)>=3, and the c^2 condition gives ord(r)>=5. If ord(q)<=3,
the term -27q^4 is uniquely leading in [c^0]Delta below order17. Hence
ord(q)>=4. If ord(r)=5, the term256r^3 uniquely leads in order15.
Thus ord(r)>=6. If ord(q)=4, -27q^4 uniquely leads in order16.
Consequently

 ord(p)>=3, ord(q)>=5, ord(r)>=6.                       (6)

Conversely every term in(3) has order>=18 under(6). The fixed nonzero
c^3 coefficient prevents an identically vanishing leading polynomial.
Therefore

 generic ord Delta>=17
 iff(6), and then generic ord Delta=18.                (7)

The comparisons cover every nonnegative higher-order coefficient and
all zero-series cases. This is an all-series valuation proof, not an
extrapolation from a truncated coefficient search.

## Two normalized branches and two simple poles

Suppose P=A/z^6 has boundary-symbol root multiplicity4 at a point away
from E, and an actual finite map supplies the companion Q with volume
unit*z^-3 dz wedge dw. PF2RULED19 gives

 local Q pole sum = I(F_c,A_w)-16>0.

The polar upper bound is18. Thus I>=17, and(7) forces I=18 and(6).
Let P3=[z^3]p and R6=[z^6]r. With z=t^2,W=t^3 X, the initial equation is

 X^4+P3*X^2+R6-c=0.                                    (8)

For generic c its four roots are distinct and nonzero: the discriminant
is16(R6-c)*(P3^2-4(R6-c))^2. The excluded parameter values are only
R6 and R6-P3^2/4. Each simple root extends uniquely as a formal power
series X(t); the involution t->-t pairs opposite initial roots. Since
W has odd order3 while z has order2, each pair gives one normalized
branch, and the two pairs give all the branches.

On each, ord(z)=2, ord(W)=3 and ord(F_W)=9: the leading derivative
coefficient2*alpha*(2*alpha^2+P3) is nonzero by the same genericity
conditions. The exact local pole formula then gives

 nu=9-(6-2)*2=1.

Thus every away common-nonic root of multiplicity2 in a finite
pole(6,9) map supplies exactly two missing cusp branches, each with a
simple Q pole. The normal form is necessary; it does not assert that
every such pencil admits a global regular Q or a finite target map.

The generic hypothesis cannot be dropped. For A=(z-w^2)^2-w^15,
the selected c=0 fibre has parameterization(w,z)=(t^2,t^4+t^15) and
polar order17. For generic c!=0, write w=t,z=t^2+t^6*X. Dividing by
t^12 gives X^2-c*(1+t^4*X)^6-t^3, with two simple initial roots;
each branch has polar order7, so the generic order is14. This pencil
fails the finite-pole condition. Thus an actual special17 fibre is
retained as a countercontrol to any stronger all-fibres claim.
For A=W^4+z^6+z^7 the generic order is18 while c=1 has order21.
The balanced W^4-c*z^6 control attains the generic two-cusp pattern.
These are local controls, not constructed global Keller maps.

## Excluding the old elliptic profile and two noncontact sectors

The [PF2RULED19 necessary elliptic/contact3 profile](POLAR_CUSP_LEADING_SECTORS.md)
has three away common-nonic roots of multiplicity2. It formerly required
polar17 at each. Equation(7) already excludes that requirement. Equally,
the two-cusp result gives at least six Q poles away from q and one more
at q, exceeding the entire finite quartic degree4. Thus the whole
elliptic(2,2,2)/contact3 family is excluded, including every leading
position and lower-pole lift. Its cover differential remains a valid
positive meromorphic period control, not a finite realization here.

For noncontact root partition(5,2,2), the two multiplicity-two roots
cost four Q poles and the other root at least one, again impossible.
For(3,2,2,2), the minimum cost is seven. Starting with the five necessary
partitions of PF2RULED19 leaves exactly(7,2),(5,4),(4,3,2).

At contact3 the finite-root types are now only a single root6 (the
cube branch), or two roots(4,2). No numerical source group is closed
in full:96/144 at contact0,80/128 at contacts3,4,5,6,7,9, and48/128
at contacts7,9 remain, with the refined root-type conditions.

## The surviving higher-polar profiles

In the noncontact(4,3,2) case the three roots cost at least1,1,2 missing
points. They therefore give exactly four points, all with simple Q poles.
At common multiplicities4,3,2 the P-symbol multiplicities are8,6,4.
The exact polar identity gives necessary local intersections

 33,25,18 respectively.

The last location has the two cusp branches just proved. The first two
must each have a single normalized branch. No realization or exclusion
of the polar33 degree-eight point is proved here.

For contact3 with away root types(4,2), the two away points cost1,2,
and q costs at least one. Again all four Q poles are simple. The away
polar intersections must be33 and18; at q the contact correction is4,
so the required polar intersection is29. This is the second surviving
profile requiring a single degree-eight boundary branch with polar33.
It remains open. The sixteen marked two-root leading candidates of
PF2RULED19 are also unchanged as leading tests, not global maps.

## Evidence and limits

The [checker](verify_quartic_polar_gap.py),
[certificates](quartic_polar_gap_certificates.json) and
[receipt](quartic_polar_gap_receipt.json) retain the universal normal-form
induction step, the entire discriminant, every valuation alternative,
the exact certificate(5), the forced-cone bounds, generic branch data,
special-fibre controls and every resulting pole cost. Formal convergence,
all-order valuation coverage and the generic fibre interpretation are
written proofs. No finite jet truncation replaces their scope.

No formal verification, independent full reproof, conditional2026 upgrade,
external review or novelty priority is claimed. The three remaining
noncontact sectors, contact single-root/cube cases, both higher-polar
profiles, full global Jacobian/cusp/finiteness, higher heights and actual
compatible embedding remain unresolved. The first reduced pairs are
still(5,9),(6,9),(7,9),(8,9); neither full quartic packet nor the original
full research goal is closed.
