# Exact first two exceptional quadratic face-forcing pairs

Put H_j(r,z)=Lambda_1((r x²+x+z)^j), T_j=H_j+jH_(j-1).
For q in{1,2}, let m=q(4q+3)², D=-3q(q+1),
L_q=Q(omega), omega²=D, and z=-8q(q+1)(2q+1)+omega.

**Theorem 1.** The polynomials T_(m-1)(r,z) and T_m(r,z)
are coprime in L_q[r] for each of q=1 and q=2.

The exact proof data are the retained [native pair certificates](certificates/q1_q2_native_pair_certificates.json),
SHA256 6b08d338b95892bb374aaf22896c8e1adefef9fe1d8555dbe7dbdb0e1547a167. The two cases use primes p109 and p487,
respectively, and m49 and m242. In each case p>2m and D is a
quadratic nonresidue modulo p. The field reduction of Z[omega]
is F_p[s]/(s²-D), an actual finite field.

## 1. Native finite-field certificate and number-field lifting

The JSON gives every coefficient, in ascending r order, of two
polynomials U,V over that finite field. Literal native evaluation is
[r^a]H_j = (j!/a!) sum_(d=0)^(j-a)
  (j+a-d)! z^d / ((j-a-d)!d!).
The checker recomputes these native coefficients, then checks exactly
U*T_(m-1)+V*T_m=1.
It checks primality by trial division, the nonresidue condition by
modular exponentiation, both full degrees, and the leading factorial
coefficients (2m-2)! and (2m)!, which are units since p>2m.
The certificate vectors are fully retained; this is not a gcd-report-only
proof input or an experiment on rounded numerical roots.

The identity implies coprimeness of the two reduced polynomials.
Equivalently their fixed-degree Sylvester resultant is nonzero in the
residue field. If the characteristic-zero resultant in Z[omega]
were zero, its reduction would be zero as well. This contradiction
proves Theorem1 over L_q. The conjugate choice of omega is covered
by the number-field automorphism.

## 2. Stronger symmetric moment prefix, using the actual compatible pair

The additional inputs in this section are exactly the necessary exceptional
family theorem in [the residual-frontier proof](QUADRATIC_FACE_FORCING_THREE_ZERO_RESIDUAL_FRONTIER.md),
and the polynomial interior recurrence in [the full symmetric quadratic proof](FULL_SYMMETRIC_BINARY_QUADRATIC_STRONG_FACTORIAL.md),
Section2. No shorter-window conjecture or old loose prefix is an input.

**Corollary 2 (prefix through2641).** Over every characteristic-zero field,
for a,e,c with aec!=0, arbitrary b, and every1<=n<=2641, the phase
f=a(x^2+y^2)+bxy+e(x+y)+c has a nonzero factorial moment among
n,...,n+4. This includes the actual support-five b=0 sector in this range.

Indeed its exact recurrence is

 M_(j+1)=(2a+b)(j+1)^2 M_j
   +j(j+1)(e^2-(2a+b)c)M_(j-1)+H_(j+1)+(j+1)eH_j,

where H_j=Lambda_1((at^2+et+c)^j). Five zero M moments starting at n
give three zero normalized forcing values at n+2,n+3,n+4 after dividing
by powers of e. Their middle is m=n+3 and the exact exceptional-family
theorem puts m=q(4q+3)^2 with q>=1,q1or2mod4 and
z=c/e=-8q(q+1)(2q+1)±sqrt(-3q(q+1)). In particular the actual adjacent
pair at m-1,m has a common finite r=a/e. Theorem1 excludes q1 andq2,
with either sign covered by conjugation. The next eligible q is5,
so m>=2645 and n>=2642. This contradiction proves the corollary.

Thus the proof uses actual interior-to-pair compatibility, not only an
inequality for previously necessary starting indices. No assertion at
all larger indices follows from this finite-prefix theorem alone.

## 3. A genuine tail reduction of the remaining universal criterion

For every eligible q let L_q and z_q be as above. Write R for the existing
assertion that T_(m-1)(r,z_q),T_m(r,z_q) are coprime in L_q[r] for ALL
q>=1,q1or2mod4. Let R5 be the same assertion restricted to q>=5.
Both are universal assertions; no proof of R5 is given here.

**Proposition 3 (full criterion iff tail criterion).** The assertions
R and R5 are equivalent, using the two proved finite pairs of Theorem1.

R implies R5 by restriction. Conversely, the only eligible positive q
less than5 are1 and2, handled by Theorem1. R5 handles every other eligible
q, proving R. Both directions use the same adjacent-pair coprimeness
assertion; no infinite claim is inferred from a finite sample.

The existing equivalence between R and the universal quadratic forcing
three-window, with its separately proved context, can therefore be
extended to include R5. This proposition is an implication context and
seeds none of those universal members. A hypothesis R5 can conditionally
recover R, the universal forcing assertion, and then the existing
symmetric support-five target through its already recorded one-way
implication. There is no converse from that target and no unconditional
closure of it here.

## 4. Portable exact replay and scope

The [read-only standard-library checker](certificates/verify_q1_q2_native_pair_certificates.py)
reconstructs every native coefficient and multiplies every supplied
cofactor. From the repository root run

    python research/factorial/image-rigidity/certificates/verify_q1_q2_native_pair_certificates.py research/factorial/image-rigidity/certificates/q1_q2_native_pair_certificates.json

The certificate bytes are pinned above. The checker does not write files
or need a network, a computer-algebra package, or any moment table.
The finite certificate is a complete proof input for exactly q1 andq2;
the extra prefix and tail implications have the written proofs above.
No universal tail, full symmetric support-five, arbitrary binary quadratic,
SFC(2), or SFC(3) theorem is concluded. All arithmetic is exact.
