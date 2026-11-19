# Exact/logarithmic separation and one-double-root odd-cubic centers

Date: 19 November 2026. The all-order arguments below are independent of the separate bounded moment controls.

## 1. A reusable exact/logarithmic separation lemma

**Lemma FunctionFieldExactLogSeparation.** Let k be an algebraically closed field of characteristic zero and L a finitely generated extension of transcendence degree one over k. For g in L, f_i in L*, and a_i in k, if

    dg + sum_i a_i df_i/f_i = 0 in Omega_(L/k),           (L1)

then g is in k and the logarithmic sum vanishes separately.

Proof. If g is nonconstant then L/k(g) is finite separable. The derivation d/dg extends uniquely to L and commutes with all k(g)-embeddings into a normal closure. Apply it to (L1) and take trace to k(g). The product-over-embeddings norm formula gives

    [L:k(g)] + sum_i a_i G_i'(g)/G_i(g) = 0,
    G_i=Norm_(L/k(g))(f_i) in k(g)*.                    (L2)

Each nonzero rational G_i has logarithmic derivative O(1/g) at infinity. The nonzero constant [L:k(g)] cannot be canceled, a contradiction. Thus g is algebraic over k and hence is in k. Equation (L1) now makes the logarithmic sum zero. This proof uses no analytic branch, valuation matrix, residue theorem, or curve compactification. The algebraic-closure hypothesis may be weakened to k relatively algebraically closed in L, but that extension is not needed here.

The norm identity follows directly by differentiating the finite product of embeddings. Uniqueness of extended derivations follows by differentiating the separable minimal polynomial. These are the same elementary inputs already audited in [ONE_POLE_ODD_CUBIC_CLASSIFICATION.md](ONE_POLE_ODD_CUBIC_CLASSIFICATION.md).

## 2. Application to every unit-center odd-cubic phase

Let p=t/x*(C+t^2 D) be a regular all-pure odd-cubic phase over a characteristic-zero field K, with C(0)=-D(0)=1. The rational-u descent in [ODD_CUBIC_HEIGHT_LOW_DEGREE_CLASSIFICATION.md](ODD_CUBIC_HEIGHT_LOW_DEGREE_CLASSIFICATION.md) gives u in K(x), u=x+O(x^2), and

    D=-C^3 u^2/x^2.

Over k=Kbar, decompose the proper rational differential (or include its polynomial part if C is constant):

    dx/C(x)=dR(x)+sum_i a_i dx/(x-r_i).                 (L3)

Here R is rational, normalized R(0)=0, and r_i run over geometric roots of C with nonzero residues. All r_i are nonzero. Let V be the normalized formal primitive of dx/C, and v(u) the inverse branch x. Purity says V(v(u))-V(v(-u))=2u.

Put x_+=v(u), x_-=v(-u), inside their common finite function field L=k(u,x_+,x_-). The rationality and nonconstancy of u make this extension finite. Differentiating the purity identity and (L3) gives

    d[R(x_+)-R(x_-)-2u]
      +sum_i a_i dlog((x_+-r_i)/(x_--r_i))=0.           (L4)

Apply Lemma FunctionFieldExactLogSeparation. The bracket is constant, and its formal value at u=0 is zero. Therefore every pure phase satisfies the TWO separate all-order conditions

    R(v(u))-R(v(-u))=2u,
    sum_i a_i dlog((v(u)-r_i)/(v(-u)-r_i))=0.           (L5)

Conversely (L5) implies the original purity criterion, since the normalized formal logarithmic primitive difference has derivative zero and value zero at u=0. Thus (L5) is an equivalence, provided rational u gives polynomial D and the stated unit-center normalization.

In particular R is nonconstant, so squarefree nonconstant C is impossible. The field E=k(u,R) is stable under the formal involution u->-u,R->R-2u. It is an intermediate rational-function field, so by Luroth E=k(W) and the involution is Mobius on W. This is a reduction, not a classification of all C. The logarithmic differential need not descend to E without a further argument; one can safely enlarge E to k(u,R,T), T=1/(C u'), which is stable since T->2-T. In the degree-one-R case below E already equals k(x), so no descent issue occurs.

## 3. A genuine wider prime-degree obstruction

**Theorem OneDoubleRootPrimeExclusion.** Suppose a regular unit-center pure odd-cubic phase has exactly one geometric double root of C, while every other geometric root of C is simple. Then its least balanced/sphere degree is composite or even; in particular it cannot be an odd prime.

Normalize C=U^2 Q(U), U=1+lambda*x, where lambda!=0, Q squarefree, Q(0)!=0, Q(1)=1, and q=deg Q. The unique repeated root descends to K (it is Galois fixed), but the degree obstruction may equivalently be proved over k=Kbar.

The rational primitive part in (L3), up to an additive constant, is

    R=-1/(lambda Q(0) U).

There is no other rational primitive part, since all other roots are simple. Thus R has degree one and E=k(u,R)=k(U). By (L5) the involution is a Mobius map sigma on U fixing U=1. Its derivative there is -1, and differentiating u=(R-Rsigma)/2 at x=0 gives 1=R'(0)=1/Q(0). Hence Q(0)=1.

All nontrivial Mobius involutions fixing 1 are

    sigma(U)=(U-alpha)/((2-alpha)U-1), alpha!=1,

with alpha=infinity interpreted as sigma(U)=1/U. For finite alpha,

    u=(U-1)*((2-alpha)U-alpha)/(2lambda U(U-alpha)),
    D=-U^4 Q(U)^3*((2-alpha)U-alpha)^2/[4(U-alpha)^2].   (L6)

If alpha is finite and nonzero, polynomial D forces Q(alpha)=0. Alpha=0 is also allowed. For alpha!=2 the numerator's linear factor has nonzero leading coefficient, so

    degD=3q+4,  leastdegree=3q+6=3(q+2).                (L7)

This is composite, including q=0.

The exceptional alpha=2 would give reflection sigma(U)=2-U. It is impossible for nonconstant Q. Indeed (L5) says that the rational logarithmic differential

    Omega=(1/Q(U)-1) dU/(lambda U^2)

is sigma-invariant. At infinity it is regular and nonzero: for q>=1 its leading term is -dU/(lambda U^2). But infinity is a fixed point of the reflection, whose derivative in coordinate z=1/U is -1. A sigma-invariant regular differential must vanish at such a fixed point (its value equals its negative), a contradiction. If q=0 then alpha=2 was already excluded by polynomiality, since Q has no root there.

It remains to treat inversion sigma(U)=1/U. Here

    u=(U-U^(-1))/(2lambda),
    D=-U^4 Q(U)^3*(U+1)^2/4,
    leastdegree=3q+8.                                  (L8)

The invariance of Omega is exactly

    1/Q(U)+U^2/Q(1/U)=1+U^2.                           (L9)

The roots of Q are simple nonzero finite poles of Omega. Infinity is regular; thus U=0 is also regular by inversion symmetry, and the root set of Q is stable under reciprocal inversion. Consequently its reciprocal polynomial Qrev(U)=U^q Q(1/U) is a scalar multiple of Q. Evaluating at U=1 and using Q(1)=1 makes that scalar 1. Therefore (L9) forces

    Q(U)=(1+U^(q+2))/(1+U^2).                          (L10)

This is a polynomial if and only if q is a multiple of 4, as seen at U=i and U=-i. Conversely for q=4s (s>=0) it is monic, squarefree, Q(0)=Q(1)=1, and satisfies (L9). The displayed u,D then satisfy both parts of (L5), so these inversion cases really are pure. Their degrees are 12s+8, all even. This completes the prime exclusion. It does not classify which finite-alpha cases are admissible beyond (L7).

## 4. Positive examples and exact scope

The finite-alpha branch is not empty. For any s>=1 put

    Q(U)=(1+(2U-1)^(2s))/2,
    C=U^2 Q(U), u=x/U, D=-U^4 Q(U)^3.

Then sigma(U)=U/(2U-1), and Phi=1/Q(U). Since Q(sigma(U))=Q(U)/(2U-1)^(2s), Phi+sigma(Phi)=2. Thus the phase is all-pure and has degree 6s+6, composite. Its single double root is U=0 and its other roots are simple. This is a test family, not a new-degree or priority claim.

More generally, whenever deg C=c and u~a*x^r at infinity, a!=0, equation D=-C^3u^2/x^2 gives

    degD=3c+2r-2,
    leastdegree=max(c,3c+2r).

In particular any u finite and nonzero at infinity forces degree 3c. The above theorem closes the only possible exceptional infinity behavior in the one-double-root class. It does not classify centers with multiple repeated roots, a triple root plus other roots, or arbitrary quintic/septic phases outside odd-cubic height.

## 5. Classical antecedents

The exact/logarithmic decomposition is classical Liouville theory. Rosenlicht, “Liouville's theorem on functions with elementary integrals,” Pacific J. Math.24 (1968),153--161, https://msp.org/pjm/1968/24-1/pjm-v24-n1-p16-p.pdf, was inspected: the proposition at printedpp154--155 gives unique algebraic extension of derivations; the theorem's algebraic-extension step at printedp159 uses trace/norm descent, and the example at printedp160 uses the impossibility of nonconstant exact derivatives being logarithmic derivatives. The short Lemma1 above is proved directly and is not a new classical theorem claim.

The retained [Milne Fields and Galois Theory v5.10](sources/milne_fields_galois_v510.pdf), Remark5.47, supplies the precise separable norm-product formula; Theorem9.19 (printedp117) supplies Luroth's theorem for the optional rational-field reduction. The one-double-root obstruction itself has R degree1, so it does not import Luroth.


The [source ledger](mixed_center_logarithmic_separation_source_ledger.json) records the exact retained descent, norm and original-contraction inputs, the optional Luroth use, and the classical Liouville antecedent. The [checker](verify_mixed_center_logarithmic_separation.py) and [receipt](mixed_center_logarithmic_separation_verification.json) pin these proof and source bytes.

```bash
python3 research/foundations/extended-geometry/verify_mixed_center_logarithmic_separation.py \
  --write-receipt research/foundations/extended-geometry/mixed_center_logarithmic_separation_verification.json
```
