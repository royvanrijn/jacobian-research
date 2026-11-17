# Gaussian radial transfer and exceptional quotient uniformity

15 November 2026. Internal written deduction from the retained
[Gaussian field theorem](FIBERWISE_POLYNOMIAL_HOM_AND_GAUSSIAN_FIELD_UNIFORMITY.md),
[coefficient-ideal algorithms](PRIMARY_THICKENING_MOMENT_IDEAL_ALGORITHM.md)
and [nilpotent lifting theorem](NILPOTENT_LIFTING_AND_FORMAL_SAFE_COMPLETION.md).
Independent specialist review remains outstanding, including review of the
Gaussian foundation. No literature-priority claim is made. The
[registry](../../MATH_STATUS.json) is authoritative.

This supplies the uniform radial bridge absent from the earlier
[Haar decision note](ENTIRE_HAAR_PERIODS_AND_EXCEPTIONAL_DECISION.md).
For each fixed matrix-entry degree there is a common linear cutoff for
every safe pure-null phase on either SL2/2O or SL2/2I, and their pure-safe
loci are closed algebraic sets. Neither whole exceptional Mathieu--Zhao
question is decided by this result.

## 1. Even representatives and homogeneous lifts

Let Gamma be either binary exceptional group 2O or 2I in SU2, and set

    B=K[x,y,z,w]/(Delta-1),  Delta=xw-yz,  A=B^Gamma.

Matrix-entry degree means the quotient filtration induced from K[x,y,z,w],
not degree in a chosen system of invariant generators. Work first over C;
algebraic embeddings and other characteristic-zero fields are handled by
a fixed field extension and descent. The finite groups act on the right.

Both groups contain -I. Reynolds averaging a representative of any P in
A of filtered degree at most d produces an invariant polynomial of degree
at most d. Its homogeneous pieces are invariant separately, and its odd
pieces vanish. Put D=2 floor(d/2), and write P=sum_(j even) P_j on Delta=1.
The polynomial

    Ptilde=sum_(j even) Delta^((D-j)/2) P_j                 (1)

is invariant and homogeneous of degree D, and restricts to P. The same
construction gives Qtilde of any even degree E at least the filtered
degree of an invariant marker Q. It is linear in all phase and marker
coefficients. For D=0 it just gives constants.

The lift at prescribed degree is independent of the representative.
A homogeneous polynomial vanishing on SU2 vanishes on its real radial
cone, after the invertible complex linear substitution below, and hence
is the zero polynomial. Equivalently one can use Zariski density and
homogeneity on Delta=1. Thus (1) identifies the filtered invariant function
space with the homogeneous invariant polynomial space of degree D.
This argument is a linear identity and extends to coefficient algebras.

## 2. The radial identity in every order

Use the quaternion coordinates

    [[x,y],[z,w]]=[[a+ib,c+id],[-c+id,a-ib]],
    Delta=a²+b²+c²+d²=r².

For four independent standard real Gaussians (variance one), write
(a,b,c,d)=r u with u uniform on S3. Its quaternion matrix is Haar SU2,
and r is independent of u. The homogeneous product has degree L=Dm+E.
Since r² has chi-square distribution with four degrees of freedom,

    E[r^L]=2^(L/2) Gamma(L/2+2)/Gamma(2)
           =2^(L/2)(L/2+1)!.

Consequently, for every integer m>=0,

    E[Qtilde Ptilde^m]
      =2^((Dm+E)/2)((Dm+E)/2+1)! h(Q P^m).               (2)

Here h is normalized Haar, equivalently quotient Haar on invariant
functions. D and E are even, so the factor is a nonzero integer and a
unit over every Q-algebra. This proves the identity in all orders;
a finite list of moment checks is only a transcription control.
Other Gaussian variance conventions change the displayed factor.

## 3. All Gaussian markers and exact cutoff ideals

The right Gamma action is orthogonal on the real quaternion space and
preserves the Gaussian. Since Ptilde is invariant, every polynomial marker
R satisfies

    E[R Ptilde^m]=E[(Reynolds_Gamma R) Ptilde^m].           (3)

Averaging preserves homogeneous degree and never increases filtered
degree. Its odd homogeneous components vanish because -I belongs to Gamma.
For an even homogeneous component of degree e, apply (2) to its restriction
on Delta=1; this restriction has filtered degree at most e. Conversely
each invariant marker of filtered degree e has the lift (1) of degree e.
Decompose arbitrary markers into homogeneous components to conclude:

* P is pure-null exactly when Ptilde is Gaussian pure-null.
* P is safe for every fixed invariant marker exactly when Ptilde is safe
  for every fixed Gaussian polynomial marker.
* For every real N>=0 the condition h(QP^m)=0 for all m>N deg(Q), m>=1,
  is equivalent to the identical Gaussian degree cutoff for Ptilde.

For qualitative safety, a finite sum of homogeneous components has a
finite maximum of their individual cutoffs. This supplies the necessary
all-marker quantifier, not just a transfer for preselected homogeneous Q.

Fix any finite invariant coefficient family over a characteristic-zero
effective field, with polynomial parameter ring R and the linear lift (1).
Let I_pure be its full ordinary pure-moment ideal and I_N the ideal of
all forbidden marked moments m>N deg(Q), including constant markers.
Then I_pure and every I_N are exactly the pullbacks of their Gaussian
ideals for the lifted family. This is equality of ordinary ideals, not
just equality of zero sets: (2)--(3) express each generator as a linear
combination of generators in the other model, multiplied by field units.
When a homogeneous ambient marker restricts to smaller filtered degree,
the inequality m>N e still implies the required quotient inequality.
In the reverse direction use the lift of the marker's filtered degree.
The fixed quaternion scalar extension is faithfully flat, so identities
and ideal equalities descend. The same conclusions hold after any base
change to a coefficient Q-algebra with the group defined over its fixed
field of definition.

## 4. Uniform field and reduced-ring consequences

Let C_G(4,D) be a common integer supplied by the retained Gaussian field
theorem. Every pure-null safe phase P of matrix-entry degree at most d
on either exceptional quotient satisfies

    h(QP^m)=0 whenever m>C_G(4,D) deg_matrix(Q).            (4)

The same constant works for both quotients and all characteristic-zero
fields over which their models are defined. The case D=0 contains only
the zero pure-null field phase and needs slope zero. The pure-null safe
locus is the inverse image of the closed Gaussian pure-safe locus under
the linear lift. Its radical ideal is theoretically computable; the
least common integer quotient slope is theoretically computable by
enumerating I_N and testing radical equality. General producers and
universal numerical values have not been computed here.

The same bound works over reduced coefficient algebras when every prime
residue-field phase is pure-null and safe. Radical cutoff ideals eventually
stabilize. This strengthens the P-dependent constants and bounded-family
decisions in the Haar note; it does not replace their explicit Haar source
or assert that all pure-null phases are safe.

## 5. Ordinary coefficients and nilpotent lifting

The exact ideal transfer gives access to the Gaussian ordinary pure-ideal
and prescribed-cutoff algorithms, retaining embedded and nilpotent
directions. They remain theoretical general algorithms. The retained
Gaussian nilpotent theorem, applied to the lifted family, also gives:

* For a pure-null phase over any coefficient Q-algebra, safety is equivalent
  to safety after removing nilpotents, or over every prime residue field.
* A fixed nilpotent thickness gives a common linear cutoff. For example,
  if Nil(A)^L=0, one sufficient slope is
  C+(CD+1)(qL-1), with C>=max(1,C_G(4,D)) and q=binomial(4+D,D).
* On a fixed ordinary pure coefficient scheme, its mixed-cutoff ideals
  are cofinal with powers of its reduced safe ideal. Thus their limit
  is the formal completion along the reduced safe subscheme.

For the last assertion one may either repeat coefficient differentiation
with the linear degree-D lift in the original regular coefficient space,
or pull back the Gaussian cofinal containments. In the latter argument
the pulled-back Gaussian safe ideal need not itself be radical: its
radical is the original safe ideal, and finite generation makes the two
power filtrations cofinal. No assertion commutes radicals with arbitrary
base change. Pure equations must be imposed before claiming cofinality.

There is no claim of a thickness-independent cutoff over all nonreduced
algebras, ordinary ideal stabilization, a sharp quotient constant, or
transfer of every ambient Gaussian failure into the invariant subspace.
In particular the retained quartic Gaussian failure is not an exceptional
counterexample merely by radial integration. Both all-degree exceptional
Mathieu--Zhao decisions remain open.

## 6. Evidence and remaining work

The user reported eight exact identities as transcription controls; their
original inputs were not supplied in this request. The separate retained
[checker](verify_gaussian_radial_transfer.py) executes eight reproducible
identities, including determinant homogenization, a zero moment and complex
quaternion entries. Those controls test the universal even sector; they do
not enumerate either exceptional group or prove the all-order statement.
The written homogeneity and Reynolds arguments above supply that statement.
The [receipt](gaussian_radial_transfer_verification.json) records this scope.

Next steps are independent review of this bridge and its Gaussian foundation,
structural classification of the invariant safe locus, actual failure
construction or an exhaustive positive theorem, and implementation of the
existing coefficient/holonomic producers when explicitly requested.
Uniform radial transfer itself is no longer an open obligation in this
internal programme. Matrix-entry bounds remain fixed throughout.
