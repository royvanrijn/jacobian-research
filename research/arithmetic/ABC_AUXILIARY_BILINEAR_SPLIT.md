# A split with an invariant additive defect

The direct use of a genuine row in a tree split fixes the
child ratio. An auxiliary bilinear transformation makes that
ratio uniform on its allowed finite-field range. Its effect
on the exact additive defect is multiplication by just one
parameter, so it preserves exact truth. The resulting
conditional laws also explain why a leaf statistic must
retain information about the original labels.

## An exact integer transformation and its unit laws

Fix positive \(A,B,C\) and positive integer parameters \(R,T\).
Write
\[
  A'=A(R+TC),\qquad B'=BR,\qquad C'=C(R+TA).
\]
The exact identity is
\[
                A'+B'-C'=R(A+B-C).
                                                               \tag{BS1}
\]
This has no height restriction and constrains the original
integer equation, not only its reduction. If the original
complete values are bounded by \(X^H\) and the parameters by
\(X^J\), the auxiliary values have polynomial height. The
three radical caps and any exact repair restriction are those
of the original row.

Work also at an odd prime \(q\) at which all three original
values are units. Write
\(\kappa=(A+B-C)/B\in\mathbb F_q\).
Choose \(R\) uniform on the units and \(s=T/R\) uniform
on the field subject to \(1+As\ne0\), \(1+Cs\ne0\),
independently of \(R\). In the units chart put
\[
                    x=B'/C',\qquad y=A'/C'.
\]
One has
\[
                x=\frac{B/C}{1+As},
             \qquad y+(1-\kappa)x=1.
\]
The first map bijects the allowed \(s\)'s with the nonzero
\(x\)'s such that \(1-(1-\kappa)x\ne0\).
In fact the missing \(s=-1/C\) corresponds to
\(x=1/(1-\kappa)\) if \(A\ne C\) in the field.
If \(A=C\) the two forbidden \(s\)'s coincide and no
nonzero \(x\) is deleted. For each fixed \(s\), the map
\(R\mapsto C'=C R(1+As)\) permutes the units.
Thus conditionally on any fixed input triple the output is
uniform on the unit triples in the pencil
\[
                           A'+(1-\kappa)B'=C'.
                                                               \tag{BS2}
\]
For \(\kappa=0\) this is exactly the uniform additive
unit law. Each transformation in this argument can
also be represented by positive integers \(R,T\le q\)
and satisfies (BS1) on the original integers.

In the true class, for any fixed nonzero parent argument
\(u\), take \(u_+=u A'/C'\) and \(u_-=-u B'/C'\).
Then \(u_+-u_-=u\) and the ratio \(u_+/u_-\) is
uniform on \(\mathbb F_q^\times\setminus\{1\}\).
Restricting it to a prescribed nonempty square coset
and retaining the valid ratios gives a uniform split in
that coset. The restriction is an ordinary conditioning
on the two auxiliary variables. Thus one can avoid
the sparse-ratio obstruction to the direct row split
at an individual node. This statement concerns a
single split at fixed parent data, rather than the law
on common slots in a pair of trees.
[Direct row split obstruction (CT3)](ABC_TWO_SIDED_AUXILIARY_TREE_CANCELLATION.md#direct-row-splits-do-not-realize-the-uniform-law)

## Independent radical products as parameters

The independent unit supply (MSU2) can provide the two
parameters at most auxiliary primes. Choose two independent
products of fixed numbers of draws from the one-role radical
dictionaries. Their reduced laws have relative densities
\(1+O(n^{-\gamma})\) against Haar measure on the units,
uniformly in the original triangle. Use one for \(R\).
For \(T\) use the other with probability \(1-1/q\),
and the positive integer \(q\) with probability \(1/q\).
Its reduction now has relative density
\(1+O(n^{-\gamma})\) against uniform field measure,
including zero. The joint law of \((R,T/R)\) is therefore
\(1+O(n^{-\gamma})\) relative to the independent field
law used above. Conditioning on the two transformed
unit conditions and pushing forward through (BS2) preserves
this relative error. The same holds if the child ratio
is further restricted to an allowed square coset.
This supplies a genuine arithmetic marginal law for one
split, with the original row fixed and its additive
defect preserved exactly by (BS1). Pairing several
splits with shared slots still requires a joint law.
[Independent support-product units (MSU1–MSU2)](ABC_MULTIPLICATIVE_SPECTATOR_UNIT_SUPPLY.md#uniform-mixing-at-auxiliary-primes)

## Price of retaining the original labels

The marginal law in (BS2) depends on the input only
through its normalized defect \(\kappa\). In particular,
on true rows a statistic depending only on the three
auxiliary values has identical expectation for every
input row. A transfer of the original row counts must
therefore retain original labels or additional records of
the transformation. For fixed \(R,T\) over the field
and given the three output values, \(B=B'/R\) is
determined and \(C-A=(C'-A')/R\) is determined.
Writing the latter value as \(D\), the original
\(A\) satisfies
\[
                        A'=A\{R+T(A+D)\}.
\]
This is a nonzero linear or quadratic equation and
has at most two solutions. Thus retaining the two
parameters makes the lift reversible at bounded fibre
cost. An unrestricted parameter catalogue, however,
has order \(q^2\). Elimination of this cost in a
signed comparison would require keeping information
about the original labels through the averaging,
as well as coordinating the law on the shared slots.
The affine multipliers \(R+TA,R+TC\) have no
automatic small-radical cap.
