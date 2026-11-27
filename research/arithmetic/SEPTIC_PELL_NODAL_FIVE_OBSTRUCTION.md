# Ordinary nodes at five do not support five septic Pell directions

**Theorem SN1.** Let \(P\in\mathbb Z_5[t]\) have degree seven
and unit leading coefficient, and suppose \(P^2+1\) is squarefree
over \(\mathbb Q_5\). If every geometric root of the reduction
\(\bar P^2+1\) has multiplicity at most two, then
\[
 \dim_{\mathbb F_3}
 \operatorname{Jac}(R^2=P^2+1)(\mathbb Q_5)[3]\le4 .
\]
Thus a degree-seven rational square-plus-one model with five
rational \(3\)-torsion directions cannot have such a reduction.
For an integral degree-preserving model at five, at least one
root of \(\bar P^2+1\) must have multiplicity at least three.
The smooth case is the obstruction of
[SP2](SEPTIC_PELL_FIVE_DIRECTION_CONSTRUCTION.md).

The last possible nodal residue class fails because one
three-torsion class on its normalization does not lift to the
generalized Jacobian. Counting the normalization, torus and
component groups separately would leave this class open.

## The local rank budget

Write \(\bar F=\bar P^2+1=g^2 F_0\) with \(g\) monic and
\(g,F_0\) squarefree and coprime. Since the derivative of
\(\bar P\) has degree six, there are at most six double roots
over the algebraic closure. In particular \(F_0\) has positive
degree. The integral double cover of the \(t\)-line, including
its infinity chart, has geometrically irreducible nodal special
fiber with normalization \(C_0:y^2=F_0\). There are two distinct
smooth points at infinity, because the leading coefficient
of \(\bar F\) is a nonzero square. The special fiber is reduced.
Each singularity is an ordinary node, and its local smoothing
has thickness \(n\ge1\) since the generic fiber is smooth.

For an irreducible factor \(h\) of \(g\) with degree \(d\), put
\(e_h=1\) if \(F_0\bmod h\) is a square in \(\mathbb F_{5^d}\)
and \(e_h=-1\) otherwise. This is the sign on the two branches
after \(d\) Frobenius steps. Let \(T\) be the torus in the
identity component of the special Néron fiber and let \(\Phi\)
be its component group. The contribution of this node orbit to
\(\#T(\mathbb F_5)\) is \(5^d-e_h\), and to
\(\dim T(\mathbb F_5)[3]\) it is one exactly when
\(5^d\equiv e_h\bmod3\).
For \(\Phi\), the geometric graph has one loop of length \(n\)
at each geometric node. The pairing on the loops is diagonal
with these lengths. Therefore the component group is a direct
sum of the \(\mathbb Z/n\) of the geometric nodes, with the
same signed permutation. An orbit contributes at most one to
\(\dim\Phi(\mathbb F_5)[3]\), and only if \(e_h=1\) and
its thickness is divisible by three. The identity component
sequence and the Frobenius-equivariant component group formula are
[DDMM, equations (2.3)–(2.7), reduction (2.18), and Lemmas 2.21–2.22](sources/dokchitser-dokchitser-maistret-morgan-local-hyperelliptic.pdf).
In the present model the graph pairing follows by resolving each
ordinary double point, which subdivides its loop into \(n\) edges.

The formal reduction kernel has no \(3\)-torsion because
the residue characteristic is five. Also the connected
group is an extension of \(\operatorname{Jac}(C_0)\) by
\(T\). We therefore have the initial upper bound
\[
 \dim J(\mathbb Q_5)[3]\le
 b_0 + b_T + b_\Phi,\qquad
 \begin{cases}
 b_T=\#\{h:5^{\deg h}\equiv e_h\bmod3\},\\
 b_\Phi=\#\{h:e_h=1\},
 \end{cases}
 \tag{1}
\]
where \(b_0\) is the minimum of the \(3\)-valuation of
\(\#\operatorname{Jac}(C_0)(\mathbb F_5)\) and the
multiplicity of \(X-1\) in its Frobenius polynomial modulo
three. This bound permits every split thickness to be divisible
by three and permits every rational normalization class of
order three to lift to the torus extension.

## Finite elimination and the remaining gluing condition

An affine substitution makes \(\bar P\) monic and kills the
degree-six coefficient over \(\mathbb F_5\), exactly as in SP2.
Enumerate \(t^7+\sum_{j=0}^5 c_jt^j\) up to
\(\bar P(t)\mapsto-\bar P(-t)\). The [checker](verify_septic_pell_nodal_five.py)
covers all \(7875\) representatives. There are \(616\)
representatives with a root of multiplicity greater than two.
Of the remaining \(7259\), exactly six have the bound (1)
at least five, and in all six the bound is five.
As a control, the computed order of the connected reduction is
divisible by seven throughout the nodal and smooth census.
Indeed the Pell function on the integral special fiber has
divisor seven times the infinity difference. This difference
cannot be principal: a function with one simple pole at a
smooth infinity would give a degree-one finite morphism from
that integral curve to \(\mathbb P^1\), forcing an isomorphism.

The six cases are in the table. Split orbits are indicated by
their degrees. In the last column for five cases, the point
\(Q\) gives a \(3\)-primary class of order at least nine in the
Jacobian of \(C_0\): take \([Q-\infty_+]\) and multiply by
the prime-to-three part of the indicated Jacobian order.
It follows that the actual dimension of its \(3\)-torsion is at
most \(v_3(\#\operatorname{Jac}(C_0)(\mathbb F_5))-1\).
In each of these five cases this improves (1) to at most four.

| Ascending coefficients \((c_0,\ldots,c_5)\) | Split orbit degrees | \(\#\operatorname{Jac}(C_0)(\mathbb F_5)\) | \(b_0\) | \(Q\) and order after projection |
| --- | --- | ---: | ---: | --- |
| \((0,0,1,4,3,2)\) | \(1\) | 4536 | 4 | \((0,1)\), 81 |
| \((0,2,0,0,0,1)\) | \(2,2\) | 21 | 1 | Gluing calculation below |
| \((0,4,2,2,0,2)\) | \(2\) | 756 | 3 | \((0,1)\), 9 |
| \((1,0,0,4,2,0)\) | \(1\) | 4536 | 4 | \((2,2)\), 27 |
| \((2,2,1,1,4,1)\) | \(1,1\) | 1512 | 3 | \((1,2)\), 27 |
| \((2,2,1,2,1,4)\) | \(1,2\) | 252 | 2 | \((2,1)\), 9 |

In the remaining row,
\[
 \begin{split}
 \bar P&=t^7+t^5+2t,\\
 g&=(t^2+2t+4)(t^2+3t+4),\\
 F_0&=(t^3+t^2+2)(t^3+4t^2+3).
 \end{split}
\]
The two irreducible quadratic factors give split nodes over
\(\mathbb F_{25}\). Put \(D_0=(0,1)-\infty_+\) on \(C_0\).
The class of \(7D_0\) generates its rational \(3\)-torsion.
There is a function \(f=A+B y\) of divisor \(21D_0\).
Its certificate is
\[
 \begin{split}
 A&=4t^{21}+3t^{20}+4t^{19}+3t^{18}+t^{16}+4t^{15}
     +t^{14}+4t^{13}+2t^{12}+t^{11}+3t^{10}+t^8\\
   &\qquad+2t^7+2t^6+4t^4+t^3+2t+1,\\
 B&=4t^{18}+3t^{17}+t^{16}+2t^{15}+4t^{14}+4t^{13}
     +t^{12}+3t^{11}+2t^{10}+t^8+4t^7\\
   &\qquad+2t^4+3t^2+3t+4,\\
 A^2-F_0 B^2&=2t^{21}.
 \end{split}
\]
The nonzero branch value at \((0,-1)\) and the pole of
order 21 at \(\infty_+\) specify the divisor from the norm
identity.

The divisor \(7D_0\) is disjoint from the nodes, so it lifts
to a divisor class \(E\) on the nodal curve. In the torus
coordinates of \(3E\), the first coordinate is the ratio
of the values of \(f\) at the two branches over
\(h=t^2+2t+4\), up to inversion. Indeed the torus in the
normalization sequence records precisely the gluing scalars
between the two fibers of a line bundle. Reducing \(f\) on
these branches yields the ratio \(a+2\) in
\(\mathbb F_{25}=\mathbb F_5[a]/(a^2+2a+4)\).
The exact identity \((a+2)^8=3a\ne1\) shows this ratio is
not a cube. Changing the lift \(E\) by a torus class changes
\(3E\) by a cube there. Thus the normalization's order-three
line has no order-three lift in the connected reduction.
The connected \(3\)-torsion has dimension two from the two
split quadratic torus orbits. The component group adds at
most two. This proves the bound in the last row and finishes
the census.

## Scope and replay

Run the [checker](verify_septic_pell_nodal_five.py) with SageMath's
Python. It checks the complete residue class census using exact
finite-field polynomial and Frobenius arithmetic, the five
explicit Jacobian order witnesses, and the last node's gluing
certificate. The passage from the residue group to the local
Jacobian uses the written semistable reduction argument.
The result concerns characteristic-zero models with a unit leading
coefficient at five and only smooth or ordinary nodal reduction
in the given equation. Multiple roots of order at least three,
failure of degree preservation, other Pell polynomial degrees,
and the five-direction construction remain open.
