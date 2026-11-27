# Closure of the primitive quartic pencil in top rank one

**Theorem QuinticRankOnePrimitiveQuarticPencilExclusion.** A
noninvertible degree-five three-variable Keller map of top rank
one cannot have primitive cancelling-pencil degree four in
[QuinticRankOneQuarticPencilResidual](QUINTIC_KELLER_REDUCED_BOUNDARY_PENCIL_RESIDUAL.md#rank-one-pencil).
Here the pencil degree is the ordinary homogeneous degree of
its coprime projective arguments.

It suffices by that result to consider
\[
              P_5=z^5,\qquad Q_4=z^4,\qquad R_4=C,
\]
where \(r=C(X,Y,1)\) has degree four and generates a relatively
algebraically closed rational subfield of \(\mathbb C(X,Y)\).
The proof uses a formal slice of \(Q\) to constrain its lower
terms. The slice is a determinant calculation on polynomial
coordinates. It does not assert the existence of a polynomial
coordinate root of \(Q\).

## Formal slice and its top binary degree

Set \(x=\tau X,\ y=\tau Y,\ z=\tau Z\) and write
\(\epsilon=\tau^{-1}\). The equation \(Q=\tau^4\) has a unique
formal solution
\[
 Z\in 1+\epsilon\mathbb C[X,Y][[\epsilon]]
\]
by the formal implicit function theorem, since its equation at
\(\epsilon=0\) is \(Z^4=1\) with derivative four at \(Z=1\).
Let \(\Pi\) be the series \(\tau^{-5}P\) on this slice and
\(\mathcal R\) the series \(\tau^{-4}R\), with
\(\mathcal R|_{\epsilon=0}=r\). The determinant identity on the
slice gives
\[
 \{P,R\}_{X,Y}=\pm\kappa\tau^2/Q_z,\qquad
 \{\Pi,\mathcal R\}_{X,Y}=O(\epsilon^{10}),
\tag{1}
\]
since \(Q_z=4\tau^3+O(\tau^2)\). For the first equality,
differentiate \(Q=\tau^4\) in \(X,Y\) at fixed \(\tau\) and
eliminate the derivatives of \(z\) in those of \(P,R\).

Write \(p_j=P_j(X,Y,1)\), \(q_j=Q_j(X,Y,1)\).
The coefficient of \(\epsilon\) in \(\Pi\) is
\(p_4-5q_3/4\). It commutes with \(r\) by (1). A binary
polynomial commuting with \(r\) belongs to \(\mathbb C[r]\),
by relative algebraic closure and a univariate Bézout
identity. Hence the coefficient is affine in \(r\), and we can
subtract a multiple of \(R\) from \(P\) to make it constant,
say \(c_0\). Now \(P_4=(5/4)zQ_3+c_0 z^4\).
The coefficient of \(\epsilon^2\) in \(\Pi\) is
\[
              p_3-\tfrac54 q_2-c_0q_3-\tfrac5{32}q_3^2 .
\tag{2}
\]
It too commutes with \(r\). For example the first two coefficients
of \(Z\) satisfy \(Z_1=-q_3/4\) and
\(4Z_2=-6Z_1^2-(Q_3)_z(X,Y,1)Z_1-q_2\). Substituting these
and the expression for \(P_4\) gives (2).
If \(q_3\) has degree three, (2) has degree six and cannot
belong to \(\mathbb C[r]\). Thus \(\deg q_3\le2\).

For use below, suppose we have a bound \(\deg Z_i\le h i\)
for each \(i\), where \(Z_i=[\epsilon^i]Z\).
Then
\[
 \deg[\epsilon^i]\mathcal R\le 4+hi-1 \qquad(i\ge1).
\tag{3}
\]
In fact a positive-order term from \(R_4(X,Y,Z)\) comes from
at least one \(Z\)-derivative, lowering its binary degree by
at least one. A term from a lower homogeneous component has
binary degree at most three before substitution. This gives
(3) using the bound on \(Z_i\).

In each of the next two cases \(\deg[\epsilon^i]\Pi\le h i\).
Use (1) successively to remove each nonconstant coefficient
through order six as follows. Once all earlier coefficients
are constant, the coefficient at order \(i\) commutes with \(r\)
and is \(H_i(r)\) for a polynomial \(H_i\). Subtract the series
\(\epsilon^i H_i(\mathcal R)\), which commutes with \(\mathcal R\).
It has the prescribed coefficient at order \(i\), and if
\(H_i\ne0\) then \(4\deg H_i\le hi\).
By (3) its coefficient at any later order \(i+j\), \(j>0\),
has degree at most \(4\deg H_i+hj-1<h(i+j)\).
Thus these formal subtractions preserve the binary degree
bounds and never change the binary term of maximal allowed
degree at a later order. They are only used to test necessary
coefficients, without any polynomial change of target
coordinates other than the already specified multiple of \(R\).

For series with coefficient degrees bounded by \(hi\),
taking at order \(i\) the homogeneous binary term of degree
\(hi\) preserves sums and products. We use it to isolate the
coefficients that the preceding subtractions cannot change.

## Excluding a quadratic binary part of \(Q_3\)

Suppose \(\deg q_3=2\) and let \(L\ne0\) be its top quadratic
part. The equation for \(Z\) gives the bound \(\deg Z_i\le2i\)
by recursion, since the coefficient at order \(i\) of the
original scaled equation has binary degree at most \(2i\).
The same bound holds for \(\Pi\), using the expression for
\(P_4\). For these bounds \(h=2\). In the highest binary
homogeneous terms the slice equation and the expression for
\(\Pi\) become
\[
 Z^4+A Z=1,\qquad
 \Pi=Z^5+\tfrac54 A Z^2=Z(1+A Z/4),
 \qquad A=\epsilon L.
\]
The appearances of \(Z,\Pi\) in this display denote their
highest binary terms at each order. All the lower \(Q_j\) and
\(P_j\) are below the indicated binary degree bounds and
contribute zero to these terms.

Solving the scalar relation at \(Z=1\) gives
\[
 Z=1-\tfrac14 A-\tfrac1{32}A^2
           +\tfrac7{2048}A^4+\tfrac1{512}A^5+O(A^6).
\]
Thus the coefficient of \(A^5\) in \(Z(1+AZ/4)\) is
\(1/256\). After the successive removals, the coefficient of
\(\epsilon^5\) must commute with \(r\) and still has nonzero
binary degree ten. This is impossible for a polynomial in a
degree-four \(r\). Therefore \(\deg q_3\le1\).

## Excluding the remaining binary terms

Under this bound (2) has degree at most three, so is constant.
In particular \(p_3\) has degree at most two and its quadratic
part is \((5/4)K+(5/32)L^2\), where \(L\) is now the linear
part of \(q_3\) and \(K\) the quadratic part of \(q_2\).
The bounds on \(Z_i,\Pi_i\) above hold with \(h=1\).
In their top binary homogeneous terms write
\(A=\epsilon L,\ B=\epsilon^2K,\ W=Z^2\). Their relations are
\[
 W^2+AW+B=1,\qquad
 \Pi=Z\left(1+AW/4+B/4+5A^2/32\right).
\tag{4}
\]
For example the top \(P_3\) term is
\(((5/4)B+(5/32)A^2)Z\), and the top \(P_4\) term is
\((5/4)AZ^3\), giving (4).

The \(\epsilon^3\) term in (4) is
\((5/128)L(L^2-8K)\). Since a coefficient of degree at
most three commuting with \(r\) is constant, this must vanish.
If \(L\ne0\), put \(K=L^2/8\). Then the scalar expansions are
\[
 \begin{aligned}
 W&=1-\tfrac12 A+\tfrac1{16}A^2-\tfrac1{512}A^4+O(A^6),\\
 Z&=1-\tfrac14 A-\tfrac1{1024}A^4-\tfrac1{4096}A^5+O(A^6).
 \end{aligned}
\]
Equation (4) has \(A^5\)-coefficient \(-1/1024\).
The coefficient at \(\epsilon^5\) must therefore have binary
degree five, again impossible. Hence \(L=0\).

Now (4) becomes \(\Pi=(1-B)^{1/4}(1+B/4)\). Its \(B^3\)
coefficient is \(-5/64\). The \(\epsilon^6\) term has
nonzero binary degree six if \(K\ne0\), which is likewise
impossible. Consequently \(K=0\).

It follows that \(Q_3\) is univariate in \(z\) and \(Q_2\)
has no term quadratic in \(x,y\). Hence
\[
                   Q=f(z)+a(z)x+b(z)y
\]
with \(f\) monic quartic and \(a,b\) affine. The polynomials
\(a,b\) cannot have a common zero, since at a common zero with
at least one nonzero slope one can solve \(Q_z=0\) for \(x,y\).
If both slopes and both values vanish, \(Q\) is univariate
and has a critical point. Thus choose constants \(r_0,s_0\)
with \(r_0a+s_0b=1\). Then \((z,Q,-s_0x+r_0y)\) is a
coordinate completion, with inverse degree at most four by
the constant determinant of its two-by-two fibre matrix.
The
[bounded coordinate-completion theorem](LOW_DEGREE_COORDINATE_KELLER_MAPS_IN_THREE_VARIABLES.md)
makes \(F\) invertible. The proof transfers to characteristic
zero by coefficient-field embedding and descent of the
polynomial inverse.
