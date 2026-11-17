"""PF111 complete peripheral modules, conditional on certified core generation."""
import math
import sympy as S
from sympy.matrices.normalforms import hermite_normal_form,smith_normal_form
import full_peripheral_linking as P


def group_of_columns(matrix):
    matrix=S.Matrix(matrix)
    if any(not x.is_Integer for x in matrix):
        raise ValueError('integral presentation required')
    D=smith_normal_form(matrix,domain=S.ZZ)
    diagonal=[abs(int(D[i,i])) for i in range(min(D.shape)) if D[i,i]!=0]
    return dict(free=matrix.rows-len(diagonal),torsion=[x for x in diagonal if x>1])


def core_generation(data):
    bs=data['boundaries'];rows=data['relation_rows']
    m=[b['meridian'] for b in bs];l=[b['longitude'] for b in bs]
    closed=P.reduced_h1(rows,m)
    if closed['free']:
        raise ValueError('rational homology sphere normal filling required')
    quotient=P.reduced_h1(rows,m+l)
    exterior=P.reduced_h1(rows)
    return dict(closed_H1=closed,exterior_H1=exterior,unseen_quotient=quotient,
                all_cores_generate=(quotient==dict(free=0,torsion=[])),
                required_outside_prime_count=len(quotient['torsion']),
                affine_realization_claimed=False)


def linking_module(matrix,*,generating_cores_certified):
    if generating_cores_certified is not True:
        raise ValueError('full integral generation certificate required; rational linking alone is insufficient')
    F=S.Matrix(matrix)
    if not F.rows or F.rows!=F.cols or F!=F.T or any(not x.is_Rational for x in F):
        raise ValueError('nonempty rational symmetric linking matrix required')
    q=math.lcm(*(int(x.q) for x in F))
    integer_columns=S.Matrix.hstack(q*S.eye(F.rows),q*F.T)
    integer_columns=S.Matrix([[int(x) for x in row] for row in integer_columns.tolist()])
    H=hermite_normal_form(integer_columns)
    basis=H/q
    meridians=basis.inv()
    longitudes=basis.inv()*F.T
    if any(not x.is_Integer for x in list(meridians)+list(longitudes)):
        raise AssertionError('rational lattice basis failed integral membership')
    return dict(F=F,basis=basis,meridians=meridians,longitudes=longitudes,
                normal_filling_H1=group_of_columns(meridians),
                row_denominators=[math.lcm(*(int(x.q) for x in F.row(i))) for i in range(F.rows)],
                affine_realization_claimed=False)


def further_filling(module,slopes):
    if len(slopes)!=module['F'].rows:
        raise ValueError('one primitive integral filling slope per boundary required')
    if any(not isinstance(a,int) or not isinstance(b,int) or math.gcd(a,b)!=1 for a,b in slopes):
        raise ValueError('primitive integer filling slopes required')
    b=S.diag(*[x[0] for x in slopes]);c=S.diag(*[x[1] for x in slopes])
    return group_of_columns(module['meridians']*b+module['longitudes']*c)


def projective_action(a,b,c,d):
    if (a*d-b*c)%5!=1:
        raise ValueError('determinant one over F5 required')
    points=list(range(5))+[None];out=[]
    for x in points:
        num,den=(a,c) if x is None else ((a*x+b)%5,(c*x+d)%5)
        y=None if den%5==0 else num*pow(den,-1,5)%5
        out.append(points.index(y))
    return tuple(out)


def six_sheet_cover(cables=(),reverse=False):
    cables=list(cables)
    if any(not isinstance(r,int) or not isinstance(k,int) or r<=1 or k<=1 or
           r%5!=1 or k%5!=0 or math.gcd(r,k)!=1 for r,k in cables):
        raise ValueError('coprime ordinary cables r=1,k=0 mod5 required')
    rho=[projective_action(0,-1,1,0),projective_action(0,-1,1,1)]
    base=P.peripheral_cover(2,3,[],rho)
    lam=base['longitude_monodromy']
    return P.peripheral_cover(2,3,cables,rho+[lam]*len(cables),reverse)


def six_sheet_transfer(cables):
    cables=list(cables)
    if any(not isinstance(r,int) or not isinstance(k,int) or r<=1 or k<=1 or
           r%5!=1 or k%5!=0 or math.gcd(r,k)!=1 for r,k in cables):
        raise ValueError('proved cable sector required')
    R=math.prod(r for r,k in cables)
    core_order=2//math.gcd(2,R);hidden_order=4//(2//math.gcd(2,R))
    return dict(winding_product=R,matrix=[[-5*R*R,R*R],[R*R,(1-R*R)//5]],
                closed_H1=dict(free=0,torsion=[4]),
                core_orders=[core_order,core_order],
                exterior_H1=dict(free=2,torsion=[hidden_order]),
                unseen_quotient=dict(free=0,torsion=[hidden_order]),
                required_outside_prime_count=1,remaining_Euler_budget=0,
                actual_Keller_family_excluded=True,realization_claimed=False)
