#!/usr/bin/env python3
"""Exact binary-frame and one-parameter radical tests; no global prefix search."""
from pathlib import Path
import hashlib,json,math
import sympy as S

Z,W,T,s,t=S.symbols('Z W T s t')


def frame_matrix(poly):
    poly=S.Poly(S.expand(poly),Z,W,T)
    assert poly.total_degree()<=3
    pieces={d:S.Add(*(c*Z**a[0]*W**a[1]*T**a[2] for a,c in poly.terms() if sum(a)==d)) for d in [2,3]}
    direction={Z:s*s,W:-t*t/2,T:s*t}
    forms=[S.expand(S.diff(pieces[3],x).subs(direction,simultaneous=True)) for x in [Z,W,T]]
    forms.append(S.expand(pieces[2].subs(direction,simultaneous=True)))
    columns=[]
    for f in forms:
        for j in range(4):
            product=S.Poly(S.expand(f*s**(3-j)*t**j),s,t)
            columns.append([product.coeff_monomial(s**(7-k)*t**k) for k in range(8)])
    return S.Matrix(8,16,lambda i,j:columns[j][i]),forms


def in_radical(g,ideal,variables):
    inverse=S.Dummy('inverse')
    basis=S.groebner(list(ideal)+[1-inverse*g],inverse,*variables,domain=S.QQ)
    return basis.reduce(S.Integer(1))[1]==0


def main():
    p=Path(__file__).resolve().parent
    v=S.Matrix([s*s,-t*t/2,s*t])
    assert S.expand(2*v[0]*v[1]+v[2]**2)==0
    coefficients=S.symbols('c0:20')
    monomials=[Z**i*W**j*T**k for d in range(4) for i in range(d+1) for j in range(d-i+1) for k in [d-i-j]]
    assert len(monomials)==20
    generic=sum(c*m for c,m in zip(coefficients,monomials))
    F3=sum(c*m for c,m in zip(coefficients,monomials) if S.Poly(m,Z,W,T).total_degree()==3)
    F2=sum(c*m for c,m in zip(coefficients,monomials) if S.Poly(m,Z,W,T).total_degree()==2)
    directional=lambda f:sum(v[i]*S.diff(f,x) for i,x in enumerate([Z,W,T]))
    expected=2*sum(S.diff(F3,x).subs({Z:v[0],W:v[1],T:v[2]},simultaneous=True)*x for x in [Z,W,T])+2*F2.subs({Z:v[0],W:v[1],T:v[2]},simultaneous=True)
    assert S.expand(directional(directional(generic))-expected)==0
    universal,forms=frame_matrix(generic)
    assert universal.shape==(8,16) and math.comb(16,8)==12870
    assert all(S.Poly(entry,*coefficients).total_degree()<=1 for entry in universal)
    examples=[('zero',S.Integer(0),True),('pure_linear',Z+W+T,True),
              ('arbitrary_quadratic',Z**2+2*W**2+T*T+3*Z*W+W*T,True),
              ('sharp_shear',W+Z*T-Z**3/2,True),
              ('homogeneous_Z_cubic',Z**3,True),
              ('nonisotropic_cube_frame_not_null',T**3,True),
              ('nonisotropic_cube_obstructed',T**3+Z**2+W**2,False),
              ('sharp_cusp',Z**2*T+W**3,False),
              ('isotropic_repeated_with_W_squared',Z**3+W**2,False),
              ('concurrent_nodes_obstructed',Z*W*T+Z**2+W**2,False)]
    rows=[]
    for name,poly,has_frame in examples:
        matrix,ff=frame_matrix(poly)
        rank=matrix.rank()
        nonzero=[S.Poly(f,s,t) for f in ff if f!=0]
        common=S.gcd_list([f.as_expr() for f in nonzero]) if nonzero else S.Integer(0)
        common_degree=None if common==0 else S.Poly(common,s,t).total_degree()
        assert (rank<8)==has_frame
        assert (common==0 or common_degree>0)==has_frame
        row={'name':name,'rank':rank,'has_frame':has_frame,'common_quartic_gcd':str(common),'gcd_degree':common_degree}
        if not has_frame:
            pivots=list(matrix.rref()[1]);assert len(pivots)==8
            determinant=S.factor(matrix[:,pivots].det());assert determinant!=0
            row.update({'maximal_minor_columns':pivots,'maximal_minor':str(determinant)})
        rows.append(row)
    # Standalone binary-quartic test retains a common point at infinity.
    assert S.gcd(s**4,s**3*t)==s**3
    # Exact completed family: no twenty-variable/global run is attempted.
    a=S.symbols('a')
    cusp=json.loads((p/'cuspidal_cubic_threshold_verification.json').read_text())
    assert cusp['completed'] and cusp['witness_M10']==1508855040000
    coefficient=S.Integer(cusp['witness_M10'])
    assert not in_radical(a**8,[],[a])
    assert in_radical(a**8,[coefficient*a**10],[a])
    assert S.rem(a**8,a**10,a)!=0
    # An explicit certificate is independently expanded, not merely inferred
    # from a Groebner verdict. Set x=u*a^8; x^2 is in (a^10).
    u=S.symbols('u')
    certificate=(1+u*a**8)*(1-u*a**8)+u*u*a**6*a**10
    assert S.expand(certificate)==1
    closure=json.loads((p/'cubic_gmc_closure_verification.json').read_text())
    assert closure['completed']
    result={'completed':True,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'prior_cubic_checker_sha256':closure['checker_sha256'],
            'prior_cusp_checker_sha256':cusp['checker_sha256'],
            'universal_matrix_shape':[8,16],'coefficient_count':20,'quadratic_cubic_coefficient_count':16,
            'maximal_minor_count':math.comb(16,8),'maximal_minor_degree':8,
            'frame_examples':rows,
            'completed_cusp_line':{'least_prefix':10,'first_nine_moments':'identically zero','M10':str(coefficient*a**10),
                                   'frame_ideal':'(a^8)','moment_ideal_at10':'(a^10)','common_radical':'(a)',
                                   'radical_membership_verified':True,'ideal_membership_false':True,
                                   'Rabinowitsch_certificate':'(1+u*a^8)*(1-u*a^8)+u^2*a^6*a^10=1'},
            'global_prefix_algorithm':'Proved terminating using the retained internal frame-plus888 theorem and exact radical-membership tests; full coefficient-space search not executed.',
            'numeric_global_prefix':None,'global_prefix_lower_bound':10,
            'written_proof_certified_by_checker':False,
            'proof_boundary':'Exact universal conic/directional derivative/matrix-size algebra, stated specialized ranks/gcds with nonzero maximal minors, and the one-parameter radical/ideal/plateau certificate. Universal common-root equivalence and global least-prefix computability are written deductions; no numerical global threshold is computed.'}
    (p/'cubic_frame_prefix_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:result[key] for key in ['completed','universal_matrix_shape','maximal_minor_count','numeric_global_prefix','global_prefix_lower_bound','written_proof_certified_by_checker']}))


if __name__=='__main__':main()
