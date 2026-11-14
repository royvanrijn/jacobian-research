"""Finite cusp-locus certificates from canonical rational norm matrices."""
from pathlib import Path
import sys,json,time,hashlib
import sympy as S
from sympy.polys.matrices import DomainMatrix
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'jc2-r7-quintic-incidence'))
from finite import Algebra,rank,modular
from norms import nonimmersion_norm
sys.path.insert(0,str(HERE.parent/'jc2-r7-complement'))
import geometry as G
A,M,s,t=S.symbols('A M s t')

def analyze(family,prime):
    start=time.monotonic()
    incidence=json.loads((HERE.parent/'jc2-r7-quintic-incidence/incidence.json').read_text())
    norm=next(v for v in incidence['nonimmersion_norms'] if v['family']==family)
    canonical=nonimmersion_norm(family)
    assert canonical['matrix']==norm['matrix'] and canonical['norm']==norm['norm']
    row=next(v for v in json.loads((HERE.parent/'jc2-r7-quintic-normalizations/normalizations.json').read_text())['families'] if v['family']==family)
    par=S.Symbol(row['parameter']);loc=dict(A=A,M=M,s=s,**{str(par):par});parse=lambda v:S.sympify(v,locals=loc)
    saved=json.loads((HERE.parent/'jc2-r7-quintic-incidence'/('nonimmersion-parameters-'+family+'.json')).read_text())
    if family=='L':
        w=S.Symbol('w');h=S.sympify(saved['primitive_polynomial'],locals={'w':w}).subs(w,par**3)
    else:h=parse(saved['primitive_polynomial'])
    assert S.rem(parse(norm['norm']),h,par)==0
    assert S.gcd(h,S.diff(h,par))==1 and h.subs(par,0)!=0
    def modpoly(value):
        return S.Poly(sum(modular(co,prime)*par**m[0] for m,co in S.Poly(value,par).terms()),par,modulus=prime)
    hp=modpoly(h).monic()
    assert hp.degree()==S.degree(h,par)
    reduce=lambda value:modpoly(value).rem(hp)
    matrix=S.Matrix([[modpoly(parse(v)).as_expr() for v in r] for r in norm['matrix']])
    domain=S.GF(prime).poly_ring(par)
    cof=[]
    for j in range(7):
        minor=DomainMatrix.from_Matrix(matrix.minor_submatrix(j,0)).convert_to(domain).det().as_expr()
        cof.append(reduce((-1)**j*minor))
    # This is the adjugate row of the characteristic-zero norm matrix,
    # reduced modulo h and the chosen prime. Verify its left-kernel identity.
    for j in range(7):
        assert sum((cof[k]*reduce(matrix[k,j]) for k in range(7)),S.Poly(0,par,modulus=prime)).rem(hp).is_zero
    assert reduce(parse(norm['norm'])).is_zero
    checkunits={}
    def unit(label,value):
        value=reduce(value);g=S.gcd(value,hp);checkunits[label]=g.degree()
        assert g.degree()==0,label
        return value
    denominator=unit('cofactor_unit',cof[0].as_expr())
    inverse=S.invert(denominator,hp)
    pos={tuple(v):i for i,v in enumerate(row['standard_monomials'])}
    av=(cof[pos[(1,0)]]*inverse).rem(hp);mv=(cof[pos[(0,1)]]*inverse).rem(hp)
    def evaluate(expression):
        out=S.Poly(0,par,modulus=prime)
        for (i,j),co in S.Poly(expression,A,M).terms():out=(out+reduce(co)*av.pow(i)*mv.pow(j)).rem(hp)
        return out
    for expression in row['necessary_sufficient_equations']:assert evaluate(parse(expression)).is_zero
    q=parse(row['normalization_q']);b,c,e,f,n=[evaluate(S.Poly(q,s).nth(j)) for j in [4,3,2,1,0]]
    linear=(e.mul_ground(2)-av*b*pow(3,-1,prime)*4).rem(hp)
    linear=unit('unique_critical_linear',linear.as_expr())
    constant=(av*av*modular(S.Rational(5,9),prime)-c*av+f).rem(hp)
    critical=(-constant*S.invert(linear,hp)).rem(hp)
    assert (critical*critical*3+av).rem(hp).is_zero
    jet=(critical.pow(3)*20+b*critical.pow(2)*6-e).rem(hp)
    clean=(e+b*critical.pow(2)*3-critical.pow(3)*2).rem(hp)
    unit('ordinary_cusp_jet',jet.as_expr());unit('clean_cusp_fibre',clean.as_expr())
    triple=(av*av+av*b*b-b*e-av*c+f).rem(hp)
    unit('no_triple_B',triple.as_expr())
    low=evaluate(parse(row['implicit_difference']['(2, 0)']))
    unit('contact_exactly38',low.as_expr())
    # Divide the collision quartic by its known diagonal root v=2*s0.
    r2=(b+2*critical).rem(hp)
    r1=(av+2*critical*r2).rem(hp)
    r0=(2*av*b-e+2*critical*r1).rem(hp)
    assert (-av*av-f+c*av+2*critical*r0).rem(hp).is_zero
    pair_variable=S.Symbol('v')
    cubic=pair_variable**3+r2.as_expr()*pair_variable**2+r1.as_expr()*pair_variable+r0.as_expr()
    pair_algebra=Algebra((par,pair_variable),(hp.degree(),3),[hp.as_expr(),cubic],prime)
    pair_matrix=pair_algebra.matrix(S.diff(cubic,pair_variable));pair_rank,_=rank(pair_matrix,prime)
    uc=(mv-critical.pow(3)*2).rem(hp)
    vc=(critical.pow(5)+b*critical.pow(4)+c*critical.pow(3)+e*critical.pow(2)+f*critical+n).rem(hp)
    pc=G.F.p.subs(G.F.line) if family=='L' else G.K.P
    qc=G.F.q.subs(G.F.line) if family=='L' else G.K.Q
    pc=pc.subs(G.F.t,t);qc=qc.subs(G.F.t,t)
    field_expression=lambda expr:S.expand(sum(modular(co,prime)*par**i*t**j for (i,j),co in S.Poly(expr,par,t).terms()))
    algebra=Algebra((par,t),(hp.degree(),6),[hp.as_expr(),field_expression(pc)-uc.as_expr()],prime)
    product=algebra.matrix(field_expression(qc)-vc.as_expr());r,pivots=rank(product,prime)
    result=dict(family=family,prime=prime,parameter=str(par),base_degree=hp.degree(),unit_gcd_degrees=checkunits,
        cusp_off_C_rank=r,cusp_off_C_dimension=len(product),three_node_rank=pair_rank,three_node_dimension=len(pair_matrix),A=str(av.as_expr()),M=str(mv.as_expr()),
        critical_parameter=str(critical.as_expr()),cusp=[str(uc.as_expr()),str(vc.as_expr())],
        h_mod_prime=str(hp.as_expr()),adjugate_row_mod_prime=[str(v.as_expr()) for v in cof],
        cusp_off_C_matrix_sha256=hashlib.sha256(json.dumps(product,separators=(',',':')).encode()).hexdigest(),
        three_node_matrix_sha256=hashlib.sha256(json.dumps(pair_matrix,separators=(',',':')).encode()).hexdigest(),
        elapsed_seconds=time.monotonic()-start,status='EXACT_MODULAR_CERTIFICATE; characteristic-zero meaning proved from the rational norm, cofactor unit and finite-free quotients.')
    return result

def main():
    family=sys.argv[1];prime=int(sys.argv[2]) if len(sys.argv)>2 else (103 if family=='L' else 101)
    result=analyze(family,prime)
    (HERE/('cusps-'+family+'.json')).write_text(json.dumps(result,indent=2)+'\n');print(family,result['cusp_off_C_rank'],result['cusp_off_C_dimension'],result['unit_gcd_degrees'],flush=True)

if __name__=='__main__':main()
