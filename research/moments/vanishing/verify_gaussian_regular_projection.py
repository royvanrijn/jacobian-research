"""Auxiliary formal Gaussian projection and stronger scalar D cofactors.

The all-order formal projection/regular connection, complementary-minor argument,
polynomial extension and complete analytic/frame cover remain written inputs.
Only retained seed columns are replayed; no previous checker main is run.
"""
from pathlib import Path
from itertools import combinations
from functools import lru_cache
import hashlib,json,time
import sympy as s
from sympy.polys.matrices import DomainMatrix

P=Path(__file__).resolve().parent
t,D,u,z,q,p,r=s.symbols('t D u z q p r')
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
tri=lambda n:n*(n-1)//2

def gaussian_identities():
    # H=h0+h1*z+h2*z^2+h3*z^3; h_i' are independent u derivatives.
    hs=s.symbols('h0:4');hp=s.symbols('hp0:4');h0,h1,h2,h3=hs
    nmax=8;E=[s.Integer(1)]
    for n in range(1,2*nmax+1):
        E.append(s.expand(sum(k*hs[k]*E[n-k]/t**2 for k in range(1,min(3,n)+1))/n))
    @lru_cache(None)
    def wick_power(a):
        return s.Poly(sum((-t*t/2)**k*s.factorial(2*k)*E[2*k-a]*D**k/s.factorial(k)
                          for k in range(nmax+1) if 2*k>=a),D)
    counts={'transverse_twisted_derivative':0,'u_connection':0,'normalized_D_connection':0}
    bound=4
    for a in range(7):
        transverse=wick_power(a+1).as_expr()+D*sum(k*hs[k]*wick_power(a+k-1).as_expr() for k in range(1,4))
        if a:transverse+=t*t*D*a*wick_power(a-1).as_expr()
        for k in range(bound+1):assert s.expand(transverse).coeff(D,k)==0
        counts['transverse_twisted_derivative']+=bound+1
        current=wick_power(a).as_expr()
        derivative=sum(hp[k]*s.diff(current,hs[k]) for k in range(1,4))
        commutator=sum(hp[k]*wick_power(a+k).as_expr() for k in range(4))-hp[0]*current-t*t*derivative
        for k in range(bound+1):assert s.expand(commutator).coeff(D,k)==0
        counts['u_connection']+=bound+1
        # Derivative of normalized exp(+z^2/(2t^2 D)) integration.
        singular=s.expand(s.diff(current,D)+current/(2*D)+wick_power(a+2).as_expr()/(2*t*t*D*D))
        for k in range(-2,bound+1):assert s.expand(singular).coeff(D,k)==0
        counts['normalized_D_connection']+=bound+3
    return {'H':'h0+h1*z+h2*z^2+h3*z^3','formal_variance':'-t^2 D',
            'identity_through_D_order':bound,'input_monomials':list(range(7)),
            'coefficient_checks':counts,'total_coefficient_checks':sum(counts.values())}

def comparison_determinants():
    A=u*u+p*u+r;v0=-A/q
    # Cusp: fast row selects v; slow rows use 1,u,u^2,u^3,u^4.
    columns=[s.Integer(1),u,u*u,D*u*v0,D*u*u*v0]
    cusp=s.Matrix([[s.expand(c).coeff(u,i) for c in columns] for i in range(5)])
    assert s.factor(cusp.det())==D**2/q**2
    # Tangent: replace the v column by v+(u^2+p*u+r)/q=z/q.
    # Formal Gaussian IBP gives pi(z)=-D*pi(H_z).  H_z at z=0
    # has leading -2*u^3/q^2, independently of the full quadratic G2.
    a,b,c,e,f,g=s.symbols('a b c e f g')
    H=u*((z-A)/q)**2+a*u*u+b*u*((z-A)/q)+c*((z-A)/q)**2+e*u+f*((z-A)/q)+g
    Hz0=s.expand(s.diff(H,z).subs(z,0))
    assert Hz0.coeff(u,3)==-2/q**2 and s.degree(Hz0,u)==3
    columns=[s.Integer(1),D*u,D*D*u*u,-D**3*Hz0/q]
    tangent=s.Matrix([[s.expand(c).coeff(u,i) for c in columns] for i in range(4)])
    assert s.factor(tangent.det())==2*D**6/q**3
    ordinary=[]
    for kind,basis,powers in [
        ('nodal',[(0,0),(1,0),(0,1),(2,0),(1,1),(0,2),(0,3)],[-3,-2,-1,0,1,2]),
        ('transverse',[(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)],[-2,-1,0,1,2]),
        ('triangle',[(0,0),(1,0),(0,1),(1,1),(0,2)],[-2,-1,0,1])]:
        cols=[s.expand(u**i*(-r/u)**k) for i,k in basis]
        fast=[s.Integer(int(i==0)) for i in range(len(cols))]
        Q=s.Matrix([fast]+[[c.coeff(u,k) for c in cols] for k in powers])
        det=s.factor(Q.det());assert det!=0 and not det.has(D)
        ordinary.append({'kind':kind,'finite_basis_powers':powers,'comparison_determinant_at_D0':str(det),
                         'axis_row_scales_are_generic_D_units':True,'D_loss':0})
    return {'ordinary':ordinary,'cusp':{'finite_rank':5,'D_loss':2,'determinant_leading_term':'D^2/q^2'},
            'tangent':{'finite_rank':4,'D_loss':6,'determinant_leading_term':'2*D^6/q^3',
                       'full_quadratic_Hz_at_z0':str(Hz0),'leading_u3_coefficient':'-2/q^2'}}

def weights(kind,n):
    return [0]+[4*j-1 if kind=='cusp' else 6*j if kind=='tangent' else 2*j for j in range(1,n+1)]

def minor_bounds(kind,n):
    w=weights(kind,n);loss={'cusp':2,'tangent':6}.get(kind,0)
    # Every transformed minor has at most one unrestricted fast row.
    projection=[max(0,sum(sorted(w[j] for j in range(n+1) if j!=k)[:n-1])-loss) for k in range(n+1)]
    old=[tri(n)]+[max(0,tri(n)-k) for k in range(1,n+1)]
    return [max(a,b) for a,b in zip(projection,old)]

def determinant_cofactors(cols,n):
    ring=s.QQ.poly_ring(t);out=[]
    for k in range(n+1):
        selected=cols[:k]+cols[k+1:]
        mat=DomainMatrix.from_list_sympy(n,n,[[x.as_expr() for x in row] for row in zip(*selected)]).convert_to(ring)
        out.append((-1)**k*s.Poly(ring.to_sympy(mat.det()),t,domain=s.QQ))
    assert all(sum((out[k]*cols[k][i] for k in range(n+1)),s.Poly(0,t)).is_zero for i in range(n))
    return out

def retained_models():
    out=[]
    for filename in ['ordinary_rank_one_lattice_certificate.json','cusp_rank_one_lattice_certificate.json','tangent_rank_one_lattice_certificate.json']:
        for m in json.loads((P/filename).read_text())['models']:
            if m['D_degree']!=1:continue
            kind=m.get('kind','cusp' if filename.startswith('cusp') else 'tangent')
            n=len(m['seed_columns']);columns=[[s.Poly(int(i==0),t,domain=s.QQ) for i in range(n)]]
            columns += [[s.Poly(s.sympify(x,locals={'t':t}),t,domain=s.QQ) for x in row['N']] for row in m['seed_columns']]
            den=s.Poly(s.sympify(m['D'],locals={'t':t}),t,domain=s.QQ)
            co=determinant_cofactors(columns,n);needs=minor_bounds(kind,n);checks=[]
            for k,C in enumerate(co):
                C.exquo(den**needs[k])
                H=s.Poly(1,t) if k==0 else s.Poly(s.sympify(m['seed_columns'][k-1]['H'],locals={'t':t}),t,domain=s.QQ)
                reduced=(H*C).exquo(den**needs[0])
                val=min(e[0] for e,c in reduced.terms() if c) if reduced else None
                if k==0:assert reduced and reduced.degree()<={'nodal':87,'transverse':97,'triangle':98,'cusp':77,'tangent':66}[kind]
                else:assert val is None or val>=3*k
                checks.append({'order':k,'minor_D_power':needs[k],'common_scalar_D_power':needs[0],
                               'reduced_scalar_degree':int(reduced.degree()) if reduced else -1,
                               'reduced_scalar_valuation':val})
            out.append({'origin':filename,'origin_sha256':sha(P/filename),'kind':kind,
                        'label':m.get('label',m.get('mode')),'D':m['D'],'rank':n,'checks':checks})
            print('Retained full cofactor relation',kind,out[-1]['label'],'passed',flush=True)
    assert len(out)==15 and sum(len(m['checks']) for m in out)==100
    return out

def arithmetic():
    grades={'nodal':[0,1,1,2,2,2,3],'transverse':[0,1,1,2,2,2],
            'triangle':[0,3,3,4,5],'cusp':[0,1,2,1,3,4],'tangent':[0,2,4,3,6]}
    # Triangle rows use the retained two-axis completed scales.
    bs=[(0,0),(1,0),(0,1),(1,1),(0,2)]
    grades['triangle']=[i+k+int(s.ceiling(i+s.Rational(k,2)))+int(s.ceiling(s.Rational(i,2)+k)) for i,k in bs]
    rows=[]
    for kind,gs in grades.items():
        slope={'nodal':5,'transverse':7,'triangle':9,'cusp':7,'tangent':9}[kind]
        offset={'transverse':2,'cusp':1}.get(kind,0);cases=[]
        for n in range(1,len(gs)+1):
            needs=minor_bounds(kind,n)
            for k in range(1,n+1):
                Hpower=4*k-1 if kind=='cusp' else 6*k if kind=='tangent' else 2*k
                assert needs[k]+Hpower>=needs[0],(kind,n,k,needs)
            for selected in combinations(range(len(gs)),n):
                if 0 not in selected:continue
                raw=slope*n*(n+1)//2-offset*n-sum(gs[k] for k in selected)
                cases.append({'rank':n,'rows':selected,'slow_seed_D_orders':weights(kind,n)[1:],
                              'minor_D_bounds':needs,'common_scalar_D_power':needs[0],
                              'raw_degree_bound':raw,'reduced_degree_bound':raw-needs[0]})
        gate=max(c['reduced_degree_bound'] for c in cases)
        assert gate=={'nodal':87,'transverse':97,'triangle':98,'cusp':77,'tangent':66}[kind]
        rows.append({'kind':kind,'linear_D_gate':gate,'all_rank_row_cases':cases})
    cover={'quadratic':7,'smooth':94,'nodal':101,'cusp':109,'transverse':106,'tangent':105,'triangle':103,'concurrent':104,'repeated':99}
    assert max(cover.values())==109
    return rows,cover

def main():
    start=time.monotonic();wick=gaussian_identities();comparison=comparison_determinants()
    print('Formal Gaussian coefficient identities and comparison determinants passed',flush=True)
    models=retained_models();cases,cover=arithmetic()
    certificate={'formal_gaussian':wick,'comparison_determinants':comparison,'retained_models':models,
                 'degree_cases':cases,'global_cover':cover}
    f=P/'gaussian_regular_projection_certificate.json';f.write_text(json.dumps(certificate,indent=2,default=int)+'\n')
    receipt={'completed':True,'checker_sha256':sha(Path(__file__)),'certificate_sha256':sha(f),
             'runtime_seconds':round(time.monotonic()-start,3),'sympy_version':s.__version__,
             'formal_gaussian_coefficient_checks':wick['total_coefficient_checks'],
             'comparison_determinants':5,'retained_full_cofactor_models':len(models),
             'full_scalar_cofactor_checks':sum(len(m['checks']) for m in models),
             'linear_D_gates':{c['kind']:c['linear_D_gate'] for c in cases},
             'constant_D_gates_retained':{'nodal':101,'transverse':106,'triangle':103,'cusp':109,'tangent':105},
             'global_cover':cover,'global_sufficient_pure_prefix':109,
             'odd_sufficient_pure_prefix_unchanged':98,'homogeneous_sufficient_pure_prefix_unchanged':96,
             'least_prefix_interval':[20,109],'odd_least_prefix_interval':[20,98],
             'homogeneous_least_prefix_interval':[14,96],
             'least_prefix_computed':False,'universal_written_proof_certified_by_checker':False,
             'independent_replay':False,'external_review':False,'formal_verification':False,
             'boundary':'Finite formal Gaussian coefficient identities, five generic comparison determinants, 15 retained original Gaussian seed/cofactor families with 100 scalar cofactor checks, all rank/row degree bounds and scalar valuations. All-order formal projection and regular finite connection, determinant-loss/complementary-minor argument, polynomial extension, full coefficient/frame and analytic exclusion proofs remain written inputs. No previous main, fresh literature review, radical search or producer package is run.'}
    (P/'gaussian_regular_projection_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))

if __name__=='__main__':main()
