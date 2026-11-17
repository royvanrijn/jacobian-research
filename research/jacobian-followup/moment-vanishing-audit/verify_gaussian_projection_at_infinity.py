"""Exact auxiliary formal Morse/projection-at-infinity degree checks.
Universal formal projection, connection, comparison-minor and analytic/frame
arguments remain written inputs. Retained original checker mains are not run.
"""
from pathlib import Path
from itertools import combinations
import hashlib,json,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
P=Path(__file__).resolve().parent
t,x,u,z,w,D,q,p,r=s.symbols('t x u z w D q p r')
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
T=lambda n:n*(n-1)//2

def truncate(f,n=5):return s.series(f,x,0,n).removeO().expand()

def formal_morse():
    h0,h1,h2,h3=s.symbols('h0 h1 h2 h3');H=lambda f:h0+h1*f+h2*f*f+h3*f**3
    phi=s.Integer(0)
    for k in range(6):phi=truncate(-x*(h1+2*h2*phi+3*h3*phi*phi),6)
    critical=truncate(phi+x*(h1+2*h2*phi+3*h3*phi*phi),6);assert critical==0
    Phi=truncate(H(phi)+phi*phi/(2*x),5);assert Phi.subs(x,0)==h0
    Hess=1+2*x*h2+6*x*h3*phi
    # Exact normal Morse coordinate y=w*sqrt(Hess+2*x*h3*w).
    kappa=truncate(Hess+2*x*h3*w-1,5)
    root=truncate(sum(s.binomial(s.Rational(1,2),j)*kappa**j for j in range(5)),5)
    y=truncate(w*root,5)
    error=truncate((phi+w)**2/2+x*H(phi+w)-x*Phi-y*y/2,5)
    assert error==0 and root.subs(x,0)==1
    # The residual phase derivative under d/dt=-x^2*d/dx is regular.
    derivative=s.expand(-x*x*s.diff(Phi/(D*x),x))
    assert s.expand(derivative-(Phi-x*s.diff(Phi,x))/D)==0
    assert s.Poly(derivative,x).degree()>=0
    return {'normal_critical_point_through_x5':str(phi),'critical_equation_zero_through_x5':True,
            'residual_phase_through_x4':str(Phi),'residual_phase_constant':'h0',
            'Morse_identity_zero_through_x4':True,'Morse_Jacobian_at_x0':1,
            'Gaussian_variance':'-D*x^2','normalized_prefactor_x_order':1,
            'time_derivation':'-x^2 partial_x','residual_phase_time_derivative':'(Phi-x*Phi_x)/D',
            'regular_time_connection_required':True}

def comparisons():
    v0=-(u*u+p*u+r)/q
    cusp=s.Matrix([[s.expand(c).coeff(u,i) for c in [1,u,u*u,x*u*v0,x*u*u*v0]] for i in range(5)])
    assert s.factor(cusp.det())==x*x/q**2
    a,b,c,e,f,g=s.symbols('a b c e f g')
    V=(z-u*u-p*u-r)/q
    G=u*V*V+a*u*u+b*u*V+c*V*V+e*u+f*V+g
    Hz0=s.expand(D*s.diff(G,z).subs(z,0));assert Hz0.coeff(u,3)==-2*D/q**2
    tangent=s.Matrix([[s.expand(v).coeff(u,i) for v in [1,x*u,x*x*u*u,-x**3*Hz0/q]] for i in range(4)])
    assert s.factor(tangent.det())==2*D*x**6/q**3
    ordinary=[]
    for kind,bs,powers in [
        ('nodal',[(0,0),(1,0),(0,1),(2,0),(1,1),(0,2),(0,3)],[-3,-2,-1,0,1,2]),
        ('transverse',[(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)],[-2,-1,0,1,2]),
        ('triangle',[(0,0),(1,0),(0,1),(1,1),(0,2)],[-2,-1,0,1])]:
        cols=[s.expand(u**i*(-r/u)**j) for i,j in bs]
        Q=s.Matrix([[int(i==0) for i in range(len(bs))]]+[[c.coeff(u,k) for c in cols] for k in powers])
        det=s.factor(Q.det());assert det and not det.has(x)
        ordinary.append({'kind':kind,'finite_basis_powers':powers,'determinant_at_x0':str(det),'loss':0})
    return {'ordinary':ordinary,'cusp':{'loss':2,'leading_determinant':'x^2/q^2'},
            'tangent':{'loss':6,'leading_determinant':'2D*x^6/q^3','leading_Hz_u3':'-2D/q^2'},
            'constant_D_row_units_omitted_only_from_comparison':True}

def cofactors(cols,n):
    ring=s.QQ.poly_ring(t);out=[]
    for k in range(n+1):
        selected=cols[:k]+cols[k+1:]
        mat=DomainMatrix.from_list_sympy(n,n,[[c.as_expr() for c in row] for row in zip(*selected)]).convert_to(ring)
        out.append((-1)**k*s.Poly(ring.to_sympy(mat.det()),t,domain=s.QQ))
    assert all(sum((out[k]*cols[k][i] for k in range(n+1)),s.Poly(0,t)).is_zero for i in range(n))
    return out

def replay():
    out=[]
    for name in ['ordinary_rank_one_lattice_certificate.json','cusp_rank_one_lattice_certificate.json','tangent_rank_one_lattice_certificate.json']:
        for m in json.loads((P/name).read_text())['models']:
            if m['D_degree']!=0:continue
            kind=m.get('kind','cusp' if name.startswith('cusp') else 'tangent');n=len(m['seed_columns'])
            columns=[[s.Poly(int(i==0),t,domain=s.QQ) for i in range(n)]]
            columns += [[s.Poly(s.sympify(c,locals={'t':t}),t,domain=s.QQ) for c in row['N']] for row in m['seed_columns']]
            co=cofactors(columns,n);checks=[]
            for k,C in enumerate(co):
                H=s.Poly(1,t) if k==0 else s.Poly(s.sympify(m['seed_columns'][k-1]['H'],locals={'t':t}),t,domain=s.QQ)
                scalar=H*C;val=min(e[0] for e,c in scalar.terms() if c) if scalar else None
                if k==0:assert scalar and scalar.degree()<={'nodal':80,'transverse':91,'triangle':93,'cusp':71,'tangent':61}[kind]
                else:assert val is None or val>=3*k
                checks.append({'order':k,'scalar_degree':int(scalar.degree()) if scalar else -1,'scalar_valuation':val})
            out.append({'origin':name,'origin_sha256':sha(P/name),'kind':kind,'label':m.get('label',m.get('mode')),
                        'rank':n,'D':m['D'],'seed_rows_bound':m['seed_columns'][-1]['row_degree_bounds'],'checks':checks})
            print('Retained constant-D model',kind,out[-1]['label'],'passed',flush=True)
    assert len(out)==8 and sum(len(m['checks']) for m in out)==53
    return out

def arithmetic():
    triangle=[(0,0),(1,0),(0,1),(1,1),(0,2)]
    gs={'nodal':[0,1,1,2,2,2,3],'transverse':[0,1,1,2,2,2],
        'triangle':[i+j+int(s.ceiling(i+s.Rational(j,2)))+int(s.ceiling(s.Rational(i,2)+j)) for i,j in triangle],
        'cusp':[0,1,2,1,3,4],'tangent':[0,2,4,3,6]}
    cases=[]
    for kind,grades in gs.items():
        slope={'nodal':4,'transverse':6,'triangle':8,'cusp':6,'tangent':8}[kind]
        offset={'transverse':2,'cusp':1}.get(kind,0);loss={'cusp':2,'tangent':6}.get(kind,0)
        weight=lambda j:3*j-1 if kind=='cusp' else 5*j if kind=='tangent' else j
        rows=[]
        for n in range(1,len(grades)+1):
            saving=max(0,sum(weight(j) for j in range(1,n))-loss)
            for selected in combinations(range(len(grades)),n):
                if 0 not in selected:continue
                raw=slope*n*(n+1)//2-offset*n-sum(grades[i] for i in selected)
                rows.append({'rank':n,'rows':selected,'raw_degree_bound':raw,'slow_seed_orders':[weight(j) for j in range(1,n+1)],
                             'comparison_loss':loss,'infinity_order_saving':saving,'degree_bound':raw-saving})
        gate=max(row['degree_bound'] for row in rows)
        assert gate=={'nodal':80,'transverse':91,'triangle':93,'cusp':71,'tangent':61}[kind]
        cases.append({'kind':kind,'constant_D_gate':gate,'all_rank_row_cases':rows})
    linear={'nodal':87,'transverse':97,'triangle':98,'cusp':77,'tangent':66}
    frame={'nodal':96,'transverse':96,'triangle':96,'cusp':95,'tangent':95}
    cover={'quadratic':7,'smooth':94,'nodal':96,'cusp':95,'transverse':97,'tangent':95,'triangle':98,'concurrent':104,'repeated':99}
    for c in cases:assert cover[c['kind']]==max(c['constant_D_gate'],linear[c['kind']],frame[c['kind']])
    assert max(cover.values())==104
    return cases,cover

def main():
    start=time.monotonic();morse=formal_morse();comparison=comparisons()
    print('Formal Morse and leading comparison identities passed',flush=True)
    models=replay();cases,cover=arithmetic()
    cert={'formal_Morse':morse,'comparisons':comparison,'retained_constant_models':models,'degree_cases':cases,'global_cover':cover}
    f=P/'gaussian_projection_at_infinity_certificate.json';f.write_text(json.dumps(cert,indent=2,default=int)+'\n')
    report={'completed':True,'checker_sha256':sha(Path(__file__)),'certificate_sha256':sha(f),'runtime_seconds':round(time.monotonic()-start,3),
            'sympy_version':s.__version__,'formal_critical_order':5,'formal_Morse_order':4,'comparison_determinants':5,
            'retained_constant_D_models':8,'full_scalar_cofactor_checks':53,
            'constant_D_gates':{c['kind']:c['constant_D_gate'] for c in cases},
            'linear_D_gates_retained':{'nodal':87,'transverse':97,'triangle':98,'cusp':77,'tangent':66},
            'global_cover':cover,'global_sufficient_pure_prefix':104,'odd_sufficient_pure_prefix_unchanged':98,'homogeneous_sufficient_pure_prefix_unchanged':96,
            'least_prefix_interval':[20,104],'odd_least_prefix_interval':[20,98],'homogeneous_least_prefix_interval':[14,96],
            'least_prefix_computed':False,'universal_written_proof_certified_by_checker':False,'independent_replay':False,'external_review':False,'formal_verification':False,
            'boundary':'Finite formal critical/Morse identities, five generic leading comparison determinants, eight retained original constant-D Gaussian families with53 scalar-cofactor checks,all rank/row degree estimates and unchanged scalar valuations. All-order formal stationary projection,regular time connection,comparison-minor degree proof,polynomial extension and complete analytic/frame cover remain written inputs. No original checker main,radical search,producer package or fresh literature review.'}
    (P/'gaussian_projection_at_infinity_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
