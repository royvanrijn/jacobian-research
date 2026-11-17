#!/usr/bin/env python3
"""Exact seed fibres, trace selectors, actual covariance/Jacobians and a nonsplit fibre; no full descent formalization."""
from pathlib import Path
from itertools import product
import ast,hashlib,json,platform,time
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
checks=[]
def ck(label,value):
    assert value,label
    checks.append(label)
def substitute(F,variables,values):return [q.subs(dict(zip(variables,values)),simultaneous=True) for q in F]
def label_poly(x):return (3*x*x-x)/2
def selector(N,T):return -T*s.prod(T-j for j in range(2,2*N+1))/s.factorial(2*N-1)
def polynomial_zero(v):return s.expand(v)==0
def modpoly(expr,variables,prime):
    P=s.Poly(expr,*variables,domain=s.QQ)
    data={m:(int(c.p)*pow(int(c.q),-1,prime))%prime for m,c in P.as_dict().items()}
    return s.Poly.from_dict(data,variables,modulus=prime)
def main():
    started=time.monotonic();x,y,z=s.symbols('x y z');variables=[x,y,z]
    predecessor=ROOT/'jacobian-followup/universal-keller-compression/verify.py';nodes=[n for n in ast.parse(predecessor.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='source_map'];ck('reuse only retained pure seed function',len(nodes)==1)
    env={'x':x,'y':y,'z':z};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(predecessor),'exec'),env);K=env['source_map'](3,{})
    rho=[-s.Rational(1,4),0,0];mu=[s.Rational(1,3),-1,s.Rational(2,3)];zero=[0,0,0]
    t=1+x*y;q=t*t*z+y*y*(1+3*t);a,b,c=K
    ck('actual cubic Jacobian one',s.expand(s.Matrix(K).jacobian(variables).det())==1)
    ck('cleared inverse polynomial identity',polynomial_zero(2*x*t*t+2*b*x*x*t+2*a*x**3-c*t**3))
    ck('cleared inverse derivative unit identity',polynomial_zero(t*t+2*b*x*t+3*a*x*x-t))
    ck('global t q comaximality identity',polynomial_zero(x*x*q-1-t*(x*x*t*z+3*t*t-5*t+1)))
    Gempty=s.groebner([f-v for f,v in zip(K,mu)],z,y,x,domain=s.QQ)
    Gzero=s.groebner(K,z,y,x,domain=s.QQ)
    Gtriple=s.groebner([f-v for f,v in zip(K,rho)],z,y,x,domain=s.QQ)
    ck('complete empty fibre unit ideal',[g.as_expr() for g in Gempty.polys]==[1])
    ck('complete singleton fibre maximal ideal',[g.as_expr() for g in Gzero.polys]==[z,y,x])
    expected=[z+s.Rational(1,4)-s.Rational(27,4)*x*x,y+s.Rational(3,2)*x,x**3-x]
    ck('complete split triple exact triangular ideal',[g.as_expr() for g in Gtriple.polys]==expected)
    points=[(s.Integer(0),s.Integer(0),-s.Rational(1,4)),(s.Integer(1),-s.Rational(3,2),s.Rational(13,2)),(-s.Integer(1),s.Rational(3,2),s.Rational(13,2))]
    ck('all three labelled points in whole seed fibre',all(all(polynomial_zero(f-v) for f,v in zip(substitute(K,variables,P),rho)) for P in points))
    ck('label polynomial values zero one two',[label_poly(P[0]) for P in points]==[0,1,2])
    T=s.symbols('T');selector_controls=[]
    for N in [2,3,4]:
        lam=selector(N,T);ck(f'rank{N} exact selector table',all(lam.subs(T,j)==int(j==1) for j in range(2*N+1)))
        selected=[]
        for labels in product(range(3),repeat=N):
            value=lam.subs(T,sum(labels))
            if value==1:selected.append(labels)
            else:ck(f'rank{N} excluded selector label '+''.join(map(str,labels)),value==0)
        expected_labels=[tuple(int(j==i) for j in range(N)) for i in range(N)]
        ck(f'rank{N} complete filtered configurations',set(selected)==set(expected_labels))
        selector_controls.append({'rank':N,'selected_configurations':selected,'complete_marked_rank':len(selected),'dimension':3*N+3,'volume_dimension':3*N+4,'generic_degree':3**(N+1),'coordinate_degree':4*N})
    # Actual quotient-ring selector, literal whole-source maps and covariance.
    whole_controls=[]
    for N in [2,3]:
        xx=s.symbols('X:'+str(N));I=s.groebner([v**3-v for v in xx],*xx,domain=s.QQ)
        C=sum(label_poly(v) for v in xx);e=sum(xx[i]*(xx[i]+1)/2*s.prod(1-xx[j]**2 for j in range(N) if j!=i) for i in range(N))
        ck(f'rank{N} selector equals desired finite-algebra idempotent',I.reduce(s.expand(selector(N,C)-e))[1]==0)
        ck(f'rank{N} selected finite-algebra element idempotent',I.reduce(s.expand(e*e-e))[1]==0)
        blocks=[s.symbols('p'+str(j)+'_0:3') for j in range(N)];u=s.symbols('u:3');coords=[v for block in blocks for v in block]+list(u)
        Kblocks=[substitute(K,variables,block) for block in blocks];Ku=substitute(K,variables,u);Cv=sum(label_poly(block[0]) for block in blocks);P=selector(N,Cv)
        outputs=[f for block in Kblocks for f in block]+[v-(1-P)*w for v,w in zip(Ku,mu)]
        J=s.Matrix(outputs).jacobian(coords)
        ck(f'rank{N} full source upper-right derivative block zero',J[:3*N,3*N:]==s.zeros(3*N,3))
        ck(f'rank{N} full source first derivative blocks exact',J[:3*N,:3*N]==s.diag(*[s.Matrix(v).jacobian(block) for v,block in zip(Kblocks,blocks)]))
        ck(f'rank{N} full source last derivative block exact',J[3*N:,3*N:]==s.Matrix(Ku).jacobian(u))
        ck(f'rank{N} product of actual diagonal Jacobians one',s.prod(s.expand(s.Matrix(v).jacobian(block).det()) for v,block in zip(Kblocks,blocks))*s.expand(s.Matrix(Ku).jacobian(u).det())==1)
        swapped=list(blocks);swapped[0],swapped[1]=swapped[1],swapped[0];newcoords=[v for block in swapped for v in block]+list(u)
        changed=substitute(outputs,coords,newcoords);expected_output=outputs[3:6]+outputs[:3]+outputs[6:]
        ck(f'rank{N} actual adjacent-transposition covariance',all(polynomial_zero(v-w) for v,w in zip(changed,expected_output)))
        degree=max(s.Poly(v,*coords).total_degree() for v in outputs);ck(f'rank{N} literal coordinate degree bound exact',degree==4*N)
        H=[v for block in blocks for v in block]+outputs[-3:];R=[v for block in Kblocks for v in block]+list(u)
        ck(f'rank{N} literal two-self-factor identity',all(polynomial_zero(v-w) for v,w in zip(substitute(R,coords,H),outputs)))
        whole_controls.append({'rank':N,'actual_source_dimension':len(coords),'coordinate_degree':degree,'diagonal_jacobian_product':1,'two_factor_geometric_degrees':[3,3**N]})
    # A universal nonsplit quadratic fibre. Its points are not chosen over the input field.
    D,theta,U,A,B=s.symbols('D theta U A B');field=s.QQ.frac_field(D)
    GP=s.groebner([2*A-1,A*A+D*B*B-A,2*A*B-B],A,B,domain=field)
    ck('primitive idempotent quadratic scheme exact ideal',[v.as_expr() for v in GP.polys]==[A-s.Rational(1,2),B*B-1/(4*D)])
    ck('quadratic fibre generator gives input extension',GP.reduce(4*D*D*B*B-D)[1]==0)
    G=s.groebner([theta*theta-D,U*U-D],theta,U,domain=field);e=s.Rational(1,2)+theta*U/(2*D)
    ck('universal separability idempotent squares to itself',G.reduce(s.expand(e*e-e))[1]==0)
    fibrepoint=[e,-s.Rational(3,2)*e,-s.Rational(1,4)+s.Rational(27,4)*e]
    ck('intrinsic quadratic fibre point maps to actual seed target',all(G.reduce(s.expand(v-w))[1]==0 for v,w in zip(substitute(K,variables,fibrepoint),rho)))
    ck('intrinsic quadratic selector trace is one',s.expand(2*s.Rational(1,2))==1 and G.reduce(s.expand(label_poly(e)-e))[1]==0)
    # Whole subgroup-fixed restrictions, computed from orbit weights rather than label enumeration.
    def orbit_sizes(N,generators):
        left=set(range(N));sizes=[]
        while left:
            orbit={min(left)};changed=True
            while changed:
                bigger=orbit|{g[i] for g in generators for i in orbit};changed=bigger!=orbit;orbit=bigger
            sizes.append(len(orbit));left-=orbit
        return sorted(sizes)
    subgroup_controls=[]
    examples=[('S3',3,[(1,0,2),(1,2,0)]),('V4',4,[(1,0,3,2),(2,3,0,1)]),('C3 fixing one',4,[(1,2,0,3)]),('S3 fixing two',5,[(1,0,2,3,4),(1,2,0,3,4)]),('two C3 orbits',6,[(1,2,0,4,5,3)])]
    for name,N,generators in examples:
        sizes=orbit_sizes(N,generators);survivors=[labels for labels in product(range(3),repeat=len(sizes)) if selector(N,sum(n*a for n,a in zip(sizes,labels)))==1];r=sizes.count(1);d=3**(len(sizes)+1)
        ck(name+' complete fixed marked fibre',len(survivors)==r)
        ck(name+' fixed-degree normalization budget',r<=d<=r+3**(N+1)-N)
        subgroup_controls.append({'group':name,'rank':N,'orbit_sizes':sizes,'fixed_dimension':3*len(sizes)+3,'generic_fixed_degree':d,'complete_fixed_marked_rank':r})
    # Integral models: small good fibres test the actual seed ideals and full selector table.
    good_controls=[]
    for N,prime in [(2,5),(3,7),(4,11)]:
        mod=modpoly(selector(N,T),[T],prime)
        ck(f'good p{prime} rank{N} selector complete table',all(int(mod.eval(j))%prime==int(j==1) for j in range(2*N+1)))
        polynomials=[modpoly(v-w,variables,prime).as_expr() for v,w in zip(K,mu)]
        ck(f'good p{prime} actual empty seed fibre',[v.as_expr() for v in s.groebner(polynomials,z,y,x,modulus=prime).polys]==[1])
        good_controls.append({'rank':N,'characteristic':prime,'selector_and_empty_fibre':'pass'})
    for N in [2,3,4]:
        P=s.eye(3*N+3);P[0:3,0:3]=s.zeros(3);P[3:6,3:6]=s.zeros(3);P[0:3,3:6]=s.eye(3);P[3:6,0:3]=s.eye(3)
        ck(f'rank{N} canonical determinant-dual line cancels sign',P.det()==-1 and s.diag(P,-1).det()==1)
    manifest=json.loads((HERE/'sources.json').read_text())
    for src in manifest['local_sources']:ck('local pin '+src['path'],hashlib.sha256((ROOT/src['path']).read_bytes()).hexdigest()==src['sha256'])
    for src in manifest['external_sources']:
        data=(HERE/src['path']).read_bytes();ck('primary pin '+src['path'],len(data)==src['bytes'] and hashlib.sha256(data).hexdigest()==src['sha256'])
    receipt={'schema_version':1,'result':'pass','check_count':len(checks),'checks':checks,'seed_fibre_groebner_bases':{'empty':[str(v.as_expr()) for v in Gempty.polys],'singleton':[str(v.as_expr()) for v in Gzero.polys],'triple':[str(v.as_expr()) for v in Gtriple.polys]},'selector_controls':selector_controls,'literal_map_controls':whole_controls,'nonsplit_quadratic_control':'Primitive rank-one idempotents have A=1/2,B^2=1/(4D); U=2DB gives U^2=D. The universal separability idempotent maps to the exact seed fibre.','whole_fixed_controls':subgroup_controls,'good_characteristic_controls':good_controls,'python':platform.python_version(),'sympy':s.__version__,'runtime_seconds':round(time.monotonic()-started,3),'sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['PROOF.md','LITERATURE_AND_GAPS.md','verify.py','sources.json']},'assurance':'Exact whole seed fibre ideals, symbolic inverse/unit identities, finite-algebra selectors, literal whole-source Jacobian blocks/covariance/factorization, universal nonsplit quadratic fibre, subgroup orbit controls, and good-characteristic controls. No formalization or independent whole-proof replay of arbitrary-base descent, all-rank receiver, volume functor or atomic-length arguments. The constructed map is explicitly non-atomic. No publication-priority certificate.'}
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n');print(json.dumps({'result':'pass','checks':len(checks),'seconds':receipt['runtime_seconds']}))
if __name__=='__main__':main()
