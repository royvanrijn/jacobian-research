#!/usr/bin/env python3
"""Full-window inverse guards and complete outward-ball boundary replay.

The analytic covering/parameter argument is in the pinned proof. These controls
certify its numerical/algebraic inputs, not a period remainder or SFC threshold.
"""
from __future__ import annotations
if not __debug__:
    raise SystemExit('Refusing -O/-OO: all exact guards must execute.')
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
import argparse, importlib.util, json, subprocess, sys, tempfile
HERE=Path(__file__).resolve().parent
ROWS=[]
def sha(p): return sha256(p.read_bytes()).hexdigest()
def st(x):
    x=F(x)
    return str(x.numerator) if x.denominator==1 else f'{x.numerator}/{x.denominator}'
def ck(name,left,right,relation='<'):
    left,right=F(left),F(right)
    assert (left<right if relation=='<' else left==right),(name,st(left),relation,st(right))
    ROWS.append(dict(name=name,left=st(left),relation=relation,right=st(right)))
def add(a,b):
    out=dict(a)
    for k,v in b.items(): out[k]=out.get(k,F(0))+v
    return {k:v for k,v in out.items() if v}
def scale(a,c): return {k:v*c for k,v in a.items() if v*c}
def mul(a,b):
    out={}
    for i,u in a.items():
        for j,v in b.items():
            k=tuple(x+y for x,y in zip(i,j));out[k]=out.get(k,F(0))+u*v
    return {k:v for k,v in out.items() if v}
def deriv(a,i):
    out={}
    for k,v in a.items():
        if k[i]:
            j=list(k);j[i]-=1;out[tuple(j)]=v*k[i]
    return out
def var(n,j): return {tuple(int(i==j) for i in range(n)):F(1)}
def tiling(leaves):
    # Assign every closed leaf to its exact original root rectangle. For every
    # open x-slab, the sorted y-intervals must partition that root with no gap
    # and no overlap. Closed faces are retained by every adjacent leaf.
    ymin,ymax=F(18999,10000),F(21001,10000);dy=(ymax-ymin)/8
    groups={(k,j):[] for j in range(8) for k in range(64)}
    for row in leaves:
        al,ah,yl,yh=map(F,row[:4]);assert al<ah and yl<yh
        k=int(al*32);j=int((yl-ymin)/dy)
        assert (k,j) in groups
        assert F(k,32)<=al<ah<=F(k+1,32)
        assert ymin+j*dy<=yl<yh<=ymin+(j+1)*dy
        groups[k,j].append((al,ah,yl,yh))
    for (k,j),boxes in groups.items():
        assert boxes
        xs=sorted({q for row in boxes for q in row[:2]})
        assert xs[0]==F(k,32) and xs[-1]==F(k+1,32)
        for xl,xh in zip(xs,xs[1:]):
            mid=(xl+xh)/2
            intervals=sorted((yl,yh) for al,ah,yl,yh in boxes if al<mid<ah)
            cursor=ymin+j*dy
            for yl,yh in intervals:
                assert yl==cursor,('leaf gap or overlap',k,j,xl,xh,yl,cursor)
                cursor=yh
            assert cursor==ymin+(j+1)*dy
    return len(groups)
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--repo-root',type=Path,required=True)
    ap.add_argument('--output',type=Path,default=HERE/'gap_two_full_endpoint_inverse_verification.json')
    args=ap.parse_args();repo=args.repo_root.resolve()
    ledgerpath=HERE/'gap_two_full_endpoint_inverse_source_ledger.json'
    ledger=json.loads(ledgerpath.read_text());assert ledger.get('schema')==1
    mandatory={('local','GAP_TWO_FULL_ENDPOINT_UNIFORM_CUBIC_INVERSE.md'),('local',Path(__file__).name),('local','gap_two_full_endpoint_inverse_inputs.json'),('local','gap_two_full_endpoint_inverse_primary_audit.md'),('local','verify_hermite_endpoint_outer_boundary_ball.py'),('local','hermite_endpoint_outer_boundary_ball_verification.json'),('repo','research/factorial/image-rigidity/QUANTITATIVE_GAP_TWO_ENDPOINT_DATA.md'),('repo','research/factorial/image-rigidity/verify_quantitative_gap_two_endpoint_data.py'),('repo','research/factorial/image-rigidity/quantitative_gap_two_endpoint_data_source_ledger.json'),('repo','research/factorial/image-rigidity/HERMITE_CRITICAL_CUBIC_INVERSE_DISK.md'),('repo','research/factorial/image-rigidity/sources/dlmf-1-10-rouche.html')}
    actual={(r['root'],r['path']) for r in ledger['files']}
    assert actual==mandatory and len(ledger['files'])==len(mandatory)
    for row in ledger['files']:
        assert row['root'] in {'local','repo'}
        path=(HERE if row['root']=='local' else repo)/row['path']
        assert path.is_file() and sha(path)==row['sha256'],row['path']
    inputs=json.loads((HERE/'gap_two_full_endpoint_inverse_inputs.json').read_text())
    assert inputs['schema']==1
    d={k:F(v) for k,v in inputs['domain'].items()};b={k:F(v) for k,v in inputs['bounds'].items()}
    assert d==dict(delta=F(1,10**8),real_z=F(1,10**5),imag_z_min=F(18999,10**4),imag_z_max=F(21001,10**4),physical_radius=F(47,50))
    assert b==dict(boundary_hermite=F(9,100),boundary_actual=F(447,5000),zeta=F(1,8),inverse_radius=F(91,200),jacobian_radius=F(2,5),jacobian_bound=F(19),divided_amplitude_radius=F(3,8),divided_amplitude_derivative_bound=F(350000))
    assert inputs['coverage']==[{'declaration':'Theorem FullActualEndpointInverse','support':'Written complex-analytic/covering proof plus exact rational margins, universal identities, complete S3 enumeration, and closed-cell boundary replay.'}]
    delta=d['delta'];r=d['physical_radius'];x=d['real_z'];Z=F(1051,1000)
    den=d['imag_z_min']/2-r
    ck('upper-half-plane log margin',F(99,10000),den)
    ck('whole matched-data disk',(F(1001,10000))**2+x*x,F(101,1000)**2)
    ck('Hermite z derivative',(r+Z)/2+1/(2*den),52)
    ck('physical v modulus',r+Z,2)
    ck('actual delta phase derivative',F(8,6)/(1-delta)**2+1/(den-delta),102)
    vc=F(3,2);dc=F(3998,10000)
    ck('critical mean delta derivative',vc**3/(6*(1-delta*vc/2)**2)+vc**2/(4*dc),2)
    ck('actual phase delta allowance',102+2,110)
    actual_boundary=b['boundary_hermite']-52*x-110*delta
    ck('actual boundary',b['boundary_actual'],actual_boundary)
    ck('cubic inverse radius',b['inverse_radius']**3/3+b['zeta']*b['inverse_radius'],b['boundary_actual'])
    ck('stationary comparison',F(11,100)+4*delta+delta*delta,r*r)
    # 1+z^2/4=i*eta+eta^2/4; |eta|<.101 on this rectangle.
    ck('Hermite stationary constant',F(101,1000)+F(101,1000)**2/4,F(11,100))
    ck('critical cubic value squared',F(4,9)*b['zeta']**3,F(3,100)**2)
    ck('critical value below boundary',F(3,100),b['boundary_actual'])
    ck('zero-curve Rouche',F(49,100)/99999+F(1,4)*F(1,10**5),F(1,10**5))
    ck('safe punctured-domain point',F(1,10**5),F(1,20))
    ck('canonical stationary t inside inverse',b['zeta'],b['inverse_radius']**2)
    ck('closed inverse disk',F(9,20),b['inverse_radius'])
    ck('actual inverse Jacobian',r/F(1,20),b['jacobian_bound'])
    ck('divided amplitude contour encloses canonical criticals',b['zeta'],b['jacobian_radius']**2)
    ck('divided amplitude t margin',b['divided_amplitude_radius'],b['jacobian_radius'])
    G=b['jacobian_radius']*b['jacobian_bound']/((b['jacobian_radius']-b['divided_amplitude_radius'])**2*(b['jacobian_radius']**2-b['zeta']))
    ck('actual divided amplitude derivative',G,b['divided_amplitude_derivative_bound'])
    # Exact polynomial identities support the analytic covering argument.
    T,B=var(2,0),var(2,1);zeta=mul(B,B)
    phi=add(scale(mul(T,mul(T,T)),F(1,3)),scale(mul(zeta,T),-1))
    assert add(phi,scale(mul(B,mul(B,B)),F(2,3)))==scale(mul(mul(add(T,scale(B,-1)),add(T,scale(B,-1))),add(T,scale(B,2))),F(1,3))
    assert add(phi,scale(mul(B,mul(B,B)),-F(2,3)))==scale(mul(mul(add(T,B),add(T,B)),add(T,scale(B,-2))),F(1,3))
    V,Z,D=var(3,0),var(3,1),var(3,2);one={(0,0,0):F(1)}
    hv=add(V,scale(mul(D,mul(V,V)),F(1,4)));hp=add(one,scale(mul(D,V),F(1,2)))
    P=add(mul(V,add(Z,hv)),scale(mul(hp,hp),-1))
    assert P=={(2,0,0):1,(1,1,0):1,(0,0,0):-1,(3,0,1):F(1,4),(1,0,1):-1,(2,0,2):-F(1,4)}
    U,Z=var(2,0),var(2,1);one={(0,0):F(1)};zhalf=scale(Z,F(1,2))
    quad=add(scale(mul(add(U,scale(zhalf,-1)),add(U,scale(zhalf,-1))),-F(1,2)),scale(mul(Z,Z),F(1,4)))
    # Clearing 2*(u+z/2) from the moving-centre matched derivative.
    lhs=add(mul(scale(add(U,zhalf),2),deriv(quad,1)),one)
    rhs=add(mul(add(U,zhalf),add(U,zhalf)),one)
    assert lhs==rhs
    # Complete S3 alternatives, rather than one favorable transposition pair.
    trans=[(1,0,2),(2,1,0),(0,2,1)];identity=(0,1,2);cases=[]
    def comp(a,b):return tuple(a[b[j]] for j in range(3))
    for a in trans:
        for q in trans:
            product=comp(a,q);fixed=[j for j in range(3) if a[j]==j and q[j]==j]
            assert (product==identity)==(a==q)==(len(fixed)==1)
            cases.append(dict(first=a,second=q,boundary=product,common_fixed_sheets=fixed))
    # The original source and executed receipt remain immutable.
    source=HERE/'verify_hermite_endpoint_outer_boundary_ball.py'
    baseline=json.loads((HERE/'hermite_endpoint_outer_boundary_ball_verification.json').read_text())
    expected=inputs['boundary_evidence']
    assert baseline['source_sha256']==sha(source)==expected['source_sha256']
    for name in ['root_boxes','attempted_boxes','certified_leaf_boxes','maximum_bisection_depth','leaf_stream_sha256']:
        assert baseline[name]==expected[name]
    leaves=baseline['leaves'];assert len(leaves)==1300
    assert sha256(json.dumps(leaves,separators=(',',':')).encode()).hexdigest()==expected['leaf_stream_sha256']
    root_boxes=tiling(leaves);assert root_boxes==512
    spec=importlib.util.spec_from_file_location('retained_boundary',source);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    for row in leaves:
        assert abs(mod.phase_box(*map(F,row[:4])))>mod.arb('9/100')
    with tempfile.TemporaryDirectory(prefix='full_endpoint_boundary_') as temp:
        out=Path(temp)/'boundary.json'
        run=subprocess.run([sys.executable,str(source),'--output',str(out)],capture_output=True,text=True)
        assert run.returncode==0,(run.stdout,run.stderr)
        fresh=json.loads(out.read_text())
    assert {k:v for k,v in fresh.items() if k!='seconds'}=={k:v for k,v in baseline.items() if k!='seconds'}
    # Strong negative controls: false phase bound, omitted closed leaf,
    # a physical log-pole domain violation, insufficient inverse radius,
    # nontrivial monodromy, and a holomorphic map with a multiple zero.
    control=mod.phase_box(F(1,2),F(1,2),F(21,10),F(21,10))
    assert abs(control)<mod.arb('1/10')
    try:tiling(leaves[1:])
    except AssertionError:pass
    else:raise AssertionError('removed leaf was not detected')
    try:mod.phase_box(F(3,2),F(3,2),F(9,5),F(9,5))
    except AssertionError:pass
    else:raise AssertionError('log-domain violation was not detected')
    assert F(23,50)**3/3+b['zeta']*F(23,50)>b['boundary_actual']
    assert comp(trans[0],trans[1])!=identity
    assert not any(trans[0][j]==j and trans[1][j]==j for j in range(3))
    assert deriv({(3,):F(1)},0).get((0,),F(0))==0
    # In this full window unmatched +-2b may be OUTSIDE the inverse disk:
    # stationary matching must use F_v*v_t=t^2-zeta, never that shortcut.
    assert 4*b['zeta']>b['inverse_radius']**2
    result=dict(schema=1,status='PASS_FULL_WINDOW_INVERSE_SUPPORT',scope=inputs['scope'],rational_guards=ROWS,universal_identity_count=4,complete_monodromy_cases=cases,boundary_replay={k:v for k,v in fresh.items() if k not in {'seconds','leaves'}},independent_closed_tiling_root_boxes=root_boxes,independently_rechecked_leaf_enclosures=len(leaves),negative_controls=7,coverage=inputs['coverage'],source_ledger_sha256=sha(ledgerpath),pins=ledger['files'])
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(dict(status='PASS',rational_guards=len(ROWS),universal_identities=4,complete_S3_cases=len(cases),boundary_leaves=len(leaves),root_boxes=root_boxes,negative_controls=7)))
if __name__=='__main__':main()
