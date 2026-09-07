#!/usr/bin/env python3
"""Independent polynomial-gcd counts, PARI Sturm topology and capacity join."""
import argparse
from pathlib import Path
import warnings
import retrospective as r
import generic_selmer_capacity as prior
import root_curve_frobenius_capacity as finite
import root_curve_real_components as real

OUTPUT=r.OUT/'rank_jump_root_curve_capacity_verification_v1.json'
COMPARISON=r.OUT/'rank_jump_root_curve_capacity_comparison_v1.json'


def check_bindings(obj):
    for p,sha in obj.get('bindings',{}).items():assert r.digest((r.ROOT/p).read_bytes())==sha


def compute():
    from sage.all import QQ,ZZ,GF,PolynomialRing,pari,prod,matrix,identity_matrix,block_diagonal_matrix
    inputs=r.read(finite.INPUT);fc=r.read(finite.OUTPUT);rc=r.read(real.OUTPUT)
    for obj in (inputs,fc,rc,r.read(prior.VERIFICATION)):check_bindings(obj)
    R=PolynomialRing(QQ,'t');t=R.gen();families=[]
    def sturm(poly,a=None,b=None):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',DeprecationWarning)
            return int(pari(poly).polsturm() if a is None else pari(poly).polsturm(a,b))
    for f in inputs['families']:
        name=f['family'];fr=next(x for x in fc['rows'] if x['family']==name);rr=next(x for x in rc['rows'] if x['family']==name)
        assert fr['status']==rr['status']=='PASS';check_bindings(fr);check_bindings(rr)
        A=R(f['A']);B=R(f['B']);D=R(f['discriminant']);assert D==-16*(4*A**3+27*B**2)
        factors=[(R(q['coefficients_ascending']),q['multiplicity']) for q in f['squarefree_factors']]
        H=prod(q for q,e in factors);g=f['root_curve_genus'];allowance=sum(q.degree() for q,e in factors if e==2)
        assert allowance in (0,1) and g==10-allowance
        counts=[]
        # Recheck selection, but count roots by gcd(f,x^p-x), not enumeration.
        selected=[]
        for p in range(5,500,2):
            if not ZZ(p).is_prime():continue
            try:Ap=A.change_ring(GF(p));Bp=B.change_ring(GF(p));Dp=D.change_ring(GF(p));qs=[(q.change_ring(GF(p)),e) for q,e in factors]
            except (ValueError,ZeroDivisionError):continue
            if Dp.degree()!=24 or Ap.gcd(Dp)!=1:continue
            if any(q.degree()!=factors[i][0].degree() or q.gcd(q.derivative())!=1 for i,(q,e) in enumerate(qs)):continue
            if any(qs[i][0].gcd(qs[j][0])!=1 for i in range(len(qs)) for j in range(i)):continue
            nodes=[];okay=True
            for q,e in qs:
                if e==1:continue
                z=-q[0]/q[1];x=-3*Bp(z)/(2*Ap(z));delta=Ap.derivative()(z)**2-6*x*(Ap.derivative(2)(z)*x+Bp.derivative(2)(z))
                if not x or not delta:okay=False;break
                nodes.append({'base':int(z),'double_root':int(x),'tangent_discriminant':int(delta),
                    'normalization_correction':1 if delta.is_square() else -1})
            if not okay:continue
            selected.append(p);stored=fr['counts'][len(selected)-1];assert stored['prime']==p and stored['nodes']==nodes
            P=PolynomialRing(GF(p),'x');x=P.gen();ffcounts=[]
            for z in range(p):
                cubic=x**3+Ap(z)*x+Bp(z)
                ffcounts.append(int(cubic.gcd(pow(x,p,cubic)-x).degree()))
            infinity=x**3+Ap[8]*x+Bp[12];infc=int(infinity.gcd(pow(x,p,infinity)-x).degree())
            assert ffcounts==stored['finite_fibre_counts'] and infc==len(stored['infinity_roots'])
            total=sum(ffcounts)+infc+sum(x['normalization_correction'] for x in nodes)
            assert total==stored['normalized_root_curve_count'] and total%2==stored['odd_count']
            assert stored['raw_root_curve_count']==sum(ffcounts)+infc
            assert stored['cubic_root_histogram']==[ffcounts.count(k) for k in range(4)]
            counts.append({'prime':p,'normalized_count':total,'trace':p+1-total})
            if len(selected)==3:break
        assert len(selected)==3
        events=rr['events'];assert sturm(H)==len(events)
        for i,event in enumerate(events):
            a,b=map(QQ,event['interval']);assert a<b and (i==0 or QQ(events[i-1]['interval'][1])<a)
            assert H(a)*H(b)<0 and sturm(H,a,b)==1 and sturm(B,a,b)==0
            assert B(a).sign()==B(b).sign()==event['B_sign']
            at_node=any(q(a)*q(b)<0 for q,e in factors if e==2)
            assert at_node==(event['type']=='I2')
            if at_node:
                z=QQ(event['node_parameter']);x0=-3*B(z)/(2*A(z));assert a<z<b and D(z)==D.derivative()(z)==0
                delta=A.derivative()(z)**2-6*x0*(A.derivative(2)(z)*x0+B.derivative(2)(z))
                assert x0==QQ(event['node_double_root']) and delta==QQ(event['node_tangent_discriminant'])
                assert delta.sign()==event['node_tangent_discriminant_sign']
        arities=rr['interval_real_root_counts'];assert len(arities)==len(events)
        for i,s in enumerate(rr['interval_samples']):
            s=QQ(s);assert s>QQ(events[i]['interval'][1])
            if i+1<len(events):assert s<QQ(events[i+1]['interval'][0])
            assert arities[i]==(3 if D(s)>0 else 1)
        assert (3 if D.leading_coefficient()>0 else 1)==arities[-1]
        # Validate all boundary connections, then count components by DFS.
        expected=[];n=len(events)
        for i,e in enumerate(events):
            l=(i-1)%n;rt=i;lo,hi=(1,2) if e['B_sign']==1 else (0,1);single=3-lo-hi
            if e['type']=='I1':
                assert sorted((arities[l],arities[rt]))==[1,3]
                triple=l if arities[l]==3 else rt;one=rt if triple==l else l
                expected.extend([[[triple,lo],[triple,hi]],[[triple,single],[one,0]]])
            elif e['node_tangent_discriminant_sign']==1:
                assert arities[l]==arities[rt]==3
                expected.extend([[[l,k],[rt,hi if k==lo else lo if k==hi else k]] for k in range(3)])
            else:
                assert arities[l]==arities[rt]==1;expected.append([[l,0],[rt,0]])
        assert rr['edges']==expected
        adj={(i,k):[] for i,a in enumerate(arities) for k in range(a)}
        for a,b in expected:a=tuple(a);b=tuple(b);adj[a].append(b);adj[b].append(a)
        unseen=set(adj);groups=[]
        while unseen:
            stack=[min(unseen)];group=set()
            while stack:
                v=stack.pop()
                if v in group:continue
                group.add(v);stack.extend(adj[v])
            unseen-=group;groups.append(sorted(group))
        assert sorted(groups)==sorted([list(map(tuple,x)) for x in rr['components']])
        s=len(groups);assert s==rr['real_components'] and rr['real_jacobian_two_torsion_dimension']==g+s-1
        families.append({'family':name,'genus':g,'node_allowance':int(allowance),'finite_counts':counts,
            'real_discriminant_places':len(events),'real_components':s,
            'real_Jacobian_two_torsion_dimension':g+s-1,'arithmetic_global_pool_dimension_upper_bound':g+s-1+int(allowance),
            'odd_frobenius_count_found':any(x['normalized_count']%2 for x in counts)})
    # Same unipotent characteristic polynomial, different fixed-space sizes.
    examples=[]
    for g in (9,10):
        for k in range(4):
            blocks=[matrix(GF(2),[[1,1 if i<k else 0],[0,1]]) for i in range(g)]
            M=block_diagonal_matrix(blocks);J=block_diagonal_matrix([matrix(GF(2),[[0,1],[1,0]])]*g)
            assert M.transpose()*J*M==J
            poly=M.charpoly();assert poly==(poly.parent().gen()+1)**(2*g)
            fixed=2*g-(M-identity_matrix(GF(2),2*g)).rank();assert fixed==2*g-k
            examples.append({'genus':g,'transvections':k,'fixed_dimension':int(fixed)})
    paths=[Path(__file__),finite.INPUT,finite.OUTPUT,real.OUTPUT,finite.PROTOCOL,real.PROTOCOL,
        Path(finite.__file__),Path(real.__file__),prior.VERIFICATION]
    return {'schema':'rank-jump.root-curve-capacity-verification.v1','status':'PASS','families':families,
        'unipotent_ambiguity_examples':examples,'bindings':prior.bind(paths),
        'boundary':'Independent finite-field gcd counts, exact PARI Sturm root isolation checks, normalization/topology graph replay and real-capacity arithmetic. Picard-descent and Frobenius/real-torsion implications are proved in the note. No exceptional points used.'}


def compare(result):
    byfamily={x['family']:x for x in result['families']};old=r.read(prior.COMPARISON);rows=[]
    for x in old['rows']:
        d=byfamily[x['family']]['arithmetic_global_pool_dimension_upper_bound'];m=x['generic_mod_two_dimension'];R=x['retained_rank_lower_bound']
        rows.append({'token':x['token'],'id':x['id'],'family':x['family'],'generic_dimension':m,
            'retained_rank_lower_bound':R,'arithmetic_global_pool_dimension_upper_bound':d,
            'extra_global_pool_capacity_upper_bound':d-m,
            'rational_fibre_dimensions_outside_pool_lower_bound':max(0,R-d),
            'additional_quotient_CT':'UNKNOWN'})
    return {'schema':'rank-jump.root-curve-capacity-comparison.v1','rows':rows,
        'bindings':prior.bind([Path(__file__),OUTPUT,prior.COMPARISON]),
        'boundary':'Labels joined after equation-only counting. These are capacity bounds and retrospective necessities, not a prospective rank predictor or an additional class basis.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();result=compute()
    if args.mode=='build':r.write_new(OUTPUT,result);r.write_new(COMPARISON,compare(result))
    else:assert result==r.read(OUTPUT);assert compare(result)==r.read(COMPARISON)
    print('PASS',[(x['family'],x['real_components'],x['arithmetic_global_pool_dimension_upper_bound']) for x in result['families']],flush=True)
