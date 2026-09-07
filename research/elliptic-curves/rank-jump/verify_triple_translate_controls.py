#!/usr/bin/env python3
"""Independent matched intersection, local obstruction and selector replay."""
import argparse
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,matrix,vector,prime_range
import retrospective as r
from verify_native_intersection_solubility import add,short_vectors
import triple_translate_controls as c

HERE=Path(__file__).resolve().parent
OUTPUT=r.OUT/'rank_jump_triple_translate_controls_verification_v1.json'


def verify():
    selection=r.read(c.SELECTION);summary=r.read(c.OUTPUT);base=r.read(c.BASE)
    for data in (selection,summary):
        for path,sha in data['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    G=matrix(ZZ,base['gram']);w=vector(ZZ,selection['trace_sum']);V=matrix(ZZ,selection['reduced_coset_basis'])
    assert abs(V.det())==2**16
    for row in V.rows():assert all(x%2==0 for x in row) or all((row[i]-w[i])%2==0 for i in range(17))
    vs,nodes=short_vectors(V*G*V.transpose(),10);actual={tuple(vector(ZZ,v)*V) for v in vs}
    assert actual=={tuple(v['vector']) for v in selection['signed_vectors']}
    def key(v):return sum(abs(int(x)) for x in v),max(abs(int(x)) for x in v),tuple(map(int,v))
    eligible=set()
    for v in actual:
        rr=vector(ZZ,v);s=(rr+w)/2
        if rr*G*rr==10 and s*G*s==10:eligible.add(tuple(map(int,s)))
    canon={min([s]+([tuple(w-vector(ZZ,s))] if tuple(w-vector(ZZ,s)) in eligible else []),key=key) for s in eligible}
    chosen=sorted(canon,key=key)[:3]
    assert list(map(list,sorted(eligible,key=key)))==selection['eligible_words']
    assert list(map(list,chosen))==selection['selected_words'] and len(canon)==selection['conjugacy_class_count']
    R=PolynomialRing(QQ,'t');K=R.fraction_field();A=R(base['A']);B=R(base['B'])
    qs=[R(x['residual_chord']['q_coefficients']) for x in base['covers']]
    tau=[vector(ZZ,x['published_basis_w']) for x in base['covers']]
    assert w==tau[0]-tau[1]+tau[2]
    delta=-16*(4*A**3+27*B**2)
    assert delta.is_squarefree() and all(q.degree()==2 and q.is_squarefree() and q.gcd(delta)==1 for q in qs)
    assert all(qs[i].gcd(qs[j])==1 for i in range(3) for j in range(i))
    def dec(x):return K(R(x['numerator']))/R(x['denominator'])
    basis=[]
    for x in base['sections']:
        xx=K(R(x['x_coefficients_low_to_high']))
        if 'y_coefficients_low_to_high' in x:yy=K(R(x['y_coefficients_low_to_high']))
        else:
            ch=x['chord'];ref=basis[ch['reference_basis_index']]
            yy=ref[1]+R(ch['slope_coefficients_low_to_high'])*(xx-ref[0])
        assert yy*yy==xx**3+A*xx+B;basis.append((xx,yy))
    rows=[];old=r.read(c.OLD)
    for i,word in enumerate(chosen):
        data=r.read(c.resultpath(i));inp=r.read(c.casepath(i))
        for source in (data,inp):
            for path,sha in source['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
        s=vector(ZZ,word);z=2*s+tau[1]-tau[2];rv=z-tau[0]
        assert s*G*s==10 and rv*G*rv==10
        aa=z*G*z/2+2;trace=2*z
        assert 8*aa-trace*G*trace==16 and 2*aa+8-tau[0]*G*trace==12
        terms=[(int(n/abs(n))*vector(ZZ,[int(j==k) for j in range(17)]),(basis[k][0],int(n/abs(n))*basis[k][1])) for k,n in enumerate(s) if n for _ in range(abs(int(n)))]
        acc=vector(ZZ,[0]*17);S=None
        while terms:
            j=min(range(len(terms)),key=lambda j:((acc+terms[j][0])*G*(acc+terms[j][0]),j));v,P=terms.pop(j)
            acc+=v;S=add(A,S,P)
        assert S==(dec(data['generic_translate']['x']),dec(data['generic_translate']['y']))
        f=R(data['intersection_polynomial']);assert f.degree()==12 and f.is_squarefree()
        factors=[R(x['coefficients']) for x in data['factorization']];product=R(1)
        for ff,row in zip(factors,data['factorization'],strict=True):assert row['multiplicity']==1;product*=ff
        assert product==f
        F=R.quotient(f,'tt')
        def red(x):x=K(x);assert x.denominator().gcd(f)==1;return F(x.numerator())/F(x.denominator())
        ps=[]
        for root,cover,point,q in zip(data['rational_roots'],inp['covers'],data['point_maps'],qs,strict=True):
            root=dec(root);assert red(root)**2==F(q)
            x0,x1,y0,y1=[R(cover['lifted_section'][k+'_coefficients']) for k in ('x0','x1','y0','y1')]
            xx=x0+root*x1;yy=y0+root*y1
            assert (xx,yy)==(dec(point['x']),dec(point['y']))
            xx,yy=red(xx),red(yy);assert yy*yy==xx**3+F(A)*xx+F(B);ps.append((xx,yy))
        assert add(F(A),add(F(A),ps[0],(ps[1][0],-ps[1][1])),ps[2])==tuple(map(red,S))
        assert R(data['excluded_polynomial']).gcd(f)==1
        # Fixed bounded prime list: reduction certificates for the finite scheme.
        integer=R(f.denominator()*f);local=[];irreducible_witness=None;nonlinear_witness=None
        nonlin=next(ff for ff in factors if ff.degree()>1);inl=R(nonlin.denominator()*nonlin)
        possible=set(range(1,int(nonlin.degree())));nonlinear_patterns=[];degree_cuts=[]
        for p in prime_range(2,132):
            Rp=PolynomialRing(GF(p),'t');fp=Rp(integer);np=Rp(inl)
            if fp.degree()!=12 or not fp.is_squarefree():continue
            degrees=sorted(int(ff.degree()) for ff,e in fp.factor() for _ in range(e))
            roots=[int(x) for x in fp.roots(multiplicities=False)]
            local.append({'prime':int(p),'factor_degrees':degrees,'rational_roots':roots})
            if degrees==[12] and irreducible_witness is None:irreducible_witness=int(p)
            if np.degree()==nonlin.degree() and np.is_squarefree():
                ds=sorted(int(ff.degree()) for ff,e in np.factor() for _ in range(e));nonlinear_patterns.append({'prime':int(p),'degrees':ds})
                sums={0}
                for d in ds:sums|={x+d for x in list(sums)}
                after=possible&sums
                if after!=possible:degree_cuts.append({'prime':int(p),'factor_degrees':ds,'remaining_proper_factor_degrees':sorted(after)})
                possible=after
                if ds==[int(nonlin.degree())] and nonlinear_witness is None:nonlinear_witness=int(p)
        no_root=next((x for x in local if not x['rational_roots']),None)
        assert not possible, 'No bounded modular irreducibility certificate'
        n=int(nonlin.degree())
        transposition=next(x for x in nonlinear_patterns if x['degrees'].count(2)==1 and all(d==2 or d%2 for d in x['degrees']))
        if ZZ(n).is_prime():primitive={'reason':'transitive action of prime degree','degree':n}
        else:
            cycle=next((x,d) for x in nonlinear_patterns for d in x['degrees'] if ZZ(d).is_prime() and d>n/2)
            primitive={'reason':'prime cycle larger than half the degree excludes every nontrivial block system',
                       'prime':cycle[0]['prime'],'factor_degrees':cycle[0]['degrees'],'cycle_length':cycle[1]}
        iscal=selection['matches_oracle_conjugacy_class'][i]
        if iscal:
            assert f==R(old['intersection_polynomial'])
            for x,y in zip(data['rational_roots'],old['rational_roots'],strict=True):assert red(dec(x)+dec(y))==0
        rows.append({'index':i,'factor_degrees':[int(ff.degree()) for ff in factors],
                     'proper_geometric_intersections_verified':12,'nonlinear_factor_irreducible_mod_prime':nonlinear_witness,
                     'irreducibility_degree_cuts':degree_cuts,
                     'nonlinear_Galois_group':f'S{n}','primitivity_certificate':primitive,
                     'transposition_certificate':transposition,
                     'polynomial_irreducible_mod_prime':irreducible_witness,
                     'first_certified_local_obstruction':no_root,'modular_factor_patterns':local,
                     'oracle_conjugacy_calibration':iscal})
    assert [row['factor_degrees'] for row in rows]==[[12],[1,11],[12]]
    assert rows[0]['first_certified_local_obstruction'] and rows[2]['first_certified_local_obstruction']
    return {'schema':'rank-jump.triple-translate-controls-verification.v1','status':'PASS','rows':rows,
        'exact_coset_enumeration_nodes':nodes,'signed_short_vectors_verified':len(actual),'eligible_translate_classes':len(canon),
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (c.PROTOCOL,c.SELECTION,c.OUTPUT,c.BASE,c.OLD,Path(__file__),HERE/'retrospective.py',HERE/'verify_native_intersection_solubility.py')},
        'boundary':'Independent rational group law and exact finite-scheme constructions; good-prime no-root reductions are local obstructions for these relation schemes, not obstructions for the whole carrier or original curve.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();data=verify()
    if a.mode=='build':r.write_new(OUTPUT,data)
    else:assert data==r.read(OUTPUT)
    for row in data['rows']:print(row['index'],row['factor_degrees'],'irreducibility prime',row['nonlinear_factor_irreducible_mod_prime'],'local obstruction',row['first_certified_local_obstruction'])
