#!/usr/bin/env sage-python
"""Replay emitted carrier identities, the strict dictionary, and local bounds."""
import importlib.machinery,importlib.util,json,hashlib
from collections import defaultdict
from pathlib import Path
from sage.all import QQ,ZZ,GF,EllipticCurve,matrix,vector,pari
path=Path(__file__).with_name('rank_triangle_geometry.sage')
loader=importlib.machinery.SourceFileLoader('triangle_geometry',str(path));spec=importlib.util.spec_from_loader(loader.name,loader);g=importlib.util.module_from_spec(spec);loader.exec_module(g)
R,F,OUT,ART,ROOT=g.R,g.F,g.OUT,g.ART,g.ROOT
def main():
    geometry=g.read(OUT/'geometry.json');native=g.read(OUT/'native-carriers.json');strict=g.read(OUT/'strict.json')
    inputs={d['id']:d for d in g.read(ART/'rank_accessibility_subsets_v1/inputs.json.gz')['curves']}
    branches=defaultdict(list);count=0;rejected=0
    for row in geometry['rows']:
        cid=row['fibre'];j=row['target'];data=inputs[cid];A,B,sections,t0=g.parent(cid)
        source=EllipticCurve(QQ,list(map(QQ,data['model'])));target=EllipticCurve(QQ,[g.at(A,t0),g.at(B,t0)])
        assert row['source_point']==data['basis_points'][16+j]
        iso=next(i for i in source.isomorphisms(target) if list(map(str,i.tuple()))==row['fibre_isomorphism'])
        P=iso(source(list(map(QQ,row['source_point']))));assert list(map(str,P.xy()))==row['parent_point']
        covers=list(row['covers'])
        assert len(covers)==(36 if cid=='302' else 35)
        assert {(c['tag']['generic_index'],c['tag']['sign']) for c in covers if c['tag']['family']=='constant_slope'}=={(i,s) for i in range(1,18) for s in [-1,1]}
        if cid=='11952':covers += [c for c in native['rows'] if c['target']==j]
        for c in covers:
            q=R(c['branch_monic_polynomial']);constant=QQ(c['constant_twist_representative']);D=constant*q
            assert q.is_monic() and q.degree()>0 and q.gcd(q.derivative())==1
            assert c['normalization_genus']==(q.degree()-1)//2 and c['branch_at_infinity']==bool(q.degree()%2)
            assert g.dec(c['raw_radical'])==D*g.dec(c['square_multiplier'])**2
            x0,x1,y0,y1=[g.dec(c['maps'][k]) for k in ['x0','x1','y0','y1']]
            assert x1 or y1  # the map retains the quadratic function field
            assert y0*y0+y1*y1*D==x0**3+3*x0*x1*x1*D+A*x0+B
            assert 2*y0*y1==3*x0*x0*x1+x1**3*D+A*x1
            base,v=map(QQ,c['lift']);assert base==t0 and v*v==D(base)
            assert g.at(x0,base)+g.at(x1,base)*v==P[0]
            assert g.at(y0,base)+g.at(y1,base)*v==P[1]
            branches[(cid,tuple(c['branch_monic_polynomial']))].append((constant,j))
            # A perturbed ordinate constant must fail a coefficient identity.
            bad=y0+1
            assert not (bad*bad+y1*y1*D==x0**3+3*x0*x1*x1*D+A*x0+B and 2*bad*y1==3*x0*x0*x1+x1**3*D+A*x1)
            rejected+=1;count+=1
        assert min(c['normalization_genus'] for c in covers)==1
        print('CARRIERS',cid,j,'PASS',flush=True)
    shared=[]
    for key,values in branches.items():
        for i,(c,j) in enumerate(values):
            for d,k in values[i+1:]:
                if j!=k and (c/d).is_square():shared.append([key[0],j,k])
    assert not shared and count==792
    B=matrix(GF(2),inputs['302']['basis_words_in_D']);C=matrix(GF(2),strict['strict_basis_in_M17_E'])
    J=matrix(GF(2),strict['strict_public_basis']);L=matrix(GF(2),strict['local_signature_matrix'])
    assert C*B==J and J*L==0 and J.rank()==10 and L.rank()==21
    assert matrix(GF(2),strict['good_character_matrix']).rank()==31
    assert (B[:17,:]*L).rank()==17
    checks=matrix(GF(2),strict['quotient_local_check_matrix']);assert checks.rank()==4 and checks*C[:,17:].transpose()==0
    assert J.row_space()==L.left_kernel()
    incidence=[]
    for row in geometry['rows']:
        if row['fibre']=='302':
            e=vector(GF(2),[int(i==row['target']-1) for i in range(14)])
            q=checks*e
            incidence.append({'target':row['target'],'quotient_local_class':list(map(int,q)),
              'generic_adjusted_strict_dimension_of_this_direction':int(q.is_zero()),
              'best_carrier_genus':1,'other_chosen_targets_sharing_any_tested_cover':[]})
    panel=g.read(OUT/'panel.json');arith={r['id']:r for r in g.read(OUT/'arithmetic.json')['results']}
    arith['302']=g.read(OUT/'arithmetic/302-cached-support/result.json');complete=0
    for row in panel['rows']:
        rec=arith[row['id']]
        if rec['status']!='PASS_EQUATION_LOCAL_ARITHMETIC':continue
        print('LOCAL_BEGIN',row['id'],flush=True)
        folder=OUT/'arithmetic'/('302-cached-support' if row['id']=='302' else row['id'])
        log=(folder/'field.log').read_text();single=dict(line.split('|',1) for line in log.splitlines() if '|' in line)
        expected=rec.get('field_log_sha256',rec.get('log_sha256'));assert g.sha(folder/'field.log')==expected
        factors=[json.loads(line.split('|')[1]) for line in log.splitlines() if line.startswith('FACTOR|')]
        E=EllipticCurve(QQ,list(map(QQ,rec['minimal_model'])));original=EllipticCurve(QQ,list(map(QQ,row['model'])))
        u=QQ(json.loads(single['TRANSFORM'])[0]);assert original.c4()==u**4*E.c4() and original.c6()==u**6*E.c6()
        assert all(ZZ(p).is_prime(proof=True) for p,e in factors)
        product=ZZ(1)
        for p,e in factors:product*=ZZ(p)**e
        assert product==abs(E.discriminant())
        f=R(rec['cubic']);assert f==R([16*E.b6(),8*E.b4(),E.b2(),1])
        assert f.discriminant()==256*E.discriminant()
        primes=[p for p,e in factors]
        # Supply the already proved discriminant factors to certification too;
        # otherwise nfcertify can unnecessarily refactor the large discriminant.
        pari.addprimes(primes)
        nf=pari.nfinit([pari(f),[2]+primes]);assert pari.nfcertify(nf)==[]
        assert str(nf.disc())==rec['field_discriminant']
        assert f.discriminant()==ZZ(rec['field_index'])**2*ZZ(rec['field_discriminant'])
        offset=2 if E.discriminant()>0 else 1;conductor=ZZ(1)
        for p,vd,ce,kod,number,c in rec['local']:
            assert E.discriminant().valuation(p)==vd and len(pari.idealprimedec(nf,p))==number
            local=pari.elllocalred(pari(E),p);assert [int(local[0]),int(local[1])]==[ce,kod]
            assert c==(int(vd%2==0) if ce==1 else number-1);offset+=c;conductor*=ZZ(p)**ce
        assert offset==rec['bk_offset'] and conductor==ZZ(rec['conductor'])
        assert int(pari.ellrootno(pari(E)))==rec['root_number']
        assert rec['g_upper'] is None and rec['rank_upper'] is None and rec['grh_g'] is None
        complete+=1
        print('LOCAL_PASS',row['id'],flush=True)
    assert complete==11
    sources=[OUT/s for s in ['protocol.json','panel.json','population.json','geometry.json','native-carriers.json','strict.json','arithmetic.json','arithmetic/302-cached-support/result.json']]
    sources += [Path(__file__),path]
    result={'status':'PASS_BOUNDED_TRIANGLE_REPLAY','bindings':{str(p):g.sha(p) for p in sources},
       'carrier_identities_checked':count,'perturbed_carriers_rejected':rejected,'targets_with_genus_one_carriers':22,
       'shared_cover_target_pairs':shared,'complete_local_arithmetic_rows':complete,'rank_upper_bounds_obtained':0,
       'strict_geometry_incidence':incidence,
       'boundary':'Emitted identity and linear-algebra replay; local number-field arithmetic still uses PARI. No independent full Selmer/class computation, exhaustive genus0 search, or basis-invariant grouping claim.'}
    g.save(OUT/'verified.json',result);print('PASS',count,'carrier identities;',complete,'local arithmetic rows; strict dictionary; no shared tested covers')
if __name__=='__main__':main()
