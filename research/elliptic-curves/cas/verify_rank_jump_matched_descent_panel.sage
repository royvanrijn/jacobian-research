#!/usr/bin/env sage-python
"""Independent full local-Hilbert-basis replay of the seven completed fibres.

Replaces the producer's local ideal-log squareclass coordinates by a complete
Hilbert-dual basis. Rechecks support, points, native17 prefix and good-prime
independence. Additional derivative tests may close the full Selmer boundary.
"""
import argparse,hashlib,json,sys
from pathlib import Path
from sage.all import AA,QQ,ZZ,GF,EllipticCurve,PolynomialRing,matrix,pari
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
PROTOCOL=ART/'rank_jump_matched_descent_panel_protocol_v1.json';INPUT=ART/'rank_jump_matched_descent_panel_v1.json'
SUPPLEMENT=ART/'rank_jump_matched_descent_dyadic_completion_v1.json'
OUTPUT=ART/'rank_jump_matched_descent_panel_verification_v1.json'
sys.path.insert(0,str(Path(__file__).parent))
from r17_60_arithmetic import native_check
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def compute():
    policy=read(PROTOCOL);data=read(INPUT);extra=read(SUPPLEMENT)
    assert data['protocol_sha256']==extra['protocol_sha256']==sha(PROTOCOL)
    for path,digest in policy['bindings'].items():assert sha(ROOT/path)==digest,path
    for path,digest in extra['bindings'].items():assert sha(ROOT/path)==digest,path
    rows={r['id']:r for r in policy['rows']};checked=[]
    for initial in data['cases']:
        d=extra if initial['id']==extra['id']else initial;r=rows[d['id']]
        assert d['rank_lower_bound']==r['rank_lower_bound']
        native_check({k:r[k]for k in ['curve','points','proof','rank_lower_bound']},r['roster_row'])
        R=PolynomialRing(QQ,'z');f=R([QQ(r['curve'][4]),QQ(r['curve'][3]),0,1])
        assert list(map(str,f.list()))==d['cubic_ascending']
        S=sorted({2,*map(int,r['retained_prime_hints'])});remain=abs(ZZ(f.discriminant()));factors=[]
        for p in S:
            assert ZZ(p).is_prime(proof=True);e=remain.valuation(p)
            if e:remain//=ZZ(p)**e;factors.append([str(p),int(e)])
        assert str(remain)==d['remaining_cofactor']and factors==d['disc_factors']
        if remain!=1:
            assert d['status']=='UNKNOWN_INCOMPLETE_SUPPORT'
            checked.append({'id':r['id'],'pair':r['pair'],'status':d['status'],'remaining_cofactor':str(remain),
                            'rank_lower_bound':r['rank_lower_bound']})
            continue
        assert d['status']=='PASS_EXACT_POINT_SUBSPACE_ANATOMY'and f.is_irreducible()
        nf=pari.nfinit([pari(f),S]);assert str(nf.disc())==d['field_discriminant']
        assert list(map(str,nf.nf_get_zk()))==d['maximal_order_basis']
        E=EllipticCurve(QQ,r['curve']);n=r['rank_lower_bound'];gs=[];polys=[]
        for point in r['points']:
            x,y=E(point).xy();den=ZZ(x.denominator()).sqrt();assert den in ZZ and den*den==x.denominator()
            a,b=ZZ(x*den**2),ZZ(y*den**3);g=R([a,-den**2]);alpha=pari.Mod(pari(g),pari(f))
            assert pari.nfeltnorm(nf,alpha)==b*b and a.gcd(den)==b.gcd(den)==1
            gs.append(alpha);polys.append(g)
        independence=[[]for _ in gs]
        for block in d['good_character_blocks']:
            p=block['prime'];assert ZZ(p).is_prime(proof=True)and 3<=p<=2000 and f.discriminant()%p
            assert sorted(block['roots'])==sorted(map(int,f.change_ring(GF(p)).roots(multiplicities=False)))
            for row,g in zip(independence,polys):
                vals=[g.change_ring(GF(p))(t)for t in block['roots']];assert all(vals)
                row.extend(int(not v.is_square())for v in vals)
        assert matrix(GF(2),independence).rank()==n
        # Append the derivative global class to every local coordinate block.
        derivative=pari.Mod(pari(-f.discriminant()*f.derivative()),pari(f));values=gs+[derivative]
        all_local=[[]for _ in values];val_constraints=[];unram_constraints=[];local_dims=[];derivative_out=[]
        for record in d['local']:
            p=record['place'];assert p in S;sig=[[]for _ in values]
            primes=pari.idealprimedec(nf,p)
            expected_dim=len(primes)-1+int(p==2)
            assert expected_dim==record['point_image_dimension']
            for P in primes:
                valuations=[int(pari.idealval(nf,g,P))%2 for g in gs]
                val_constraints.append(valuations);unram_constraints.append(valuations)
                pi=pari.nfbasistoalg(nf,pari.idealappr(nf,P));assert pari.idealval(nf,pi,P)==1
                exponent=2*int(P[2])+1 if p==2 else 1
                bid=pari.idealstar(nf,pari.idealpow(nf,P,exponent),2)
                units=[]
                for order,u in zip(bid.bid_get_cyc(),bid.bid_get_gen()):
                    if int(order)%2:continue
                    u=pari.nfbasistoalg(nf,u)if u.type()=='t_COL'else u
                    assert pari.idealval(nf,u,P)==0;units.append(u)
                basis=[pi]+units
                hilbert=matrix(GF(2),[[int(pari.nfhilbert(nf,x,y,P)==-1)for y in basis]for x in basis])
                expected_square_dim=(int(P[2])*int(P[3])+2)if p==2 else 2
                assert len(basis)==hilbert.rank()==expected_square_dim
                for row,g in zip(sig,values):row.extend(int(pari.nfhilbert(nf,g,y,P)==-1)for y in basis)
                if p==2:
                    for u in units:unram_constraints.append([int(pari.nfhilbert(nf,g,u,P)==-1)for g in gs])
            for row,part in zip(all_local,sig):row.extend(part)
            H=matrix(GF(2),sig);rank=int(H[:n,:].rank());grank=int(H[:17,:].rank())
            assert rank==record['known_image_rank']and grank==record['generic_image_rank']
            # Equality of local point span and full local Kummer dimension is
            # required before derivative nonmembership is used globally.
            if rank==expected_dim and H.rank()>rank:derivative_out.append(p)
            local_dims.append(expected_dim)
        assert len(d['local'])==len(S)and set(t['place']for t in d['local'])==set(S)
        roots=f.roots(AA,multiplicities=False);real=[[int(g(t)<0)for t in roots]for g in polys]
        assert real==d['real_signatures']
        derivative_poly=-f.discriminant()*f.derivative();derivative_real=[int(derivative_poly(t)<0)for t in roots]
        for row,bits in zip(all_local,real+[derivative_real]):row.extend(bits)
        unram_constraints.extend([list(v)for v in zip(*real)])
        realrank=int(matrix(GF(2),real).rank());realexpected=int(len(roots)==3)
        if realrank==realexpected and matrix(GF(2),real+[derivative_real]).rank()>realrank:derivative_out.append('infinity')
        H=matrix(GF(2),all_local[:n]);V=matrix(GF(2),val_constraints);U=matrix(GF(2),unram_constraints)
        assert val_constraints==d['valuation_matrix']
        assert U.right_kernel()==matrix(GF(2),d['unramified_constraint_matrix']).right_kernel()
        assert [list(map(int,w))for w in H.left_kernel().basis()]==d['strict_words']
        assert [list(map(int,w))for w in H[:17,:].left_kernel().basis()]==d['generic_strict_words']
        expected={'generic':17,'known':n,'generic_local':int(H[:17,:].rank()),'known_local':int(H.rank()),
                  'generic_strict':17-int(H[:17,:].rank()),'known_strict':n-int(H.rank()),
                  'generic_even':17-int(V[:,:17].rank()),'known_even':n-int(V.rank()),
                  'generic_unramified':17-int(U[:,:17].rank()),'known_unramified':n-int(U.rank()),
                  'local_product':sum(local_dims)+realexpected}
        unit_bound=(3+len(roots))//2-1
        expected['relative_ideal_lower_bound']=max(0,expected['known_even']-expected['generic_even']-unit_bound)
        assert expected==d['dimensions']and unit_bound==d['unit_dimension_upper_bound']
        boundary_closed=bool(derivative_out)and expected['local_product']==expected['known_local']+1
        checked.append({'id':r['id'],'pair':r['pair'],'status':'PASS_INDEPENDENT_LOCAL_HILBERT_REPLAY',
                        'rank_lower_bound':n,'dimensions':expected,'signature':d['signature'],
                        'strict_relative_dimension':expected['known_strict']-expected['generic_strict'],
                        'ordinary_unramified_relative_dimension':expected['known_unramified']-expected['generic_unramified'],
                        'derivative_nonmembership_places':derivative_out,'full_Selmer_boundary_closed':boundary_closed,
                        'full_Selmer_dimension':str(expected['known_local'])+'+c_S'if boundary_closed else 'UNKNOWN',
                        'c_S_lower_bound':expected['known_strict']})
        print('VERIFIED',r['id'],expected['known_strict'],expected['known_unramified'],'BOUNDARY',boundary_closed,flush=True)
    return {'schema':'rank-jump.matched-descent-panel-verification.v1','status':'PASS',
            'bindings':{str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),PROTOCOL,INPUT,SUPPLEMENT]},
            'cases':checked,'completed_cases':sum(x['status']=='PASS_INDEPENDENT_LOCAL_HILBERT_REPLAY'for x in checked),
            'boundary':'Seven exact known-point subspace calculations and three retained support gaps. These are not complete class groups, exact elliptic ranks, or prospective predictor validation.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args();result=compute()
    if args.check:assert read(OUTPUT)==result
    else:
        with OUTPUT.open('x')as stream:json.dump(result,stream,indent=2,sort_keys=True);stream.write('\n')
    print('PASS matched arithmetic panel; completed',result['completed_cases'])
