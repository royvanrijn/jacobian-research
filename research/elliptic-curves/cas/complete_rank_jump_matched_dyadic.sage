#!/usr/bin/env sage-python
"""Complete the one retained dyadic residue-degree-three case via full unit bid."""
import hashlib,json,sys
from pathlib import Path
from sage.all import AA,QQ,ZZ,GF,EllipticCurve,PolynomialRing,matrix,pari,prime_range
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
PROTOCOL=ART/'rank_jump_matched_descent_panel_protocol_v1.json';INITIAL=ART/'rank_jump_matched_descent_panel_v1.json'
OUTPUT=ART/'rank_jump_matched_descent_dyadic_completion_v1.json'
sys.path.insert(0,str(Path(__file__).parent))
from research_runtime.local_kummer import LocalSquareclasses
from r17_60_arithmetic import native_check
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def compute():
    case='11952-high-02';p=read(PROTOCOL);r=next(r for r in p['rows']if r['id']==case)
    initial=next(x for x in read(INITIAL)['cases']if x['id']==case)
    assert initial['status']=='UNKNOWN_DYADIC_RESIDUE_DEGREE'and initial['remaining_cofactor']=='1'
    native_check({k:r[k]for k in ['curve','points','proof','rank_lower_bound']},r['roster_row'])
    R=PolynomialRing(QQ,'z');f=R(initial['cubic_ascending']);S=sorted({2,*map(int,r['retained_prime_hints'])})
    remaining=abs(ZZ(f.discriminant()))
    for q in S:
        assert ZZ(q).is_prime(proof=True);remaining//=ZZ(q)**remaining.valuation(q)
    assert remaining==1 and f.is_irreducible()
    nf=pari.nfinit([pari(f),S]);polys=[];gs=[];E=EllipticCurve(QQ,r['curve']);n=r['rank_lower_bound']
    for point in r['points']:
        Q=E(point);x,y=Q.xy();d=ZZ(x.denominator()).sqrt();assert d in ZZ and d*d==x.denominator()
        aa,bb=ZZ(x*d*d),ZZ(y*d**3);g=R([aa,-d*d]);beta=pari.Mod(pari(g),pari(f))
        assert pari.nfeltnorm(nf,beta)==bb*bb and aa.gcd(d)==bb.gcd(d)==1
        polys.append(g);gs.append(beta)
    chars=[[]for _ in gs];blocks=[]
    for q0 in prime_range(3,2001):
        q=int(q0)
        if f.discriminant()%q==0:continue
        roots=f.change_ring(GF(q)).roots(multiplicities=False)
        if not roots:continue
        values=[[int(g.change_ring(GF(q))(t))for t in roots]for g in polys]
        if any(v==0 for row in values for v in row):continue
        for row,vs in zip(chars,values):row.extend(int(pow(v,(q-1)//2,q)==q-1)for v in vs)
        blocks.append({'prime':q,'roots':list(map(int,roots))})
        if matrix(GF(2),chars).rank()==n:break
    assert matrix(GF(2),chars).rank()==n
    joint=[[]for _ in gs];local=[];val=[];unram=[];dyadic=[]
    for q in S:
        L=LocalSquareclasses(nf,q);signatures=[list(map(int,L.signature(g)))for g in gs]
        for row,sig in zip(joint,signatures):row.extend(sig)
        local.append({'place':q,'point_image_dimension':int(L.point_kummer_dimension),'signatures':signatures,
                      'generic_image_rank':int(matrix(GF(2),signatures[:17]).rank()),'known_image_rank':int(matrix(GF(2),signatures).rank())})
        for j,P in enumerate(pari.idealprimedec(nf,q)):
            vals=[int(pari.idealval(nf,g,P))%2 for g in gs];val.append(vals);unram.append(vals)
            if q!=2:continue
            bid=pari.idealstar(nf,pari.idealpow(nf,P,2*int(P[2])+1),2)
            cyc=bid.bid_get_cyc();gens=bid.bid_get_gen();records=[]
            for order,g in zip(cyc,gens):
                if int(order)%2:continue
                unit=pari.nfbasistoalg(nf,g)if g.type()=='t_COL'else g
                assert pari.idealval(nf,unit,P)==0
                bits=[int(pari.nfhilbert(nf,beta,unit,P)==-1)for beta in gs]
                unram.append(bits);records.append({'order':str(order),'unit':str(unit),'bits':bits})
            dyadic.append({'prime_index':j,'ramification':int(P[2]),'residue_degree':int(P[3]),'units':records})
    roots=f.roots(AA,multiplicities=False);real=[[int(g(root)<0)for root in roots]for g in polys]
    for row,sig in zip(joint,real):row.extend(sig)
    unram.extend([list(v)for v in zip(*real)])
    J,V,U=matrix(GF(2),joint),matrix(GF(2),val),matrix(GF(2),unram)
    unitbound=(3+len(roots))//2-1
    dims={'generic':17,'known':n,'generic_local':int(J[:17,:].rank()),'known_local':int(J.rank()),
          'generic_strict':17-int(J[:17,:].rank()),'known_strict':n-int(J.rank()),
          'generic_even':17-int(V[:,:17].rank()),'known_even':n-int(V.rank()),
          'generic_unramified':17-int(U[:,:17].rank()),'known_unramified':n-int(U.rank()),
          'local_product':sum(x['point_image_dimension']for x in local)+(1 if len(roots)==3 else 0)}
    dims['relative_ideal_lower_bound']=max(0,dims['known_even']-dims['generic_even']-unitbound)
    result={**initial,'schema':'rank-jump.matched-descent-dyadic-completion.v1','status':'PASS_EXACT_POINT_SUBSPACE_ANATOMY','stage':'complete',
            'bindings':{str(path.relative_to(ROOT)):sha(path)for path in [Path(__file__),PROTOCOL,INITIAL]},
            'field_discriminant':str(nf.disc()),'maximal_order_basis':list(map(str,nf.nf_get_zk())),
            'signature':[len(roots),(3-len(roots))//2],'dimensions':dims,'good_character_blocks':blocks,
            'local':local,'real_signatures':real,'valuation_matrix':val,'unramified_constraint_matrix':unram,
            'dyadic_full_units':dyadic,'unit_dimension_upper_bound':unitbound,
            'strict_words':[list(map(int,w))for w in J.left_kernel().basis()],
            'generic_strict_words':[list(map(int,w))for w in J[:17,:].left_kernel().basis()],
            'full_class_group':'UNKNOWN','full_Selmer':'UNKNOWN',
            'bounds':{'wall_seconds':60,'cases':1,'point_searches':0,'factorization':0},
            'boundary':'Exact completion of the same frozen residue-degree-three field, using all even cyclic factors of(O_K/P^(2e+1))*. No parameter replacement.'}
    if OUTPUT.exists():assert read(OUTPUT)==result
    else:
        with OUTPUT.open('x')as out:json.dump(result,out,indent=2,sort_keys=True);out.write('\n')
    print('PASS',case,dims,flush=True)
if __name__=='__main__':compute()
