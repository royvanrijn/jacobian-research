#!/usr/bin/env sage-python
"""Nine fixed half ideals against six generic plus one seed-derived character.

No factoring: remove only fixed bad-prime ideals, then use a rational Jacobi
symbol only when the remaining quotient is odd cyclic and every input is a
unit. Otherwise record UNKNOWN without changing the ideal. <=25 seconds.
"""
import hashlib,json
from pathlib import Path
from sage.all import QQ,ZZ,GF,AA,PolynomialRing,NumberField,pari,matrix,vector
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_seed_artin_v2'
ARITH=ART/'rank_jump_curve302_strict_constructor_arithmetic_v1.json'
ANCHORS=ART/'rank_jump_generic_only_class_anchors_v1.json'
SEED=ART/'det1092_seed_half_ideal_v1/302-first-unlock.json'
VIRTUAL=ART/'det1092_seed_half_ideal_v1/virtual-unit.json'
GEN=[ART/('det1092_generic_virtual_units_v1/basis-%02d.json'%i) for i in range(8)]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,data):
    path=OUT/name
    if path.exists():assert json.loads(path.read_text())==data
    else:
        with path.open('x') as stream:json.dump(data,stream,indent=2,sort_keys=True);stream.write('\n')
def rows(A):return [[str(x) for x in row] for row in A.rows()]
OUT.mkdir(exist_ok=True)
save('protocol.json',{'classification':'retrospective class-image test; not prospective construction',
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [ARITH,ANCHORS,SEED,VIRTUAL,*GEN,Path(__file__)]},
    'limits':{'seconds':25,'ideal_columns':9,'generic_characters':6,
              'seed_derived_characters_max':1,'factorizations':0,'class_groups':0,
              'unit_groups':0,'point_searches':0,'ideal_retries':0},
    'method':'Reuse six frozen generic unramified characters; append a canonical generic correction of the seed only if it is unramified. Use local Hilbert symbols at fixed bad primes and a factor-free Jacobi symbol on an odd cyclic good quotient. A failed quotient/unit gate stays UNKNOWN.'})
d=json.loads(ARITH.read_text());anc=json.loads(ANCHORS.read_text())['cases'][1]
s=json.loads(SEED.read_text());sv=json.loads(VIRTUAL.read_text())
R=PolynomialRing(QQ,'z');f=R(d['cubic_ascending']);K=NumberField(f,'theta')
nf=pari.nfinit([pari(f),d['S_finite']]);assert ZZ(nf.disc())==ZZ(d['field_discriminant'])
def pa(a):return pari.Mod(pari(R(a.list())),pari(f))
def polynomial(a):return R([QQ(pari.lift(a).polcoef(i)) for i in range(3)])
alphas=[K(R(x['beta_ascending'])) for x in d['generic_classes']]+[K(R(s['alpha']))]
gammas=[pa(a) for a in alphas];roots=f.roots(AA,multiplicities=False)
places=[(p,i,P) for p in d['S_finite'] for i,P in enumerate(pari.idealprimedec(nf,p))]
constraints=[];dyadic=[]
for p,i,P in places:
    constraints.append([int(pari.idealval(nf,a,P))%2 for a in gammas])
    if p!=2:continue
    e=int(P[2]);assert int(P[3])==1
    pi=pari(2) if e==1 else pari.nfbasistoalg(nf,pari.idealappr(nf,P))
    assert int(pari.idealval(nf,pi,P))==1
    for k in range(1,2*e+1):constraints.append([int(pari.nfhilbert(nf,a,1+pi**k,P)==-1) for a in gammas])
    dyadic.append(P)
constraints += [[int(R(a.list())(root)<0) for a in alphas] for root in roots]
A=matrix(GF(2),constraints);Ag=A[:,:17]
assert Ag.rank()==11 and Ag.right_kernel().dimension()==6
words=[vector(GF(2),[(int(mask)>>i)&1 for i in range(17)]+[0]) for mask in anc['generic_coefficient_masks']]
assert all(A*w==0 for w in words)
if A.rank()==Ag.rank():
    correction=Ag.solve_right(A.column(17));words.append(vector(GF(2),list(correction)+[1]))
else:correction=None
assert matrix(GF(2),words).rank()==len(words)
char_checks=[]
for word in words:
    alpha=pari.Mod(1,pari(f))
    for bit,a in zip(word,gammas):
        if bit:alpha*=a
    checks=[]
    for P in dyadic:
        square=bool(pari.nfislocalpower(nf,P,alpha,2))
        five_square=bool(pari.nfislocalpower(nf,P,5*alpha,2))
        assert square or five_square;checks.append([square,five_square])
    assert all(polynomial(alpha)(root)>0 for root in roots)
    char_checks.append({'word':list(map(int,word)),'alpha':list(map(str,polynomial(alpha).list())),
                        'dyadic_square_or_5square':checks})
save('characters.json',{'status':'PASS_UNRAMIFIED_CHARACTER_EXTENSION',
    'generic_constraint_rank':int(Ag.rank()),'augmented_constraint_rank':int(A.rank()),
    'constraints':rows(A),'characters':char_checks,
    'classification':'six inherited characters and at most one retrospective seed-derived character'})
coords=[list(map(ZZ,pari.nfalgtobasis(nf,a))) for a in gammas]
localcache={};outputs=[]
sources=[json.loads(p.read_text()) for p in GEN]+[{'ideal':sv['reduced_square_root_ideal'],'word':sv['generic_word']+[1]}]
for index,source in enumerate(sources):
    I=pari(matrix(QQ,source['ideal']));good=I;bad=[];bad_bits=vector(GF(2),[0]*len(words))
    for p,i,P in places:
        e=int(pari.idealval(nf,I,P))
        if e:
            good=pari.idealmul(nf,good,pari.idealpow(nf,P,-e))
            bad.append({'p':p,'index':i,'valuation':e})
        if e%2:
            key=(p,i)
            if key not in localcache:
                pi=pari.nfbasistoalg(nf,pari.idealappr(nf,P));assert pari.idealval(nf,pi,P)==1
                raw=vector(GF(2),[int(pari.nfhilbert(nf,a,pi,P)==-1) for a in gammas])
                localcache[key]={'p':p,'index':i,'uniformizer':list(map(str,polynomial(pi).list())),
                                 'raw_bits':list(map(int,raw)),'bits':[int(raw*w) for w in words]}
            bad_bits+=vector(GF(2),localcache[key]['bits'])
    H=matrix(QQ,good.sage());N=ZZ(pari.idealnorm(nf,good))
    assert N>0 and all(N%p for p in d['S_finite'])
    cyclic=N%2==1 and H[0,0]==N and H[1,1]==H[2,2]==1
    record={'index':index,'ideal':source['ideal'],'bad_parts':bad,'good_ideal':rows(H),
            'good_norm':str(N),'odd_cyclic':bool(cyclic),'bad_bits':list(map(int,bad_bits)),
            'status':'UNKNOWN','reason':'noncyclic good quotient'}
    if cyclic:
        residues=[ZZ(c[0]-H[0,1]*c[1]-H[0,2]*c[2])%N for c in coords]
        gcds=[v.gcd(N) for v in residues]
        record.update(residues=list(map(str,residues)),gcds=list(map(str,gcds)))
        if all(g==1 for g in gcds):
            raw=vector(GF(2),[int(pari.kronecker(v,N)==-1) for v in residues])
            good_bits=vector(GF(2),[raw*w for w in words]);bits=bad_bits+good_bits
            record.update(status='PASS',reason=None,good_raw_bits=list(map(int,raw)),
                          good_bits=list(map(int,good_bits)),artin_bits=list(map(int,bits)))
        else:record['reason']='a marked factor is not a unit on the good quotient'
    save('column-%02d.json'%index,record);outputs.append(record)
    print('column',index,record['status'],record.get('artin_bits',record['reason']),flush=True)
complete=all(r['status']=='PASS' for r in outputs)
summary={'status':'PASS' if complete else 'PARTIAL','classification':'retrospective exact Artin projection',
    'character_count':len(words),'complete':complete,'local_symbols':list(localcache.values()),
    'boundary':'A rank increase proves nongeneric ideal-class information. No increase is inconclusive for principal ideals or units, since the character roster can be incomplete and quadratic characters annihilate2Cl.'}
if complete:
    M=matrix(GF(2),[r['artin_bits'] for r in outputs]);g=M[:8,:].rank();full=M.rank()
    summary.update(generic_image_rank=int(g),augmented_image_rank=int(full),relative_image_rank=int(full-g))
    if full>g:
        witness=M.solve_right(vector(GF(2),[0]*8+[1]));summary['separating_character_word']=list(map(int,witness))
save('summary.json',summary)
print('summary',summary.get('generic_image_rank'),summary.get('augmented_image_rank'),flush=True)
