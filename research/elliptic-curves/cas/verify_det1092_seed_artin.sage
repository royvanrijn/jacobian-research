#!/usr/bin/env sage-python
"""Independent Artin replay: local square tests, integer Jacobi, cubic HNF.

No Hilbert symbols, factorization, or PARI ideal multiplication. The two
relative classes are also exported as explicit integral cubic norm equations.
Run under timeout25s.
"""
import hashlib,json
from pathlib import Path
from sage.all import QQ,ZZ,GF,AA,PolynomialRing,NumberField,pari,matrix,vector,lcm
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_seed_artin_v2';ARITH=ART/'rank_jump_curve302_strict_constructor_arithmetic_v1.json'
SEED=ART/'det1092_seed_half_ideal_v1/302-first-unlock.json'
VIRTUAL=ART/'det1092_seed_half_ideal_v1/virtual-unit.json'
GEN=[ART/('det1092_generic_virtual_units_v1/basis-%02d.json'%i) for i in range(8)]
INPUTS=[ARITH,SEED,VIRTUAL,*GEN,OUT/'characters.json',OUT/'summary.json',
        OUT/'relative-summary.json',*[OUT/('column-%02d.json'%i) for i in range(9)],
        *[OUT/('relative-%02d.json'%i) for i in range(2)]]
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,data):
    path=OUT/name
    if path.exists():assert read(path)==data
    else:
        with path.open('x') as stream:json.dump(data,stream,indent=2,sort_keys=True);stream.write('\n')
def entries(A):return [[str(x) for x in row] for row in A.rows()]
d=read(ARITH);seed=read(SEED);sv=read(VIRTUAL);gen=[read(p) for p in GEN]
chars=read(OUT/'characters.json')['characters'];summary=read(OUT/'summary.json')
R=PolynomialRing(QQ,'z');f=R(d['cubic_ascending']);K=NumberField(f,'theta')
nf=pari.nfinit([pari(f),d['S_finite']]);assert ZZ(nf.disc())==ZZ(d['field_discriminant'])
zk=[K(R(x)) for x in d['maximal_order_basis']]
Z=matrix(QQ,[a.vector() for a in zk]).transpose();Zi=Z.inverse();O=matrix.identity(QQ,3)
def pa(a):return pari.Mod(pari(R(a.list())),pari(f))
def coordinates(a):return Zi*a.vector()
def element(v):return sum((v[i]*zk[i] for i in range(3)),K.zero())
def canonical(A):
    scale=lcm([q.denominator() for q in A.list()])
    return matrix(ZZ,(scale*A).transpose()).hermite_form(include_zero_rows=False).transpose()/scale
def product(A,B):
    return canonical(matrix(QQ,[coordinates(element(a)*element(b)) for a in A.columns() for b in B.columns()]).transpose())
def principal(a):return canonical(matrix(QQ,[coordinates(a*b) for b in zk]).transpose())
def parse(text):return matrix(QQ,[[QQ(x) for x in row.split(',')] for row in text[1:-1].split(';')])
def jacobi(a,n):
    a=int(a)%int(n);n=int(n);assert n>0 and n%2;sign=1
    while a:
        while not a%2:
            a//=2
            if n%8 in [3,5]:sign=-sign
        a,n=n,a
        if a%4==n%4==3:sign=-sign
        a%=n
    return sign if n==1 else 0
alphas=[K(R(a['beta_ascending'])) for a in d['generic_classes']]+[K(R(seed['alpha']))]
generic_ideals=[]
for g in gen:
    a=K.one()
    for bit,x in zip(g['word'],alphas[:17]):
        if bit:a*=x
    beta=K(R(g['beta']));gamma=K(R(g['gamma']));I=matrix(QQ,g['ideal'])
    assert a==K(R(g['alpha'])) and beta==a*gamma**2 and product(I,I)==principal(beta)
    generic_ideals.append(I)
words=matrix(GF(2),[c['word'] for c in chars]);assert words.rank()==7
places={}
for block in d['prime_decomposition']:
    p=block['p'];primes=pari.idealprimedec(nf,p)
    for i,row in enumerate(block['primes']):
        B=parse(row['hnf']);P=primes[i]
        assert matrix(QQ,pari.idealhnf(nf,P).sage())==B
        places[p,i]=(B,P)
roots=f.roots(AA,multiplicities=False);char_elements=[];local={}
for c in chars:
    a=K.one()
    for bit,x in zip(c['word'],alphas):
        if bit:a*=x
    assert a==K(R(c['alpha'])) and all(R(a.list())(root)>0 for root in roots)
    for (p,i),(B,P) in places.items():
        assert int(pari.idealval(nf,pa(a),P))%2==0
        if p==2:assert pari.nfislocalpower(nf,P,pa(a),2) or pari.nfislocalpower(nf,P,pa(5*a),2)
    char_elements.append(a)
for item in summary['local_symbols']:
    key=item['p'],item['index'];P=places[key][1]
    bits=[int(not pari.nfislocalpower(nf,P,pa(a),2)) for a in char_elements]
    assert bits==item['bits'];local[key]=bits
columns=[]
for i in range(9):
    c=read(OUT/('column-%02d.json'%i));assert c['status']=='PASS'
    I=matrix(QQ,c['ideal']);H=matrix(QQ,c['good_ideal']);reconstructed=H
    bad=vector(GF(2),[0]*7)
    for part in c['bad_parts']:
        key=part['p'],part['index'];B=places[key][0]
        for _ in range(part['valuation']):reconstructed=product(reconstructed,B)
        if part['valuation']%2:bad+=vector(GF(2),local[key])
    assert canonical(I)==canonical(reconstructed)
    N=ZZ(c['good_norm']);assert abs(H.det())==N and N%2 and H[0,0]==N and H[1,1]==H[2,2]==1
    def residue(a):
        v=coordinates(a);assert all(x.denominator()==1 for x in v)
        return ZZ(v[0]-H[0,1]*v[1]-H[0,2]*v[2])%N
    # Direct ring homomorphism check on an integral basis.
    for a in zk:
        for b in zk:assert (residue(a*b)-residue(a)*residue(b))%N==0
    residues=[residue(a) for a in alphas];assert list(map(str,residues))==c['residues']
    symbols=[jacobi(a,N) for a in residues];assert all(x in [-1,1] for x in symbols)
    raw=vector(GF(2),[int(x==-1) for x in symbols]);good=words*raw;bits=bad+good
    assert list(map(int,bits))==c['artin_bits'];columns.append(bits)
M=matrix(GF(2),columns);assert M[:8,:6].rank()==6 and M[:8,:].rank()==M.rank()==7
relative=read(OUT/'relative-summary.json');A=M[:8,:];w0=vector(GF(2),relative['compatible_words'][0])
ker=A.left_kernel().basis_matrix();assert ker.nrows()==1 and list(map(int,ker.row(0)))==relative['generic_kernel']
assert all(vector(GF(2),w)*A==M.row(8) for w in relative['compatible_words'])
assert vector(GF(2),relative['compatible_words'][1])==w0+ker.row(0)
# Seven class characters are independent on the eight inherited ideals.
minor_rows=list(A.transpose().pivots());assert A.matrix_from_rows(minor_rows).det()==1
# Encode the two remaining principality questions as integral cubic forms.
T=PolynomialRing(QQ,['u','v','w']);variables=T.gens();norm_forms=[]
for index in range(2):
    r=read(OUT/('relative-%02d.json'%index));alpha=K(R(sv['beta']))
    I=matrix(QQ,r['unreduced_ideal']);J=matrix(QQ,r['ideal']);gamma=K(R(r['gamma']));beta=K(R(r['beta']))
    for bit,g in zip(r['generic_virtual_basis_word'],gen):
        if bit:alpha/=K(R(g['beta']))
    assert alpha==K(R(r['alpha'])) and beta==alpha*gamma**2
    assert product(I,principal(gamma))==canonical(J) and product(J,J)==principal(beta)
    # Independently verify the relative class relation without ideal division.
    reconstruct=I
    for bit,B in zip(r['generic_virtual_basis_word'],generic_ideals):
        if bit:reconstruct=product(reconstruct,B)
    assert canonical(reconstruct)==canonical(matrix(QQ,sv['reduced_square_root_ideal']))
    N=abs(J.det());assert N==QQ(r['norm']) and beta.norm()==N*N
    basis=[element(c) for c in J.columns()]
    multiplication=matrix(T,3,3,lambda i,j:sum(variables[k]*(basis[k]*K.gen()**j)[i] for k in range(3)))
    F=T(multiplication.det()/N);assert all(q.denominator()==1 for q in F.coefficients())
    # A rational point exists automatically; principality requires integers.
    point=(Z*J).solve_right((K(N)/beta).vector());assert F(*point)==1
    form={'index':index,'normalized_norm_form':str(F),
          'coefficients':[{'exponents':list(map(int,e)),'coefficient':str(c)} for e,c in F.dict().items()],
          'ideal_basis':entries(J),'ideal_norm':str(N),'rational_point':list(map(str,point)),
          'rational_point_identity':'a=N(J)/beta has field norm N(J), since N(beta)=N(J)^2',
          'principal_ideal_criterion':'J principal iff F(u,v,w)=1 has a solution in Z^3',
          'local_statement':'The integral equation is soluble over every Z_p and over R; proof is in the canonical note.'}
    save('norm-form-%02d.json'%index,form);norm_forms.append(form)
report={'status':'PASS_INDEPENDENT_ARTIN_AND_TWO_NORM_FORM_REPLAY',
    'classification':'verified application and new exact arithmetic reduction',
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [*INPUTS,Path(__file__)]},
    'generic_character_rank':6,'with_seed_character_rank':7,
    'generic_ideal_artin_rank':7,'with_seed_ideal_artin_rank':7,
    'generic_minor_rows':minor_rows,'compatible_relative_ideals':2,
    'generic_unit_intersection_dimension_upper_bound':1,
    'norm_positive_unit_classes_outside_generic_exist':True,
    'norm_forms_have_explicit_rational_points':True,
    'boundary':'Neither principal-ideal question is solved. Artin-kernel membership does not prove a unit or generic ideal class. The new class character uses the old seed, so no prospective seed producer is certified.'}
save('independent-replay.json',report);print(report['status'],flush=True)
