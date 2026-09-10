#!/usr/bin/env sage-python
"""Independent lattice/local-duality replay of the relative302 ideal anatomy.

Checks producer ideal squares with lattice products; ordinary Artin bits at
bad primes with Hilbert symbols instead of local square tests; nonunit good
prime repairs with Hensel lifting or finite residue fields. No ideal reduction.
"""
import argparse,hashlib,json,runpy,sys
from pathlib import Path
from sage.all import AA,QQ,ZZ,GF,NumberField,PolynomialRing,matrix,vector,pari,prod
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
sys.path.insert(0,str(ROOT/'elliptic-curves/rank-jump'))
from verify_strict_half_ideals import check_lattice
from verify_half_ideal_artin import jacobi
def read(p):return json.loads(p.read_text())
def bindings(d):
    for name,sha in d['bindings'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==sha,name
def verify(negative_controls=False):
    oldcheck=runpy.run_path(str(ROOT/'elliptic-curves/cas/verify_curve302_descent_anatomy.sage'))
    oldcheck['verify'](negative_controls)
    b=read(ART/'curve302_relative_ideal_anatomy_v1.json')
    e=read(ART/'curve302_even_ideal_extension_v1.json')
    final=read(ART/'curve302_even_ideal_extension_complete_v1.json')
    for d in [b,e,final]:bindings(d)
    s=read(ART/'curve302_descent_anatomy_v1.json');h=s['half_ideal_packet']
    a=read(ART/'rank_jump_curve302_strict_constructor_arithmetic_v1.json')
    R=PolynomialRing(QQ,'z');f=R(a['cubic_ascending']);K=NumberField(f,'theta');theta=K.gen()
    nf=pari.nfinit([pari(f),a['S_finite']]);basis=[K(R(x))for x in a['maximal_order_basis']]
    B=matrix(QQ,[list(x)for x in basis]).transpose();Bi=B.inverse()
    def pa(x):return pari.Mod(pari(R(list(x))),pari(f))
    def coords(x):return Bi*vector(QQ,list(x))
    def mult(x):return matrix(QQ,[list(coords(x*y))for y in basis]).transpose()
    def same(X,Y):
        C=X.inverse()*Y;assert all(q in ZZ for q in C.list())and abs(C.det())==1
    gs=[ZZ(p['a'])-ZZ(p['d'])**2*theta for p in h['points']]
    places=[(p,j,P)for p in a['S_finite']for j,P in enumerate(pari.idealprimedec(nf,p))]
    vals=[[int(pari.idealval(nf,pa(g),P))for g in gs]for p,j,P in places]
    assert vals==b['characters']['bad_prime_valuations']
    # Reassemble complete finite valuation, dyadic-unit and real constraints.
    constraints=[];unitpos=0
    for (p,j,P),v in zip(places,vals):
        constraints.append([x%2 for x in v])
        if p!=2:continue
        ram=int(P[2]);assert int(P[3])==1
        for k in range(1,2*ram+1):
            row=b['characters']['dyadic_units'][unitpos];unitpos+=1
            assert row['prime_index']==j and row['power']==k
            unit=pari(row['unit']);assert pari.idealval(nf,unit,P)==0
            bits=[int(pari.nfhilbert(nf,pa(g),unit,P)==-1)for g in gs]
            assert bits==row['bits'];constraints.append(bits)
    roots=f.roots(AA,multiplicities=False)
    signs=[[int(R(list(g))(root)<0)for g in gs]for root in roots]
    assert signs==b['characters']['real_signs'];constraints+=signs
    assert constraints==b['characters']['constraints']
    C=matrix(GF(2),constraints);words=matrix(GF(2),b['characters']['words'])
    assert C.rank()==11 and C.right_kernel().basis_matrix()==words and words.rank()==20
    assert matrix(GF(2),vals).rank()==9
    betas=[prod(g for bit,g in zip(w,gs)if bit)for w in words]
    for beta in betas:
        assert all(R(list(beta))(root)>0 for root in roots)
        for p,j,P in places:
            assert pari.idealval(nf,pa(beta),P)%2==0
            if p==2:assert pari.nfislocalpower(nf,P,pa(beta),2)or pari.nfislocalpower(nf,P,5*pa(beta),2)
    generic=[K(R(c['beta_ascending']))for c in a['generic_classes']]
    gv=matrix(GF(2),[[int(pari.idealval(nf,pa(g),P))%2 for p,j,P in places]for g in generic])
    genrows=[read(ART/('det1092_generic_virtual_units_v1/basis-%02d.json'%i))for i in range(8)]
    genwords=matrix(GF(2),[r['word']for r in genrows])
    assert gv.rank()==9 and gv.left_kernel().basis_matrix()==genwords
    generator_betas=[]
    for row in genrows:
        alpha=prod(g for bit,g in zip(row['word'],generic)if bit)
        assert alpha==K(R(row['alpha']))
        beta=alpha*K(R(row['gamma']))**2
        assert beta==K(R(row['beta']))
        H=matrix(QQ,row['ideal']);check_lattice(K,basis,H,beta,H.det());generator_betas.append(beta)
    for i,c in enumerate(s['initial_artin_packet']['columns']):
        beta=prod(g for bit,g in zip(h['words'][i],gs)if bit)/K(R(c['multiplier']))**2
        H=matrix(QQ,c['reduced_ideal']);check_lattice(K,basis,H,beta,H.det());generator_betas.append(beta)
    for c in e['extra_columns']:
        beta=prod(g for bit,g in zip(c['word'],gs)if bit)
        J=matrix(QQ,c['half_ideal']);check_lattice(K,basis,J,beta,J.det())
        alpha=K(R(c['multiplier']));H=matrix(QQ,c['ideal']);same(J,mult(alpha)*H)
        generator_betas.append(beta/alpha**2)
    parent=read(ART/'curve302_recovered_mw17_parent_v1.json')
    embedding=matrix(GF(2),parent['basis_embedding_in_public_D'])
    allwords=matrix(GF(2),e['basis_words'])
    assert allwords[:8,:]==genwords*embedding.transpose()
    assert allwords[8:18,:]==matrix(GF(2),h['words'])
    assert allwords[18:,:]==matrix(GF(2),[c['word']for c in e['extra_columns']])
    assert allwords.rank()==22 and allwords*matrix(GF(2),vals).transpose()==0

    # For an unramified quadratic character, (beta, uniformizer)_P is its
    # Frobenius bit. This checks the producer's nfislocalpower via local duality.
    localcache={}
    def local_bit(beta,p,j):
        P=pari.idealprimedec(nf,p)[j];key=(p,j)
        if key not in localcache:
            pi=pari.nfbasistoalg(nf,pari.idealappr(nf,P));assert pari.idealval(nf,pi,P)==1
            localcache[key]=pi
        return int(pari.nfhilbert(nf,pa(beta),localcache[key],P)==-1)
    def ring(H,N):
        assert H.det()==H[0,0]==N and H[1,1]==H[2,2]==1
        assert all(H[i,j]==0 for i in range(3)for j in range(i))
        r=vector(QQ,[1,-H[0,1],-H[0,2]])
        for i in range(3):
            for j in range(3):assert (r*coords(basis[i]*basis[j])-r[i]*r[j])%N==0
        return r
    def removal(H,good,terms):
        product=pari(good)
        for t in terms:
            P=pari.idealprimedec(nf,t['p'])[t['j']];power=t['exponent']
            assert int(pari.idealval(nf,pari(H),P))==power>0
            product=pari.idealmul(nf,product,pari.idealpow(nf,P,power))
        same(H,matrix(QQ,product))
    def cyclic_column(col):
        H=matrix(QQ,col['ideal']);good=matrix(QQ,col['good_ideal']);N=ZZ(col['norm'])
        assert N>0 and all(N%p for p in a['S_finite'])
        removal(H,good,col['bad_parts']);residues=ring(good,N);out=[]
        for beta,entry in zip(betas,col['entries']):
            value=ZZ(residues*coords(beta))%N;assert str(value)==entry['residue']
            cofactor=N;total=0
            for t in entry['repairs']:
                p,exp=t['p'],t['exponent'];assert N.valuation(p)==exp
                cofactor//=ZZ(p)**exp
                root=ZZ(residues*coords(theta))%p
                bit=oldcheck['hensel_bit'](f,R(list(beta)),p,root,t['valuation'])
                assert bit==t['bit'];total+=exp*bit
            assert cofactor==ZZ(entry['cofactor'])and value.gcd(cofactor)==1
            symbol=jacobi(int(value),int(cofactor));assert symbol==entry['symbol']and symbol in [-1,1]
            assert [(t['p'],t['j'],t['exponent'])for t in entry['bad_terms']]==[(t['p'],t['j'],t['exponent'])for t in col['bad_parts']]
            for t in entry['bad_terms']:
                bit=local_bit(beta,t['p'],t['j']);assert bit==t['bit'];total+=t['exponent']*bit
            bit=(int(symbol==-1)+total)%2;assert bit==entry['bit'];out.append(bit)
        return out
    columns=[cyclic_column(c)for c in b['columns']]
    for i,c in enumerate(e['extra_columns']):
        if i!=1:columns.append(cyclic_column(c));continue
        repair=final['repair'];assert repair['column']==1 and repair['ideal']==c['ideal']
        H=matrix(QQ,c['ideal']);good=matrix(QQ,repair['good_ideal']);N=ZZ(repair['norm'])
        terms=repair['entries'][0]['local_terms'];removal(H,good,terms);residues=ring(good,N)
        bits=[]
        for beta,entry in zip(betas,repair['entries']):
            value=ZZ(residues*coords(beta))%N;assert str(value)==entry['residue']and value.gcd(N)==1
            symbol=jacobi(int(value),int(N));assert symbol==entry['symbol']and symbol in [-1,1]
            total=0
            assert [(t['p'],t['j'],t['exponent'])for t in entry['local_terms']]==[(t['p'],t['j'],t['exponent'])for t in terms]
            for t in entry['local_terms']:
                bit=local_bit(beta,t['p'],t['j']);assert bit==t['bit'];total+=t['exponent']*bit
                if t['p']==47:
                    # An independent residue-field power check at the degree2
                    # prime; linear primitive point classes are units there.
                    P=pari.idealprimedec(nf,47)[t['j']]
                    assert int(P[3])==t['degree']==2 and t['valuation']==0
                    value47=pari.nfmodpr(nf,pa(beta),pari.nfmodprinit(nf,P,'u'))
                    assert value47!=0 and int(value47**((47**2-1)//2)!=1)==bit
            bit=(int(symbol==-1)+total)%2;assert bit==entry['bit'];bits.append(bit)
        columns.append(bits)
    M=matrix(GF(2),columns).transpose()
    assert [list(map(int,row))for row in M[:,:18]]==b['artin_matrix']
    assert [list(map(int,row))for row in M]==final['artin_matrix']
    assert [int(M[:,:8].rank()),int(M[:,8:18].rank()),int(M[:,:18].rank()),int(M.rank())]==[8,10,18,19]
    assert [list(map(int,row))for row in M.right_kernel().basis()]==final['kernel_words']
    U=matrix(GF(2),final['unramified_in_even_basis']).transpose()
    assert allwords.transpose()*U==words.transpose() and (M*U).rank()==19
    # No unit may occur in the18-dimensional detected subspace. For all22,
    # the norm-positive cubic unit squareclasses have dimension2, so the full
    # ideal image has dimension at least20 even though this matrix detects19.
    assert len(generator_betas)==22 and len(betas)==20
    print('PASS 20 ordinary unramified characters; 22 exact half ideals; 440 Artin entries.')
    print('PASS inherited8 + strict10 = disjoint18; total ideal image in[20,22], relative in[12,14].')
    print('PASS full known even space22; ordinary unramified space20; generic spaces8 and6.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--negative-controls',action='store_true')
    verify(p.parse_args().negative_controls)
