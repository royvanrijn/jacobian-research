#!/usr/bin/env sage-python
"""Four canonical remaining virtual-unit directions; fixed Artin comparison.

Extend the certified18-dimensional G_even+strict subspace to all22 everywhere-
even known classes. One reduction per new ideal; no BNF/point search. Any
undetected kernel is UNKNOWN. Run under timeout120s.
"""
import hashlib
import json
from pathlib import Path
from sage.all import QQ, ZZ, GF, PolynomialRing, matrix, vector, pari, prod, prime_range

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
BASE = ART/'curve302_relative_ideal_anatomy_v1.json'
STRICT = ART/'curve302_descent_anatomy_v1.json'
ARITH = ART/'rank_jump_curve302_strict_constructor_arithmetic_v1.json'
PARENT = ART/'curve302_recovered_mw17_parent_v1.json'
GEN = [ART/('det1092_generic_virtual_units_v1/basis-%02d.json' % i) for i in range(8)]
WORK = ROOT/'artifacts/local/elliptic-curves/curve302-even-ideal-extension-v1'
OUTPUT = ART/'curve302_even_ideal_extension_v1.json'


def read(p): return json.loads(p.read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def rows(M): return [[str(x) for x in row] for row in matrix(QQ, M).rows()]
def save(path, data):
    if path.exists(): assert read(path) == data
    else:
        with path.open('x') as out:
            json.dump(data, out, indent=2, sort_keys=True); out.write('\n')


def compute():
    WORK.mkdir(parents=True, exist_ok=True)
    bindings = {str(p.relative_to(ROOT)): sha(p) for p in [Path(__file__), BASE, STRICT, ARITH, PARENT, *GEN]}
    protocol = {'schema': 'curve302.even-ideal-extension.protocol.v1', 'bindings': bindings,
                'selection': 'Append the first independent vectors of the canonical full valuation kernel to G_even+V; no Artin data select them.',
                'bounds': {'wall_seconds': 120, 'workers': 1, 'new_ideal_columns': 4,
                           'nonunit_gcd_trial_prime_bound': 1009, 'point_searches': 0, 'bnf_calls': 0},
                'failure_semantics': 'Nonunit support or noncyclic ring remains UNKNOWN; a character kernel is not a unit relation.'}
    save(WORK/'protocol.json', protocol)
    b, s, a, parent = map(read, [BASE, STRICT, ARITH, PARENT])
    assert b['union_detected_rank'] == 18
    R = PolynomialRing(QQ, 'z'); f = R(a['cubic_ascending'])
    nf = pari.nfinit([pari(f), a['S_finite']])
    h = s['half_ideal_packet']
    gs = [pari.Mod(pari(R([QQ(p['a']), -QQ(p['d'])**2])), pari(f)) for p in h['points']]
    embedding = matrix(GF(2), parent['basis_embedding_in_public_D'])
    first = [embedding*vector(GF(2), read(p)['word']) for p in GEN]
    first += [vector(GF(2), w) for w in h['words']]
    valuation = matrix(GF(2), b['characters']['bad_prime_valuations'])
    E = valuation.right_kernel().basis_matrix()
    assert E.nrows() == 22 and matrix(GF(2), first).rank() == 18
    assert all(valuation*w == 0 for w in first)
    allwords = first[:]
    for w in E.rows():
        if matrix(GF(2), allwords+[w]).rank() > len(allwords): allwords.append(w)
    assert len(allwords) == 22
    extra = allwords[18:]
    charwords = b['characters']['words']
    chars = [prod(g for bit, g in zip(w, gs) if bit) for w in charwords]
    places = [(p,j,P) for p in a['S_finite'] for j,P in enumerate(pari.idealprimedec(nf,p))]
    gcds = [pari.idealadd(nf, QQ(p['b']), g) for p,g in zip(h['points'],gs)]
    columns = []
    for index, word in enumerate(extra):
        path = WORK/('column-%02d.json' % index)
        if path.exists(): columns.append(read(path)); continue
        beta = prod(g for bit,g in zip(word,gs) if bit)
        J = pari.idealhnf(nf,1)
        for bit,I in zip(word,gcds):
            if bit: J=pari.idealmul(nf,J,I)
        for p,j,P in places:
            v=int(pari.idealval(nf,beta,P)); assert v%2==0
            correction=v//2-int(pari.idealval(nf,J,P))
            if correction: J=pari.idealmul(nf,J,pari.idealpow(nf,P,correction))
        assert pari.idealpow(nf,J,2)==pari.idealhnf(nf,beta)
        H,alpha=pari.idealred(nf,[J,1]); assert pari.idealmul(nf,H,alpha)==J
        good=H; bad=[]
        for p,j,P in places:
            e=int(pari.idealval(nf,H,P));assert e>=0
            if e:
                bad.append((p,j,P,e));good=pari.idealmul(nf,good,pari.idealpow(nf,P,-e))
        N=ZZ(pari.idealnorm(nf,good))
        assert N>0 and N%2 and all(N%p for p in a['S_finite'])
        cyclic=good[0,0]==N and good[1,1]==good[2,2]==1
        entries=[]
        for chi in chars:
            if not cyclic: entries.append({'status':'UNKNOWN_NONCYCLIC'});continue
            cs=pari.nfalgtobasis(nf,chi);res=ZZ(cs[0]-good[0,1]*cs[1]-good[0,2]*cs[2])%N
            missing=res.gcd(N);ps=[]
            for p0 in prime_range(3,1010):
                if missing%p0:continue
                while missing%p0==0:missing//=p0
                ps.append(int(p0))
            if missing!=1:entries.append({'status':'UNKNOWN_NONUNIT','remaining':str(missing)});continue
            cofactor=N;repairs=[]
            for p in ps:
                exponent=int(N.valuation(p));cofactor//=ZZ(p)**exponent
                Ps=[P for P in pari.idealprimedec(nf,p) if pari.idealval(nf,good,P)>0]
                assert len(Ps)==1;P=Ps[0]
                assert P[2]==P[3]==1 and pari.idealval(nf,good,P)==exponent
                repairs.append({'p':p,'exponent':exponent,'valuation':int(pari.idealval(nf,chi,P)),
                                'bit':int(not pari.nfislocalpower(nf,P,chi,2))})
            symbol=int(pari.kronecker(res%cofactor,cofactor));assert symbol in [-1,1]
            terms=[{'p':p,'j':j,'exponent':e,'bit':int(not pari.nfislocalpower(nf,P,chi,2))} for p,j,P,e in bad]
            bit=(int(symbol==-1)+sum(t['bit']*t['exponent'] for t in repairs+terms))%2
            entries.append({'residue':str(res),'cofactor':str(cofactor),'symbol':symbol,
                            'bad_terms':terms,'repairs':repairs,'bit':bit})
        row={'index':index,'word':list(map(int,word)),'half_ideal':rows(J),'ideal':rows(H),
             'multiplier':[str(pari.lift(pari.nfbasistoalg(nf,alpha)).polcoef(i)) for i in range(3)],
             'good_ideal':rows(good),'norm':str(N),'entries':entries,
             'bad_parts':[{'p':p,'j':j,'exponent':e} for p,j,P,e in bad]}
        save(path,row);columns.append(row)
        print('EXTRA',index,'COMPLETE',sum('bit'in e for e in entries),flush=True)
    result={'schema':'curve302.even-ideal-extension.v1','bindings':bindings,'protocol':protocol,
            'basis_words':[list(map(int,w)) for w in allwords],'extra_columns':columns,
            'complete':all(all('bit'in e for e in c['entries']) for c in columns)}
    if result['complete']:
        M=matrix(GF(2),b['artin_matrix']).augment(matrix(GF(2),[[e['bit'] for e in c['entries']] for c in columns]).transpose())
        rank=int(M.rank())
        result.update(artin_matrix=[list(map(int,row)) for row in M],detected_rank=rank,
                      kernel_words=[list(map(int,row)) for row in M.right_kernel().basis()],
                      full_even_ideal_rank_interval=[rank,22],
                      full_even_unit_kernel_interval=[0,22-rank],
                      relative_ideal_rank_interval=[rank-8,14])
        # Arithmetic unramified words in the all-even basis also have square
        # half ideals. Record their exact image rank without an upper claim.
        T=matrix(GF(2),allwords).transpose()
        U=T.solve_right(matrix(GF(2),charwords).transpose())
        result['unramified_in_even_basis']=[list(map(int,row)) for row in U.transpose()]
        result['unramified_half_ideal_detected_rank']=int((M*U).rank())
        print('FULL EVEN IDEAL',rank,'UNRAMIFIED IDEAL',result['unramified_half_ideal_detected_rank'],
              'RELATIVE',result['relative_ideal_rank_interval'],'KERNEL',result['kernel_words'],flush=True)
    save(OUTPUT,result)


if __name__=='__main__':compute()
