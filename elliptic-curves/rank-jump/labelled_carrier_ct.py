#!/usr/bin/env python3
"""Observe a labelled Cassels pairing in an isolated, arithmetic-unchanged PARI copy."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import re
import subprocess
import tarfile
import urllib.request
import retrospective as r

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'LABELLED_CARRIER_CT_PROTOCOL.json'
TARGET=r.OUT/'rank_jump_disjoint_soluble_carriers_v1.json'
CLASS=r.OUT/'rank_jump_carrier_sha_class_v1.json'
INPUT=r.OUT/'rank_jump_labelled_carrier_ct_inputs_v2.json'
OUTPUT=r.OUT/'rank_jump_labelled_carrier_ct_v2.json'
WORK=r.ROOT/'artifacts/local/rank-jump-labelled-carrier-ct-v1'
RUNTIME=Path('/home/royvanrijn/.local/share/jacobian-sage-10.9')

EXPORTS=['hyperell_locally_soluble','nf_hyperell_locally_soluble','nfhilbert','nfhilbert0',
         'ellrankinit','ell2selmer_basis','ellrank','ell2cover']
MARKER='  if (DEBUGLEVEL) err_printf("Selmer rank: %ld\\n", dim);'
OBSERVATIONS='''
  err_printf("RJ_MODEL %Ps\\n", mkvec5(gen_0,ell_get_a2(ell),gen_0,ell_get_a4(ell),ell_get_a6(ell)));
  err_printf("RJ_CUBIC %Ps\\n", pol);
  for (i = 1; i <= dim; i++)
    err_printf("RJ_BASIS %ld %Ps\\n", i, RgXQV_factorback(LS2,gel(selmer,i),pol));
'''
MAIN='''
int main(int argc, char **argv)
{
  GEN E, answer;
  if (argc != 2) return 2;
  pari_init(67108864, 500000);
  DEBUGLEVEL_ellrank = 2;
  setrand(gen_1);
  E = ellinit(gp_read_str(argv[1]), NULL, DEFAULTPREC);
  answer = ellrank(E, 0, NULL, DEFAULTPREC);
  pari_printf("RJ_RESULT %Ps\\n", answer);
  pari_close();
  return 0;
}
'''


def new_bytes(path,data):
    with path.open('xb') as f:f.write(data)


def capture():
    policy=r.read(PROTOCOL);WORK.mkdir(parents=True,exist_ok=True)
    source_path=WORK/'upstream_ellrank.c';copyright_path=WORK/'COPYING'
    if not source_path.exists():
        archive=urllib.request.urlopen(policy['upstream']['url'],timeout=20).read()
        assert r.digest(archive)==policy['upstream']['archive_sha256']
        tar=tarfile.open(fileobj=io.BytesIO(archive),mode='r:gz')
        source=tar.extractfile(next(m for m in tar.getmembers() if m.name.endswith('/src/basemath/ellrank.c'))).read()
        license=tar.extractfile(next(m for m in tar.getmembers() if m.name.endswith('/COPYING') and m.name.count('/')==1)).read()
        new_bytes(source_path,source);new_bytes(copyright_path,license)
    raw=source_path.read_bytes();assert r.digest(raw)==policy['upstream']['ellrank_c_sha256']
    text=raw.decode();assert text.count(MARKER)==1
    prefix=''.join('#define '+name+' rank_jump_'+name+'\n' for name in EXPORTS)
    instrumented=prefix+text.replace(MARKER,MARKER+'\n'+OBSERVATIONS)+MAIN
    modified=WORK/'observed_ellrank.c';binary=WORK/'observed_ellrank'
    if not modified.exists():new_bytes(modified,instrumented.encode())
    assert modified.read_text()==instrumented
    command=['cc','-O2','-I'+str(RUNTIME/'include/pari'),str(modified),'-L'+str(RUNTIME/'lib'),
             '-Wl,-rpath,'+str(RUNTIME/'lib'),'-lpari','-lm','-o',str(binary)]
    if not (WORK/'compile.json').exists():
        try:
            p=subprocess.run(command,capture_output=True,text=True,timeout=30)
            result={'returncode':p.returncode,'stdout':p.stdout,'stderr':p.stderr,'command':command}
        except subprocess.TimeoutExpired:result={'status':'TIMEOUT','command':command}
        r.write_new(WORK/'compile.json',result)
    compiled=r.read(WORK/'compile.json');assert compiled.get('returncode')==0,compiled
    model=r.read(TARGET)['rows'][0]['geometry']['minimal_Jacobian_model']
    if not (WORK/'descent.json').exists():
        try:
            p=subprocess.run([str(binary),'['+','.join(model)+']'],capture_output=True,text=True,timeout=30)
            result={'returncode':p.returncode,'stdout':p.stdout,'stderr':p.stderr}
        except subprocess.TimeoutExpired:result={'status':'TIMEOUT'}
        r.write_new(WORK/'descent.json',result)
    descent=r.read(WORK/'descent.json');assert descent.get('returncode')==0,descent
    stderr=descent['stderr'];lines=stderr.splitlines()
    def one(prefix):
        hits=[line[len(prefix):] for line in lines if line.startswith(prefix)]
        assert len(hits)==1,(prefix,len(hits));return hits[0]
    basis=[line[len('RJ_BASIS '):].split(' ',1) for line in lines if line.startswith('RJ_BASIS ')]
    assert [int(i) for i,x in basis]==list(range(1,6))
    matrix=one('Cassels Pairing: ')
    resultline=next(x[len('RJ_RESULT '):] for x in descent['stdout'].splitlines() if x.startswith('RJ_RESULT '))
    r.write_new(INPUT,{'schema':'rank-jump.labelled-carrier-ct-inputs.v1',
                       'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (PROTOCOL,TARGET,CLASS,Path(__file__))},
                       'upstream':policy['upstream'],'instrumented_source_sha256':r.digest(modified.read_bytes()),
                       'instrumentation':{'renamed_exports':EXPORTS,'insertion_marker':MARKER,'observations':OBSERVATIONS,'main':MAIN},
                       'compile':compiled,'linked_runtime_sha256':r.digest((RUNTIME/'lib/libpari-gmp-tls.so.2.17.3').read_bytes()),
                       'original_target_model':model,'internal_model_GP':one('RJ_MODEL '),'internal_cubic_GP':one('RJ_CUBIC '),
                       'basis_GP':[x for i,x in basis],'Cassels_matrix_GP':matrix,'rank_result_GP':resultline,
                       'descent_transcript':descent,
                       'boundary':'Only observational output and exported-name isolation are added to the pinned PARI algorithm.'})
    print('Captured labelled basis and Cassels matrix',matrix,flush=True)


def build():
    from sage.all import QQ,PolynomialRing,EllipticCurve,pari,matrix,GF,vector
    from cover_experiment import sqrtq
    inp=r.read(INPUT);old=r.read(CLASS);mapping=old['mapping']
    internal=[str(x) for x in pari(inp['internal_model_GP'])]
    polynomial=PolynomialRing(QQ,'x');x=polynomial.gen()
    def convert(text):
        f=pari(text)
        return polynomial([QQ(f.polcoef(i)) for i in range(int(f.poldegree())+1)])
    cubic=convert(inp['internal_cubic_GP'])
    basis=[convert(b) for b in inp['basis_GP']]
    assert all(b.degree()<=2 and b.gcd(cubic)==1 for b in basis)
    E=EllipticCurve(QQ,list(map(QQ,internal)));Er=EllipticCurve(QQ,list(map(QQ,mapping['Jacobian_model'])))
    assert E.is_isomorphic(Er) and E.a1()==E.a3()==0
    assert cubic==x**3+E.a2()*x*x+E.a4()*x+E.a6()
    short,_=r.short(internal,[])
    lam2=sqrtq(r.F(mapping['Jacobian_model'][3])/r.F(short[3]));assert lam2 is not None
    lam=sqrtq(lam2);assert lam is not None
    assert r.F(mapping['Jacobian_model'][4])==lam**6*r.F(short[4])
    b0,b1=map(QQ,mapping['beta'])
    beta=polynomial(b0+b1*QQ(str(lam2))*(x+E.a2()/3))
    polys=basis+[beta];sigs=[0]*6;primes=[];col=0
    for p in r.primes(503):
        if cubic.discriminant().denominator()%p==0 or r.mod(str(cubic.discriminant()),p)==0:continue
        if any(c.denominator()%p==0 for f in polys+[cubic] for c in f):continue
        coeff=[r.mod(str(c),p) for c in cubic.list()]
        roots=[t for t in range(p) if sum(c*pow(t,i,p) for i,c in enumerate(coeff))%p==0]
        if len(roots)!=3:continue
        vals=[[sum(r.mod(str(c),p)*pow(t,i,p) for i,c in enumerate(f.list()))%p for t in roots] for f in polys]
        if any(v==0 for vs in vals for v in vs):continue
        primes.append(p)
        for j,vs in enumerate(vals):
            for k,val in enumerate(vs):sigs[j]|=int(pow(val,(p-1)//2,p)==p-1)<<(col+k)
        col+=3
    assert r.rank(sigs[:5])==5 and r.rank(sigs)==5
    coordinates=[]
    for mask in range(32):
        value=0
        for i in range(5):
            if mask>>i&1:value^=sigs[i]
        if value==sigs[5]:coordinates.append(mask)
    assert len(coordinates)==1;mask=coordinates[0]
    C=matrix(GF(2),pari(inp['Cassels_matrix_GP']).sage());assert C.nrows()==C.ncols()==5
    assert C==C.transpose() and all(C[i,i]==0 for i in range(5)) and C.rank()==2
    coeff=vector(GF(2),[(mask>>i)&1 for i in range(5)]);pairings=C*coeff
    rank=pari(inp['rank_result_GP']);assert tuple(int(rank[i]) for i in range(3))==(2,2,2)
    rational=not any(pairings)
    result={'schema':'rank-jump.labelled-carrier-ct.v1','status':'PASS',
            'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (INPUT,PROTOCOL,TARGET,CLASS,Path(__file__),HERE/'retrospective.py')},
            'internal_model':internal,'cubic_coefficients':list(map(str,cubic.list())),
            'Selmer_basis_coefficients':[list(map(str,f.list())) for f in basis],
            'carrier_beta_coefficients':list(map(str,beta.list())),'raw_to_internal_short_scaling':str(lam),
            'proof_primes':primes,'Selmer_fingerprints':sigs[:5],'carrier_fingerprint':sigs[5],
            'Selmer_character_rank':5,'carrier_coordinates':[int(c) for c in coeff],
            'Cassels_matrix':[[int(c) for c in row] for row in C.rows()],
            'carrier_pairing_row':[int(c) for c in pairings],
            'Cassels_rank':2,'radical_dimension':3,'full_rational_Kummer_dimension':3,
            'carrier_global_solubility':'YES' if rational else 'NO',
            'carrier_nonzero_Sha_class':not rational,
            'boundary':'Complete local solubility and the exact rank/CT dimensions come from the pinned preceding descent. The labelled pairing identifies this particular carrier class; it does not classify other covers or original fibres.'}
    r.write_new(OUTPUT,result)
    print('carrier coordinates',result['carrier_coordinates'],'pairing',result['carrier_pairing_row'],'global',result['carrier_global_solubility'])


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','build']);a=p.parse_args()
    capture() if a.mode=='capture' else build()
