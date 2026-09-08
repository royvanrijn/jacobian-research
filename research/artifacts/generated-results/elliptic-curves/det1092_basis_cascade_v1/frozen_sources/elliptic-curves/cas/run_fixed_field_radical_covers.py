#!/usr/bin/env sage
"""Three fixed-field radical covers: minimal models, lattice search and descent.

The default is a cheap, offline exact replay. --prepare builds just the six
original/Q-translated presentations and standalone Magma inputs. Online jobs
are opt-in, bounded by the public calculator, with exact input/output retained.
Neither a search miss nor a resource failure is an obstruction certificate.
"""
from __future__ import annotations

import argparse
import ast
import gzip
import hashlib
import json
from pathlib import Path
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

from sage.all import QQ, ZZ, matrix, vector
from sage.version import version as sage_version
import run_fixed_field_point_realization as base

ROOT = base.ROOT
CT = ROOT / 'artifacts/generated-results/elliptic-curves/fixed_cubic_u_minus1_cassels_tate_v1.json'
WORK = ROOT / 'artifacts/local/fixed-field-radical-solve'
SUMMARY = ROOT / 'artifacts/generated-results/elliptic-curves/fixed_field_radical_models_v1.json'
EVIDENCE = SUMMARY.with_name('fixed_field_radical_models_evidence_v1.json.gz')
MASKS = (1047173, 596921, 450876)
LABELS = {596921: 'r_1', 450876: 'r_2', 1047173: 'r_1+r_2'}
BOUND = 10**7


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def approved_masks():
    ct = json.loads(CT.read_text())['arithmetic']
    assert ct['restricted_radical_dimension'] == 2
    masks = [z['anchor_mask'] for z in ct['restricted_radical_basis']]
    assert set(masks + [masks[0] ^ masks[1]]) == set(MASKS)
    return MASKS


def initial_quadrics(coefficients):
    """Twice the symmetric matrices of XZ-Y^2 and W^2-F(X,Y,Z)."""
    e, d, c, b, a = map(QQ, coefficients)
    q1 = matrix(QQ, 4)
    q1[0, 2] = q1[2, 0] = 1
    q1[1, 1] = -2
    q2 = matrix(QQ, [[-2*a, -b, 0, 0], [-b, -2*c, -d, 0],
                     [0, -d, -2*e, 0], [0, 0, 0, 2]])
    return [q1, q2]


def invariants(qs):
    """Integral degree-four invariants from det(x A+B), using Hessians."""
    f = (base.x * qs[0] + qs[1]).det()
    e, d, c, b, a = [f[i] for i in range(5)]
    I = 12*a*e - 3*b*d + c*c
    J = 72*a*c*e - 27*a*d*d - 27*b*b*e + 9*b*c*d - 2*c**3
    return I, J/2, (I**3-(J/2)**2)/1728


def verify_quadric_model(row, model, E):
    assert not any(map(QQ,row['reduced_quartic_model'][1])), \
        'complete the square before forming this quadric intersection'
    qs = initial_quadrics(row['reduced_quartic_model'][0])
    eq = matrix(QQ, 2, model['equation_transform'])
    var = matrix(QQ, 4, model['variable_transform'])
    new = [matrix(QQ, 4, z) for z in model['quadric_matrices']]
    assert eq.det() and var.det()
    for i in range(2):
        assert new[i] == sum((eq[i, j] * var * qs[j] * var.transpose()
                              for j in range(2)), matrix(QQ, 4))
        assert new[i].is_symmetric()
        assert all(z in ZZ for z in new[i].list())
        assert all(new[i][j,j] % 2 == 0 for j in range(4))
    c4, c6, disc = invariants(new)
    # A nonsingular integral locally soluble genus-one model has nonnegative
    # level. Equality with the minimal Jacobian discriminant proves level zero.
    Emin = E.minimal_model()
    assert disc == Emin.discriminant()
    assert c4 == Emin.c4() and c6 == Emin.c6()
    return new, var


def quadric_to_quartic(coordinates):
    X, Y, Z, W = map(QQ, coordinates)
    assert X*Z == Y*Y
    if Z:
        s, t, v = Y/Z, QQ(1), W/Z
    else:
        assert X and not Y
        s, t, v = QQ(1), QQ(0), W/X
    return s,t,v


def map_point(row, model, coordinates, data, A, B, E):
    qs, var = verify_quadric_model(row, model, E)
    pt = vector(QQ, coordinates)
    assert len(pt) == 4 and any(pt)
    assert all(pt*q*pt == 0 for q in qs)
    s,t,v = quadric_to_quartic(pt * var)
    assert v*v == base.homogeneous(row['reduced_quartic_model'][0], s,t,4)
    reduced_point = [str(z) for z in (s,t,v)]
    s,t,v = base.inverse_change(base.inverse_change((s,t,v), row['reduction_change']),
                                row['minimization_change'])
    gamma = matrix(QQ, row['conic_parameter_matrix']) * vector([s*s,s*t,t*t])
    d = v / QQ(row['quartic_ordinate_scale'])
    beta, root = tuple(map(QQ, row['beta'])), QQ(row['norm_square_root'])
    result = base.verify_point(beta,root,gamma,d,A,B,-1,E)
    result['reduced_quartic_point'] = reduced_point
    result['quadric_point'] = list(map(str,pt))
    if row['translated_by_universal_point']:
        beta0,root0 = base.class_input(data,row['mask'],A,B)
        result['original_W_class_realization'] = base.translate_back(
            result,beta0,root0,beta,A,B,E)
    return result


def magma_input(row):
    assert not any(map(QQ,row['reduced_quartic_model'][1]))
    e,d,c,b,a = row['reduced_quartic_model'][0]
    return f'''SetColumns(0); SetSeed(1); GetVersion(); Q:=Rationals();
S<X,Y,Z,W>:=PolynomialRing(Q,4);
g4:=GenusOneModel([X*Z-Y^2,W^2-(({a})*X^2+({b})*X*Y+({c})*Y^2+({d})*Y*Z+({e})*Z^2)]);
gm4,tr4,lev4:=Minimise(g4); gr4,rr4:=Reduce(gm4);
assert (rr4*tr4)*g4 eq gr4;
print "QI_MATRICES",[Eltseq(m):m in Matrices(gr4)]; print "EQ_TRANS",Eltseq(Tuple(rr4*tr4)[1]);print "VAR_TRANS",Eltseq(Tuple(rr4*tr4)[2]);
SetVerbose("QISearch",1); print "SEARCH_START"; pts:=PointsQI(Curve(gr4),10^7 : OnlyOne:=true); print "POINTS",[Eltseq(pt):pt in pts];print "DONE";
'''


def four_descent_input(coefficients, primes):
    return ('SetColumns(0); SetSeed(1); GetVersion();\n'
            'StoreFactor([' + ','.join(map(str, primes)) + ']);\n'
            'P<x>:=PolynomialRing(Integers());\n'
            'f:=P![' + ','.join(coefficients) + '];\n'
            'SetVerbose("FourDescent",3); SetVerbose("Factorization",1); '
            'SetVerbose("ClassGroup",1);\n'
            'print "START_FOUR_DESCENT";\ncs:=FourDescent(f);\n'
            'print "FOUR_DESCENT_COUNT",#cs;\n'
            'for c in cs do print DefiningEquations(c); end for;\nprint "DONE";\n')


def run_online(path):
    """Preserve server warnings too: a DONE after an error is not success."""
    request = urllib.request.Request('https://magma.maths.usyd.edu.au/xml/calculator.xml',
        data=urllib.parse.urlencode({'input':path.read_text()}).encode(),
        headers={'Referer':'https://magma.maths.usyd.edu.au/calc/',
                 'Content-Type':'application/x-www-form-urlencoded'})
    started = time.monotonic()
    try:
        raw = urllib.request.urlopen(request,timeout=80).read()
        path.with_suffix('.xml').write_bytes(raw)
        _,text = parse_xml(raw.decode())
        path.with_suffix('.out').write_text(text)
    except Exception as exc:
        path.with_suffix('.error').write_text(repr(exc))
        raise
    print(path.name, round(time.monotonic()-started,2), 'seconds', flush=True)


def parse_xml(raw):
    root = ET.fromstring(raw)
    text = '\n'.join(''.join(z.itertext()) for z in root.findall('.//results/line'))
    error = ('Runtime error','User error','System Error','Syntax error',
             'memory limit','Internal error')
    good = not root.findall('.//warning') and not any(s in text for s in error)
    return good and 'DONE' in text, text


def literal(text, label):
    start = text.index(label) + len(label)
    start = text.index('[',start)
    depth = 0
    for end in range(start,len(text)):
        depth += (text[end] == '[') - (text[end] == ']')
        if depth == 0:
            # Only numbers, brackets, commas and whitespace are accepted.
            expr = text[start:end+1]
            assert re.fullmatch(r'[\[\],\s0-9+\-/]+',expr)
            # ast parses structure; fractions are handled without eval.
            tree = ast.parse(expr,mode='eval').body
            def read(node):
                if isinstance(node,ast.List): return [read(z) for z in node.elts]
                if isinstance(node,ast.Constant) and type(node.value) is int:
                    return QQ(node.value)
                if isinstance(node,ast.UnaryOp) and isinstance(node.op,ast.USub):
                    return -read(node.operand)
                if isinstance(node,ast.BinOp) and isinstance(node.op,ast.Div):
                    return read(node.left)/read(node.right)
                raise ValueError('unexpected output syntax')
            return read(tree)
    raise ValueError('unterminated output list')


def parse_model(raw):
    ok,text = parse_xml(raw)
    assert ok, 'incomplete or failed Magma job'
    return {'quadric_matrices':[[str(z) for z in r] for r in literal(text,'QI_MATRICES')],
            'equation_transform':list(map(str,literal(text,'EQ_TRANS'))),
            'variable_transform':list(map(str,literal(text,'VAR_TRANS'))),
            'returned_points':[[str(z) for z in r] for r in literal(text,'POINTS')]}


def descent_status(raw):
    if raw is None:
        return 'TRANSPORT_ERROR_NO_DESCENT_RESULT'
    good,text = parse_xml(raw)
    if good and re.search(r'FOUR_DESCENT_COUNT\s+\d+',text):
        return 'COMPLETED_REQUIRES_ARITHMETIC_REPLAY'
    if 'memory limit' in text:
        return 'MEMORY_LIMIT_NO_DESCENT_RESULT'
    if 'time limit' in raw:
        return 'TIME_LIMIT_NO_DESCENT_RESULT'
    return 'ERROR_NO_DESCENT_RESULT'


def prepare(work):
    work.mkdir(parents=True,exist_ok=True)
    _,run,_,_,_ = base.context(base.SOURCE,-1)
    for mask in approved_masks():
        for tr in (0,1):
            path=work/f'cover-{mask}-{tr}.json'
            if not path.exists():
                # Height one is a smoke/infinity check, not the point method.
                base.worker(base.SOURCE,-1,mask,1,path,bool(tr))
            row=json.loads(path.read_text())
            (work/f'qi-{mask}-{tr}.m').write_text(magma_input(row))
            (work/f'descent-{mask}-{tr}.m').write_text(four_descent_input(
                row['reduced_quartic_model'][0],run['complete_finite_place_support']))


def audit(evidence):
    assert evidence['source_sha256'] == sha(base.SOURCE)
    assert evidence['pairing_summary_sha256'] == sha(CT)
    data,run,A,B,E = base.context(base.SOURCE,-1)
    assert len(evidence['covers']) == 6
    assert {(r['cover']['mask'],r['cover']['translated_by_universal_point'])
            for r in evidence['covers']} == {(m,bool(t)) for m in approved_masks() for t in (0,1)}
    realized=[]
    for record in evidence['covers']:
        row = record['cover']
        assert row['source_sha256'] == sha(base.SOURCE)
        base.replay_row(row,data,A,B,-1,E)
        assert record['magma_input'] == magma_input(row)
        model = parse_model(record['magma_xml'])
        assert model == record['quadric_model']
        verify_quadric_model(row,model,E)
        mapped=[map_point(row,model,pt,data,A,B,E) for pt in model['returned_points']]
        assert mapped == record['verified_points']
        if mapped: realized.append(row['mask'])
    attempts = evidence['descent_attempts']
    assert len(attempts) == 6
    assert {(a['mask'],a['translated']) for a in attempts} == {(m,t) for m in MASKS for t in (0,1)}
    for attempt in attempts:
        row=next(r['cover'] for r in evidence['covers'] if
                 r['cover']['mask']==attempt['mask'] and
                 r['cover']['translated_by_universal_point']==attempt['translated'])
        assert attempt['magma_input']==four_descent_input(
            row['reduced_quartic_model'][0],run['complete_finite_place_support'])
        assert attempt['status']==descent_status(attempt['magma_xml'])
        assert attempt['magma_xml'] is not None or attempt.get('transport_error')
        assert attempt['status']!='COMPLETED_REQUIRES_ARITHMETIC_REPLAY', \
            'a successful descent must be examined and certified before publishing'
    # Every point is certified in the original anchor squareclass basis;
    # Q is independently separated by valuation parity above 19.
    qcert=base.universal_point_certificate(data,A,B,E)
    rank=base.f2_rank([[m>>j&1 for j in range(20)] for m in realized])
    return {'exact_minimal_quadric_models':6,'jacobian_minimal_discriminant':str(E.minimal_model().discriminant()),
            'verified_points':sum(len(r['verified_points']) for r in evidence['covers']),
            'realized_anchor_masks':sorted(set(realized)),
            'certified_realized_subspace_dimension':rank,
            'certified_curve_rank_lower_bound':rank+qcert['certified_rank_lower_bound'],
            'incomplete_four_descent_attempts':len(attempts),
            'point_or_sha':{str(m):('RATIONAL_POINT' if m in realized else 'UNKNOWN') for m in MASKS}}


def collect(work):
    data,run,A,B,E=base.context(base.SOURCE,-1)
    evidence={'source_sha256':sha(base.SOURCE),'pairing_summary_sha256':sha(CT),
              'covers':[],'descent_attempts':[]}
    if (work/'quartic_bnf.py').exists() and (work/'quartic_bnf.log').exists():
        evidence['local_field_pilot']={
            'input':(work/'quartic_bnf.py').read_text(),
            'output':(work/'quartic_bnf.log').read_text(),
            'command':'timeout 300 sage -python -u artifacts/local/fixed-field-radical-solve/quartic_bnf.py',
            'wall_seconds_limit':300,'pari_stack_bytes':1024_000_000,
            'status':'NO_BNF_RESULT_NOT_AN_OBSTRUCTION',
            'note':'bnfinit would require certification before any unconditional use; it did not return.'}
    for mask in approved_masks():
        for tr in (0,1):
            row=json.loads((work/f'cover-{mask}-{tr}.json').read_text())
            raw=(work/f'qi-{mask}-{tr}.xml').read_text()
            model=parse_model(raw)
            evidence['covers'].append({'cover':row,'magma_input':magma_input(row),
                'magma_xml':raw,'quadric_model':model,
                'verified_points':[map_point(row,model,pt,data,A,B,E) for pt in model['returned_points']]})
            descent_path=work/f'descent-{mask}-{tr}.xml'
            descent_raw=descent_path.read_text() if descent_path.exists() else None
            evidence['descent_attempts'].append({'mask':mask,'translated':tr,
                'magma_input':four_descent_input(row['reduced_quartic_model'][0],run['complete_finite_place_support']),
                'magma_xml':descent_raw,'status':descent_status(descent_raw),
                'transport_error':None if descent_raw is not None else descent_path.with_suffix('.error').read_text()})
    result=audit(evidence)
    summary={'schema':'elliptic-curves.fixed-field-radical-models.v1',
             'status':'EXACT_MINIMAL_MODELS_WITH_BOUNDED_LATTICE_SEARCH',
             'arithmetic':result,'limits':{'classes':3,'Q_translations_per_class':2,
                'pointsqi_height_argument':BOUND,'only_one':True,'calculator_seconds_per_job':60},
             'method':'Magma degree-four Minimise/Reduce and Elkies p-adic lattice PointsQI',
             'software':{'sage':sage_version,'magma':'2.29-9'},
             'claim_boundary':['No bounded miss or failed descent is a Sha certificate.',
                'No complete Mordell-Weil basis, full Selmer group or exact rank is asserted.'],
             'evidence':str(EVIDENCE.relative_to(ROOT)),
             'source_hashes':{str(Path(__file__).resolve().relative_to(ROOT)):sha(__file__),
                 str(Path(base.__file__).resolve().relative_to(ROOT)):sha(base.__file__)}}
    packed=gzip.compress((json.dumps(evidence,sort_keys=True)+'\n').encode(),mtime=0)
    EVIDENCE.write_bytes(packed)
    summary['evidence_sha256']=sha(EVIDENCE)
    base.save(SUMMARY,summary)
    print(json.dumps(result,indent=2))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--workdir',type=Path,default=WORK)
    ap.add_argument('--prepare',action='store_true')
    ap.add_argument('--online-qi',action='store_true')
    ap.add_argument('--online-descent',action='store_true')
    ap.add_argument('--collect',action='store_true')
    args=ap.parse_args()
    if args.prepare: prepare(args.workdir)
    if args.online_qi or args.online_descent:
        for m in approved_masks():
            for tr in (0,1):
                for kind,active in [('qi',args.online_qi),('descent',args.online_descent)]:
                    if active:
                        path=args.workdir/f'{kind}-{m}-{tr}.m'
                        if not path.with_suffix('.xml').exists() and not path.with_suffix('.error').exists():
                            try:run_online(path)
                            except Exception as exc:print(path.name,'NO_RESULT',repr(exc),flush=True)
    if args.collect: collect(args.workdir)
    if not any([args.prepare,args.online_qi,args.online_descent,args.collect]):
        summary=json.loads(SUMMARY.read_text())
        assert sha(EVIDENCE)==summary['evidence_sha256']
        for path,digest in summary['source_hashes'].items():assert sha(ROOT/path)==digest
        result=audit(json.loads(gzip.decompress(EVIDENCE.read_bytes())))
        assert result==summary['arithmetic']
        print('PASS_EXACT_SIX_RADICAL_MODELS',json.dumps(result))


if __name__=='__main__':main()
