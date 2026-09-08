#!/usr/bin/env sage-python
"""Export and replay exact compact equations and all generic sections of six R17 families."""
import argparse
from importlib.machinery import SourceFileLoader
from pathlib import Path
import sys
from sage.all import PolynomialRing,QQ,ZZ,Matrix

ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import certify_compact_r17_candidates as cert
helpers=SourceFileLoader('atlas_base_helpers',str(CAS/'reduce_r17_family_base.sage')).load_module()
R=PolynomialRing(QQ,'t');t=R.gen();K=R.fraction_field()


def rational(record):
    return K(R(list(map(QQ,record['numerator_coefficients_low_to_high']))))/R(list(map(QQ,record['denominator_coefficients_low_to_high'])))


def encode(value):
    value=K(value)
    return {'numerator_coefficients_low_to_high':list(map(str,value.numerator().list())),
            'denominator_coefficients_low_to_high':list(map(str,value.denominator().list()))}


def native(source,family):
    if family=='074d9':
        model=source['representative']
        sections=[(K(R(list(map(QQ,r['x_coefficients_low_to_high'])))),K(R(list(map(QQ,r['y_coefficients_low_to_high'])))))
                  for r in source['native_chart_sections']['norm12-orbit-074d9']]
        W=Matrix(ZZ,source['generic_basis']['word_matrix_rows']);G=Matrix(ZZ,source['generic_basis']['height_gram'])
        if abs(W.det())!=1:raise ArithmeticError('native word change is not unimodular')
        gram=W.inverse()*G*W.inverse().transpose()
    else:
        model=source['weierstrass_model'];sections=[(rational(r['X']),rational(r['Y'])) for r in source['sections']['records']]
        gram=Matrix(ZZ,source['sections']['height_gram'])
    if len(sections)!=17 or gram.det()!=948 or any(q not in ZZ for q in gram.list()):raise ArithmeticError('native MW17 certificate changed')
    return model,sections,gram


def transport(row):
    source=ROOT/row['source']
    if cert.hashed(source)!=row['source_sha256']:raise ArithmeticError('literal source changed')
    model,sections,gram=native(cert.read(source),row['family'])
    A,B=(R(list(map(QQ,model[k]))) for k in ('A_coefficients_low_to_high','B_coefficients_low_to_high'))
    new_A,new_B=(R(list(map(QQ,row[k]))) for k in ('A_coefficients_low_to_high','B_coefficients_low_to_high'))
    a,b,c,d=map(QQ,row['base_matrix_a_b_c_d']);u=QQ(row['total_scale_from_literal_source'])
    if not u or a*d==b*c:raise ArithmeticError('singular map')
    if helpers.homogeneous(A,8,a,b,c,d)!=u**4*new_A or helpers.homogeneous(B,12,a,b,c,d)!=u**6*new_B:raise ArithmeticError('coefficient identity failed')
    base=K(a*t+b)/(c*t+d);answer=[]
    for i,(x,y) in enumerate(sections):
        if y**2!=x**3+A*x+B:raise ArithmeticError('native section identity failed')
        xx=K((c*t+d)**4/u**2)*x(base);yy=K((c*t+d)**6/u**3)*y(base)
        if yy**2!=xx**3+new_A*xx+new_B:raise ArithmeticError('transported section identity failed')
        answer.append({'basis_index':i,'X':encode(xx),'Y':encode(yy)})
    return answer,[[int(q) for q in r] for r in gram.rows()]


def build(output):
    if output.exists():raise FileExistsError('use a new immutable export')
    families=[];bindings={}
    for f in ('103b2','11952','074d9','07ca9','08234','08f72'):
        dirname='r17-compactification-160bit-v1' if f in ('11952','08234') else 'r17-compactification-v1'
        directory=ROOT/'artifacts/local/elliptic-curves'/dirname;path=directory/(f+'.json');row=cert.read(path)
        if row['status']!='PASS_EXACT_BASE_CHANGE':raise ArithmeticError('compactification unfinished')
        sections,gram=transport(row)
        families.append({**row,'sections':sections,'generic_height_gram':gram,'generic_height_gram_determinant':948,
            'parameter_convention':'t is the new compact coordinate; native parameter is (a*t+b)/(c*t+d)',
            'base_change_degree':1,'generic_rank_lower_bound':17})
        for q in (path,directory/'protocol.json'):
            bindings[str(q.relative_to(ROOT))]=cert.hashed(q)
        print('TRANSPORTED COMPACT',f,'sections',len(sections),'bits',row['after_bits'],flush=True)
    source_paths=(Path(__file__).resolve(),CAS/'reduce_r17_family_base.sage',CAS/'audit_r17_constant_scaling.py')
    cert.write(output,{'schema':'elliptic-curves.compact-six-r17-atlas.v1','status':'PASS_EXACT_EQUATIONS_AND_102_SECTION_TRANSPORTS',
        'families':families,'discovery_bindings':bindings,'checker_sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in source_paths},
        'claim_boundary':'Six existing generic-rank-17 fibrations in new compact coordinates. No new fibration, generic-rank increase, globally minimal model, or new specialized curve is claimed. Generic independence is inherited through degree-one base/Weierstrass isomorphisms from the pinned native Gram certificates.'})


def check(path):
    data=cert.read(path)
    for name,h in data['checker_sources'].items():
        if cert.hashed(ROOT/name)!=h:raise ArithmeticError('checker changed')
    for row in data['families']:
        sections,gram=transport(row)
        if sections!=row['sections'] or gram!=row['generic_height_gram']:raise ArithmeticError('section or Gram transport changed')
        print('REPLAYED COMPACT FAMILY',row['family'],'all 17 sections',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);g=p.add_mutually_exclusive_group(required=True);g.add_argument('--output',type=Path);g.add_argument('--check',type=Path);a=p.parse_args()
    check(a.check) if a.check else build(a.output)
