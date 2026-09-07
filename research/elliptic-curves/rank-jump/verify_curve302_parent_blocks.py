#!/usr/bin/env python3
"""Independent PARI group-law, quotient and all-prime octic replay (60s cap)."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'artifacts/generated-results/elliptic-curves'
INPUT = OUT/'curve302_parent_blocks_inputs_v1.json'
RESULT = OUT/'curve302_parent_blocks_v1.json'
OUTPUT = OUT/'curve302_parent_blocks_pari_replay_v1.json'
LOCAL = ROOT/'artifacts/local/curve302-parent-block'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mat(rows):
    return '['+';'.join(','.join(map(str, row)) for row in rows)+']'


def vec(values):
    return '['+','.join(map(str, values))+']'


def compute():
    inputs = json.loads(INPUT.read_text()); result = json.loads(RESULT.read_text())
    for data in (inputs, result):
        for path, sha in data['bindings'].items():
            assert digest(ROOT/path) == sha
    # Parse the original public list through Python's literal AST, without
    # importing the generator or its elliptic arithmetic.
    import ast
    tree = ast.parse((ROOT/'elliptic-curves/cas/icarm_curve302.py').read_text())
    literal = next(node.value for node in tree.body if isinstance(node, ast.Assign)
                   and any(isinstance(t, ast.Name) and t.id == 'POINTS' for t in node.targets))
    public = [[ast.literal_eval(coordinate.args[0]) for coordinate in point.elts]
              for point in literal.elts]
    program = [f'E=ellinit({vec(inputs["model"])});',
               'must(x)={if(!x,error("replay failed"))};',
               f'P={vec([vec(p) for p in public])};',
               'P=vector(#P,i,[36*P[i][1]+15,108*(2*P[i][2]+P[i][1]+1)]);',
               'for(i=1,#P,must(ellisoncurve(E,P[i])));', 'print(version());']
    summaries = []
    for parent, row in zip(inputs['parents'], result['parents']):
        assert parent['id'] == row['id']
        m = row['parent_rank']; q = row['residual_free_rank']
        program += [f'C={mat(parent["columns"])}; W={mat(parent["residual_lift_columns"])};',
                    f'B=matconcat([C,W]); K={mat(row["full_inverse_rows"])};',
                    'must(abs(matdet(B))==1);must(K*B==matid(31));',
                    f'Proj={mat(row["projection_rows"])};must(Proj*C==matrix({q},{m}));must(Proj*W==matid({q}));',
                    f'Q=vector({q});']
        for j, point in enumerate(row['residual_points_short_model'], 1):
            program += [f'R=[0];for(i=1,31,R=elladd(E,R,ellmul(E,P[i],W[i,{j}])));',
                        f'must(R=={vec(point)});Q[{j}]=R;']
        program += [f'H=vector({row["pair_count"]});Sc=vector(#H);']
        for k, pair in enumerate(row['cochains'], 1):
            i, j = [x+1 for x in pair['indices_zero_based']]
            program += [f'xp=Q[{i}][1];yp=Q[{i}][2];xq=Q[{j}][1];yq=Q[{j}][2];c=xq-xp;Sc[{k}]=c;',
                        f'H[{k}]=x^8-4*(yp+yq)*x^6+6*c^2*(xp+xq)*x^4-4*c^3*(yq-yp)*x^2+c^6;',
                        f'must(Vecrev(H[{k}])=={vec(pair["octic_ascending"])});']
        scalar_matches = row['scalar_norm_dictionary']['equal_rational_squareclass_pairs_zero_based']
        program += ['matches=List();for(i=1,#Sc,for(j=i+1,#Sc,if(issquare(Sc[i]/Sc[j]),listput(matches,[i-1,j-1]))));',
                    f'must(Vec(matches)=={vec([vec(pair) for pair in scalar_matches])});',
                    f'must(matrank(Mod({mat([p["psi"] for p in row["inert_prime_table"]])},2))=={row["sample_bit_matrix_rank"]});']
        for prime_row in row['inert_prime_table']:
            p = prime_row['prime']
            program.append(f'p={p}; expected={vec(prime_row["psi"])};')
            program.append('for(k=1,#H,F=factormod(H[k],p);must(vecmax(F[,2])==1);d=vecsort(vector(matsize(F)[1],i,poldegree(F[i,1])));must(d==if(expected[k],[2,6],[1,1,3,3])));')
        c = row['bad_components']
        program += [f'D=matdiagonal({vec(c["moduli"])}); J={mat(c["parent_image_columns"])};',
                    f'must(vecsort(Vec(matsnf(matconcat([D,J]))))=={vec(sorted(c["product_mod_parent_smith"]))});']
        summaries.append({'id': parent['id'], 'exact_lift_identities': q,
                          'exact_octic_formulas': row['pair_count'],
                          'independent_prime_factorizations': row['pair_count']*len(row['inert_prime_table']),
                          'independent_scalar_squareclass_comparisons': row['pair_count']*(row['pair_count']-1)//2,
                          'sample_binary_rank': row['sample_bit_matrix_rank'],
                          'integral_projection_and_component_product': 'PASS'})
    program.append('print("PASS");')
    source = '\n'.join(program)+'\n'
    LOCAL.mkdir(parents=True, exist_ok=True)
    (LOCAL/'independent-replay.gp').write_text(source)
    run = subprocess.run(['gp', '-q'], input=source, text=True, capture_output=True, timeout=60, check=True)
    (LOCAL/'independent-replay.out').write_text(run.stdout+run.stderr)
    assert run.stdout.strip().splitlines()[-1] == 'PASS' and not run.stderr.strip(), run.stderr
    return {'schema': 'curve302.parent-blocks-pari-replay.v1', 'status': 'PASS',
            'pari_version': run.stdout.strip().splitlines()[0], 'worker_timeout_seconds': 60,
            'bindings': {str(p.relative_to(ROOT)): digest(p) for p in (INPUT, RESULT, Path(__file__))},
            'rows': summaries,
            'boundary': 'Independent arithmetic replay; governing-field and carrier implications require the written proofs.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('mode', choices=['build', 'check'])
    args = parser.parse_args(); result = compute()
    if args.mode == 'build':
        with OUTPUT.open('x') as f:
            json.dump(result, f, sort_keys=True, indent=2); f.write('\n')
    else:
        assert json.loads(OUTPUT.read_text()) == result
    print('PASS: independent PARI quotient, point, octic and Frobenius replay')
