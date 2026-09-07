#!/usr/bin/env python3
"""Bounded oracle-assisted metric translations in the fixed old27 subgroup."""
from importlib.machinery import SourceFileLoader
from sage.all import EllipticCurve, QQ, ZZ, matrix, vector, pari
from fpylll import GSO, Enumeration, IntegerMatrix
import inventory188_nearest_translate_control_v2 as control
from research_runtime.store import checkpoint

geometry = SourceFileLoader('nearest_geometry', str(control.CAS/'prospective_half_lattice_v2.sage')).load_module()


def main():
    protocol = control.protocol()
    if control.OUT.exists():
        raise FileExistsError('preserve nearest-translate result')
    public = control.read(control.PUBLIC)
    E = EllipticCurve(QQ, public['curve'])
    old = [E(list(map(QQ, p))) for p in public['points'][:27]]
    charts = control.charts()
    output = {'status': 'RUNNING', 'protocol_sha256': control.hashed(control.D/'protocol.json'),
              'charts': charts, 'rows': []}
    checkpoint(control.D/'checkpoint.json', output)
    for index in protocol['public_indices']:
        P = E(list(map(QQ, public['transported_public_points'][index])))
        points = old+[P]
        gram, asym = geometry.canonical_height_gram(tuple(map(control.F, public['curve'])),
            [tuple(map(control.F, map(str, p.xy()))) for p in points])
        full = matrix(ZZ, geometry.rounded_gram(gram, 1000000))
        assert full.is_positive_definite()
        G = full[:27, :27]
        U = matrix(ZZ, pari(G).qflllgram()).transpose()
        assert abs(U.det()) == 1
        reduced = U*G*U.transpose()
        cross = U*full.column(27)[:27]
        target = -reduced.change_ring(QQ).solve_right(cross)
        gso = GSO.Mat(IntegerMatrix.from_matrix([list(map(int,r)) for r in reduced.rows()]),
                      gram=True, float_type='dd', update=True)
        transformed_target = [float(target[i])+sum(float(target[j])*gso.get_mu(j,i)
                                                  for j in range(i+1,27)) for i in range(27)]
        radius = float(sum(abs(x) for x in reduced.list())/4)+1.0
        solutions = Enumeration(gso).enumerate(0, 27, radius, 0, target=transformed_target)
        assert solutions
        reported, coeff = solutions[0]
        nearest = vector(ZZ, [round(x) for x in coeff])
        assert all(abs(float(x)-int(y)) < 1e-7 for x,y in zip(coeff, nearest))
        offset = nearest-target
        distance = offset*reduced*offset
        assert abs(float(distance)-reported) < 0.001
        choices = [vector(ZZ, [0]*27), nearest]
        for i in range(27):
            for sign in (-1,1):
                choice = vector(ZZ, list(nearest)); choice[i] += sign; choices.append(choice)
        choices = list(dict.fromkeys(tuple(map(int,x)) for x in choices))
        assert len(choices) in (55,56), 'complete fixed neighbour roster required'
        row = {'public_index': index, 'metric_gram': [[str(x) for x in r] for r in gram],
               'rounded_full_gram': [list(map(int,r)) for r in full.rows()],
               'change_of_basis': [list(map(int,r)) for r in U.rows()],
               'reduced_gram': [list(map(int,r)) for r in reduced.rows()],
               'target_reduced_coefficients': list(map(str,target)),
               'nearest_reduced_word': list(map(int,nearest)), 'exact_rounded_distance': str(distance),
               'floating_distance_error': abs(float(distance)-reported), 'representatives': []}
        for reduced_word in choices:
            word = vector(ZZ, reduced_word)*U
            Q = P+sum((c*B for c,B in zip(word,old)), E(0))
            assert Q and Q in E
            observations = []
            for chart in charts:
                x,y = map(QQ, chart['base_point']); a,b,c,d = map(QQ,chart['matrix'])
                assert Q[0] != x
                for sign in (-1,1):
                    t = (sign*Q[1]+y)/(Q[0]-x)
                    num,den = d*t-b,a-c*t
                    if den:
                        s = num/den; coordinate = [str(s.numerator()),str(s.denominator())]
                        height = max(abs(s.numerator()),s.denominator()); infinity = False
                    else:
                        coordinate = ['1','0']; height = ZZ(1); infinity = True
                    observations.append({'arm':chart['arm'],'chart':chart['index'],'sign':sign,
                                         'coordinate':coordinate,'height':int(height),'infinity':infinity,
                                         'within_height_or_infinity':bool(infinity or height <= protocol['height'])})
            row['representatives'].append({'reduced_word':list(reduced_word), 'old_word':list(map(int,word)),
                'point':list(map(str,Q.xy())), 'observations':observations,
                'best':min(observations,key=lambda x:(x['height'],x['arm'],x['chart'],x['sign']))})
        row['best'] = min(({'representative_index':i,**r['best']} for i,r in enumerate(row['representatives'])),
                          key=lambda x:(x['height'],x['representative_index'],x['arm'],x['chart'],x['sign']))
        output['rows'].append(row)
        checkpoint(control.D/'checkpoint.json',output)
        print('WITNESS',index,'REPRESENTATIVES',len(choices),'BEST',row['best'],flush=True)
    output['status'] = 'PASS_EXACT_TRANSLATE_COORDINATES'
    output['observation_count'] = sum(len(r['observations']) for row in output['rows'] for r in row['representatives'])
    output['visible_count'] = sum(v['within_height_or_infinity'] for row in output['rows'] for r in row['representatives'] for v in r['observations'])
    output['claim_boundary'] = protocol['scope']
    checkpoint(control.OUT,output)
    checkpoint(control.D/'checkpoint.json',output)
    print(output['status'],output['observation_count'],output['visible_count'],flush=True)


if __name__ == '__main__':
    main()
