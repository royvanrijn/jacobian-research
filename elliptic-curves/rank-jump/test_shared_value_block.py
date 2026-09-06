from fractions import Fraction as F
import unittest
import retrospective as r
import shared_value_soluble_block as first
import shared_value_soluble_block_completion as completed


def add(P,Q,a2,a4):
    if P is None:return Q
    if Q is None:return P
    x,y=map(F,P);u,v=map(F,Q)
    if x==u and y==-v:return None
    slope=(3*x*x+2*a2*x+a4)/(2*y) if P==Q else (v-y)/(u-x)
    z=slope*slope-a2-x-u
    return z,-y+slope*(x-z)


class SharedValueTests(unittest.TestCase):
    def test_cyclic_anchor_collapse_independent_group_law(self):
        R=(F(0),F(1));twice=add(R,R,-1,-2)
        self.assertEqual(twice,(2,1))
        self.assertEqual(add(twice,R,-1,-2),(-1,-1))

    def test_S3_anchor_dependence_independent_group_law(self):
        left=add((0,3),(-1,5),7,-10)
        twice=add((2,5),(2,5),7,-10)
        self.assertEqual(left,(twice[0],-twice[1]))

    def test_retained_finite_fingerprints_from_integer_arithmetic(self):
        for row in r.read(completed.OUTPUT)['rows']:
            anchor=row['anchor'];m=F(anchor['m']);c=row['c']
            for record in anchor['finite_fingerprints']:
                p=record['prime']
                roots=[x for x in range(p) if (x**3+int(m)*x*x-(int(m)+3)*x+c*c)%p==0]
                self.assertEqual(roots,record['roots'])
                expected=[]
                for raw in anchor['rational_points']:
                    x=F(raw[0]);bits=[]
                    for root in roots:
                        value=(r.mod(x,p)-root)%p
                        if not value:value=(3*root*root+2*int(m)*root-int(m)-3)%p
                        bits.append(int(pow(value,(p-1)//2,p)==p-1))
                    expected.append(r.pack(bits))
                self.assertEqual(expected,record['point_signatures'])

    def test_no_extra_rational_Galois_direction(self):
        for row in r.read(completed.OUTPUT)['rows']:
            fixed=[j for j in range(5) if all(g[j]==1 for g in row['constant_galois_diagonals'])]
            self.assertEqual(fixed,[0,1,2])
            self.assertEqual(row['base_geometric_generic_rank'],5)
            self.assertEqual(row['new_arithmetic_generic_rank'],2)


if __name__=='__main__':unittest.main()
