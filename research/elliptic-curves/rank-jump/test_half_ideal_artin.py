import unittest

import retrospective as r
import half_ideal_artin as art
import half_ideal_artin_completion as complete
import half_ideal_class_blocks as blocks
from verify_half_ideal_artin import jacobi


class HalfIdealArtin(unittest.TestCase):
    def test_integer_jacobi_against_prime_factor_definition(self):
        for n in range(1,200,2):
            factors=[];m=n;p=3
            while p*p<=m:
                exponent=0
                while m%p==0:
                    m//=p;exponent+=1
                if exponent:factors.append((p,exponent))
                p+=2
            if m>1:factors.append((m,1))
            for a in range(-6,15):
                expected=1
                for p,e in factors:
                    v=pow(a%p,(p-1)//2,p)
                    expected*=(-1 if v==p-1 else v)**e
                self.assertEqual(jacobi(a,n),expected)

    def test_completion_changes_exactly_the_frozen_missing_entries(self):
        raw={x["case_index"]:x for x in r.read(art.INPUT)["cases"]}
        repairs={x["case_index"]:x for x in r.read(complete.INPUT)["cases"]}
        rows=r.read(complete.OUTPUT)["rows"];missing=0
        for row in rows:
            index=row["case_index"]
            known={(x["column"],x["character"]):x for x in repairs[index]["repairs"]}
            for j,column in enumerate(raw[index]["columns"]):
                for i,old in enumerate(column["evaluations"]):
                    if "artin_bit" in old:
                        self.assertNotIn((j,i),known)
                        bit=old["artin_bit"]
                    else:
                        missing+=1
                        bit=known[j,i]["artin_bit"]
                    self.assertEqual(row["matrix_rows"][i][j],bit)
        self.assertEqual(missing,19)

    def test_selected_factors_contain_the_generic_half_ideal_images(self):
        rows=r.read(blocks.OUTPUT)["rows"]
        self.assertEqual([x["split_elementary_S_class_factor_dimension"] for x in rows],[9,8,6])
        self.assertEqual([x["relative_dimension_inside_selected_factor"] for x in rows],[8,6,0])
        for row in rows:
            g=row["generic_strict_dimension"]
            self.assertEqual(row["selected_half_ideal_indices"][:g],list(range(g)))
            self.assertIn("epsilon",row["remaining_S_class_two_rank"])

    def test_dual_words_split_selected_characters(self):
        matrices={x["case_index"]:x["matrix_rows"] for x in r.read(complete.OUTPUT)["rows"]}
        for row in r.read(blocks.OUTPUT)["rows"]:
            A=matrices[row["case_index"]]
            for k,word in enumerate(row["dual_half_ideal_words"]):
                actual=[sum(A[i][j] for j in range(len(A)) if word>>j&1)%2
                        for i in row["selected_character_indices"]]
                self.assertEqual(actual,[int(i==k) for i in range(len(actual))])

    def test_distinct_left_and_right_kernels_are_not_conflated(self):
        row=r.read(blocks.OUTPUT)["rows"][0]
        A=r.read(complete.OUTPUT)["rows"][0]["matrix_rows"]
        left=row["left_character_kernel_masks"][0]
        right=row["right_half_ideal_kernel_masks"][0]
        self.assertEqual((left,right),(896,637))
        self.assertNotEqual(left,right)
        self.assertTrue(all(sum(A[i][j] for i in range(10) if left>>i&1)%2==0 for j in range(10)))
        self.assertTrue(all(sum(A[i][j] for j in range(10) if right>>j&1)%2==0 for i in range(10)))

    def test_invisible_two_torsion_need_not_be_zero(self):
        # In C=Z/4 x Z/2 the nonzero element (2,0) is killed by every quadratic character.
        elements=[(x,y) for x in range(4) for y in range(2)]
        h=(2,0)
        self.assertNotEqual(h,(0,0))
        self.assertEqual((2*h[0]%4,2*h[1]%2),(0,0))
        self.assertTrue(all((a*h[0]+b*h[1])%2==0 for a in range(2) for b in range(2)))
        self.assertTrue(any((2*x%4,2*y%2)==h for x,y in elements))


if __name__=="__main__":
    unittest.main()
