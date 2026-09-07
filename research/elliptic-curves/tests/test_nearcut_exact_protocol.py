import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
from nearcut60v2_mw16_pari_batch import EXACT_ROW_FIELDS,project_rows
from research_runtime.store import digest
from audit_nearcut60v3_mw16_accounting import check_exact_roster


class ExactProtocolTests(unittest.TestCase):
    def test_observation_timings_do_not_enter_arithmetic_protocol(self):
        row={k:1 for k in EXACT_ROW_FIELDS}
        row.update(id='fixture',model=['0','0','0','-1','1'],wall_seconds=0.125,scores={'old_cpu_seconds':0.25})
        with self.assertRaises(TypeError):digest([row])
        result=project_rows([row])
        self.assertTrue(digest(result))
        self.assertEqual(result[0]['model'],row['model'])
        self.assertEqual(row['wall_seconds'],0.125)
        self.assertNotIn('wall_seconds',result[0])
        check_exact_roster({'rows':result},[row])
        result[0]['model']=['0','0','0','-2','1']
        with self.assertRaises(ArithmeticError):check_exact_roster({'rows':result},[row])

    def test_inexact_arithmetic_field_is_rejected_before_maps(self):
        row={k:1 for k in EXACT_ROW_FIELDS};row['combined_late_units']=1.5
        with self.assertRaises(TypeError):project_rows([row])
        with self.assertRaises(ArithmeticError):check_exact_roster({'rows':[row]},[row])


if __name__=='__main__':unittest.main()
