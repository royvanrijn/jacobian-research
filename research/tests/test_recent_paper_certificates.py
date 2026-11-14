"""Fail-closed controls for the newly adopted finite certificates."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
PILOT = ROOT/"artifacts/generated-results/hill1152-transfer-20261113-v2"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT/path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


X = load("recent_xns121", "elliptic-curves/cas/verify_xns121_localization.py")
H = load("recent_hill1152", "elkies-k3/scripts/verify_hill1152_transfer.py")
E = load("recent_icarm4239", "elliptic-curves/cas/verify_icarm_curve4239_rank32.py")


class RecentCertificateTamperTests(unittest.TestCase):
    def test_local_formal_direction_in_selmer_span_is_rejected(self):
        data = X.parse_assignments((X.DATA/"even_selmer_0.m").read_text())
        data["FormalRows"][0] = data["SelmerLocalRows"][0]
        with self.assertRaises(ValueError):
            X.check_even(data, 5)

    def test_data_parser_does_not_execute_source(self):
        with self.assertRaises((ValueError, SyntaxError)):
            X.parse_assignments("Rows := __import__('os').getcwd();")

    def tampered_pilot(self, filename, mutate):
        with tempfile.TemporaryDirectory(prefix="hill1152-tamper-") as tmp:
            tmp = Path(tmp)
            for name in ("input.json", "pencils.json", "result.json", "mestre.json"):
                data = json.loads((PILOT/name).read_text())
                if name == filename:
                    mutate(data)
                # The input binding is exact bytes; preserve its original bytes.
                if name == "input.json":
                    (tmp/name).write_bytes((PILOT/name).read_bytes())
                else:
                    (tmp/name).write_text(json.dumps(data))
            with self.assertRaises(ValueError):
                H.verify(tmp)

    def test_branch_matrix_tamper_is_rejected(self):
        def mutate(data):
            data[0]["branch_matrix"][0][0] = str(int(data[0]["branch_matrix"][0][0].split("/")[0])+1)
        self.tampered_pilot("pencils.json", mutate)

    def test_false_conic_witness_is_rejected(self):
        self.tampered_pilot("result.json",
            lambda d: d["quadratic_covers"][0].update(conic_point=["0", "1", "0"]))

    def test_false_mestre_local_root_is_rejected(self):
        self.tampered_pilot("mestre.json",
            lambda d: d["fields"]["A"]["local_place"].update(root=0))

    def test_rank32_finite_row_tamper_is_rejected(self):
        data = json.loads(E.DEFAULT.read_text())
        data["signatures"][0]["rows"][0][0] ^= 1
        with self.assertRaises(AssertionError):
            E.verify(data)


if __name__ == "__main__":
    unittest.main()
