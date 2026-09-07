"""Checks of quotient invariance and fail-closed boundary prerequisites."""
import copy
import unittest
import affine_selmer as af
import retrospective as r
from selmer_comparison import OUTPUT, analyze, quotient_rows

try:
    from sage.all import GF, VectorSpace, matrix
except ImportError:
    VectorSpace = None


class SelmerComparisonTests(unittest.TestCase):
    def test_quotient_is_linear_and_has_exact_kernel(self):
        # Include a non-reduced basis and two different generating sets.
        for generators in ([0b0011, 0b0110], [0b0101, 0b0110]):
            span = {0, generators[0], generators[1],
                    generators[0] ^ generators[1]}
            images = quotient_rows(range(16), generators)
            self.assertEqual({i for i, v in enumerate(images) if v == 0}, span)
            for x in range(16):
                for y in range(16):
                    self.assertEqual(images[x ^ y], images[x] ^ images[y])

    def test_reject_incomplete_local_point_image(self):
        case = copy.deepcopy(r.read(af.INPUT)["cases"][0])
        local = next(x for x in case["local"] if x["point_dimension"])
        local["point_signature_rows"].pop()
        with self.assertRaises(AssertionError):
            analyze(case)

    def test_reject_lost_transporter_ramification(self):
        case = copy.deepcopy(r.read(af.INPUT)["cases"][0])
        for local in case["local"]:
            local["class_signature_rows"][20] = [
                0 for _ in local["class_signature_rows"][20]]
        with self.assertRaises(AssertionError):
            analyze(case)

    @unittest.skipIf(VectorSpace is None, "optional independent Sage algebra")
    def test_independent_sage_quotients_and_strict_kernels(self):
        report = r.read(OUTPUT)
        for case in r.read(af.INPUT)["cases"]:
            rows = [[] for _ in range(21)]
            half_dimension = 0
            for local in case["local"]:
                V = VectorSpace(GF(2), len(local["class_signature_rows"][0]))
                old = V.subspace(local["class_signature_rows"][:20])
                new = V.subspace(local["point_signature_rows"])
                common, total = old.intersection(new), old + new
                half_dimension += old.dimension() - common.dimension()
                Q = V.quotient(common)
                for i, signature in enumerate(local["class_signature_rows"]):
                    self.assertIn(V(signature), total)
                    rows[i].extend(Q(V(signature)))
            M = matrix(GF(2), rows)
            result = next(x for x in report["cases"] if x["u"] == case["u"])
            self.assertEqual(M.rank(), half_dimension)
            self.assertEqual(M.rank(), result["complete_global_boundary_dimension"])
            kernel = M.left_kernel()
            self.assertEqual(kernel.dimension(), len(result["strict_kernel_anchor_masks"]))
            for mask in result["strict_kernel_anchor_masks"]:
                self.assertIn(kernel.ambient_vector_space()(
                    [(mask >> j) & 1 for j in range(21)]), kernel)


if __name__ == "__main__":
    unittest.main()
