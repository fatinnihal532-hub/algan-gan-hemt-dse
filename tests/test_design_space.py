import unittest
from src.design_space import run_sweep, pareto_front, ron_bv_pareto_front


class TestDesignSpace(unittest.TestCase):
    def test_sweep_produces_expected_row_count(self):
        rows = run_sweep(
            gate_lengths_m=[0.15e-6, 0.5e-6],
            gate_drain_spacings_m=[1e-6, 2e-6, 4e-6],
            thicknesses_m=[15e-9, 25e-9],
            al_fractions=[0.2, 0.3],
        )
        self.assertEqual(len(rows), 2 * 3 * 2 * 2)

    def test_pareto_front_is_subset_of_input(self):
        rows = run_sweep(
            gate_lengths_m=[0.15e-6, 0.25e-6, 0.5e-6],
            gate_drain_spacings_m=[1e-6, 2e-6, 4e-6, 8e-6],
            thicknesses_m=[15e-9, 20e-9],
            al_fractions=[0.2, 0.3],
        )
        front = ron_bv_pareto_front(rows)
        self.assertGreaterEqual(len(front), 1)
        self.assertLessEqual(len(front), len(rows))

    def test_pareto_front_hand_checked_example(self):
        # (1,5), (3,3) and (5,1) each beat every other point on at least
        # one objective, so none is dominated. (2,2) is dominated by (3,3),
        # which is at least as good on both and strictly better on both.
        rows = [
            {"a": 1, "b": 5},
            {"a": 3, "b": 3},
            {"a": 5, "b": 1},
            {"a": 2, "b": 2},
        ]
        front = pareto_front(rows, "a", "b")
        pairs = sorted((r["a"], r["b"]) for r in front)
        self.assertEqual(pairs, [(1, 5), (3, 3), (5, 1)])


if __name__ == "__main__":
    unittest.main()
