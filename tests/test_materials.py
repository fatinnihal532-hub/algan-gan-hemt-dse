import unittest
from src.materials import AlGaN


class TestAlGaN(unittest.TestCase):
    def test_endpoints_recover_binary_values(self):
        gan = AlGaN(0.0)
        aln = AlGaN(1.0)
        self.assertAlmostEqual(gan.bandgap_eV, 3.42, places=6)
        self.assertAlmostEqual(aln.bandgap_eV, 6.20, places=6)
        self.assertAlmostEqual(gan.eps_r, 8.9, places=6)
        self.assertAlmostEqual(aln.eps_r, 8.5, places=6)

    def test_bandgap_bowing_is_downward(self):
        # Bowing should pull the mid-composition bandgap below the naive
        # linear interpolation of the endpoints.
        mid = AlGaN(0.5)
        linear = 0.5 * 3.42 + 0.5 * 6.20
        self.assertLess(mid.bandgap_eV, linear)

    def test_polarization_charge_increases_with_al_fraction(self):
        low = AlGaN(0.10).total_polarization_charge_Cm2
        high = AlGaN(0.40).total_polarization_charge_Cm2
        self.assertGreater(high, low)

    def test_delta_ec_is_positive_and_increases_with_x(self):
        low = AlGaN(0.10).delta_ec_eV
        high = AlGaN(0.40).delta_ec_eV
        self.assertGreater(low, 0.0)
        self.assertGreater(high, low)

    def test_rejects_out_of_range_fraction(self):
        with self.assertRaises(ValueError):
            AlGaN(1.5)
        with self.assertRaises(ValueError):
            AlGaN(-0.1)


if __name__ == "__main__":
    unittest.main()
