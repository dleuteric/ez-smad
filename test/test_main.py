# tests/test_main.py

import unittest
from libraries.physicslib import BlackbodyRadiation, DetectionPhysics
from libraries.mathlib_0 import NumericalIntegration
import main


class TestMainFunctions(unittest.TestCase):

    def test_input_validation(self):
        # Test with invalid inputs
        invalid_inputs = {
            'temperature': -100,  # Invalid temperature
            'emissivity': 1.5,  # Invalid emissivity
            'aperture_size': -0.1,  # Invalid aperture size
            'wavelength': 0,  # Invalid wavelength
            'snr': -5,  # Invalid SNR
            'nefd': -1e-12  # Invalid NEFD
        }

        result = main.validate_inputs(invalid_inputs)
        self.assertFalse(result)  # Should return False for invalid inputs

    def test_input_validation_valid(self):
        # Test with valid inputs
        valid_inputs = {
            'temperature': 500,  # Valid temperature
            'emissivity': 0.9,  # Valid emissivity
            'aperture_size': 0.3,  # Valid aperture size
            'wavelength': 4e-6,  # Valid wavelength
            'snr': 5,  # Valid SNR
            'nefd': 1e-12  # Valid NEFD
        }

        result = main.validate_inputs(valid_inputs)
        self.assertTrue(result)  # Should return True for valid inputs


if __name__ == '__main__':
    unittest.main()
