# Test functions to verify the propagation functionalities of OrbitLib
from poliastro.bodies import Earth
from astropy import units as u
from astropy.time import Time
from libraries.orbitlib import OrbitLib


def test_propagate_orbit():
    """
    Test the propagation function of OrbitLib.
    """
    # Create an orbit using classical orbital elements
    test_orbit = OrbitLib.from_classical(
        name="Test Orbit",
        attractor=Earth,
        a=7000,  # Semi-major axis in km
        ecc=0.01,  # Eccentricity
        inc=98.5,  # Inclination in degrees
        raan=25.0,  # RAAN in degrees
        argp=40.0,  # Argument of pericenter in degrees
        nu=0.0,  # True anomaly in degrees
        epoch=Time("J2000")
    )

    # Propagate the orbit by 3600 seconds (1 hour)
    propagated_orbit = test_orbit.propagate(3600)

    # Print the propagated position and velocity
    print(f"Propagated position: {propagated_orbit.r}")
    print(f"Propagated velocity: {propagated_orbit.v}")

    # Validate by checking against expected values (if known)
    # This part would involve asserting values for automated testing
    print("Propagation test complete.")


# Example usage of the test function
if __name__ == "__main__":
    test_propagate_orbit()
