# Test functions to verify the plot functionalities of OrbitPlotLib
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from astropy.time import Time
from astropy import units as u
from libraries.plotlib import OrbitPlotLib

def test_plot_orbit_2d():
    """
    Test the 2D plotting function of OrbitPlotLib.
    """
    # Create an example orbit from classical elements
    orbit = Orbit.from_classical(
        Earth,
        7000 * u.km,  # Semi-major axis
        0.01 * u.one,  # Eccentricity
        98.5 * u.deg,  # Inclination
        25.0 * u.deg,  # RAAN
        40.0 * u.deg,  # Argument of pericenter
        0.0 * u.deg,   # True anomaly
        epoch=Time("J2000")
    )
    # Call the plot function
    OrbitPlotLib.plot_orbit_2d(orbit, "Test Orbit 2D")

def test_plot_orbit_3d():
    """
    Test the 3D plotting function of OrbitPlotLib.
    """
    # Create an example orbit from classical elements
    orbit = Orbit.from_classical(
        Earth,
        7000 * u.km,  # Semi-major axis
        0.01 * u.one,  # Eccentricity
        98.5 * u.deg,  # Inclination
        25.0 * u.deg,  # RAAN
        40.0 * u.deg,  # Argument of pericenter
        0.0 * u.deg,   # True anomaly
        epoch=Time("J2000")
    )
    # Call the plot function
    OrbitPlotLib.plot_orbit_3d(orbit, "Test Orbit 3D")

def test_plot_multiple_orbits_2d():
    """
    Test the 2D plotting function for multiple orbits in OrbitPlotLib.
    """
    # Create multiple example orbits
    orbit1 = Orbit.from_classical(
        Earth,
        7000 * u.km,
        0.01 * u.one,
        98.5 * u.deg,
        25.0 * u.deg,
        40.0 * u.deg,
        0.0 * u.deg,
        epoch=Time("J2000")
    )
    orbit2 = Orbit.from_classical(
        Earth,
        7200 * u.km,
        0.02 * u.one,
        55.0 * u.deg,
        30.0 * u.deg,
        45.0 * u.deg,
        10.0 * u.deg,
        epoch=Time("J2000")
    )
    # Call the plot function for multiple orbits
    OrbitPlotLib.plot_multiple_orbits_2d([orbit1, orbit2], ["Orbit 1", "Orbit 2"])

# Example usage of the test functions
if __name__ == "__main__":
    test_plot_orbit_2d()
    test_plot_orbit_3d()
    test_plot_multiple_orbits_2d()
