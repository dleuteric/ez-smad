import numpy as np
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from astropy import units as u
from astropy.time import Time
from poliastro.twobody.propagation import FarnocchiaPropagator, CowellPropagator
from scipy.integrate import solve_ivp

from libraries.mathlib import MathLib
from libraries.physicslib import PhysicsLib


class OrbitLib:
    def __init__(self, name, initial_orbit=None, epoch=None, frame="inertial",
                 initial_conditions=None, mass=None, thrust=None, I_sp=None,
                 C_D=None, S=None):
        self.name = name
        self.initial_orbit = initial_orbit
        self.epoch = epoch
        self.frame = frame
        self.propagator = FarnocchiaPropagator() if initial_orbit else None
        self.initial_conditions = initial_conditions
        self.mass = mass
        self.thrust = thrust  # Thrust vector [Tx, Ty, Tz]
        self.I_sp = I_sp      # Specific impulse
        self.C_D = C_D        # Drag coefficient
        self.S = S            # Cross-sectional area


    @classmethod
    def from_classical(cls, name, attractor, a, ecc, inc, raan, argp, nu, epoch=Time("J2000"), plane="inertial"):
        """
        Create an orbit from classical orbital elements.
        """
        orbit = Orbit.from_classical(attractor, a * u.km, ecc * u.one, inc * u.deg,
                                     raan * u.deg, argp * u.deg, nu * u.deg, epoch)
        return cls(name, orbit, epoch, plane)

    @classmethod
    def from_vectors(cls, name, attractor, r, v, epoch=Time("J2000"), plane="inertial"):
        """
        Create an orbit from position and velocity vectors.
        """
        orbit = Orbit.from_vectors(attractor, r * u.km, v * (u.km / u.s), epoch)
        return cls(name, orbit, epoch, plane)

    @classmethod
    def from_ballistic(cls, name, launch_lat, launch_lon, impact_lat, impact_lon, apogee, epoch, mass):
        """
        Create a ballistic trajectory from launch and impact points.
        """
        # Convert lat/lon to Cartesian for the launch and impact points
        radius_earth = 6371.0  # Earth radius in km
        x_launch, y_launch, z_launch = MathLib.latlon_to_cartesian(launch_lat, launch_lon, radius_earth)
        x_impact, y_impact, z_impact = MathLib.latlon_to_cartesian(impact_lat, impact_lon, radius_earth)

        # Initial position is the launch point
        initial_position = np.array([x_launch, y_launch, z_launch])  # in km

        # Estimate an initial velocity vector towards apogee (simplified)
        velocity_magnitude = np.sqrt(2 * 9.81 * (apogee * 1000)) / 1000  # Approx initial velocity in km/s
        initial_velocity = np.array([0, velocity_magnitude, 0])  # Placeholder velocity vector in km/s

        initial_conditions = np.concatenate((initial_position, initial_velocity))

        return cls(name=name, epoch=epoch, frame="inertial", initial_conditions=initial_conditions, mass=mass)

    def set_propagator(self, propagator_type="farnocchia"):
        """
        Set the propagator type for orbit propagation.
        """
        if propagator_type.lower() == "farnocchia":
            self.propagator = FarnocchiaPropagator()
        elif propagator_type.lower() == "cowell":
            self.propagator = CowellPropagator()
        else:
            raise ValueError("Unsupported propagator type. Choose from: 'farnocchia', 'cowell'.")

    def propagate(self, time_delta):
        """
        Propagate the orbit by a given time delta.
        """
        return self.initial_orbit.propagate(time_delta * u.s, method=self.propagator)

    def get_position(self, epoch):
        """
        Get the position of the orbit at a specific epoch.
        """
        propagated_orbit = self.initial_orbit.propagate(epoch - self.epoch, method=self.propagator)
        return propagated_orbit.r

    def to_string(self):
        """
        Get a string representation of the orbit.
        """
        return f"Orbit: {self.name}, Epoch: {self.epoch}, Frame: {self.frame}"

    def boost_equations_of_motion(self, t, y):
        x, y_pos, z, vx, vy, vz, m = y
        r = np.array([x, y_pos, z])
        v = np.array([vx, vy, vz])

        # Compute accelerations
        a_r = PhysicsLib.thrust_acceleration(self.thrust, m)
        h = np.linalg.norm(r) - 6371000  # Altitude in meters
        a_d = PhysicsLib.drag_acceleration(v, h, self.C_D, self.S, m)
        g = PhysicsLib.gravity_acceleration(r)
        a_ni = PhysicsLib.non_inertial_acceleration(v, r)

        # Total acceleration
        a_total = a_r + a_d + g + a_ni

        # Debug prints
        print(f"Time: {t}, Altitude: {h}, a_r: {a_r}, a_d: {a_d}, g: {g}, a_total: {a_total}, m: {m}")

        # Mass rate change
        m_dot = PhysicsLib.mass_rate_change(self.thrust, self.I_sp)

        return [vx, vy, vz, a_total[0], a_total[1], a_total[2], m_dot]

    def ballistic_equations_of_motion(self, t, y):
        """
        Equations of motion for the ballistic phase.
        """
        x, y_pos, z, vx, vy, vz, m = y  # Include mass (constant during this phase)
        r = np.array([x, y_pos, z])
        v = np.array([vx, vy, vz])

        # Compute accelerations
        g = PhysicsLib.gravity_acceleration(r)
        a_ni = PhysicsLib.non_inertial_acceleration(v, r)

        # Total acceleration
        a_total = g + a_ni

        return [vx, vy, vz, a_total[0], a_total[1], a_total[2], 0]  # Mass rate is zero

    def reentry_equations_of_motion(self, t, y):
        """
        Equations of motion for the reentry phase.
        """
        x, y_pos, z, vx, vy, vz, m = y  # Include mass (constant during this phase)
        r = np.array([x, y_pos, z])
        v = np.array([vx, vy, vz])

        # Compute accelerations
        h = np.linalg.norm(r) - 6371000  # Altitude in meters
        a_d = PhysicsLib.drag_acceleration(v, h, self.C_D, self.S, m)
        g = PhysicsLib.gravity_acceleration(r)
        a_ni = PhysicsLib.non_inertial_acceleration(v, r)

        # Total acceleration
        a_total = a_d + g + a_ni

        return [vx, vy, vz, a_total[0], a_total[1], a_total[2], 0]  # Mass rate is zero

    def propagate_boost(self, t_span, callback=None):
        """
        Propagate the boost phase.
        """
        if self.initial_conditions is None or self.mass is None:
            raise ValueError("Initial conditions and mass must be provided.")

        y0 = np.concatenate((self.initial_conditions, [self.mass]))

        def boost_wrapper(t, y):
            return self.boost_equations_of_motion(t, y)

        solution = solve_ivp(boost_wrapper, t_span, y0, method='RK45', max_step=1.0)
        return solution

    def propagate_ballistic(self, t_span, callback=None):
        """
        Propagate the ballistic phase.
        """
        if self.initial_conditions is None:
            raise ValueError("Initial conditions must be provided.")

        def ballistic_wrapper(t, y):
            return self.ballistic_equations_of_motion(t, y)

        solution = solve_ivp(ballistic_wrapper, t_span, self.initial_conditions, method='RK45')
        return solution

    def propagate_reentry(self, t_span, callback=None):
        """
        Propagate the reentry phase.
        """
        if self.initial_conditions is None:
            raise ValueError("Initial conditions must be provided.")

        def reentry_wrapper(t, y):
            return self.reentry_equations_of_motion(t, y)

        solution = solve_ivp(reentry_wrapper, t_span, self.initial_conditions, method='RK45')
        return solution