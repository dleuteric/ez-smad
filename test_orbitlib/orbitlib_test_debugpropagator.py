# Test script for ballistic trajectory phases
from libraries.orbitlib import OrbitLib
from libraries.mathlib import MathLib
from astropy.time import Time
import numpy as np
import matplotlib.pyplot as plt

# Example parameters
launch_lat = 28.5  # degrees
launch_lon = -80.6  # degrees
apogee = 300e3  # Apogee in meters
epoch = Time("2023-01-01 00:00:00")
mass = 50000  # kg
thrust = [2e6, 0, 0]  # Thrust vector in Newtons
I_sp = 300  # Specific impulse in seconds
C_D = 0.3  # Drag coefficient
S = 10.0  # Cross-sectional area in m^2

# Convert launch coordinates to Cartesian
x0, y0, z0 = MathLib.latlon_to_cartesian(launch_lat, launch_lon, 6371000)
initial_position = np.array([x0, y0, z0])

# Initial velocity (simplified)
initial_velocity = np.array([0, 0, 500])  # Start from rest for testing

# Combine initial conditions
initial_conditions = np.concatenate((initial_position, initial_velocity))

# Create OrbitLib instance
missile_orbit = OrbitLib(
    name="Test Missile",
    epoch=epoch,
    initial_conditions=initial_conditions,
    mass=mass,
    thrust=thrust,
    I_sp=I_sp,
    C_D=C_D,
    S=S
)

# Define time spans
t_boost = (0, 150)  # Boost phase
t_ballistic = (150, 600)  # Ballistic phase
t_reentry = (600, 750)  # Reentry phase

# Propagate each phase
boost_solution = missile_orbit.propagate_boost(t_boost)
missile_orbit.initial_conditions = boost_solution.y[:, -1]

# Propagate ballistic phase
ballistic_solution = missile_orbit.propagate_ballistic(t_ballistic)

# Update initial conditions for reentry phase (retain mass)
missile_orbit.initial_conditions = ballistic_solution.y[:, -1]

# Propagate reentry phase
reentry_solution = missile_orbit.propagate_reentry(t_reentry)

# Combine results
time_total = np.concatenate((boost_solution.t, ballistic_solution.t, reentry_solution.t))
state_total = np.hstack((boost_solution.y, ballistic_solution.y, reentry_solution.y))


# Plot results
plt.figure(figsize=(10, 6))
plt.plot(state_total[0, :], state_total[1, :], label="Trajectory")
plt.xlabel("X Position (m)")
plt.ylabel("Y Position (m)")
plt.title("Ballistic Missile Trajectory")
plt.legend()
plt.grid()
plt.show()

# Calculate altitude over time
altitudes = np.linalg.norm(state_total[:3, :], axis=0) - 6371000  # Earth radius subtracted
plt.figure(figsize=(10, 6))
plt.plot(time_total, altitudes)
plt.xlabel("Time (s)")
plt.ylabel("Altitude (m)")
plt.title("Altitude vs Time")
plt.grid()
plt.show()