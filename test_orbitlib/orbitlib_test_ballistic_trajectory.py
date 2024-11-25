import numpy as np
from astropy.time import Time
from libraries.orbitlib import OrbitLib
import matplotlib.pyplot as plt
from tqdm import tqdm

# Separate validation for each phase of the ballistic trajectory

def validate_boost_phase():
    # Define parameters for the ballistic trajectory
    launch_lat = 46.303  # Launch latitude in degrees
    launch_lon = 48.092  # Launch longitude in degrees
    impact_lat = 48.25  # Impact latitude in degrees
    impact_lon = 34.9  # Impact longitude in degrees
    apogee = 600  # Apogee in km
    epoch = Time("2025-01-01 12:00:00")
    mass = 500.0  # Mass of the missile in kg

    # Create ballistic trajectory object using OrbitLib
    missile_trajectory = OrbitLib.from_ballistic(
        name="Missile",
        launch_lat=launch_lat,
        launch_lon=launch_lon,
        impact_lat=impact_lat,
        impact_lon=impact_lon,
        apogee=apogee,
        epoch=epoch,
        mass=mass
    )

    # Set time span for boost phase
    t_boost_span = (0, 120)  # Boost phase duration in seconds

    # Progress bar setup
    with tqdm(total=t_boost_span[1] - t_boost_span[0], desc="Boost Phase Progress") as pbar:
        def boost_callback(t, y):
            pbar.update(1)

        # Propagate boost phase
        print("Starting boost phase propagation...")
        boost_phase_solution = missile_trajectory.propagate_boost(t_boost_span, callback=boost_callback)
        print("Boost phase propagation completed.")

    # Plot the boost phase
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(boost_phase_solution[0], boost_phase_solution[1], boost_phase_solution[2], label="Boost Phase", color="red")
    ax.set_title("3D Boost Phase Trajectory")
    ax.set_xlabel("X Position (km)")
    ax.set_ylabel("Y Position (km)")
    ax.set_zlabel("Z Position (km)")
    ax.legend()
    plt.show()


def validate_ballistic_phase():
    # Define parameters for the ballistic trajectory
    launch_lat = 46.303  # Launch latitude in degrees
    launch_lon = 48.092  # Launch longitude in degrees
    impact_lat = 48.25  # Impact latitude in degrees
    impact_lon = 34.9  # Impact longitude in degrees
    apogee = 600  # Apogee in km
    epoch = Time("2025-01-01 12:00:00")
    mass = 500.0  # Mass of the missile in kg

    # Create ballistic trajectory object using OrbitLib
    missile_trajectory = OrbitLib.from_ballistic(
        name="Missile",
        launch_lat=launch_lat,
        launch_lon=launch_lon,
        impact_lat=impact_lat,
        impact_lon=impact_lon,
        apogee=apogee,
        epoch=epoch,
        mass=mass
    )

    # Set time span for ballistic phase
    t_ballistic_span = (120, 600)  # Ballistic phase duration in seconds

    # Progress bar setup
    with tqdm(total=t_ballistic_span[1] - t_ballistic_span[0], desc="Ballistic Phase Progress") as pbar:
        def ballistic_callback(t, y):
            pbar.update(1)

        # Propagate ballistic phase
        print("Starting ballistic phase propagation...")
        ballistic_phase_solution = missile_trajectory.propagate_ballistic(t_ballistic_span, callback=ballistic_callback)
        print("Ballistic phase propagation completed.")

    # Plot the ballistic phase
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(ballistic_phase_solution[0], ballistic_phase_solution[1], ballistic_phase_solution[2], label="Ballistic Phase", color="green")
    ax.set_title("3D Ballistic Phase Trajectory")
    ax.set_xlabel("X Position (km)")
    ax.set_ylabel("Y Position (km)")
    ax.set_zlabel("Z Position (km)")
    ax.legend()
    plt.show()


def validate_reentry_phase():
    # Define parameters for the ballistic trajectory
    launch_lat = 46.303  # Launch latitude in degrees
    launch_lon = 48.092  # Launch longitude in degrees
    impact_lat = 48.25  # Impact latitude in degrees
    impact_lon = 34.9  # Impact longitude in degrees
    apogee = 600  # Apogee in km
    epoch = Time("2025-01-01 12:00:00")
    mass = 500.0  # Mass of the missile in kg

    # Create ballistic trajectory object using OrbitLib
    missile_trajectory = OrbitLib.from_ballistic(
        name="Missile",
        launch_lat=launch_lat,
        launch_lon=launch_lon,
        impact_lat=impact_lat,
        impact_lon=impact_lon,
        apogee=apogee,
        epoch=epoch,
        mass=mass
    )

    # Set time span for reentry phase
    t_reentry_span = (600, 900)  # Reentry phase duration in seconds

    # Progress bar setup
    with tqdm(total=t_reentry_span[1] - t_reentry_span[0], desc="Reentry Phase Progress") as pbar:
        def reentry_callback(t, y):
            pbar.update(1)

        # Propagate reentry phase
        print("Starting reentry phase propagation...")
        reentry_phase_solution = missile_trajectory.propagate_reentry(t_reentry_span, callback=reentry_callback)
        print("Reentry phase propagation completed.")

    # Plot the reentry phase
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(reentry_phase_solution[0], reentry_phase_solution[1], reentry_phase_solution[2], label="Reentry Phase", color="blue")
    ax.set_title("3D Reentry Phase Trajectory")
    ax.set_xlabel("X Position (km)")
    ax.set_ylabel("Y Position (km)")
    ax.set_zlabel("Z Position (km)")
    ax.legend()
    plt.show()


if __name__ == "__main__":
    validate_boost_phase()
    validate_ballistic_phase()
    validate_reentry_phase()
