# Imports
import math
from os.path import Path


class Threat:
    """
    Class representing a threat object with its properties.
    """

    def __init__(self, name: str, temperatureK: float, emissivity: float, area: float, trajectory: str):
        """
        Initialize a threat object.

        :param name: Name of the threat object
        :param temperatureK: Temperature in Kelvin [K]
        :param emissivity: Emissivity (dimensionless)
        :param area: Emitting area in square meters [m^2]
        :param trajectory: Trajectory description (e.g., 'LEO', 'Suborbital')
        """
        self.name = name
        self.temperatureK = temperatureK
        self.emissivity = emissivity
        self.area = area
        self.trajectory = trajectory