import numpy as np

class NumericalIntegration:

    @staticmethod
    def simpson_rule(y, x):
        """
        Perform numerical integration using Simpson's rule.

        :param y: Array of function values at sample points
        :param x: Array of sample points
        :return: Approximate integral
        """
        if len(x) % 2 == 0:
            raise ValueError("Number of samples must be odd for Simpson's rule.")
        h = (x[-1] - x[0]) / (len(x) - 1)
        S = y[0] + y[-1] + 4 * np.sum(y[1:-1:2]) + 2 * np.sum(y[2:-2:2])

        return h * S / 3

class Interpolation:
    @staticmethod
    def linear_interpolate(x, xp, yp):
        """
        Perform linear interpolation.

        :param x: The x-coordinate at which to interpolate
        :param xp: Array of x-coordinates of the data points
        :param yp: Array of y-coordinates of the data points
        :return: Interpolated y-coordinate
        """
        return np.interp(x, xp, yp)

class Statistics:
    @staticmethod
    def mean(data):
        return sum(data) / len(data)

    @staticmethod
    def variance(data):
        mu = Statistics.mean(data)
        return sum((x - mu) ** 2 for x in data) / len(data)

class RootFinding:
    @staticmethod
    def bisection_method(f, a, b, tol=1e-6, max_iter=100):
        """
        Find a root of function f in interval [a, b] using the bisection method.

        :param f: Function for which the root is sought
        :param a: Lower bound of the interval
        :param b: Upper bound of the interval
        :param tol: Tolerance for convergence
        :param max_iter: Maximum number of iterations
        :return: Approximate root value
        """
        for _ in range(max_iter):
            c = (a + b) / 2
            if f(c) == 0 or (b - a) / 2 < tol:
                return c
            if np.sign(f(c)) == np.sign(f(a)):
                a = c
            else:
                b = c
        raise ValueError("Root not found within the given tolerance and maximum iterations.")

class VectorOperations:
    @staticmethod
    def add(v1, v2):
        """
        Add two vectors.

        :param v1: First vector (array-like)
        :param v2: Second vector (array-like)
        :return: Sum of vectors v1 and v2
        """
        return np.add(v1, v2)

    @staticmethod
    def subtract(v1, v2):
        """
        Subtract one vector from another.

        :param v1: First vector (array-like)
        :param v2: Second vector (array-like)
        :return: Difference of vectors v1 and v2
        """
        return np.subtract(v1, v2)

    @staticmethod
    def dot_product(v1, v2):
        """
        Calculate the dot product of two vectors.

        :param v1: First vector (array-like)
        :param v2: Second vector (array-like)
        :return: Dot product of v1 and v2
        """
        return np.dot(v1, v2)

    @staticmethod
    def cross_product(v1, v2):
        """
        Calculate the cross product of two vectors.

        :param v1: First vector (array-like)
        :param v2: Second vector (array-like)
        :return: Cross product of v1 and v2
        """
        return np.cross(v1, v2)

    @staticmethod
    def magnitude(v):
        """
        Calculate the magnitude of a vector.

        :param v: Vector (array-like)
        :return: Magnitude of vector v
        """
        return np.linalg.norm(v)

    @staticmethod
    def distance(v1, v2):
        """
        Calculate the Euclidean distance between two vectors.

        :param v1: First vector (array-like)
        :param v2: Second vector (array-like)
        :return: Euclidean distance between v1 and v2
        """
        return np.linalg.norm(np.subtract(v1, v2))

class Trigonometry:
    @staticmethod
    def sin(angle_deg):
        """
        Calculate the sine of an angle in degrees.

        :param angle_deg: Angle in degrees
        :return: Sine of the angle
        """
        return np.sin(np.radians(angle_deg))

    @staticmethod
    def cos(angle_deg):
        """
        Calculate the cosine of an angle in degrees.

        :param angle_deg: Angle in degrees
        :return: Cosine of the angle
        """
        return np.cos(np.radians(angle_deg))

    @staticmethod
    def tan(angle_deg):
        """
        Calculate the tangent of an angle in degrees.

        :param angle_deg: Angle in degrees
        :return: Tangent of the angle
        """
        return np.tan(np.radians(angle_deg))

    @staticmethod
    def atan(value):
        """
        Calculate the arctangent of a value and return the angle in degrees.

        :param value: Value to calculate the arctangent for
        :return: Arctangent of the value in degrees
        """
        return np.degrees(np.arctan(value))


class MathLib:
    @staticmethod
    def latlon_to_cartesian(lat, lon, radius):
        """
        Convert latitude and longitude to Cartesian coordinates.

        Parameters:
        - lat: Latitude in degrees
        - lon: Longitude in degrees
        - radius: Radius of the body (e.g., Earth)

        Returns:
        - Cartesian coordinates as a tuple (x, y, z) in km
        """
        lat = np.radians(lat)
        lon = np.radians(lon)

        x = radius * np.cos(lat) * np.cos(lon)
        y = radius * np.cos(lat) * np.sin(lon)
        z = radius * np.sin(lat)

        return x, y, z

    @staticmethod
    def cartesian_to_latlon(x, y, z):
        """
        Convert Cartesian coordinates to latitude and longitude.

        Parameters:
        - x, y, z: Cartesian coordinates in km

        Returns:
        - Latitude and longitude in degrees as a tuple (lat, lon)
        """
        radius = np.sqrt(x ** 2 + y ** 2 + z ** 2)
        lat = np.degrees(np.arcsin(z / radius))
        lon = np.degrees(np.arctan2(y, x))

        return lat, lon
