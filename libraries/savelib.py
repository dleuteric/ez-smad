import json
import os

class SaveLib:
    """
    Library for saving and exporting results in different formats.
    Currently supports saving data to JSON.
    """

    @staticmethod
    def save_to_json(data, file_name, directory="/Users/daniele_leuteri/PycharmProjects/thirdeye/results"):
        """
        Save data to a JSON file.

        :param data: Dictionary or list to save to a JSON file.
        :param file_name: The name of the JSON file (without extension).
        :param directory: Directory where the JSON file should be saved. Defaults to 'results'.
        """
        # Ensure the directory exists
        if not os.path.exists(directory):
            os.makedirs(directory)

        # File path
        file_path = os.path.join(directory, f"{file_name}.json")

        # Write the data to JSON
        with open(file_path, "w") as json_file:
            json.dump(data, json_file, indent=4)

        print(f"Data saved to {file_path}")

    @staticmethod
    def load_from_json(file_path):
        """
        Load data from a JSON file.

        :param file_path: Path to the JSON file to be loaded.
        :return: Data loaded from the JSON file (usually a dictionary or list).
        """
        # Check if file exists
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"The file {file_path} does not exist.")

        # Load the data from JSON
        with open(file_path, "r") as json_file:
            data = json.load(json_file)

        return data

# Example usage:
# SaveLib.save_to_json(wavelength_band_detection_ranges, "detection_analysis")
# loaded_data = SaveLib.load_from_json("results/detection_analysis.json")
