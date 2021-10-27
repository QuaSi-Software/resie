import os
from json import loads
from datetime import datetime

DATE_FORMAT = "%Y-%m-%d %H:%M:%S"

def discover_profiles(root_dir: str) -> dict:
    profiles = {}

    for dirpath, _, filenames in os.walk(root_dir):
        for filename in filenames:
            filepath = os.path.join(dirpath, filename)

            if filename[-5:].lower() == ".json":
                with open(filepath, "r") as file:
                    profile = Profile(loads(file.read()))
                    profiles[profile.name] = profile

    return profiles

class Profile():
    def __init__(self, config: dict) -> None:
        self.name: str = str(config["name"])
        self.description: str = str(config["description"])

        type: str = str(config["time_series_type"])
        if type == "csv":
            raise NotImplementedError("CSV profiles are not implemented yet")
        elif type == "json":
            self.time_series_data: dict = config["time_series_data"]

    def get_value(self, time: datetime) -> float:
        timestamp = time.strftime(DATE_FORMAT)
        if timestamp not in self.time_series_data:
            return 0 # todo: log warning instead
        return self.time_series_data[timestamp]
