from collections import defaultdict
import matplotlib.pyplot as plt
from project import find_unit_by_bas

class VarLogger():
    """Provides utility to log the values of variables during a run."""
    logged_values = defaultdict(lambda: {})

    def __init__(self, vars_to_log, district) -> None:
        self.variables = []
        for bas_key, var_name in vars_to_log:
            self.variables.append(
                (bas_key, var_name, find_unit_by_bas(district, bas_key))
            )

    def log_vars(self, time) -> None:
        for bas_key, var_name, unit in self.variables:
            name = bas_key + "." + var_name
            self.logged_values[time][name] = getattr(unit, var_name)

    def plot(self) -> None:
        time_series = [t for t in self.logged_values]
        figure, axis = plt.subplots()

        for bas_key, var_name, _ in self.variables:
            series_name = bas_key + "." + var_name
            axis.plot(
                time_series,
                [self.logged_values[t][series_name] for t in time_series], # data series
                label=series_name
            )

        axis.set_xlabel("Time")
        axis.set_ylabel("Energy [Wh]")
        axis.set_title("Energy use in district")
        axis.legend()
        plt.show()
