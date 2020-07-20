import pandas as pd

def __get_neighbors(val, array):
    # returns array index value of neighbors
    if val < array[0]:
        return (0, array[1])
    elif val > array[-1]:
        return (array[-2], array[-1])
    for index, i in enumerate(array):
        if index < len(array) + 1:
            if i < val <= array[index + 1]:
                return (index + 1, index)

def interpolate1d(val, val_array, interp_array):
    sorted_interp_array = [x for _, x in sorted(zip(val_array, interp_array))]
    sorted_val_array = list(reversed(val_array))
    print()

    neighbor_indices = __get_neighbors(val=val, array=sorted_val_array)
    val_neighbors = (sorted_val_array[neighbor_indices[0]], sorted_val_array[neighbor_indices[1]])
    interp_neighbors = (sorted_interp_array[neighbor_indices[0]], sorted_interp_array[neighbor_indices[1]])
    print(val, val_neighbors)
    print(interp_neighbors)
    return (((val_neighbors[1] - val) / (val_neighbors[1] - val_neighbors[0])) * interp_neighbors[0]) + \
           (((val_neighbors[0] - val) / (val_neighbors[1] - val_neighbors[0])) * interp_neighbors[1])

phase_df = pd.read_fwf("./phase_data/duniteS_vapour_curve.txt", skiprows=1,
                                        names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                               "entropy_sol_liq", "entropy_vap"])

i = interpolate1d(val=7.9993E+03, val_array=phase_df["temperature"], interp_array=phase_df["entropy_vap"])
print(i)
