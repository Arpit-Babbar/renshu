coordinates_min = -1.0 # minimum coordinate
coordinates_max = 1.0  # maximum coordinate

initial_condition_sine_wave(x) = 1.0 + 0.5 * sin(pi * x)

n_elements = 16 # number of elements

dx = (coordinates_max - coordinates_min) / n_elements # length of one element



