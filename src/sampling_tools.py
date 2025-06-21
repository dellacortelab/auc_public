import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import cumulative_trapezoid

def find_representative_points(x, y, n):
    if n < 3:
        raise ValueError("n must be at least 3 to satisfy inclusion constraints")

    # Normalize cumulative integral
    cum_integral = cumulative_trapezoid(y, x, initial=0)
    cum_integral /= cum_integral[-1]

    # Interpolate cumulative function
    integral_interp = interp1d(cum_integral, x, bounds_error=False)

    x_01 = float(integral_interp(0.01 * max(y)))
    x_99 = float(integral_interp(1 - 0.01 * max(y)))
    x_peak = float(x[np.argmax(y)])

    # Start building the set
    key_points = [x_01, x_peak, x_99]
    key_points = sorted(key_points)

    x_left, x_mid, x_right = key_points
    remaining_points = n - 3

    # Distance of intervals
    len_left = x_mid - x_left
    len_right = x_right - x_mid

    if remaining_points == 0:
        return np.array(key_points)

    if remaining_points % 2 == 0:
        n_left = n_right = remaining_points // 2
    else:
        if len_left >= len_right:
            n_left = remaining_points // 2 + 1
            n_right = remaining_points // 2
        else:
            n_left = remaining_points // 2
            n_right = remaining_points // 2 + 1

    # Generate evenly spaced points (excluding endpoints)
    left_points = np.linspace(x_left, x_mid, n_left + 2)[1:-1]
    right_points = np.linspace(x_mid, x_right, n_right + 2)[1:-1]

    all_points = np.sort(np.concatenate(([x_01, x_peak, x_99], left_points, right_points)))
    return np.array(all_points)


# 4. Generate Random Values
def extract_data(gridpoints, n, m, random_values):
    """
    Extract data by randomly sampling m participant rows from random_values
    and calculating the mean and standard deviation in n equally spaced columns.

    Parameters:
    n (int): Number of equally spaced columns - gridpoints.
    m (int): Number of participant rows to sample - participants.
    random_values (numpy.ndarray): Array of random values with shape (participants, timepoints).

    Returns:
    tuple: A tuple containing:
        - means (numpy.ndarray): Mean values for n equally spaced columns.
        - stds (numpy.ndarray): Standard deviation values for n equally spaced columns.
    """
    # Randomly sample m participant rows
    sampled_rows = random_values[np.random.choice(random_values.shape[0], m, replace=False), :]

    # Determine the indices for n equally spaced columns
    indices = np.linspace(0, sampled_rows.shape[1] - 1, n, dtype=int)

    # Calculate the mean and standard deviation for the selected columns
    means = sampled_rows[:, indices].mean(axis=0)
    stds = sampled_rows[:, indices].std(axis=0)
    current_grid = gridpoints[indices]

    return current_grid, means, stds

