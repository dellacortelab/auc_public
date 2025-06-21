# Model Zoo Generator

import numpy as np
from scipy import stats
from scipy.special import gamma

def generate_response_curve(x, curve_type='skewed_gaussian', params=None, threshold=0.01):
    """
    Generate various normalized response curves commonly used in pharmacokinetics and glycemic responses.
    
    Parameters:
    - x (array-like): Input time values (should be non-negative).
    - curve_type (str): Type of response curve to generate. Options:
        - 'skewed_gaussian': Skewed Gaussian curve (default)
        - 'biexponential': Two-compartment model (absorption/elimination)
        - 'gamma': Gamma distribution
        - 'lognormal': Log-normal distribution
        - 'weibull': Weibull function
        - 'bateman': Bateman function (absorption/elimination)
        - 'inverse_gaussian': Inverse Gaussian distribution
    - params (dict): Parameters for the specified curve type. If None, default values are used.
    - threshold (float): Fraction of peak value to define meaningful support region.
    
    Returns:
    - response (array-like): Normalized response values.
    - support_start (float): Start of region where response > threshold * peak.
    - support_end (float): End of region where response > threshold * peak.
    """
    x = np.asarray(x)
    
    # Ensure x is non-negative for biological plausibility
    x_valid = x.copy()
    x_valid[x < 0] = 0
    
    # Default parameters for each curve type
    default_params = {
        'skewed_gaussian': {'A': 100, 'mu': 15, 'sigma': 10, 'skew': 2},
        'biexponential': {'A': 1, 'alpha': 0.1, 'beta': 0.01, 'ka': 0.2},
        'gamma': {'alpha': 3, 'beta': 0.1},
        'lognormal': {'mu': 1.5, 'sigma': 0.5},
        'weibull': {'k': 2, 'lambda': 30},
        'bateman': {'ka': 0.2, 'ke': 0.1, 'F': 1, 'D': 1, 'V': 1},
        'inverse_gaussian': {'mu': 30, 'lambda': 100}
    }
    
    # Use default parameters if none provided
    if params is None:
        params = default_params[curve_type]
    else:
        # Merge with defaults for any missing parameters
        for key, value in default_params[curve_type].items():
            if key not in params:
                params[key] = value
    
    # Generate the appropriate curve
    if curve_type == 'skewed_gaussian':
        A, mu, sigma, skew = params['A'], params['mu'], params['sigma'], params['skew']
        gaussian = A * np.exp(-0.5 * ((x_valid - mu) / sigma) ** 2)
        skew_factor = 1 + skew * (x_valid - mu) / sigma
        response = gaussian * np.clip(skew_factor, 0, None)
    
    elif curve_type == 'biexponential':
        A, alpha, beta, ka = params['A'], params['alpha'], params['beta'], params['ka']
        response = A * (np.exp(-beta * x_valid) - np.exp(-ka * x_valid))
        response[x_valid <= 0] = 0  # Set response to 0 for x <= 0
    
    elif curve_type == 'gamma':
        alpha, beta = params['alpha'], params['beta']
        # Using gamma PDF directly
        response = stats.gamma.pdf(x_valid, a=alpha, scale=1/beta)
    
    elif curve_type == 'lognormal':
        mu, sigma = params['mu'], params['sigma']
        # Using lognormal PDF directly
        response = stats.lognorm.pdf(x_valid, s=sigma, scale=np.exp(mu))
        response[x_valid <= 0] = 0  # Set response to 0 for x <= 0
    
    elif curve_type == 'weibull':
        k, lam = params['k'], params['lambda']
        # Using Weibull PDF directly
        response = stats.weibull_min.pdf(x_valid, c=k, scale=lam)
        response[x_valid <= 0] = 0  # Set response to 0 for x <= 0
    
    elif curve_type == 'bateman':
        ka, ke = params['ka'], params['ke']
        F, D, V = params['F'], params['D'], params['V']
        
        # Bateman function (oral administration model)
        response = (F * D * ka) / (V * (ka - ke)) * (np.exp(-ke * x_valid) - np.exp(-ka * x_valid))
        response[x_valid <= 0] = 0  # Set response to 0 for x <= 0
    
    elif curve_type == 'inverse_gaussian':
        mu, lam = params['mu'], params['lambda']
        # Using inverse Gaussian PDF directly
        response = stats.invgauss.pdf(x_valid, mu=mu, scale=lam)
        response[x_valid <= 0] = 0  # Set response to 0 for x <= 0
    
    else:
        raise ValueError(f"Unknown curve type: {curve_type}")
    
    # Normalize the curve to have integral = 1
    integral = np.trapz(response, x)
    if integral > 0:  # Avoid division by zero
        response /= integral
    
    # Determine region of significant response
    peak = np.max(response)
    active_mask = response > (threshold * peak)
    if np.any(active_mask):
        support_start = x[active_mask][0]
        support_end = x[active_mask][-1]
    else:
        support_start, support_end = x[0], x[-1]  # fallback if response is too weak

    return response, support_start, support_end

parameter_bounds = {
    'skewed_gaussian': {
        'A': (10, 500),       # Amplitude
        'mu': (10, 100),      # Mean (time of peak response)
        'sigma': (5, 40),     # Standard deviation (width of response)
        'skew': (-2, 5)       # Skewness parameter
    },
    'biexponential': {
        'A': (0.5, 5),        # Amplitude factor
        'alpha': (0.05, 0.5), # Elimination rate constant
        'beta': (0.005, 0.05),# Terminal elimination rate constant
        'ka': (0.1, 0.4)      # Absorption rate constant
    },
    'gamma': {
        'alpha': (1, 5),      # Shape parameter
        'beta': (0.05, 0.3)   # Rate parameter
    },
    'lognormal': {
        'mu': (1, 3),         # Mean of logarithm
        'sigma': (0.2, 0.8)   # Standard deviation of logarithm
    },
    'weibull': {
        'k': (1.5, 4),        # Shape parameter
        'lambda': (15, 60)    # Scale parameter
    },
    'bateman': {
        'ka': (0.1, 0.5),     # Absorption rate constant
        'ke': (0.01, 0.2),    # Elimination rate constant
        'F': (0.5, 1),        # Bioavailability fraction
        'D': (1, 1),          # Dose (fixed for normalization)
        'V': (1, 1)           # Volume of distribution (fixed for normalization)
    },
    'inverse_gaussian': {
        'mu': (15, 60),       # Mean (location parameter)
        'lambda': (50, 200)   # Shape parameter
    }
}

def sample_curves(parameter_bounds, num_samples=10):
    """
    Samples curves for each response type using the provided parameter bounds.

    Args:
        parameter_bounds (dict): Dictionary containing parameter bounds for each response type.
        num_samples (int): Number of curves to sample for each response type.

    Returns:
        dict: A dictionary where keys are response types and values are lists of sampled parameter dictionaries.
    """
    sampled_curves = {}

    for curve_type, bounds in parameter_bounds.items():
        sampled_curves[curve_type] = []
        for _ in range(num_samples):
            params = {}
            for param, (min_val, max_val) in bounds.items():
                params[param] = np.random.uniform(min_val, max_val)
            sampled_curves[curve_type].append(params)

    return sampled_curves