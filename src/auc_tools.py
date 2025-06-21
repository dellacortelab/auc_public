import numpy as np
from numpy import trapz
from scipy.interpolate import interp1d, UnivariateSpline
from scipy.stats import norm

# 1. Standard AUC Calculation
def standard_auc(x_vals, y_vals, stds, ns=5):
    """Calculate AUC using the trapezoidal rule."""
    lower, upper = estimate_auc_bounds(x_vals, y_vals, stds)
    return np.trapz(y_vals, x=x_vals), lower, upper

def estimate_auc_bounds(x_vals, means, stds):
    """
    Estimate rough lower and upper bounds for AUC using trapezoidal rule on (means ± stds).
    
    Parameters:
    - x_vals: Timepoints
    - means: Mean values at each timepoint
    - stds: Standard deviations at each timepoint

    Returns:
    - (lower_bound_auc, upper_bound_auc)
    """
    lower_vals = means - stds
    upper_vals = means + stds

    auc_lower = np.trapz(lower_vals, x_vals)
    auc_upper = np.trapz(upper_vals, x_vals)

    return auc_lower, auc_upper

# 2. Monte Carlo AUC Calculation
def monte_carlo_auc(x_vals, means, stds, n_sim=1000):
    """Calculate AUC using Monte Carlo simulation."""
    aucs = []
    for _ in range(n_sim):
        sampled_means = np.random.normal(means, stds)
        f_interp = UnivariateSpline(x_vals, sampled_means, s=0)
        fine_x = np.linspace(x_vals[0], x_vals[-1], 1000)
        fine_y = f_interp(fine_x)
        iauc = trapz(fine_y, fine_x)
        aucs.append(iauc)
    return np.mean(aucs), np.std(aucs)

def calculate_performance_metrics(estimates, true_values, ses=None, confidence_level=0.95):
    """
    Calculate performance metrics for AUC estimates.
    
    Parameters:
    -----------
    estimates : array-like
        The estimated AUC values
    true_values : array-like
        The true AUC values
    ses : array-like, optional
        Standard errors for the estimates
    confidence_level : float, optional
        Confidence level for interval coverage
        
    Returns:
    --------
    dict
        Dictionary containing performance metrics
    """
    estimates = np.array(estimates)
    true_values = np.array(true_values)
    n_sims = len(estimates)
    
    # Bias
    bias = np.mean(estimates - true_values)
    bias_mcse = np.std(estimates - true_values) / np.sqrt(n_sims)
    
    # Empirical Standard Error
    empirical_se = np.std(estimates)
    emp_se_mcse = empirical_se / np.sqrt(2 * (n_sims - 1))
    
    # Root Mean Squared Error
    rmse = np.sqrt(np.mean((estimates - true_values)**2))
    
    results = {
        'bias': bias,
        'empirical_se': empirical_se,
        'rmse': rmse,
        'bias_mcse': bias_mcse,
        'emp_se_mcse': emp_se_mcse
    }
    
    # Coverage (if standard errors provided)
    if ses is not None:
        ses = np.array(ses)
        z_value = norm.ppf((1 + confidence_level) / 2)
        lower = estimates - z_value * ses
        upper = estimates + z_value * ses
        coverage = np.mean((lower <= true_values) & (true_values <= upper))
        coverage_mcse = np.sqrt(coverage * (1 - coverage) / n_sims)
        
        results.update({
            'coverage': coverage,
            'coverage_mcse': coverage_mcse
        })
    
    return results

def format_performance_metrics(metrics, decimals=3):
    """Format performance metrics for display."""
    formatted = {}
    for k, v in metrics.items():
        formatted[k] = round(v, decimals)
    return formatted

def compare_methods_performance(true_values, methods_dict):
    """Compare performance of different AUC methods."""
    results = {}
    for method_name, estimates in methods_dict.items():
        metrics = calculate_performance_metrics(estimates, true_values)
        results[method_name] = format_performance_metrics(metrics)
    return results
