# AUC Analysis Tools

This repository contains code for calculating Area Under the Curve (AUC) metrics using standard and Monte Carlo methods. It implements reliable tools for AUC analysis from graph-derived summary data.

## Repository Structure

```
auc_public/
├── src/                    # Source code
│   ├── auc_tools.py       # Core functionality for AUC calculations
│   ├── response_curves.py # Implementation of response curve models
│   └── sampling_tools.py  # Utilities for sampling and analysis
├── notebooks/             # Jupyter notebooks
│   ├── clean_protocol.ipynb      # Main analysis notebook
│   ├── create_curves_and_fits.ipynb
│   ├── plot_results.ipynb
│   ├── post_calc_iauc.ipynb
│   └── use_case_figs.ipynb
├── results/               # Generated figures and data
│   └── response_curves_zoo.pdf
├── docs/                  # Documentation
├── requirements.txt       # Python dependencies
├── LICENSE               # MIT License
└── README.md            # This file
```

## Installation

To run this code, you need Python 3.8 or higher. Clone this repository and install the required packages:

```bash
git clone [repository-url]
cd auc_public
pip install -r requirements.txt
```

## Usage

The repository provides two main methods for AUC calculation:

1. **Standard AUC**: A traditional method using the trapezoidal rule with error bounds
2. **Monte Carlo AUC**: A simulation-based approach that provides robust estimates with uncertainty quantification

The main analysis can be reproduced by running the Jupyter notebooks in the `notebooks/` directory.

## Citation

If you use this code in your research, please cite our article:

> Titensor, Sean, Joshua Ebbert, Karen Della Corte, and Dennis Della Corte. "Estimating area under the curve from graph-derived summary data: a systematic comparison of standard and Monte Carlo approaches." BMC Medical Research Methodology 25, no. 1 (2025): 197.

## License

This project is licensed under the MIT License - see the LICENSE file for details.
