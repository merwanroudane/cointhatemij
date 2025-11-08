# Hatemi-J Cointegration Test with Two Unknown Regime Shifts

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

Python implementation of the Hatemi-J (2008) cointegration test that accounts for two unknown structural breaks in the cointegrating relationship.

## Overview

This package implements three residual-based test statistics (ADF, Zt, Za) for testing cointegration between time series variables when the long-run relationship may experience two regime shifts during the sample period. The timing of these shifts is determined endogenously by the data.

## Features

- ✅ **Three Model Specifications:**
  - Model 1 (C): Level shift - changes in intercept only
  - Model 2 (C/T): Level shift with trend
  - Model 3 (C/S): Regime shift - changes in both intercept and slopes

- ✅ **Three Test Statistics:**
  - ADF: Augmented Dickey-Fuller test
  - Zt: Phillips-Perron Zt test
  - Za: Phillips-Perron Zα test

- ✅ **Flexible Options:**
  - Multiple information criteria for lag selection (AIC, BIC, t-stat)
  - Seven variance estimation methods (iid, Bartlett, QS, SPC, Kurozumi)
  - Customizable trimming rates and bandwidth

- ✅ **Accurate Implementation:**
  - Follows Hatemi-J (2008) methodology exactly
  - Critical values from Monte Carlo simulations in the paper
  - Proper computation of Zt and Za statistics (equations 3-6)
  - Compatible with GAUSS code structure

## Installation

```bash
pip install cointhatemij
```

Or install from source:
```bash
git clone https://github.com/merwanroudane/cointhatemij.git
cd cointhatemij
pip install -e .
```

## Quick Start

```python
import numpy as np
from cointhatemij import coint_hatemi_j

# Generate sample data
np.random.seed(42)
n = 200
x = np.cumsum(np.random.randn(n, 1), axis=0)  # I(1) process
y = 0.5 + 0.3 * x[:, 0] + np.random.randn(n) * 0.1

# Test for cointegration with two regime shifts
results = coint_hatemi_j(y, x, model=3)

# Access results
print(f"ADF statistic: {results['ADF_min']:.3f}")
print(f"First break: {results['TB1_adf']}")
print(f"Second break: {results['TB2_adf']}")
```

## Models

### Model 1 (C): Level Shift
Changes in intercept only:
```
y_t = α0 + α1*D1t + α2*D2t + β'x_t + u_t
```

**Use when:** Only the level of the cointegrating relationship changes, but the slopes remain constant.

### Model 2 (C/T): Level Shift with Trend
Changes in intercept with time trend:
```
y_t = α0 + α1*D1t + α2*D2t + γ*t + β'x_t + u_t
```

**Use when:** The cointegrating relationship has a time trend and the level shifts.

### Model 3 (C/S): Regime Shift
Changes in both intercept and slopes:
```
y_t = α0 + α1*D1t + α2*D2t + β0'x_t + β1'(D1t*x_t) + β2'(D2t*x_t) + u_t
```

**Use when:** Both the level and the cointegrating coefficients change (most general case).

## Detailed Usage

### Basic Test with All Models

```python
from cointhatemij import coint_hatemi_j

# Model 1: Level shift only
results1 = coint_hatemi_j(y, x, model=1)

# Model 2: Level shift with trend
results2 = coint_hatemi_j(y, x, model=2)

# Model 3: Regime shift (default)
results3 = coint_hatemi_j(y, x, model=3)
```

### Custom Options

```python
results = coint_hatemi_j(
    y, x,
    model=3,           # Regime shift model
    bwl=5,            # Bandwidth for long-run variance
    ic=1,             # 1=AIC, 2=BIC, 3=t-stat
    pmax=12,          # Maximum lags for ADF test
    varm=2,           # 1=iid, 2=Bartlett, 3=QS, ...
    trimm=0.15,       # Trimming rate (15%)
    verbose=True      # Print detailed results
)
```

### Using the Class Interface

```python
from cointhatemij import HatemiJTest

# Initialize test
test = HatemiJTest(y, x, model=3, ic=1, pmax=8)

# Run test
results = test.fit()

# Print summary
test.summary()

# Access individual components
print(f"Number of observations: {test.n}")
print(f"Number of regressors: {test.k}")
print(f"Bandwidth: {test.bwl}")
```

## Interpreting Results

### Test Statistics

All three statistics (ADF, Zt, Za) should be **negative**. More negative values provide stronger evidence against the null hypothesis of no cointegration.

```python
results = coint_hatemi_j(y, x, model=3, verbose=False)

# Test statistics
adf_stat = results['ADF_min']  # Should be < 0
zt_stat = results['Zt_min']    # Should be < 0
za_stat = results['Za_min']    # Should be < 0

# Critical values
cv = results['cv_adf_zt']  # For ADF and Zt
print(f"10% critical value: {cv['10%']}")  # e.g., -5.653
print(f"5% critical value: {cv['5%']}")    # e.g., -6.015
print(f"1% critical value: {cv['1%']}")    # e.g., -6.503
```

### Decision Rule

Reject H0 (no cointegration) if:
```
test_statistic < critical_value
```

Example:
```python
if results['ADF_min'] < results['cv_adf_zt']['10%']:
    print("Reject H0: Cointegration detected (10% level)")
else:
    print("Fail to reject H0: No cointegration")
```

### Break Dates

The test estimates two break points for each statistic:

```python
# Break points for ADF test
tb1_adf = results['TB1_adf']  # First break (observation number)
tb2_adf = results['TB2_adf']  # Second break (observation number)

# Express as fractions of sample
frac1 = tb1_adf / results['n']  # e.g., 0.30 (30% through sample)
frac2 = tb2_adf / results['n']  # e.g., 0.70 (70% through sample)
```

## Examples

### Example 1: Financial Market Integration

Testing integration between US and UK stock markets (similar to paper's application):

```python
import numpy as np
from cointhatemij import coint_hatemi_j

# Load your data (here we simulate)
np.random.seed(1991)
n = 550  # Weekly data over ~10 years

# Simulate log stock indices
us_index = 100 + np.cumsum(np.random.randn(n) * 0.02)
uk_index = 10 + 0.8 * us_index + np.random.randn(n) * 2

# Test for cointegration with regime shifts
results = coint_hatemi_j(uk_index, us_index.reshape(-1, 1), model=3)

# The paper found significant cointegration with breaks in
# early 1991 (Gulf War) and end 1992 (exchange rate crisis)
```

### Example 2: Multiple Regressors

```python
# Generate data with 3 regressors
n = 200
x = np.cumsum(np.random.randn(n, 3), axis=0)
y = 1 + 0.5*x[:,0] + 0.3*x[:,1] - 0.2*x[:,2] + np.random.randn(n)*0.5

# Test with multiple regressors
results = coint_hatemi_j(y, x, model=3)

# Critical values automatically adjust for k=3
print(f"Critical value (10%): {results['cv_adf_zt']['10%']}")  # -7.118
```

### Example 3: Comparing Models

```python
models = [1, 2, 3]
model_names = ['Level shift', 'Level shift + trend', 'Regime shift']

for model, name in zip(models, model_names):
    print(f"\n{name}:")
    results = coint_hatemi_j(y, x, model=model, verbose=False)
    print(f"  ADF: {results['ADF_min']:.3f}")
    print(f"  Breaks: {results['TB1_adf']}, {results['TB2_adf']}")
```

## Critical Values

Critical values from Hatemi-J (2008, Table 1) are built into the package:

| Regressors (k) | Test | 1% | 5% | 10% |
|----------------|------|----|----|-----|
| 1 | ADF, Zt | -6.503 | -6.015 | -5.653 |
| 1 | Za | -90.704 | -76.003 | -52.232 |
| 2 | ADF, Zt | -6.928 | -6.458 | -6.224 |
| 2 | Za | -99.458 | -83.644 | -76.806 |
| 3 | ADF, Zt | -7.833 | -7.352 | -7.118 |
| 3 | Za | -118.577 | -104.860 | -97.749 |
| 4 | ADF, Zt | -8.353 | -7.903 | -7.705 |
| 4 | Za | -140.135 | -123.870 | -116.169 |

## Parameters

### Main Parameters

- **y**: `np.ndarray` - Dependent variable (n×1)
- **x**: `np.ndarray` - Independent variables (n×k), k ≤ 4
- **model**: `int` - Model specification (1, 2, or 3)

### Optional Parameters

- **bwl**: `int` - Bandwidth for long-run variance. Default: `round(4*(n/100)^(2/9))`
- **ic**: `int` - Information criterion for ADF lag selection:
  - `1` = AIC (Akaike)
  - `2` = BIC (Schwarz)
  - `3` = t-stat significance (default)
- **pmax**: `int` - Maximum lags for ADF test. Default: `8`
- **varm**: `int` - Variance estimation method:
  - `1` = iid (default)
  - `2` = Bartlett kernel
  - `3` = Quadratic Spectral (QS)
  - `4` = SPC with Bartlett
  - `5` = SPC with QS
  - `6` = Kurozumi with Bartlett
  - `7` = Kurozumi with QS
- **trimm**: `float` - Trimming rate. Default: `0.10` (10%)
- **verbose**: `bool` - Print detailed results. Default: `True`

## Return Values

The `coint_hatemi_j()` function returns a dictionary with:

```python
{
    'ADF_min': -7.123,           # ADF test statistic
    'TB1_adf': 45,               # First break (ADF)
    'TB2_adf': 135,              # Second break (ADF)
    'Zt_min': -7.234,            # Zt test statistic
    'TB1_zt': 44,                # First break (Zt)
    'TB2_zt': 136,               # Second break (Zt)
    'Za_min': -89.123,           # Za test statistic
    'TB1_za': 44,                # First break (Za)
    'TB2_za': 136,               # Second break (Za)
    'cv_adf_zt': {...},          # Critical values for ADF/Zt
    'cv_za': {...},              # Critical values for Za
    'n': 200,                    # Sample size
    'k': 2,                      # Number of regressors
    'model': 3,                  # Model used
    'model_name': '...',         # Model description
    'trimm': 0.10,              # Trimming rate
    'bwl': 4,                   # Bandwidth
    'ic': 3,                    # IC method
    'pmax': 8,                  # Max lags
    'varm': 1                   # Variance method
}
```

## Important Notes

### Data Requirements

1. **Stationarity**: Variables should be integrated of the same order (typically I(1))
2. **Sample size**: Minimum n ≈ 50-100 for reliable inference
3. **Breaks**: At least 15% (default trimming) of observations between breaks

### When to Use Each Model

- **Model 1 (C)**: When you suspect only the constant term changes
  - Example: Regulatory change that shifts price levels
  
- **Model 2 (C/T)**: When there's a time trend and level shifts
  - Example: Technological progress with structural breaks
  
- **Model 3 (C/S)**: When the entire relationship changes
  - Example: Regime changes in market integration
  - **Most general and commonly used**

### Computation Time

The test searches over all possible break combinations:
- For n=100: ~2,000 combinations
- For n=200: ~11,000 combinations  
- For n=500: ~70,000 combinations

Typical runtime: 1-30 seconds depending on sample size.

## Mathematical Details

### Test Statistics

**ADF Test:** t-statistic on lagged residual level

**Zt Test (Equation 6):**
```
Zt = (ρ̂* - 1) × √(Σû²ₜ / lrv)
```

**Za Test (Equation 5):**
```
Za = n(ρ̂* - 1)
```

Where ρ̂* is the bias-corrected first-order autocorrelation (Equation 3):
```
ρ̂* = [Σûₜûₜ₊₁ - Σw(j/B)γ̂(j)] / Σû²ₜ
```

### Long-Run Variance

Computed as:
```
lrv = γ̂(0) + 2Σⱼ₌₁ᴮ w(j/B)γ̂(j)
```

Where γ̂(j) is the autocovariance function (Equation 4).

## Validation

This implementation has been thoroughly tested against:
- ✅ Original GAUSS code structure
- ✅ Hatemi-J (2008) paper methodology
- ✅ Critical values from Table 1
- ✅ Multiple test scenarios

All test statistics produce:
- ✅ Negative values (as expected)
- ✅ Reasonable magnitudes
- ✅ Correct rejection rates

## References

1. **Hatemi-J, A. (2008).** Tests for cointegration with two unknown regime shifts with an application to financial market integration. *Empirical Economics*, 35, 497-505. [DOI: 10.1007/s00181-007-0175-9](https://doi.org/10.1007/s00181-007-0175-9)

2. **Gregory, A. W., & Hansen, B. E. (1996).** Residual-based tests for cointegration in models with regime shifts. *Journal of Econometrics*, 70, 99-126.

3. **Phillips, P. C. B. (1987).** Time series regression with a unit root. *Econometrica*, 55, 277-301.

4. **Engle, R., & Granger, C. (1987).** Cointegration and error correction: representation, estimation and testing. *Econometrica*, 35, 251-276.

## Citation

If you use this package in your research, please cite:

```bibtex
@software{cointhatemij2025,
  author = {Roudane, Merwan},
  title = {cointhatemij: Hatemi-J Cointegration Test for Python},
  year = {2025},
  url = {https://github.com/merwanroudane/cointhatemij}
}

@article{hatemi2008tests,
  title={Tests for cointegration with two unknown regime shifts with an application to financial market integration},
  author={Hatemi-J, Abdulnasser},
  journal={Empirical Economics},
  volume={35},
  pages={497--505},
  year={2008},
  publisher={Springer}
}
```

## License

MIT License - see LICENSE file for details

## Author

**Dr. Merwan Roudane**  
Email: merwanroudane920@gmail.com  
GitHub: [github.com/merwanroudane](https://github.com/merwanroudane)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## Support

If you encounter any issues or have questions:
- Open an issue on [GitHub](https://github.com/merwanroudane/cointhatemij/issues)
- Email: merwanroudane920@gmail.com

## Changelog

### Version 1.1.0 (2025-11-08)
- ✅ Added Model 1 (C - Level shift)
- ✅ Added Model 2 (C/T - Level shift with trend)
- ✅ Fixed Zt and Za computation (now correctly negative)
- ✅ Implemented proper long-run variance estimation
- ✅ Added comprehensive examples
- ✅ Improved documentation

