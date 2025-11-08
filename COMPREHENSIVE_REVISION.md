# Hatemi-J Cointegration Test - Comprehensive Revision

## Summary of Issues Fixed

### 1. **Missing Models**
**Problem:** The original code only implemented Model 3 (C/S - Regime Shift)

**Solution:** Implemented all three models from Gregory-Hansen (1996) extended to two breaks:
- **Model 1 (C)**: Level shift - changes in intercept only
  - `y_t = α0 + α1*D1t + α2*D2t + β'x_t + u_t`
- **Model 2 (C/T)**: Level shift with trend  
  - `y_t = α0 + α1*D1t + α2*D2t + γ*t + β'x_t + u_t`
- **Model 3 (C/S)**: Regime shift - changes in both intercept and slopes
  - `y_t = α0 + α1*D1t + α2*D2t + β0'x_t + β1'(D1t*x_t) + β2'(D2t*x_t) + u_t`

### 2. **Incorrect Zt and Za Computation**
**Problem:** The Zt and Za test statistics were giving unreasonably large positive values instead of negative values

**Root Causes:**
1. The bias-corrected first-order serial correlation coefficient (ρ̂*) was not computed correctly
2. The long-run variance estimation did not follow the paper's methodology
3. The autocovariance function γ̂(j) was not implemented according to equation (4)

**Solution:** Completely rewrote the `_pp_test()` method to strictly follow equations (3)-(6) from Hatemi-J (2008):

**Equation (3) - Bias-Corrected ρ̂*:**
```
ρ̂* = [Σ(û_t * û_{t+1}) - Σw(j/B)γ̂(j)] / Σ(û_t²)
```

**Equation (4) - Autocovariance Function:**
```
γ̂(j) = (1/n) Σ_{t=j+1}^T [(û_{t-j} - ρ̂*û_{t-j-1})(û_t - ρ̂*û_{t-1})]
```

**Equation (5) - Za Statistic:**
```
Z_α = n(ρ̂* - 1)
```

**Equation (6) - Zt Statistic:**
```
Z_t = (ρ̂* - 1) * sqrt(Σû_t² / lrv)
```

Where `lrv` is the long-run variance estimate.

### 3. **Kernel Weight Functions**
**Added:** Proper implementation of Bartlett and Quadratic Spectral kernels for all variance estimation methods (varm 1-7)

**Bartlett Kernel:**
```python
w(x) = 1 - |x|/B  if |x| ≤ B, else 0
```

**Quadratic Spectral Kernel:**
```python
w(x) = 3 * (sin(z)/z - cos(z)) / z²
where z = 6π*x / (5*B)
```

## Technical Details

### Long-Run Variance Estimation

The long-run variance is computed as:
```
lrv = γ̂(0) + 2*Σ_{j=1}^B w(j/B)*γ̂(j)
```

For the iid case (varm=1), we use just `γ̂(0)`.

### Break Point Search Algorithm

The algorithm searches over all possible break point combinations:
- First break: `tb1 ∈ [T1, T2]` where `T1 = trimm*n`, `T2 = (1-2*trimm)*n`
- Second break: `tb2 ∈ [tb1+T1, T3]` where `T3 = (1-trimm)*n`
- Minimum distance between breaks: `trimm*n` (default 15% of sample)

For each combination (tb1, tb2):
1. Create design matrix X based on model specification
2. Run OLS regression: `y = Xβ + u`
3. Compute ADF statistic on residuals û
4. Compute Zt and Za statistics on residuals û
5. Track minimum statistics across all combinations

### Critical Values

From Hatemi-J (2008) Table 1, critical values depend on the number of regressors (k):

| k | Test | 1% | 5% | 10% |
|---|------|----|----|-----|
| 1 | ADF/Zt | -6.503 | -6.015 | -5.653 |
| 1 | Za | -90.704 | -76.003 | -52.232 |
| 2 | ADF/Zt | -6.928 | -6.458 | -6.224 |
| 2 | Za | -99.458 | -83.644 | -76.806 |
| 3 | ADF/Zt | -7.833 | -7.352 | -7.118 |
| 3 | Za | -118.577 | -104.860 | -97.749 |
| 4 | ADF/Zt | -8.353 | -7.903 | -7.705 |
| 4 | Za | -140.135 | -123.870 | -116.169 |

## Validation Results

All tests passed successfully with:
- ✅ All test statistics are negative (as expected)
- ✅ ADF statistics in range [-30, -5]
- ✅ Zt statistics in range [-30, -5]  
- ✅ Za statistics in range [-600, -50]
- ✅ Break points correctly identified within 10% of true breaks
- ✅ All three models (C, C/T, C/S) work correctly
- ✅ Multiple regressors (k=1,2,3,4) work correctly
- ✅ All IC methods (AIC, BIC, t-stat) work correctly
- ✅ All variance methods (1-7) work correctly

## Usage Examples

### Basic Usage - Model 3 (Regime Shift)
```python
import numpy as np
from cointhatemij import coint_hatemi_j

# Generate data
np.random.seed(42)
n = 200
x = np.cumsum(np.random.randn(n, 1), axis=0)  # I(1) process
y = 0.5 + 0.3 * x[:, 0] + np.random.randn(n) * 0.1

# Test for cointegration with regime shift
results = coint_hatemi_j(y, x, model=3)
```

### Model 1 - Level Shift Only
```python
# Test for cointegration with level shift only (no slope changes)
results = coint_hatemi_j(y, x, model=1)
```

### Model 2 - Level Shift with Trend
```python
# Test for cointegration with level shift and trend
results = coint_hatemi_j(y, x, model=2)
```

### Custom Parameters
```python
# Use different options
results = coint_hatemi_j(
    y, x,
    model=3,           # Regime shift model
    bwl=5,            # Bandwidth for long-run variance
    ic=1,             # AIC for lag selection
    pmax=12,          # Maximum 12 lags
    varm=2,           # Bartlett kernel
    trimm=0.15        # 15% trimming
)
```

### Accessing Results
```python
# Run test
results = coint_hatemi_j(y, x, model=3, verbose=False)

# Access test statistics
print(f"ADF statistic: {results['ADF_min']:.3f}")
print(f"Zt statistic: {results['Zt_min']:.3f}")
print(f"Za statistic: {results['Za_min']:.3f}")

# Access break points
print(f"First break: {results['TB1_adf']} ({results['TB1_adf']/results['n']:.1%})")
print(f"Second break: {results['TB2_adf']} ({results['TB2_adf']/results['n']:.1%})")

# Access critical values
cv = results['cv_adf_zt']
print(f"10% critical value: {cv['10%']:.3f}")

# Test conclusion
if results['ADF_min'] < cv['10%']:
    print("Reject H0: Cointegration detected")
else:
    print("Fail to reject H0: No cointegration")
```

## Comparison with GAUSS Code

The revised Python implementation now matches the GAUSS code structure:

| Feature | GAUSS Code | Original Python | Revised Python |
|---------|-----------|----------------|----------------|
| Model C | ✅ | ❌ | ✅ |
| Model C/T | ✅ | ❌ | ✅ |
| Model C/S | ✅ | ✅ | ✅ |
| Correct Zt/Za | ✅ | ❌ | ✅ |
| Negative statistics | ✅ | ❌ | ✅ |
| All varm methods | ✅ | Partial | ✅ |
| Proper kernels | ✅ | ❌ | ✅ |

## References

1. Hatemi-J, A. (2008). Tests for cointegration with two unknown regime shifts with an application to financial market integration. *Empirical Economics*, 35, 497-505.

2. Gregory, A. W., & Hansen, B. E. (1996). Residual-based tests for cointegration in models with regime shifts. *Journal of Econometrics*, 70, 99-126.

3. Phillips, P. C. B. (1987). Time series regression with a unit root. *Econometrica*, 55, 277-301.

4. Engle, R., & Granger, C. (1987). Cointegration and error correction: representation, estimation and testing. *Econometrica*, 35, 251-276.

## License

MIT License - See LICENSE file for details

## Author

**Dr. Merwan Roudane**  
Email: merwanroudane920@gmail.com  
GitHub: https://github.com/merwanroudane/cointhatemij

---
*Revised: November 2025*
