---
editor_options: 
  markdown: 
    wrap: 72
---

# African Farm Size Distribution Analysis

## Overview

This repository contains a comprehensive analysis of African farm size
distributions from 2000-2060, combining spatial farm structure data with
demographic projections to create corrected and harmonized datasets for
agricultural research and policy applications.

## Mathematical Methodology

### 1. Data Integration Framework

The analysis employs a three-case processing framework:

-   **Case A**: Countries with both spatial farm data (fsfix) and
    demographic projections (SSP2T)
-   **Case B**: Countries with demographic projections but limited/no
    spatial farm data\
-   **Case C**: Countries with only spatial farm data

### 2. Farm Size Distribution Sampling

For each farm size class *i* with *N*<sub>i</sub> farms, the algorithm
generates representative farm sizes:

#### Uniform Sampling (Classes 0-20 ha)

```         
s_i ~ U(min_i, max_i)
```

where *s*<sub>i</sub> are sampled farm sizes and *U* denotes uniform
distribution.

#### Log-normal Sampling (Class 20+ ha)

```         
s_i ~ LogNormal(μ = log(50), σ = 0.8)
```

constrained to *s*<sub>i</sub> ∈ [20, 1000] hectares.

#### Sample Weighting

Each sample receives weight:

```         
w_i = N_i / n_samples_i
```

where *n_samples*<sub>i</sub> = min(*N*<sub>i</sub>, MAX_SAMPLE_SIZE).

### 3. Mean-Shift Correction

To ensure consistency with agricultural area constraints, the algorithm
applies mean-shift correction:

#### Target Mean Calculation

```         
A_m = (Agricultural_Area_kha × 1000) / Total_Farms
```

#### Observed Mean Calculation

```         
A_o = Σ(s_i × w_i) / Σ(w_i)
```

#### Scaling Factor

```         
α = A_m / A_o
```

#### Corrected Farm Sizes

```         
s'_i = s_i × α
```

### 4. Rebinning Algorithm

Corrected farm sizes are rebinned into standardized classes:

```         
N'_0-1 = Σ w_i for s'_i ∈ (0, 1]
N'_1-2 = Σ w_i for s'_i ∈ (1, 2]  
N'_2-5 = Σ w_i for s'_i ∈ (2, 5]
N'_5-10 = Σ w_i for s'_i ∈ (5, 10]
N'_10-20 = Σ w_i for s'_i ∈ (10, 20]
N'_20+ = Σ w_i for s'_i > 20
```

### 5. Country-Specific Processing

#### Case A: Full Data Countries

```         
N'_ij(t) = (N_ij / Σ_j N_ij) × P_i(t)
```

where *N'*<sub>ij</sub>(t) is corrected farms in size class *j* for
admin unit *i* at time *t*, and *P*<sub>i</sub>(t) is projected total
farms.

#### Case B: SSP2T Only Countries

```         
N'_ij(t) = p_j^continental × P_i(t)
```

where *p*<sub>j</sub><sup>continental</sup> is the continental average
proportion for size class *j*.

#### Case C: Spatial Data Only Countries

```         
N'_ij(t) = N_ij(2000) × γ_ETH(t)
```

where *γ*<sub>ETH</sub>(t) is Ethiopia's growth factor relative to base
year 2000.

## Data Sources and Corrections

### Input Datasets

1.  **LUGE Smallholder Map (2020)**
    -   Source: SmallholderMap_20201202.shp
    -   Spatial farm size distribution data
    -   Admin-level disaggregation for African countries
2.  **SSP2T Demographic Projections**
    -   Source: ssp2t.rds
    -   Farm number projections 2000-2060
    -   Agricultural area projections
3.  **FAO World Census of Agriculture (2010)**
    -   Source: WCA_2010.csv
    -   Validation and correction data
    -   URL: <https://www.fao.org/world-census-agriculture/en>

### Key Data Corrections Applied

#### South Africa Correction

Farm structure adjusted using Aliber study estimates for non-commercial
farms: - Target: 2.6 million total farms - Source: Aliber & Hart (2009)
smallholder farmer support study - Applied area-weighted scaling using
Ricciardi et al. (2018) gross area estimates

#### South Sudan Demographic Constraint

Farm totals capped at 1 million farms based on: - Population (2000):
6,032,267 (World Bank) - Average household size: 5.9 - Estimated
households: 1,022,418 - Agricultural participation rate: 95% - **Maximum
farms**: 971,297 ≈ 1,000,000

Mathematical constraint:

```         
N_SSD(t) ≤ 1,000,000 ∀t
```

#### FAO Category Harmonization

Original FAO categories mapped to standardized bins: - Holdings 0-\<1,
1-\<2 → N0_1, N1_2\
- Holdings 2-\<3, 3-\<4, 4-\<5 → N2_5 - Holdings 5-\<10 → N5_10 -
Holdings 10-\<20 → N10_20 - Holdings ≥20 → N20\_

### Missing Data Handling

#### Administrative Unit Gaps

```         
x̄_missing = (1/n) × Σ_k x_k^observed
```

where missing values are filled with country-level means.

#### Continental Averages (Case B)

```         
p_j^continental = (1/|C|) × Σ_{c∈C} (N_cj / Σ_j N_cj)
```

where *C* is the set of countries with complete data.

## Output Dataset

### Primary Output: `fsfix_corrected_2000_2060.shp`

**Spatial Coverage**: All African countries with available data
**Temporal Range**: 2000, 2010, 2020, 2030, 2040, 2050, 2060
**Resolution**: Administrative level 1 (states/provinces)

#### Variables

| Variable         | Description                 | Units             |
|------------------|-----------------------------|-------------------|
| `GID_0`          | ISO3 country code           | \-                |
| `NAME_0`         | Country name                | \-                |
| `NAME_1`         | Administrative unit name    | \-                |
| `Year`           | Projection year             | \-                |
| `N0_1`           | Farms 0-1 hectares          | count             |
| `N1_2`           | Farms 1-2 hectares          | count             |
| `N2_5`           | Farms 2-5 hectares          | count             |
| `N5_10`          | Farms 5-10 hectares         | count             |
| `N10_20`         | Farms 10-20 hectares        | count             |
| `N20_`           | Farms 20+ hectares          | count             |
| `treatment_flag` | Processing method applied   | categorical       |
| `Median.smooth`  | SSP2T projected total farms | count             |
| `Agarea.kha`     | Agricultural area           | thousand hectares |

#### Treatment Flags

-   `"Case A: Data from both fsfix and SSP2T"` - Full methodology
    applied
-   `"Case B: SSP2T only - continental average disaggregation"` -
    Continental averages used
-   `"Case C: fsfix only - Ethiopia growth trend applied"` - Growth
    trends extrapolated
-   `"Case A/C Hybrid: fsfix data exists, but SSP2T missing for this year. Used ETH growth."` -
    Mixed approach
-   `"No farm data available"` - No processing possible

### Additional Outputs

1.  **`full_farm_distributions_all_countries.rds`**
    -   Individual farm size samples before mean-shift correction
    -   Detailed distribution data for statistical analysis
2.  **`full_farm_distributions_shifted.rds`**
    -   Individual farm size samples after mean-shift correction
    -   Used for validation and distribution analysis
3.  **Validation Plots** (in `figs/` directory)
    -   `validation_totals_check.png` - SSP2T vs output totals
    -   `validation_areas_check.png` - Agricultural area consistency
    -   `africa_farms_2000_2050.png` - Spatial distribution maps
    -   `africa_farm_size_distribution_evolution.png` - Distribution
        evolution

## Validation and Quality Control

### 1. Farm Total Consistency

Validates that output totals match SSP2T projections within tolerance:

```         
|N'_total - N_SSP2T| / N_SSP2T < ε
```

### 2. Agricultural Area Consistency

Validates mean farm sizes align with area constraints:

```         
|A_implied - A_target| / A_target < ε
```

### 3. Distribution Preservation

Ensures farm size distribution shapes remain realistic after
corrections.

## Usage

### Requirements

``` r
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
```

### Running the Analysis

``` r
# Data preprocessing and harmonization
source("fsdistfix.R")

# Full analysis with projections  
source("test.R")
```

### Loading Output Data

``` r
# Load spatial results
farm_data <- st_read("output/fsfix_corrected_2000_2060.shp")

# Load distribution data
distributions <- readRDS("output/full_farm_distributions_shifted.rds")
```

## References

-   Aliber, M., & Hart, T. G. (2009). Should subsistence agriculture be
    supported as a strategy to address rural food insecurity? *Agrekon*,
    48(4), 434-458.

-   FAO. (2010). World Census of Agriculture. Rome: Food and Agriculture
    Organization. Available:
    <https://www.fao.org/world-census-agriculture/en>

-   Ricciardi, V., et al. (2018). How much of the world's food do
    smallholders produce? *Global Food Security*, 17, 64-72.

-   World Bank. Population data. Available:
    <https://data.worldbank.org/>

## Citation

If you use this dataset, please cite:

```         
Mehrabi, Zia (2025). African Farm Size Distribution Analysis. Better Planet Laboratory.
```

## License

MIT

## Contact
