# Data Description

## Processed Dataset

The file `Empirical_Data_3DCCE_Paper.xlsx` is included in this folder.
It contains the fully processed bilateral trade panel used for empirical
estimation. A researcher can replicate the empirical results by running
the four scripts in `code/empirical/` directly from this file — no
additional data download is required.

For researchers who wish to extend or modify the panel (different country
coverage, different years, or updated BACI releases), the raw data sources
and processing instructions are described below.

---

## Workbook Structure

| Sheet | Contents |
|-------|----------|
| `README` | Panel description and data source summary |
| `1.BilateralPanel` | Full balanced panel — one row per origin × destination × year (9,660 rows) |
| `2.GDP` | Origin and destination GDP series in long format |
| `3.EstimationResults` | Mean group estimates (populated after running E03) |
| `4.CD_Diagnostics` | Cross-sectional dependence diagnostics (populated after running E02) |
| `5.PairSummary` | Per-pair summary statistics |

---

## Variable Definitions — Sheet `1.BilateralPanel`

| Column | Variable | Unit | Source |
|--------|----------|------|--------|
| `Origin (ISO3)` | Origin economy — 3-letter ISO alpha code | — | CEPII BACI |
| `Destination (ISO3)` | Destination economy — 3-letter ISO alpha code | — | CEPII BACI |
| `Year` | Calendar year | — | — |
| `i (origin index)` | Integer index for origin, 1 to N | — | Constructed |
| `j (dest index)` | Integer index for destination, 1 to M | — | Constructed |
| `t (time index)` | Integer index for year, 1 to T | — | Constructed |
| `Trade value (USD thousands)` | Bilateral exports summed over all HS product codes | USD thousands | CEPII BACI |
| `ln(Trade+1)` | ln(trade value + 1) — outcome variable | Log USD thousands | Constructed |
| `ln(GDP_origin)` | ln(GDP, origin economy) | Log current USD | World Bank / IMF |
| `ln(GDP_dest)` | ln(GDP, destination economy) | Log current USD | World Bank / IMF |

---

## Panel Coverage

| Dimension | Value |
|-----------|-------|
| N — origin economies | 18 |
| M — destination economies | 19 |
| T — years | 30 (1995–2024) |
| Balanced pairs | 322 out of 342 possible (N × M − N, excluding self-flows) |
| Total observations | 9,660 |

**Origin economies (N = 18):**
AUS, CAN, CHE, CHN, DEU, ESP, FRA, GBR, IND, ITA, JPN, KOR, MEX, MYS,
NLD, RUS, TAP, USA

**Destination economies (M = 19):**
CAN, CHE, CHN, DEU, ESP, FRA, GBR, HKG, IND, ITA, JPN, KOR, MEX, NLD,
POL, SGP, TAP, USA, VNM

---

## Data Sources and Download Instructions

### Source 1 — CEPII BACI Bilateral Trade Database

**URL:** https://www.cepii.fr/CEPII/en/bdd_modele/bdd_modele_item.asp?id=37

**Files used in the processed dataset:**
- HS92 classification, years 1995–2016, version V202601
  (`BACI_HS92_Y{YYYY}_V202601.csv` for each year)
- HS17 classification, years 2017–2024, version V202601
  (`BACI_HS17_Y{YYYY}_V202601.csv` for each year)

**Raw file structure:**

| Column | Description |
|--------|-------------|
| `t` | Year |
| `i` | Exporter (UN numeric country code) |
| `j` | Importer (UN numeric country code) |
| `k` | HS product code |
| `v` | Trade value (USD thousands) |
| `q` | Quantity |

**Processing applied:** Within each year, rows are grouped by (`t`, `i`, `j`)
and values are summed over all product codes `k` to produce bilateral totals.
UN numeric codes are mapped to ISO3 alpha codes. Chinese Taipei is assigned
code 490 in BACI and mapped to `TAP` in the processed data.

**Note on redistribution:** The raw BACI files are not redistributed in
this repository. They must be downloaded from the CEPII website. The
processed `Empirical_Data_3DCCE_Paper.xlsx` derived from these files is
included and is sufficient for replication.

---

### Source 2 — World Bank World Development Indicators

**Indicator:** NY.GDP.MKTP.CD (GDP, current US dollars)

**URL:** https://databank.worldbank.org/source/world-development-indicators

**Direct API:**
```
https://api.worldbank.org/v2/country/{ISO3}/indicator/NY.GDP.MKTP.CD?format=json&date=1995:2024
```

**Coverage:** All economies in the panel except Chinese Taipei (not
available in World Bank databases).

---

### Source 3 — IMF World Economic Outlook (Chinese Taipei only)

**Indicator:** NGDPD (GDP, current prices, billions of USD)

**URL:** https://www.imf.org/external/datamapper/NGDPD/TWN

**Direct API:**
```
https://www.imf.org/external/datamapper/api/v1/NGDPD/TWN
```

**Note:** IMF reports in billions of USD. The processing scripts multiply
by 1,000,000,000 before taking the logarithm to maintain consistent units
with the World Bank series.

---

## Notes on Data Construction

- The panel is balanced: every origin-destination pair appears in all 30 years.
- Zero-trade flows are included. They enter the model as ln(0 + 1) = 0.
- The integer indices `i`, `j`, `t` are assigned by sorting economies
  alphabetically and years chronologically. These indices are used by
  `E01_data_from_excel.R` to construct the three-dimensional arrays.
- Researchers extending the panel to additional years should append new
  BACI year files and re-run the full assembly pipeline. The top-20
  origins and destinations are selected by average trade value over the
  full sample; extending the sample may change the composition marginally.
