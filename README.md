## Overview

This repository contains the code developed for my Bachelor's thesis in the Nolte Group.

The core script, **`Bindfitnplot.py`**, automates the calculation of association constants from triplicate fluorescence titrations. The program evaluates three different 1:1 binding models:

- **Stern–Volmer**
- **Tsukube**
- **Connors**

For each model, the script fits the corresponding binding isotherms and determines parameters such as:

- **Ka** — association constant  
- **F₀** — fluorescence of the host in the absence of guest  
- Additional proportionality constants depending on the chosen model  

The script also calculates the initial concentrations \([H]_0\) and \([G]_0\), which are required as input for all three binding equations.

The fitting routine has been validated against previous experimental datasets and existing binding-analysis programs. These comparisons show that the program performs reliably for **1:1 host–guest complexes**.

The main features and required data specifiers of the fitting program are detailed in the sections below.

## Publication
https://www.nature.com/articles/s41467-020-18596-1

---

## Fitting Methods

The script provides three fitting functions corresponding to three different binding equations. Each function computes intermediate concentrations (such as \([G]\) or \([HG]\)) using the appropriate analytical expressions and then fits the model parameters to the experimental data.

### 1. `Volmer_Stern`
This function fits the data using Equation 26. The free guest concentration \([G]\) is calculated using Equation 9.  
The fitted parameters are:

- **Ka** — association constant  
- **F₀** — fluorescence of the host before any guest is added  

Although the script allows the user to fix \(F₀\) at the intensity of the initial peaks, an internal comparison check is implemented. For best performance, it is recommended to **let the program optimize \(F₀\)**.

### 2. `Tsukube`
This function fits the data using Equation 24. The concentration of the host–guest complex \([HG]\) is computed via Equation 12.  
The fitted parameters are:

- **Ka**
- **F₀**
- **kΔHG** — proportionality constant \(k_{\Delta HG}\)

### 3. `Connors`
This is the most flexible (and most parameter-rich) fitting function. It fits the data using Equation 25, again obtaining \([HG]\) from Equation 12.  
The parameters fitted are:

- **kH** — proportionality constant of free host in presence of guest  
- **kH0** — proportionality constant for the initial free host  
- **kHG** — proportionality constant of the complex  
- **Ka**  
- **F₀**

Because this model includes **five fitting parameters**, it is easier to obtain seemingly perfect fits, but this can reduce the accuracy of the resulting \(K_a\). Only **ratios** such as \(k_H/k_H^0\) and \(k_{HG}/k_H^0\) are physically meaningful.  
If \(k_H/k_H^0 \approx 1\), the assumption of **no dynamic quenching** of the host is justified.

## Fitting Program: Data Specifiers

The fitting program uses a set of data specifiers that must be provided for each experiment.  
Below is an overview of all required fields in the order they appear in the script.

### Experiment Metadata
| Specifier | Description |
|----------|-------------|
| **Expnr** | Experiment name (e.g., `PT002`) |
| **Host** | Host name (e.g., `(+)-H23`) |
| **Guest** | Guest name (e.g., `(S,S)-V2`) |

### Host Stock Solution
| Specifier | Description |
|----------|-------------|
| **M_Host** | Host molecular weight (g/mol) |
| **m_Hstock** | Mass of dissolved Host (g) |
| **V_solv_Hstock** | Volume of host stock solution (mL) |
| **m_solv_Hstock** | Mass of host stock solution (g) |

### Guest Stock Solution
| Specifier | Description |
|----------|-------------|
| **M_Guest** | Guest molecular weight (g/mol) |
| **m_Gstock** | Mass of dissolved Guest (g) |
| **V_solv_Gstock** | Volume of guest stock solution (mL) |
| **m_solv_Gstock** | Mass of guest stock solution (g) |

### Measurement Solutions
| Specifier | Description |
|----------|-------------|
| **V_solv_Hmeasure** | Total volume of host measurement solution (mL) |
| **m_Hstock_diluted** | Mass of diluted host stock solution (g) |
| **V_H_titrated** | Volume of host measurement solution added in first step (typically 2.0 mL) |
| **V_solv_Gmeasure** | Total volume of guest measurement solution (mL) |
| **m_Gstock_diluted** | Mass of diluted guest stock solution (g) |
| **m_Hstock_diluted_forGmeas** | Mass of diluted host stock added to the guest measurement solution for each duplicate (g) |

### Measurement Settings
| Specifier | Description |
|----------|-------------|
| **T** | Temperature in degrC |
| **n** | Number of measurement sets (triplicate: `n = 3`) |
| **m** | Number of spectra measured per set |
| **lambda1** | First wavelength monitored (e.g., 650 nm) |
| **lambda2** | Second wavelength monitored (e.g., 715 nm) |
| **F0_free** | Boolean: if `True`, the program fits \( F_0 \); otherwise it is fixed |
| **smoothening** | Boolean: whether to smooth spectra before extracting data |
| **formula** | Fitting formula: `Stern_Volmer`, `Tsukube`, or `Connors` |
| **method** | Measurement mode: `"single_wl"` (single wavelength) or `"max_wl"` (peak maxima) |
| **Guest_add_i** | Volume of guest added per addition and per spectrum \(i\) (µL). |


---

## Fitting Program: Fitting Parameters

The script `bindfitnplot.py` fits a selection of parameters depending on the chosen model.  
These fitted parameters are written to the output file.  
All formulas assume 1:1 binding.

### 1. Stern–Volmer
| Parameter | Meaning |
|----------|----------|
| **Ka** | Association constant \( K_a \) |
| **F0** | Fluorescence of host solution before guest addition \( F_0 \) |

### 2. Tsukube
| Parameter | Meaning |
|----------|----------|
| **Ka** | Association constant \( K_a \) |
| **F0** | Initial fluorescence \( F_0 \) |
| **kDHG** | Proportionality constant \( k_{\Delta HG} = k_{HG} - k_{H} \) |

### 3. Connors
| Parameter | Meaning |
|----------|----------|
| **Ka** | Association constant \( K_a \) |
| **F0** | Initial fluorescence \( F_0 \) |
| **kH** | Proportionality constant of free host in presence of guest \( k_H \) |
| **kH0** | Proportionality constant for initial free host \( k_H^0 \) |
| **kHG** | Proportionality constant of the complex \( k_{HG} \) |


## Additional Settings

Besides selecting the fitting formula and choosing whether to optimize \(F₀\), the program provides two more user-adjustable settings:

### Smoothing
- **Enabled:** Noise is reduced before extracting intensities, preventing artificially high or low values.  
- **Disabled:** Raw peak intensities are used.

### Intensity Extraction Method
- **`single_wl`** — Extracts intensity at a fixed wavelength for each peak.  
  - Prevents errors if significant red/blue shifts occur.  
- **`max_wl`** — Automatically identifies the maximum intensity near the target wavelength.  
  - Corrects for small wavelength shifts that would otherwise cause systematically low intensities.

## Recommended Settings

Comparison tests show that the following combination provides the most reliable results:

- **Connors** model  
- **Smoothing** enabled  
- **Optimize \(F₀\)**  
- **`single_wl`** method  

Under these conditions, the program can follow changes in fluorescence at any wavelength without requiring quenching behaviour, as the silent-complex assumption is not applied.

## Output

The fitting program generates three sets of figures along with a detailed output file.

### 1. Spectra Overview
For each titration set, all recorded spectra are plotted together (intensity vs. wavelength).  
These superimposed plots make it easy to visualise how the fluorescence signal changes as increasing amounts of guest are added.

### 2. Binding Isotherms
For each titration and for each monitored wavelength, the program plots the binding isotherm (intensity vs. equivalents \([G]_0/[H]_0\)).  
These graphs show how the signal evolves throughout the titration and form the basis of the fitting procedure.

### 3. Fitted Curves
The third output figure overlays the fitted curves on top of the experimental binding isotherms.  
This comparison allows the user to verify the quality of the fit and determine whether the resulting association constants are reliable.

## Output File

In addition to the figures, the script produces an output file.  
This file contains:

- All settings used during the run  
- Input values  
- Initial parameter guesses  
- Calculated concentrations  
- Guest/host equivalents  
- Fitted parameters **per titration**  
- **Average** fitted parameters with **standard deviations**

The output file is automatically saved in the same directory that contains the fluorescence data.
