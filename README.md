# Bindfitnplot
This repository contains the code written for my Bachelor's Thesis at the Nolte Group

# Publication
https://www.nature.com/articles/s41467-020-18596-1

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
| **m_Hstock** | Mass of host dissolved (g) |
| **V_solv_Hstock** | Volume of host stock solution (mL) |
| **m_solv_Hstock** | Mass of host stock solution (g) |

### Guest Stock Solution
| Specifier | Description |
|----------|-------------|
| **M_Guest** | Guest molecular weight (g/mol) |
| **m_Gstock** | Mass of guest dissolved (g) |
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
| **n** | Number of measurement sets (triplicate: `n = 3`) |
| **m** | Number of spectra measured per set |
| **lambda1** | First wavelength monitored (e.g., 650 nm) |
| **lambda2** | Second wavelength monitored (e.g., 715 nm) |
| **F0_free** | Boolean: if `True`, the program fits \( F_0 \); otherwise it is fixed |
| **smoothening** | Boolean: whether to smooth spectra before extracting data |
| **formula** | Fitting formula: `Stern_Volmer`, `Tsukube`, or `Connors` |
| **method** | Measurement mode: `"single_wl"` (single wavelength) or `"max_wl"` (peak maxima) |
| **Guest_add_i** | Volume of guest added per addition and per spectrum \(i\) (µL) |

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
