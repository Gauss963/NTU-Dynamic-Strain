# Data Analysis for Cohesive Crack Modeling

This repository contains Python scripts to perform data analysis related to cohesive crack modeling and experimental strain measurements. Below is a brief description of the provided scripts, their functionalities, and how to use them.

---

## Table of Contents
1. [Introduction](#introduction)
2. [File Structure](#file-structure)
   - [CohesiveCrack.py](#cohesivecrackpy)
   - [DataProcessor.py](#dataprocessorpy)
   - [FolderActions.py](#folderactionspy)
   - [Analyzing Code (Example)](#analyzing-code-example)
3. [Getting Started](#getting-started)
4. [Usage](#usage)
5. [Dependencies](#dependencies)

---

## Introduction

The scripts in this repository focus on:
- Modeling cohesive crack behavior using analytical formulas (`CohesiveCrack.py`).
- Converting experimental strain gauge data from voltage to stress, and then applying signal processing operations like high-pass filtering (`DataProcessor.py`).
- Loading and plotting experimental strain data (`FolderActions.py`).
- An example analysis routine demonstrating how to fit model parameters to experimental data.

---

## File Structure

### CohesiveCrack.py
This file contains functions and definitions for modeling cohesive cracks in a solid under stress. Key functions include:
- **`alpha_s(C_f, C_s)` / `alpha_d(C_f, C_d)`**: Calculate the dimensionless parameters related to wave speeds.
- **`D(alpha_s, alpha_d)`**: Core function for the cohesive crack formulation.
- **`delta_sigma_xy` / `delta_sigma_xx`**: Compute shear stress fluctuations (`\Delta \sigma_{xy}` or `\Delta \sigma_{xx}`) at positions \((x, y)\).
- **`main()`**: Example driver that uses the above functions to plot the stress fluctuations for various \(y\)-offset values.

### DataProcessor.py
Contains utility functions for processing experimental data:
- **`voltage_to_strain(raw_voltage)`**: Converts voltage measurements to strain using a Wheatstone bridge configuration.
- **`shear_strain_to_stress(E, poisson_ratio, strain)`**: Converts shear strain to shear stress based on the shear modulus.
- **`highpass_filter(data, cutoff, fs, order=4)`**: Applies a high-pass Butterworth filter to remove low-frequency components from a signal.
- **`fitting_function(X_c, C_f, Gamma, x, y)`** / **`chi_square(X_c, Gamma, C_f, X, Y)`**: Demonstrates how to integrate the cohesive crack modeling function in a curve-fitting or parameter estimation routine.

### FolderActions.py
- Demonstrates loading experimental data from `.npz` files.
- Shows basic plotting of raw strain measurements over time.
- Includes references to the `DataProcessor.py` or `CohesiveCrack.py` modules for more advanced data analysis (though not explicitly shown in the snippet).

### Analyzing Code (Example)
Shows a practical workflow for:
1. Loading and preparing experimental data (strain measurements).
2. Converting from voltage to strain, then to stress.
3. Filtering the data with a high-pass filter to remove noise or baseline drift.
4. Identifying key points of interest (peaks, rupture speed, etc.).
5. Fitting the cohesive crack model to the experimental data using **`iminuit`** or a similar optimization library.
6. Visualizing both the raw and fitted data results.

---

## Getting Started

1. **Clone or download** this repository to your local machine.
2. **Ensure you have the required Python environment** set up (see [Dependencies](#dependencies)).
3. Run or modify any of the scripts according to your needs.

---

## Usage

1. **CohesiveCrack**  
   - To run the `main()` function that generates a plot of stress fluctuation:
     ```bash
     python CohesiveCrack.py
     ```
   - This will create a PDF file in the `./Plot/` directory by default.

2. **DataProcessor**  
   - Import the functions into your own script:
     ```python
     import DataProcessor

     strain_data = DataProcessor.voltage_to_strain(voltage_readings)
     stress_data = DataProcessor.shear_strain_to_stress(E, nu, strain_data)
     filtered_data = DataProcessor.highpass_filter(stress_data, cutoff_freq, sample_rate)
     ```
   - Use `fitting_function` and `chi_square` as needed to integrate modeling with experimental data.

3. **FolderActions**  
   - Edit `filename` and other parameters as necessary to load your `.npz` data files:
     ```python
     python FolderActions.py
     ```
   - The script will produce plots and potentially save them to a `../Plot/` directory.

4. **Analyzing Code (Example)**  
   - This script is a workflow demonstration of the end-to-end process: from loading data, filtering, parameter fitting (with `iminuit`), and plotting results.

---

## Dependencies

Make sure you have the following libraries installed:
- Python 3.8+  
- `numpy`  
- `matplotlib`  
- `scipy`  
- `iminuit` (for parameter fitting)  

Install dependencies using:
```bash
pip install numpy matplotlib scipy iminuit
