# SIRD Model Calibration for 90-Day Epidemic Data

This repository implements a SIRD (Susceptible, Infectious, Recovered, Deceased) model to simulate and calibrate epidemic dynamics over a 90-day period. Using Euler integration and a grid search method, the project calibrates the model parameters to best fit real-world epidemic data.

## Overview

- **Model:** SIRD epidemic model with differential equations:
  - **Susceptibles (S):**  \( \frac{dS}{dt} = -\beta S I \)
  - **Infectés (I):**      \( \frac{dI}{dt} = \beta S I - \gamma I - \mu I \)
  - **Rétablis (R):**     \( \frac{dR}{dt} = \gamma I \)
  - **Décès (D):**        \( \frac{dD}{dt} = \mu I \)
- **Numerical Method:** Euler integration.
- **Parameter Calibration:** Grid search over ranges of \(\beta\), \(\gamma\), and \(\mu\) using the RMSE (Root Mean Square Error) metric.
- **Data:** A CSV file containing 90 days of epidemic data with columns: `Jour`, `Susceptibles`, `Infectés`, `Rétablis`, and `Décès`.

## Repository Contents

- **`sird_model.ipynb`**  
  A Jupyter Notebook that contains the complete implementation, including:
  - Data import and preprocessing.
  - Implementation of the SIRD model.
  - Parameter grid search.
  - Visualization of simulation results versus ground truth data.

- **`population.csv`** (or similar)  
  Your dataset file containing 90 days of epidemic data.  
  **Note:** Make sure the CSV file is formatted correctly and placed in the appropriate directory or update the file path in the notebook.

- **`README.md`**  
  This file, providing an overview of the project and instructions for use.

## Setup and Installation

1. **Clone the Repository:**

   ```bash
   git clone https://github.com/yourusername/your-repo-name.git
   cd your-repo-name
