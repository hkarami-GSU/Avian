# Avian Influenza Disease Control Optimization

Computational models for optimizing avian influenza control strategies across multiple poultry farms and in wild bird populations using optimal control theory and network disease dynamics.

## Project Structure

- **poultry_control/** - Optimal control strategies for disease mitigation across multiple **farm units**
  - Models multiple poultry farms as interconnected network nodes
  - CA/ - California multi-farm system
  - MN/ - Minnesota multi-farm system

- **poultry_inverse/** - Inverse modeling to infer transmission parameters from **aggregated farm-level** infection data
  - CA/ - California case study
  - MN/ - Minnesota case study

- **wildbird/** - Disease surveillance and forecasting for **wild bird infection incidence**
  - CA/ - California wild bird infection incidence
  - MN/ - Minnesota wild bird infection incidence

## Requirements

- MATLAB R2020a or later
- Optimization Toolbox
- ODE solver functions (ode15s)

## Data Inputs

Each regional folder requires:
- `observation_data.xlsx` - Observed poultry (farms) infection cases
- `transmission_rates.xlsx` - Time-dependent transmission rate data by group
- `wildbird_Iw_*.xlsx` - Wild bird infected cases I_w
- `wbird.xlsx` - Wild bird incidence data 

## Main Scripts

### poultry_control/*/control_*_beta_tdep.m
Optimizes poultry disease control strategies under different weight parameters (lambda values).

**Outputs:**
- SIC dynamics (Susceptible, Infected, Carrier) for each control strategy
- Control effort trajectories
- Incidence and cumulative case curves
- Cost functional values for each scenario
- Results saved to `Results/` directory

### poultry_inverse/*/Poultry_Inverse_*.m
Estimates transmission parameters from observed farm-level data.

### wildbird/*/wildbird_Inverse_*.m
Analyzes disease dynamics from wild bird infection incidence data.

## Usage

1. Navigate to desired regional folder (e.g., `poultry_control/CA/`)
2. Ensure data files are present in the same directory
3. Run the main MATLAB script:
   ```matlab
   control_CA_beta_tdep
   ```
4. Results will be saved to `Results/` subdirectory

