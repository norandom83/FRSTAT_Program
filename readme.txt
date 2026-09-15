@ FRSTAT REPLICATION CODE

Requirements

MATLAB R2024a or later
R for the supplementary finite-sample simulations

Files

Main_FRSTAT_vF.m: Main MATLAB program
FRSTAT_Data.xlsx: Processed analysis inputs
functions/: Functions called by the main MATLAB program
simulation_betweencointegration.r: R script for the simulations in Section D.4, “Finite-sample performance,” of the Supplementary Material

Run

Open Main_FRSTAT_vF.m in MATLAB and select Run. A full run generates Figures 1–3, displays Tables 1–3 and Supplementary Table 3, and prints the associated numerical diagnostics.

After a full run, individual Figure or Table sections at the end of the main program can be run separately using the variables retained in the workspace.

Run simulation_betweencointegration.r in R to reproduce the finite-sample simulations in Section D.4 of the Supplementary Material.

Data
  Sample_Info              Sample years and grid-cell counts
  Forcing_Components       Annual effective radiative forcing components
  Baseline_Density         Annual baseline anomaly densities
  Hemispheric_Density      Equal-weight hemispheric densities
  Fixed_Grid_Density       Annual densities for the fixed 190-cell grid
  Breitung_CV              Simulated null draws for the within-forcing test
  Coverage_Annual          Annual coverage counts

  FRSTAT_Data.xlsx contains the processed density inputs used in the paper.
  The main program reproduces all reported downstream figures, tables, and
  numerical diagnostics from these inputs.
