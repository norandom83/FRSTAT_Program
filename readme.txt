@ FRSTAT REPLICATION CODE

Requirements
  MATLAB R2024a or later

Files
  Main_FRSTAT_vF.m   Main program
  FRSTAT_Data.xlsx   Processed analysis inputs
  functions\          Functions called by the main program

Run
  Open Main_FRSTAT_vF.m in MATLAB and select Run. A full run opens
  Figures 1--3, displays Tables 1--3 and Supplement Table 1, and prints
  the reported diagnostics. 

  After a full run, any named Figure or Table section at the end of the
  main program can be run separately using the variables in the workspace.

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
