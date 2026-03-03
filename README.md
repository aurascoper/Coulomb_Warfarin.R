# Coulomb_Warfarin.R 

An interactive computational framework for simulating and visualizing Warfarin-Human Serum Albumin (HSA) molecular interactions and system dynamics.

# Overview

Coulomb_Warfarin.R is a molecular dynamics and pharmacokinetic simulation tool built in R. It provides a real-time, interactive environment to explore how the anticoagulant drug Warfarin binds to Human Serum Albumin (HSA). By combining spatial distribution analysis with system dynamics, this framework bridges the gap between atomic-level binding visualization and macroscopic pharmacological effects.
Key Features

3D Molecular Visualization: Real-time, interactive rendering of the HSA protein structure (PDB: 1AO6) and a simplified Warfarin ligand using a Cartesian coordinate system.

Electron Density & Binding Analysis: * Identifies and clusters binding pocket coordinates specifically targeting the Sudlow Site I (residues 195-198).

Calculates and visualizes spatial electron density using distance-weighted Gaussian distributions.

Pharmacokinetic System Dynamics: Models temporal evolution of state variables, including free warfarin concentration, antibody dynamics, and erythrocyte metrics using threshold-based state transitions and deterministic decay.

Interactive Dashboard: Built with Shiny and Plotly to allow users to manipulate atomic visualization parameters (opacity, radius) alongside multi-metric time-series plots.

Tech Stack & Dependencies

This project relies on the following R ecosystem tools:

R - Core statistical computing environment

Shiny - Web application framework for the interactive UI

Bio3D - Molecular structure analysis and trajectory processing

RGL - 3D visualization system

Tidyverse - Data manipulation and processing

Plotly - Interactive time-series data visualization

Installation & Usage

Clone the repository:


    git clone https://github.com/aurascoper/Coulomb_Warfarin.R.git
    cd Coulomb_Warfarin.R

Install required R packages:
Open your R console or RStudio and run:
R

    install.packages(c("shiny", "rgl", "tidyverse", "plotly"))
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("bio3d")

Run the Application:


    # Example command to run the main simulation (adjust based on your exact script execution)
    shiny::runApp("Warfarin_Coulomb.R")
    # or
    shiny::runApp("Warfarin_Maxwell.R")

Architecture Under the Hood

The framework is divided into four main architectural pillars:

Molecular Structure Representation: Maintains spatial coordinates in a crystallographic reference frame, utilizing an eight-atom simplified representation for Warfarin.

Spatial Distribution: Implements a 3D grid system for distance-weighted density calculations.

Temporal Evolution: Handles the mathematical iteration of system states, computing linear response functions and deterministic decay.

Asynchronous Rendering Pipeline: Ensures the 3D molecular views and 2D pharmacokinetic plots update smoothly without blocking computational threads.

Future Roadmap

Expansion to model multiple ligand interactions and secondary binding sites.

Implementation of complex binding kinetics beyond linear decay approximations.

Performance optimizations including GPU acceleration and parallel computation for higher-resolution grid densities.

Contributing

Contributions, issues, and feature requests are welcome! If you're interested in molecular dynamics or R-based data visualization, feel free to fork the repository and submit a pull request.
