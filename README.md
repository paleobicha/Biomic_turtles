# Macroevolutionary Processes in Turtles: Biomic Specialization and Historical Climatic Changes #

This repository contains the R scripts and data used to perform the analyses presented in the manuscript:

**"Macroevolutionary processes in Turtles (Testudines): a view from biomic specialization and historical climatic changes"**
Authors: Juan S. Thomas, Sara Gamboa, Manuel Hernández Fernández, Oscar Murillo, Jonathan S. Pelegrin

## Overview ##
The analyses in this study explore the macroevolutionary patterns in turtles (Testudines) with a focus on biomic specialization and the impact of historical climatic shifts on diversification. This repository provides all data and scripts necessary to replicate the analyses and figures presented in the manuscript.

## Contents ##
### Data Files ###
- matrix_complete.txt: Contains biomic occupation data for all species of Testudines included in the study. Each row corresponds to a species, and columns represent different biomes, detailing the occupancy for each species across biomic categories.

- matrix_phylogeny.txt: Includes biomic data for the turtle species included in the phylogeny established by Thomson, Spinks, and Shaffer (2021).

## R Scripts ##
The R scripts in this repository are organized to guide users through each major analysis step, including:

- **BSIfunction.R:** This script performs Monte Carlo analyses and generates tables and figures of the results. It calculates the Biomic Specialization Index (BSI) across turtle species, providing insights into biomic specialization patterns.
- **DRtortugas.R** Contains code for calculating species diversification rates (Following Jetz et al., 2012), generating figures of Diversification Rate (DR) vs. BSI and DR by biome. This script also includes the statistical analyses necessary for testing the significance of diversification rate differences across biomes and levels of biomic specialization.
