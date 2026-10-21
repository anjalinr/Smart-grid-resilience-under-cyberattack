# Smart Grid Resilience under Cyberattack

## Project Overview

This project investigates the detection and mitigation of False Data Injection Attacks (FDIA) in a smart grid using a Cyber-Physical System (CPS)-based monitoring approach.

The simulation is developed in MATLAB using the IEEE 14-bus power system. Coordinated attacks on demand measurements are studied to evaluate their effect on state estimation and the performance of residual-based attack detection.

## Objectives

- Model coordinated False Data Injection Attacks on selected demand measurements.
- Implement residual-based attack detection using a Kalman-filter approach.
- Evaluate meter isolation and measurement-based mitigation strategies.
- Compare estimation performance before and after mitigation.

## Methodology

1. Model the IEEE 14-bus power system and its measurements.
2. Introduce coordinated false data injection attacks.
3. Calculate residuals between measured and estimated values.
4. Apply a detection threshold to identify suspicious measurements.
5. Evaluate mitigation using meter isolation and trusted measurements.
6. Analyse the resulting estimation errors and detection performance.

## Tools and Technologies

- MATLAB
- IEEE 14-bus test system
- Kalman-filter-based residual detection
- False Data Injection Attack (FDIA) modelling

## Results

The simulation evaluates attack detection, false alarms, state-estimation error and demand-estimation error under the tested attack scenarios. Numerical results and plots are provided in the results folder.

## How to Run

1. Download or clone this repository.
2. Open MATLAB.
3. Add the project code folder and required subfolders to the MATLAB path.
4. Run the main simulation script identified in the project files.
5. Review the generated plots and output files.

Refer to the source code for the exact script names, input requirements and execution sequence.

## Limitations

The conclusions are based on the implemented simulation model and tested attack scenarios. The results should not be interpreted as proof of performance for every real-world power grid or cyberattack.



## Academic Use

This repository documents the simulation, methodology and results of an academic project on smart grid monitoring and resilience.

