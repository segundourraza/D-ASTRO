# D-ASTRO: Determination of Aerocapture Successful Trajectories and Robust Optimization

Welcome to the D-ASTRO algorithm repository! 
<!---
This repo contains the algorithms developed in ["A New Algorithm for Aerocapture Mission Design and Multi-disciplinary Design Optimisation (2024)"](https://arc.aiaa.org/loi/aiaaj) (NOT YET PUBLISHED).
--->
## Description

D-ASTRO is designed to address the complex challenges in the field of aerocapture maneuvers. This innovative tool redefines the optimization framework for these maneuvers, filling a significant gap in current research. Developed in MATLAB, D-ASTRO provides an efficient and robust solution for the rapid analysis and optimization of generic aerocapture missions.

Key features of D-ASTRO include its unique capability to quickly optimize non-thrusting, fixed angle-of-attack spacecraft and trajectories. It requires minimal input, needing only basic aeroshell information and mission parameters. Utilizing a nonlinear, single-objective constrained optimization strategy, D-ASTRO identifies the optimal flight path angle (γ) at the atmospheric interface (AI) and the fixed ballistic coefficient (β) of the aeroshell.

## Getting Started

### Dependencies
1. MATLAB 2023b, or newer
2. MATLAB Optimization Toolbox
3. MATLAB Parallel Computing Toolbox (Optional)

### Installing

Simply clone this repo!

### Executing program

* D-ASTRO is controlled by two main files:
  * **inputSim.m** is used to modified and control general parameters of the algorithm. These are:
    * Planetray constants
    * Atmospheric model
    * Spacecraft characteristics
    * Insertion trajectory parameters
    * Optimizer and Ordinary Differential Equation (ODE) solver tolerances and options
    * Binary search algorithms tolerances and options
  * **Run_DASTRO.m** defines problem specific constants, such as:
    * Saving and plotting options
    * Name of output files 
    * β range of interets
    * Array of γ (If D-ASTRO is being used for only trajectory anaysis).
    * Optimization weights
    * Target operational orbit
* D-ASTRO is executed by running **Run_DASTRO.m** script
* After completion, a new folder named "Output files" is added to the working directory where a ".mat" file will be saved. If individual tajectory simulations are performed, these will be saved to separate ".csv" files with the follwoing column order:

| Time | Radial distance | Velocity | Flight path angle | Longitude | Latitude | Heading angle | Heat Transfer | Heat load |
| --- |:---:|:---:|:---:|:---:|:---:|:---:|:---:|---:|
| s | m | m/s | radians | radians | radians | radians | W/m^2 | J/m^2 |

## Authors
Segundo Urraza Atue

Linkedin: [@segundourrazaatue](https://www.linkedin.com/in/segundourrazaatue/)

Email: segundourraza@gmail.com

## Version History
* 0.1
  * Initial Release

## License

This project is licensed under the MIT License - see the LICENSE.md file for details
