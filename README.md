# Statistical Distribution Simulation and MLE Performance Analysis in R

This repository contains R scripts for conducting Monte Carlo simulations to analyze the performance of Maximum Likelihood Estimators (MLEs) for two specific statistical distributions:

1.  **Exponential Power Skew-Normal (EPS) Distribution:** A flexible distribution for modeling positive-valued data.
2.  **Log-Exponential Power Skew-Normal (LEPS) Distribution:** A distribution often used in survival analysis due to its properties, here demonstrated in a regression context.

The simulations are designed to assess the bias, Mean Squared Error (MSE), and Absolute Error (AE) of the estimated parameters under various sample sizes.

---

## Table of Contents

* [Project Overview](#project-overview)
* [Repository Structure](#repository-structure)
* [Getting Started](#getting-started)
    * [Prerequisites](#prerequisites)
    * [Installation](#installation)
    * [Running the Simulations](#running-the-simulations)
* [Key Functions and Logic](#key-functions-and-logic)
    * [EPS Simulation (`eps_simulation.R`)](#eps-simulation-eps_simulationr)
    * [LEPS Simulation (`leps_simulation.R`)](#leps-simulation-leps_simulationr)
* [Simulation Results](#simulation-results)
* [Contributing](#contributing)
* [License](#license)
* [Contact](#contact)

---

## Project Overview

This project aims to provide reproducible R code for statistical simulation studies. Specifically, it focuses on:

* **Data Generation:** Implementing custom functions to generate random variates from the EPS and LEPS distributions using acceptance-rejection sampling (for EPS) and inverse transform sampling (for LEPS).
* **Maximum Likelihood Estimation:** Defining negative log-likelihood functions for both distributions to facilitate parameter estimation using numerical optimization (`optim` function in R).
* **Performance Evaluation:** Conducting Monte Carlo simulations to evaluate the statistical properties (Bias, MSE, AE) of the MLEs across different sample sizes, providing insights into their finite-sample behavior.

This work is particularly relevant for researchers and practitioners working with non-standard distributions or exploring the efficiency of estimation methods.

---

## Repository Structure
