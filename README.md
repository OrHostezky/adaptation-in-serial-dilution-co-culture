# Adaptation in a Serial-Dilution Co-culture

A consumer-resource modeling framework in a serial-dilution setup, where $m$ species compete for $p$ nutrients in a series of batches. At each batch, a nutrient mixture is supplied with a fixed composition $\lbrace c_i(t=0)\rbrace_{i=1}^p$ and total amount $c_0=\sum_i c_i(0)$. Similarly, a mixture of species $\lbrace\rho_\sigma(t=0)\rbrace_{\sigma=1}^m$ is added, with a fixed total amount $\rho_0=\sum_\sigma \rho_\sigma(0)$. Species grow until the depletion of nutrients, the culture is diluted, and a new batch starts, while keeping relative species populations.
Each species has a *metbolic strategy*, $\vec{\alpha}_ \sigma$, which is a vector of enzyme levels, while given a fixed enzyme budget $E_\sigma=\sum_i \alpha_{\sigma,i}$, leading to a metabolic trade-off. Nutrient consumption rates are given by $j_{\sigma,i}=\frac{c_i}{K+c_i}\alpha_{\sigma,i}$.

Dynamics within each batch are given by [^1]:
```math
\frac{d\rho_\sigma}{dt}=\rho_\sigma \sum_{i=1}^p j_{\sigma,i}   \quad,
```
```math
\frac{dc_i}{dt}=-\sum_{\sigma=1}^m \rho_\sigma j_{\sigma,i}   \quad.
```

The main feature of this framework is the inclusion of species *adaptation* to changing nutrient levels throughout the batch, by allowing the dynamics of the metabolic strategies, given by [^2]:
```math
\frac{d\alpha_{\sigma,i}}{dt}=(\mathbb{P}_{\sigma,i} E_\sigma-\alpha_{\sigma,i})\sum_{i'=1}^p j_{\sigma,i'}   \quad,
```
where $\mathbb{P}_{\sigma,i}$ is an indicator function which is 1 whenever species $\sigma$ produces enzyme $i$, and 0 otherwise. An adaptor population can only produce a single enzyme-type at a time.
The framework includes a few adaptation models (see [`app__simulations.m`](Code/app__simulations.m)) that operate under these guidelines, and is receptible to the addition of different models that may work in this context. Note that this adaptation feature is based on the *2-nutrient* case and thus is limited to $p=2$ (however, $p=1$ and $p>2$ dynmics can be simulated with no adaptation).

The dynamics are numerically solved for using *MATLAB*'s built-in `ode89` Runga-Kutta solver with adaptive step size. Simulations can be applied in a specific, instant manner, or by executing and collecting data in parallel from large sets of simulations, using the Split-Apply-Combine method (see [Script Index](#script-index) for protocol).



## Script Index

> [!NOTE]
> **Here, the general structure of this repository's [*Code*](Code/) section and workflow are described. For a more specific script description, look at each script specifically.**

* [`odefun.m`](Code/odefun.m) and the [`eventfun*.m`](Code/) functions (one for each adaptation model) are the most basic functions, used by `ode89` to solve for the dynamics. [`odefun.m`](Code/odefun.m) contains the actual dynamics, while the [`eventfun*.m`](Code/) functions are used to track events througout the dynamics.

* The [`sim__*.m`](Code/) functions carry out the actual simulations (for simulation types, see [`app__simulations.m`](Code/app__simulations.m)). [`sim__batch.m`](Code/sim__batch.m) is the most basic of these, which solves for the dynamics within each batch, and is executed by the other, higher [`sim__*.m`](Code/) functions (except for [`sim__invasibility_map.m`](Code/sim__invasibility_map.m) which uses it indirectly). All functions may optionally plot the results, and (except for [`sim__batch.m`](Code/sim__batch.m)) save the data. The high-end functions [`sim__serial__interbatch.m`](Code/sim__serial__interbatch.m), [`sim__serial__full_dynamics.m`](Code/sim__serial__full_dynamics.m), and [`sim__invasibility_map.m`](Code/sim__invasibility_map.m) are executed using the apply scripts (see below).

* The [`app__*.m`](Code/) scripts are the main scripts that apply the simulations, and are directly used by the user for setting parameter values and execute the simulations. [`app__simulations.m`](Code/app__simulations.m) is used for running specific simulations, whereas the [`app__slurm__*.cmd`](Code/) files are used for implementing large simulation sets in parallel (see [below](#split-apply-combine)).

* ff



### Split-Apply-Combine

### General repository scheme:
  $\quad$  **Apply scripts**  $\quad \longmapsto \quad$  **Simulation functions**  $\quad \longmapsto \quad$  **Step functions**  $\quad \longmapsto \quad$  **Basic interaction functions**



## Data and Plots

Data and corresponding figures are automatically saved with easily identifiable, corresponding file names (including the simulation type and important parameter values), in the [*Data*](Data/) and [*Plots*](Plots/) directories, respectively.




## References

[^1]: Amir Erez, Jaime G. Lopez, Benjamin G. Weiner, Yigal Meir, and Ned S. Wingreen. *Nutrient levels and trade-offs control diversity in a serial dilution ecosystem*. **eLife**, September 2020 [(Go to paper)](https://doi.org/10.7554/eLife.57790).

[^2]: Amir Erez, Jaime G. Lopez, Yigal Meir, and Ned S. Wingreen. *Enzyme regulation and mutation in a model serial-dilution ecosystem*. **Physical Review E**, October 2021 [(Go to papar)](https://pubmed.ncbi.nlm.nih.gov/34781576/).
