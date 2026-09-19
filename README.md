# Simulation-based comparison between Machine learning, classical and hybrid methods

This study investigates the prediction performance and fixed-effect estimation accuracy of classical statistical, machine learning and hybrid methods for non-normal longitudinal data with correlated observations. A comprehensive simulation study varying distributional assumptions, sample sizes and data dimensionality was carried. Results show that classical methods, particularly the Linear Mixed-effects Model (LMM) and Linear Quantile Mixed Model (LQMM), consistently achieve lower prediction error across most scenarios by appropriately accounting for within subject correlation. Notably, the hybrid Mixed-Effects Random Forest (MERF) emerges as a competitive alternative, occasionally outperforming pure machine learning approaches such as Quantile Random Forest (QRF) and Quantile Extreme Gradient Boost (QXG), suggesting that incorporating correlation structure into machine learning frameworks meaningfully improves predictive performance. Findings were validated on a real health dataset (LISS), confirming simulation conclusions. These results offer applied researchers clear guidance: classical mixed-effects models remain preferable for non-normal longitudinal data, while hybrid approaches represent a promising direction when machine learning flexibility is desired.

## Code availability

The simulation code and analysis scripts used in this project are available in this repository. The repository provides the full reproducible workflow used to generate the simulation study results.

## Reproducibility

The code in this repository can be used to reproduce the simulation scenarios and evaluate the performance of the competing methods under the study settings described in the manuscript.

## Citation

If you use the code or adapt the simulation framework for your own work, please cite this repository using the metadata in [CITATION.cff](CITATION.cff).

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.