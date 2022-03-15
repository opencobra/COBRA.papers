Scripts to reproduce the figures and tables within the following paper:

Fleming RMT, Haraldsdottir HS, Le HM, Vuong PT, Hankemeier T, Thiele I. Cardinality optimisation in constraint-based modelling: Application to human metabolism, 2022 (submitted). 

The interface to each optimisation problem, as well as the the Difference of Convex Cardinality Optimisation (DCCO) algorithm were implemented in MATLAB (Mathworks Inc) and integrated with the COBRA toolbox, a comprehensive, cross-platform software suite for constraint-based modelling of biochemical networks, available at: https://github.com/opencobra/cobratoolbox

Within the COBRA Toolbox, the optimizeCardinality.m function is the general purpose interface to solve cardinality optimisation problems using our implementation of the Difference of Convex Cardinality Optimisation (DCCO) Algorithm described in the paper.

Stoichiometric consistency testing is implemented by the function findStoichConsistentSubset.m

Leak and siphon testing is implemented by the function findMassLeaksAndSiphons.m

Fux consistency testing is implemented by the function findFluxConsistentSubset.m

Thermodynamic flux consistency testing is implemented by the function findThermoConsistentFluxSubset.m

Sparse flux balance analysis is implemented by the function sparseFBA.m 

Relaxed flux balance analysis is implemented by the function relaxedFBA.m. 

A set of tutorials illustrating each cardinality optimisation application are also available within the COBRA toolbox tutorial collection https://github.com/opencobra/COBRA.tutorials. 

The function checkModelProperties.m provides a single interface that analyses leaks/siphons as well as stoichiometric, flux, thermodynamic flux, and a previously developed implementation that analyses a form of kinetic consistency (Fleming RMT, Vlassis N, Thiele I, Saunders MA. Conditions for duality between fluxes and concentrations in biochemical networks. Journal of Theoretical Biology. 2016 Nov 21;409:1–10. Available from: http://www.sciencedirect.com/science/article/pii/S0022519316301679)
