This repository contains the code used in the manuscript "Measures of genetic diversification in somatic tissues at bulk and single cell resolution suggest sources of unknown stochasticity". It is organized in the following manner:

1. The code for running Gillespie simulations of neutral mutation accumulation in the haematopoietic stem cell compartment is found in `src/Vafsim/vafsim.jl`, and is written in Julia.
2. The code for numerically evolving the partial differential equation for the site frequency spectrum of a haematopoietic stem cell compartment is found in `src/Vafdyn/vafdyn.jl`, and is written in Julia.
3. The code for performing the parameter fitting based on the site frequency spectrum PDE is found in `src/Vafdyn/inferencePipeline.jl`.
4. The code for simulating the mutational burden accumulation is found in `src/BurdenSimulation CPP/` and is written in C++.