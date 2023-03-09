# BarrierRound

BarrierRound is a rounding algorithm implemented in `Matlab` that preserves sparsity. It is a significantly fast algorithm that maintains similar sampling time compared to existing rounding algorithms. 

BarrierRound takes as input a constrained-based model `P' = {x' | A' * x' = b' and l' <= x' <= u'}` and outputs a rounded constrained-based model `P = {x | A * x = b and l <= x <= u}`.

## 1. Preliminaries

### 1-1. Rounding algorithms

These files can reproduce the benchmark result in `BarrierRound: Log-barrier Rounding for Constraint-based Models`, 2023 (submitted). The paper compares BarrierRound, PolyRound, and IsoRound.
1. BarrierRound: our rounding algorithm based on a self-concordant barrier of a polytope.
2. PolyRound: an improved version of rounding via John's ellipsoid, presented in `PolyRound: polytope rounding for random sampling in metabolic networks`.
3. IsoRound: a rounding algorithm for isotropy. It first draws uniform samples from a constrained-model and then computes the empirical covariance, which will be applied to the model to obtain a rounded model that is close to isotropy.


### 1-2. Directory

1. `Instances`: constrained-based models. The folder `0raw` contains original constrained-based models, while `1testbed-polyRound`, `2testbed-logBarrier`, and `3testbed-isotropic` contain constrained-based models rounded by the three algorithms.
2. `1polyRound`: files for `PolyRound` that are downloaded from gitlab (https://gitlab.com/csb.ethz/PolyRound).
3. `2logBarrier`: files for `BarrierRound`.
4. `3isotropic`: files for `IsoRound`. This algorithm relies on Constrained Riemannian Hamiltonian Monte Carlo (CRHMC), a sampling algorithm downloaded from github (https://github.com/ConstrainedSampler/PolytopeSamplerMatlab).
5. `CHAR-test`: files for coordinate Hit-and-Run (CHAR), a sampling algorithm running on full-dimensional constrained-based models. This algorithm is used to measure the sampling time when the benchmark compares the effectiveness of each rounding algorithm.


## 2. Reproduction of the benchmark

### 2-1. `BarrierRound` in `2logBarrier`

Run `logBarrierRound_main.m` for rounding. It rounds each constrained-based model in `Instance/0raw`, finds the full-dimensional representations of rounded models, and save them with rounding time in `Instance/2testbed-logBarrier`.

Run `testpaper_logBarrier.m` for sampling via CHAR from the rounded models in `Instance/2testbed-logBarrier`. CHAR draws uniform samples from the rounded models  and then saves results (on the sampling time) in `2logBarrier/logBarrier_test`.


### 2-2. `IsoRound` in `3isotropic`

First of all, remove all files in `2logBarrier` from matlab path.

Run `isotropicRound_main.m` for rounding. It rounds each constrained-based model in `Instance/0raw`, finds the full-dimensional representations of rounded models, and save them with rounding time in `Instance/3testbed-isotropic`.

Run `testpaper_isotropic.m` for sampling via CHAR from the rounded models in `Instance/3testbed-isotropic`. CHAR draws uniform samples from the rounded models  and then saves results (on the sampling time) in `3isotropic/isotropic_test`.

### 2-3. `PolyRound` in `1polyRound`

For rounding, run the Jupyter notebook `polyRound_main.ipynb` in `1polyRound/PolyRound`. Unlike the previous two algorithms, it takes as input xml files (for constrained-based models) stored in `PolyRound/PolyRound/models` and `PolyRound/PolyRound/modelsMat`, and then outputs rounded models and rounding times in csv formats in `PolyRound/PolyRound/output`. One should run `change_csv_to_mat.m` in `1polyRound/PolyRound` in order to change the format from csv to mat. These mat files are saved in `Instances/1testbed-polyRound`.

Run `testpaper_polyRound.m` for sampling via CHAR from the rounded models in `Instance/1testbed-polyRound`. CHAR draws uniform samples from the rounded models  and then saves results (on the sampling time) in `1isotropic/PolyRound_test`.

### 2-4. Plots

Run `testpaper_plot.m` to draw figurs comparing the three algorithms in terms of rounding time and mixing time.

