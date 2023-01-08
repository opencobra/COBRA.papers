# PolyRound
[![PyPI version](https://badge.fury.io/py/PolyRound.svg)](https://badge.fury.io/py/PolyRound)

Efficient random sampling in convex polytopes relies on a 'rounding' preprocessing step, in which the polytope is rescaled so that the width is as uniform as possible across different dimensions.
PolyRound rounds polytopes on the general form:

![equation](https://latex.codecogs.com/gif.latex?P&space;:=&space;\{x&space;\in&space;\mathcal{R}^n:&space;A_{eq}x&space;=&space;b_{eq},&space;A_{ineq}x&space;\leq&space;b_{ineq}\}) with matrices ![equation](https://latex.codecogs.com/gif.latex?A_{eq}&space;\in&space;\mathcal{R}^{m,n}) and ![equation](https://latex.codecogs.com/gif.latex?A_{ineq}\in&space;\mathcal{R}^{k,n}) and vectors ![equation](https://latex.codecogs.com/gif.latex?b_{eq}&space;\in&space;\mathcal{R}^{m}) and ![equation](https://latex.codecogs.com/gif.latex?b_{ineq}\in&space;\mathcal{R}^{k}).

This formulation often arises in Systems Biology as the flux space of a metabolic network.

As output, PolyRound produces a polytope on the form ![equation](https://latex.codecogs.com/gif.latex?P^{r}&space;:=&space;\{v&space;\in&space;\mathcal{R}^l:&space;A^{r}_{ineq}v&space;\leq&space;b^{r}_{ineq}\}) where ![equation](https://latex.codecogs.com/gif.latex?l&space;\leq&space;n) and the zero vector is a stricly interior point. For transforming points back to the original space, it also provides a matrix ![equation](https://latex.codecogs.com/gif.latex?S&space;\in&space;\mathcal{R}^{n,l}) and a vector ![equation](https://latex.codecogs.com/gif.latex?t&space;\in&space;\mathcal{R}^{n}), so that ![equation](https://latex.codecogs.com/gif.latex?x&space;=&space;Sv&space;&plus;&space;t).

Currently, PolyRound is supported for python 3.7 and 3.8.

PolyRound no longer depends on a Gurobi installation and uses optlang (https://github.com/opencobra/optlang) to delegate linear programs to GLPK in case Gurobi is not installed. However, PolyRound is more reliable with Gurobi. Free Gurobi licenses for academic use can be obtained at https://www.gurobi.com/. Once the license is installed, the easiest way to get gurobi to work in python is through Anaconda https://www.anaconda.com/. Installation of gurobi in a conda environment is done with "conda install -c gurobi gurobi".

An easy example of how to get started is presented in the jupyter notebook cells below.


They show how to: <br>
1) create a polytope object from a file path <br>
2) simplify, reduce, and round a polytope in separate steps, togehter with some printed checks <br>
3) simplify, reduce and round a polytope in one step <br>
4) save the rounded polytope in various formats

``` python
import os
from PolyRound.api import PolyRoundApi
from PolyRound.static_classes.lp_utils import ChebyshevFinder
from PolyRound.settings import PolyRoundSettings
from pathlib import Path
model_path = os.path.join("PolyRound", "models", "e_coli_core.xml")
```

``` python
# Create a settings object with the default settings.
settings = PolyRoundSettings()
```

``` python
# Import model and create Polytope object
polytope = PolyRoundApi.sbml_to_polytope(model_path)
```

``` python
# Remove redundant constraints and refunction inequality constraints that are de-facto equalities.
# Due to these inequalities, the polytope is empty (distance from chebyshev center to boundary is zero)
x, dist = ChebyshevFinder.chebyshev_center(polytope, settings)
print(dist)
simplified_polytope = PolyRoundApi.simplify_polytope(polytope)
# The simplified polytope has non-zero border distance
x, dist = ChebyshevFinder.chebyshev_center(simplified_polytope, settings)
print(dist)
```

``` python
transformed_polytope = PolyRoundApi.transform_polytope(simplified_polytope)
# The distance from the chebyshev center to the boundary changes in the new coordinate system
x, dist = ChebyshevFinder.chebyshev_center(transformed_polytope, settings)
print(dist)
```

``` python
rounded_polytope = PolyRoundApi.round_polytope(transformed_polytope)
# After rounding, the distance from the chebyshev center to the boundary is set to be close to 1
x, dist = ChebyshevFinder.chebyshev_center(rounded_polytope, settings)
print(dist)

# The chebyshev center can be back transformed into an interior point in the simplified space.
print(simplified_polytope.border_distance(rounded_polytope.back_transform(x)))

```

``` python
# simplify, transform and round in one call
one_step_rounded_polytope = PolyRoundApi.simplify_transform_and_round(polytope)
```

``` python
#save to hdf5
out_hdf5 = os.path.join("PolyRound", 'output', 'rounded_e_coli_core.hdf5')
PolyRoundApi.polytope_to_hdf5(one_step_rounded_polytope, out_hdf5)
#save to csv
out_csv_dir = os.path.join("PolyRound", 'output', 'e_coli_core')
Path(out_csv_dir).mkdir(parents=True, exist_ok=True)
PolyRoundApi.polytope_to_csvs(one_step_rounded_polytope, out_csv_dir)
```

``` python
# Special use case: remove redundant constraints without removing zero facettes. This will leave th polytope with its original border distance.
x, dist = ChebyshevFinder.chebyshev_center(polytope, settings)
print(dist)
settings.simplify_only = True
simplified_polytope = PolyRoundApi.simplify_polytope(polytope, settings=settings)
# The simplified polytope still has zero border distance
x, dist = ChebyshevFinder.chebyshev_center(simplified_polytope, settings)
print(dist)
```
