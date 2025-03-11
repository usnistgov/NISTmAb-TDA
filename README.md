# NISTmAb-TDA

Code and data to accompany the manuscript 

**Spatial and Sequential Topological Analysis of Molecular Dynamics Simulations of IgG1 Fc Domains**

Melinda Kleczynski
melinda.kleczynski@nist.gov

Christina Bergonzo
christina.bergonzo@nist.gov

Anthony J. Kearsley
anthony.kearsley@nist.gov

## Software requirements

Certain software are identified in this paper in order to specify the data analysis procedure adequately. Such 
identification is not intended to imply recommendation or endorsement of any product or service by NIST, 
nor is it intended to imply that the software identified are necessarily the best available for the purpose.

### Python implementation 

* [Python programming language](https://www.python.org/)
* [Rust](https://www.rust-lang.org/)

#### Minimal dependencies (GCCD matrix construction)

* [NumPy](https://numpy.org/doc/stable/index.html)
* [Open Applied Topology (OAT)](https://github.com/OpenAppliedTopology/oat_python)
* [pandas](https://pandas.pydata.org/docs/index.html)
* [SciPy](https://docs.scipy.org/doc/scipy/index.html)

#### Needed for additional analysis (classification and visualization)

* [Matplotlib](https://matplotlib.org/stable/)
* [random](https://docs.python.org/3/library/random.html)
* [scikit-learn](https://scikit-learn.org/stable/index.html)

### Julia implementation 

* [Julia programming language](https://julialang.org/)
* [CSV](https://csv.juliadata.org/stable/index.html)
* [DataFrames](https://dataframes.juliadata.org/stable/)
* [Distances](https://github.com/JuliaStats/Distances.jl)
* [Distributions](https://juliastats.org/Distributions.jl/stable/)
* [Eirene](https://github.com/henselman-petrusek/Eirene.jl)
