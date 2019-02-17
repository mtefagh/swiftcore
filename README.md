# swiftcc
[swiftcc](https://github.com/mtefagh/swiftcore/blob/master/swiftcc.m) is a tool for finding the largest flux consistent subnetwork of the original metabolic network. [swiftccTest](https://github.com/mtefagh/swiftcore/blob/master/swiftccTest.m) provides a benchmark to compare its performance against [fastcc](https://wwwen.uni.lu/research/fstc/life_sciences_research_unit/research_areas/systems_biology/software).

## Usage
### `consistent = swiftcc(S, rev [, solver])`

### Inputs:
1. `S`: the associated sparse **stoichiometric matrix**
2. `rev`: the 0-1 vector with 1's corresponding to the **reversible reactions**

### Optional inputs
* `solver`: the **LP solver** to be used; the currently available options are _gurobi_, _linprog_, and _cplex_ with the default value of _linprog_

### Outputs:
* `consistent`: the 0-1 indicator vector of the reactions constituting the maximum **flux consistent metabolic subnetwork**

***

# swiftcore
[swiftcore](https://github.com/mtefagh/swiftcore/blob/master/swiftcore.m) is a tool for the context-specific reconstruction of genome-scale metabolic networks. [swiftcoreTest](https://github.com/mtefagh/swiftcore/blob/master/swiftcoreTest.m) provides a benchmark to compare its performance against [fastcore](https://wwwen.uni.lu/research/fstc/life_sciences_research_unit/research_areas/systems_biology/software).

## Usage
### `reconstruction = swiftcore(S, rev, coreInd, weights, reduction [, solver])`

### Inputs
1. `S`: the associated sparse **stoichiometric matrix**
2. `rev`: the 0-1 vector with 1's corresponding to the **reversible reactions**
3. `coreInd`: the set of indices corresponding to the **core reactions**
4. `weights`: the **weight vector** for the penalties associated with each reaction
5. `reduction`: boolean enabling the **metabolic network reduction** preprocess

### Optional inputs
* `solver`: the **LP solver** to be used; the currently available options are _gurobi_, _linprog_, and _cplex_ with the default value of _linprog_

### Outputs
* `reconstruction`: the 0-1 indicator vector of the reactions constituting the **flux consistent metabolic network** reconstructed from the core reactions

# License
The *swift* family is distributed under the [GNU General Public License v3.0](http://www.gnu.org/copyleft/gpl.html).
