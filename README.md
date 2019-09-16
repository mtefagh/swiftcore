# swiftcc
[swiftcc](https://github.com/mtefagh/swiftcore/blob/master/src/swiftcc.m) is a tool for finding the largest flux consistent subnetwork of the original metabolic network. [swiftccTest](https://github.com/mtefagh/swiftcore/blob/master/test/swiftccTest.m) provides a benchmark to compare its performance against [fastcc](https://wwwen.uni.lu/research/fstc/life_sciences_research_unit/research_areas/systems_biology/software).

## Usage
### `consistent = swiftcc(S, rev [, solver])`

### Inputs:
* `S`: the associated sparse **stoichiometric matrix**
* `rev`: the 0-1 vector with 1's corresponding to the **reversible reactions**

### Optional inputs
* `solver`: the **LP solver** to be used; the currently available options are _gurobi_, _linprog_, and _cplex_ with the default value of _linprog_. It fallbacks to the COBRA LP solver interface if another supported solver is called.

### Outputs:
* `consistent`: the 0-1 indicator vector of the reactions constituting the maximum **flux consistent metabolic subnetwork**

***

# swiftcore
[swiftcore](https://github.com/mtefagh/swiftcore/blob/master/src/swiftcore.m) is a tool for the context-specific reconstruction of genome-scale metabolic networks. [swiftcoreTest](https://github.com/mtefagh/swiftcore/blob/master/test/swiftcoreTest.m) and [weightedTest](https://github.com/mtefagh/swiftcore/blob/master/test/weightedTest.m) provide a benchmark to compare its performance against [fastcore](https://wwwen.uni.lu/research/fstc/life_sciences_research_unit/research_areas/systems_biology/software).

## Usage
### `[reconstruction, reconInd, LP] = swiftcore(model, coreInd, weights, tol, reduction [, solver])`

### Inputs
1. `model`: the **metabolic network** with fields:
* `.S` - the associated sparse **stoichiometric matrix**
* `.lb` - **lower bounds** on reaction rates
* `.ub` - **upper bounds** on reaction rates
* `.rxns` - the cell array of **reaction abbreviations**
* `.mets` - the cell array of **metabolite abbreviations**
2. `coreInd`: the set of indices corresponding to the **core reactions**
3. `weights`: the **weight vector** for the penalties associated with each reaction
4. `tol`: **zero-tolerance**, i.e., the smallest flux value considered nonzero
5. `reduction`: boolean enabling the **metabolic network reduction** preprocess

### Optional inputs
6. `solver`: the **LP solver** to be used; the currently available options are _gurobi_, _linprog_, and _cplex_ with the default value of _linprog_. It fallbacks to the COBRA LP solver interface if another supported solver is called.

### Outputs
1. `reconstruction`: the **flux consistent metabolic network** reconstructed from the core reactions
2. `reconInd`: the 0-1 **indicator vector** of the reactions constituting the reconstruction
3. `LP`: the number of **solved LPs**

# Requirements
To run the test files in the `test` folder, one needs to download the [FASTCORE](http://wwwen.uni.lu/recherche/fstc/life_sciences_research_unit/research_areas/systems_biology/software) package and the [Recon3D](https://www.vmh.life/files/reconstructions/Recon/3D.01/Recon3D_301.zip) model and add them to the Matlab path.

# License
The *swift* family is distributed under the [GNU General Public License v3.0](http://www.gnu.org/copyleft/gpl.html).
