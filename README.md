# swiftcore
*swiftcore* is a tool for the context-specific reconstruction of genome-scale metabolic networks.

## Usage
### `reconstruction = swiftcore(S, rev, coreInd, weights [, solver])`

### Inputs
1. `S`: the associated sparse **stoichiometric matrix**
2. `rev`: the 0-1 vector with 1's corresponding to the **reversible reactions**
3. `coreInd`: the set of indices corresponding to the **core reactions**
4. `weights`: the **weight vector** for the penalties associated with each reaction

### Optional inputs
* `solver`: the **LP solver** to be used; the currently available options are either _gurobi_ or _linprog_ with the default value of _linprog_

***

### Outputs
* `reconstruction`: the 0-1 indicator vector of the reactions constituting the **consistent metabolic network** reconstructed from the core reactions

## License
swiftcore is distributed under the [GNU General Public License v3.0](http://www.gnu.org/copyleft/gpl.html).
