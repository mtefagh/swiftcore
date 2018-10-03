# swiftcore
## Usage
### `reconstruction = swiftCore(S, rev, coreInd [, solver])`
- **Inputs**
1. `S`: the associated sparse stoichiometric matrix
2. `rev`: the 0-1 vector with 1's corresponding to the reversible reactions
3. `coreInd`: the set of indices corresponding to the core reactions
4. `solver`: the LP solver to be used; the currently available options are either _gurobi_ or _linprog_ with the default value of _linprog_
- **Outputs**
5. `reconstruction`: the 0-1 indicator vector of the reactions constituting the consistent metabolic network reconstructed from the core reactions