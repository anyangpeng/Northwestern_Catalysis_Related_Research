# This folder contains the DFT calculations scripts to construct a Au[111] surface and the corresponding convergence test

### Objectives
The goal of this calculation is to get familiar with the DFT code, ASE package and VASP package specifically. 
### Method
To creat a Au surface, one must know the lattice constant for the bulk material. To get the optimzed lattice constant, a bulk Au structure that extends infinitely in the 3D space was build and allowed to relax (both positions and the cell volume were allowedd to relax).

The Au[111] surface was cut out from the bulk Au material with the optimized lattice constant. The convergence test was performed on static structures to shorten the computation time.

### Resutls
The optimized bulk Au has a lattice constant of 4.15402 Å.

  <img src="https://user-images.githubusercontent.com/66216181/109046243-230c9380-769a-11eb-9b49-58306968afc4.png" width="700" height="500">

The cutoff energy converged at 400 eV

<img src="https://user-images.githubusercontent.com/66216181/109047788-e477d880-769b-11eb-9885-9ec1eb4d1b40.png" width="440" height="300">

The number of K-points converged at 6

<img src="https://user-images.githubusercontent.com/66216181/109046663-9f06db80-769a-11eb-9912-bf5dfe088651.png" width="400" height="260">

Unfortunately, the number of layers did not converged, another test should be conducted with respect to the property of interest.

<img src="https://user-images.githubusercontent.com/66216181/109048314-92838280-769c-11eb-8b47-b549302ac07b.png" width="400" height="260">

