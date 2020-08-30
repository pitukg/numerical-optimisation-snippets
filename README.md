# numerical-optimisation-snippets
Numerical Snippets for Gradient Descent in OCaml

#### Vectors

Stored as built-in linked-lists.

#### Matrices

Stored as a pair of vector of vectors: one storing the rows and one storing the columns.
This approach is taken so that matrix multiplication can be done in O(n^3) time.

#### Linear Systems

Gradient Descent requires a linear system to be solved.
This is achieved by calculating the __Cholesky-decomposition__ of the coefficient matrix and solving the system with the simpler upper and lower triangular matrices
(by forward or backward substitution).

#### Multidimensional Calculus

Newton's method and Gradient Descent require partial & total derivatives and the Hessian of the matrix.
The derivative approximation is based on __functional streams__ to calculate the limit of the sequence by stopping when the subsequent elements are less than _epsilon_ apart.

#### Gradient Descent and Newton's Method

Given the derivatives, Hessian and Cholesky-decomposition methods we can apply these algorithms to functions of arbitrary dimensions, thus finding a local minimum.
