# GPE SOLVER

GPE Solver is C++ software for solving GPE (Gross - Pitaevski Equation).
Solver is meant to be modular with the KISS principle taken as main design goal.

## GPE 

**Time dependent Grass Pitaevski Eqaution**:
$$ i\hbar \partial_t \psi(\mathbf{r}, t) = \left[ \frac{-\hbar^2}{2m} + V(\mathbf{r}) + g |\psi|^2 + \psi_{dd} \right]\psi(\mathbf{r}, t) $$

## Dependencies

FFTW3 library is required for fast fourier transforms - [http://www.fftw.org/](http://www.fftw.org/)
Plotting is done with Julia code using excellent Makie library.

## TODOS:
- [*] make fft calculator
- [*] write visualiser / plotter (meaby live) (3d isolines?)
- [ ] Code refactor - managing object
- [ ] Fix gaussians
- [ ] Write initializator class
- [ ] Move iteration numbers to parameter file
- [ ] Write references
