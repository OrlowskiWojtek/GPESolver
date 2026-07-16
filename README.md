# GPE SOLVER

GPE Solver is C++ software for solving GPE (Gross - Pitaevski Equation).
Solver is meant to be modular with the KISS principle taken as main design goal.

## GPE 

**Time dependent Grass Pitaevski Eqaution**:

$$
i\hbar \frac{\partial \psi(\mathbf{r}, t)}{\partial t} 
    = H_{\mathrm{eGP}} \psi(\mathbf{r}, t)
$$

with hamiltonian in form:

$$
H_{\mathrm{eGP}} =
    -\frac{\hbar^2}{2m}\nabla^2
    + V_{\mathrm{ext}}({\bf r},t) 
    + g |\Psi({\bf r},t)|^2 
    + \gamma(\varepsilon_{dd}) |\Psi({\bf r},t)|^3 
    + V_{dd}({\bf r},t), 
$$

## Dependencies

- **FFTW3** - Required for CPU builds
- **CUDA** - Optional, for GPU acceleration
- **JSON** - For parameter files
- **Julia with Makie** - For visualization/plotting (optional)

In debian-like system FFTW3 and JSON libraries can be simply installed by running:
```bash
sudo apt install libfftw3-dev libffw3-thread nlohmann-json3-dev
```

## Building the Project

```bash
# Ubuntu/Debian besides needed libraries
sudo apt install cmake build-essential

# For CUDA support, CUDA compatibility required
sudo apt install nvidia-cuda-toolkit
```

Tested on cards:
- NVIDIA GeForce RTX 3050
- Quadro P6000
- NVIDIA GeForce GTX 1050 

### CPU only build:
**CPU-only build (default):**
```bash
mkdir build && cd build
cmake ..
make
```
**With CUDA support:**
```bash
mkdir build && cd build
cmake -DUSE_CUDA=TRUE ..
make
```

## Running

At the start program reads input parameters file `gpe_params.json`.
Example input file has been shown below:

```gpe_params.json
{
    "dd": 1500.0,
    "dx": 120.0,
    "dy": 120.0,
    "dz": 600.0,
    "edd": 1.50,
    "load_filename": "initial_state",
    "iter_imag": 1000,
    "iter_real": 300000,
    "m": 163.929,
    "n_atoms": 20000.0,
    "nx": 128,
    "ny": 128,
    "nz": 32,
    "fftw_n_threads": 4,
    "calc_strategy": "IT",
    "init_strategy": "MULTIPLE_GAUSS",
    "pote_strategy": "MEXICAN",
    "initial_maximas": 2,
    "bec_droplets_x": 3,
    "bec_droplets_y": 3,
    "bec_droplets_z": 2,
    "omega_x": 60,
    "omega_y": 60,
    "omega_z": 120
}
```

Parameters can be described as:
| Parameter | Description | Unit |
|----------|------|-----------|
| `dx, dy, dz` | spatial grid step | nm |
| `nx, ny, nz` | spatial grid size | - |
| `n_atoms` | number of atoms in the system | - |
| `m` | mass of single atom | Da |
| `omega_x, y, z` | trap frequencies (used in generating all potentials) | Hz |
| `edd` | $\varepsilon_{dd}$ - parameter from Gross-Pitaevski Equation  | - |
| `iter_imag` | number of imaginary time iterations | - |
| `iter_real` | number of real time iterations | - |
| `calc_strategy` | type of calculations - options below. | - |
| `pote_strategy` | type of potential - options below. | - |
| `init_strategy` | type of initialization - options below. | - |
| `fftw_threads` | number of threads used for FFTW in CPU version. | - |
| `load_filename` | filename without extension in case of initializing with file (both binary and text). | - |
| `dd` | distance between double-well minimas. Used in MEXICAN hat potential | nm |

Potential options:
- HARMONIC (harmonic in every direction)
- MEXICAN (mexican hat with two minimas separated by dd)
- MEXICAN_FREE (mexican hat after removal of a barrier)
- CYLINDRICAL (harmonic in form of V(r, z), where $\omega_r$ = $\frac{\omega_x + \omega_y}{2}$)
- FREE (droplets bound only in z plane)
- MEXICAN_ASYMETRIC (MEXICAN but added slight linear asymetry)
- MEXICAN_ASYMETRIC_FREE (MEXICAN_FREE but added slight linear asymetry)

Calculations options:
- IT (imaginary time evolution)
- RT (real time evolution)
- FS (full simulation IT -> RT)

Initialization options:
- TEXT_FILE (loads from text file ended with .gpe.dat),
- BINARY_FILE (loads from binary file with extension .gpe.bin),
- SETUP_GAUSS (Many Gaussian-distributed condensate centers, evenly distributed in the x y z directions), uses bec_droplets_x, y, z parameters,
- MULTIPLE_GAUSS (Many Gaussian-distributed condensate centers, distributed in xy plane configurations), uses initial_maxias parameter,
- GAUSS (Single Gaussian density initialization in the middle) 
- COS (Starting from cosine function zeroed at the boundaries of the box) 

## Project Structure

Project tree structure can be presented as:

```
.
├── CMakeLists.txt          # Main build configuration
├── cmake/                  # CMake modules
├── src/                    # C++ source code
│   ├── solver/             # GPE solvers (CPU & GPU)
│   │   ├── fft_solvers/    # FFT-based solvers
│   │   ├── cpu_solver/     # CPU implementation
│   │   └── cuda_solver/    # GPU implementation
│   ├── parameters/         # Parameter handling
│   ├── initializer/      # Initial conditions
│   ├── context/          # Simulation context
│   └── ...
├── lib/                    # External libraries
├── tests/                  # Unit tests
├── plotters/               # Visualization utilities
└── main.cpp                # Entry point
```

Whole computational code has been enclosed within solver/ directory.

## Plotting

Plotting utilities are done in Julia and gnuplot.
All data analysis is done using Julia scripts, but final plots are often done using gnuplot.
