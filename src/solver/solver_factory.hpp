#ifndef SOLVER_FACTORY_HPP
#define SOLVER_FACTORY_HPP

#include "solver/solver.hpp"

#ifdef USE_CUDA
#include "solver/cuda_solver/gpu_solver.hpp"
#else
#include "solver/cpu_solver/cpu_solver.hpp"
#endif

class SolverFactory {
public:
    static std::unique_ptr<AbstractGrossPitaevskiSolver> create(AbstractSimulationMediator* mediator) {      
#ifdef USE_CUDA
        return std::make_unique<GpuGrossPitaevskiSolver>(mediator);
#else
        return std::make_unique<CpuGrossPitaevskiSolver>(mediator);
#endif
    }
};

#endif
