#ifndef SOLVER_FACTORY_HPP
#define SOLVER_FACTORY_HPP

#include "parameters/parameters.hpp"
#include "solver/solver.hpp"

#ifdef USE_CUDA
#include "solver/cuda_solver/gpu_it_solver.hpp"
#include "solver/cuda_solver/gpu_rt_solver.hpp"
#else
#include "solver/cpu_solver/cpu_it_solver.hpp"
#include "solver/cpu_solver/cpu_rt_solver.hpp"
#endif

class SolverFactory {
public:
    static std::unique_ptr<AbstractGrossPitaevskiSolver>
    create(AbstractSimulationMediator *mediator) {
        auto params = PhysicalParameters::getInstance();
#ifdef USE_CUDA
        switch (params->calc_strategy.type) {
        case CalcStrategy::Type::IMAGINARY_TIME:
            return std::make_unique<GpuITGrossPitaevskiSolver>(mediator);
            break;
        case CalcStrategy::Type::REAL_TIME:
            return std::make_unique<GpuRTGrossPitaevskiSolver>(mediator);
            break;
        }

        throw std::runtime_error("CAN'T PRODUCE GPU SOLVER");
        return std::make_unique<GpuITGrossPitaevskiSolver>(mediator);
#else
        switch (params->calc_strategy.type) {
        case CalcStrategy::Type::IMAGINARY_TIME:
            return std::make_unique<CpuITGrossPitaevskiSolver>(mediator);
            break;
        case CalcStrategy::Type::REAL_TIME:
            return std::make_unique<CpuRTGrossPitaevskiSolver>(mediator);
            break;
        }

        throw std::runtime_error("CAN'T PRODUCE CPU SOLVER");
        // TODO default to fix to throw
        return std::make_unique<CpuITGrossPitaevskiSolver>(mediator);
#endif
    }
};

#endif
