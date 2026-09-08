#include "initializer/initializer.hpp"
#include "output.hpp"
#include "parameters/potentials.hpp"
#include "parameters/wavefunctions.hpp"

DataInitializer::DataInitializer(AbstractSimulationMediator *_mediator)
    : params(PhysicalParameters::getInstance())
    , p_mediator(_mediator)
    , p_sctx(SimulationContext::getInstance()) {
}

void DataInitializer::initialize_wavefunction() {
    OutputFormatter::printInfo("Initializing wavefunction");

    if (params->wvf_key == "BINARY_FILE") {
        init_from_binary_file();
    } else if (params->wvf_key == "TEXT_FILE") {
        init_from_text_file();
    } else {
        init_wavefunction();
    }

    p_mediator->on_data_initialized(_data);
}

void DataInitializer::initialize_potential() {
    OutputFormatter::printInfo("Initializing potential");

    init_pote();
    p_mediator->on_pote_initialized(_pote);
}

void DataInitializer::init_from_text_file() {
    p_mediator->request_load_from_text(_data);
}

void DataInitializer::init_from_binary_file() {
    p_mediator->request_load_from_binary(_data);
}

void DataInitializer::init_pote() {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    _pote_func = PotentialRegistry::instance().get_function(params->pote_key);
    _pote.resize(nx, ny, nz);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                double x = p_sctx->get_x(i);
                double y = p_sctx->get_y(j);
                double z = p_sctx->get_z(k);

                _pote(i, j, k) = _pote_func(x, y, z);
            }
        }
    }
}

void DataInitializer::init_wavefunction() {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    _data_func = InitializerRegistry::instance().get_function(params->wvf_key);
    _data.resize(nx, ny, nz);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                double x = p_sctx->get_x(i);
                double y = p_sctx->get_y(j);
                double z = p_sctx->get_z(k);

                _data(i, j, k) = _data_func(x, y, z);
            }
        }
    }
}

void DataInitializer::change_potential(std::string pote_key) {
    params->pote_key = pote_key;

    init_pote();
    p_mediator->on_pote_initialized(_pote);
}
