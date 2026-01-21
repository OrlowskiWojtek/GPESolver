#include "initializer/initializer.hpp"
#include "output.hpp"
#include "units.hpp"

DataInitializer::DataInitializer(AbstractSimulationMediator *_mediator)
    : params(PhysicalParameters::getInstance())
    , p_mediator(_mediator)
    , p_sctx(SimulationContext::getInstance()) {
}

void DataInitializer::initialize_wavefunction() {
    OutputFormatter::printInfo("Initializing wavefunction");

    switch (params->init_strategy.type) {
    case InitializationOption::Type::COS:
        init_with_cos();
        break;
    case InitializationOption::Type::GAUSS:
        init_with_gaussian();
        break;
    case InitializationOption::Type::MULTIPLE_GAUSS:
        init_with_multiple_gaussian();
        break;
    case InitializationOption::Type::FROM_BINARY_FILE:
        init_from_binary_file();
        break;
    case InitializationOption::Type::FROM_TEXT_FILE:
        init_from_text_file();
        break;
    }

    p_mediator->on_data_initialized(_data);
}

void DataInitializer::init_with_cos() {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;
    _data.resize(nx, ny, nz);

    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int k = 1; k < nz - 1; k++) {
                double x = p_sctx->get_x(i) - params->dd;
                double y = p_sctx->get_y(j);
                double z = p_sctx->get_z(k);

                double rrr = (static_cast<int>(params->nx / 2) * params->dx);

                double cos_x = std::cos(2 * M_PI * x / rrr);
                double cos_y = std::cos(2 * M_PI * y / rrr);
                double cos_z = std::cos(1 * M_PI * z / rrr);

                double val                = cos_x * cos_y * cos_z;
                std::complex<double> cval = std::complex<double>(val, 0.);
                _data(i, j, k)            = cval;
            }
        }
    }
}

void DataInitializer::init_with_gaussian() {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;
    _data.resize(nx, ny, nz);

    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int k = 1; k < nz - 1; k++) {
                double x = p_sctx->get_x(i);
                double y = p_sctx->get_y(j);
                double z = p_sctx->get_z(k);

                double sigma_x = (params->nx * params->dx) / 10.;
                double sigma_y = (params->ny * params->dy) / 10.;
                double sigma_z = (params->nz * params->dz) / 10.;

                double val = std::exp(-0.5 * (x * x) / (sigma_x * sigma_x)) *
                             std::exp(-0.5 * (y * y) / (sigma_y * sigma_y)) *
                             std::exp(-0.5 * (z * z) / (sigma_z * sigma_z));

                std::complex<double> cval = std::complex<double>(val, 0.);
                _data(i, j, k)            = cval;
            }
        }
    }
}

void DataInitializer::init_with_multiple_gaussian() {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;
    _data.resize(nx, ny, nz);

    int n_maximas = params->n_gauss_max;
    std::vector<double> centers_x(n_maximas);
    std::vector<double> centers_y(n_maximas);

    for (int idx = 0; idx < n_maximas; idx++) {
        centers_x[idx] = idx % 2 ? params->dd : -params->dd;

        if (n_maximas % 2 == 1) {
            double y_offset = (idx + 1) * (params->ny * params->dy) / (n_maximas + 1);

            centers_y[idx] = p_sctx->get_y(0) + y_offset;
        }

        if (n_maximas % 2 == 0) {
            int y_idx       = (idx / 2 + 1);
            int y_maximas   = (n_maximas / 2 + 1);
            double y_offset = y_idx * (params->ny * params->dy) / y_maximas;

            centers_y[idx] = p_sctx->get_y(0) + y_offset;
        }
    }

    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int k = 1; k < nz - 1; k++) {
                double x = p_sctx->get_x(i);
                double y = p_sctx->get_y(j);
                double z = p_sctx->get_z(k);

                double sigma_x = (params->nx * params->dx) / 10.;
                double sigma_y = (params->ny * params->dy) / 10.;
                double sigma_z = (params->nz * params->dz) / 10.;

                double val = 0.0;
                for (int m = 0; m < n_maximas; m++) {
                    double xc = x - centers_x[m];
                    double yc = y - centers_y[m];
                    double zc = z;

                    val += std::exp(-0.5 * (xc * xc) / (sigma_x * sigma_x)) *
                           std::exp(-0.5 * (yc * yc) / (sigma_y * sigma_y)) *
                           std::exp(-0.5 * (zc * zc) / (sigma_z * sigma_z));
                }

                std::complex<double> cval = std::complex<double>(val, 0.);
                _data(i, j, k)            = cval;
            }
        }
    }
}

void DataInitializer::init_from_text_file() {
    p_mediator->request_load_from_text(_data);
}

void DataInitializer::init_from_binary_file() {
    p_mediator->request_load_from_binary(_data);
}
