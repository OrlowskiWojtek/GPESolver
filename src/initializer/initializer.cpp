#include "initializer/initializer.hpp"
#include "output.hpp"
#include <tuple>

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
    case InitializationOption::Type::SETUP_GAUSS:
        init_with_setup_gaussian();
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

void DataInitializer::initialize_potential() {
    OutputFormatter::printInfo("Initializing potential");

    switch (params->pote_strategy.type) {
    case PotentialType::Type::REGULAR:
        set_pote_regular();
        break;
    case PotentialType::Type::MEXICAN:
        set_pote_mexican();
        break;
    case PotentialType::Type::CRADLE:
        set_pote_cradle();
        break;
    }

    init_pote();
    p_mediator->on_pote_initialized(_pote);
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

                double sigma_x = (params->nx * params->dx) / 20.;
                double sigma_y = (params->ny * params->dy) / 20.;
                double sigma_z = (params->nz * params->dz) / 20.;

                double val = std::exp(-0.5 * (x * x) / (sigma_x * sigma_x)) *
                             std::exp(-0.5 * (y * y) / (sigma_y * sigma_y)) *
                             std::exp(-0.5 * (z * z) / (sigma_z * sigma_z));

                std::complex<double> cval = std::complex<double>(val, 0.);
                _data(i, j, k)            = cval;
            }
        }
    }
}

void DataInitializer::init_with_setup_gaussian() {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;
    _data.resize(nx, ny, nz);

    int n_maximas = params->bec_droplets_x * params->bec_droplets_y * params->bec_droplets_z;
    std::vector<std::tuple<double, double, double>> centers(n_maximas);

    double x_tot = params->nx * params->dx;
    double y_tot = params->ny * params->dy;
    double z_tot = params->nz * params->dz;

    for (int i = 0; i < params->bec_droplets_x; i++) {
        for (int j = 0; j < params->bec_droplets_y; j++) {
            for (int k = 0; k < params->bec_droplets_z; k++) {
                int idx = i * params->bec_droplets_y * params->bec_droplets_z +
                          j * params->bec_droplets_z + k;

                double spacing_x = x_tot / static_cast<double>(params->bec_droplets_x + 1);
                double spacing_y = y_tot / static_cast<double>(params->bec_droplets_y + 1);
                double spacing_z = z_tot / static_cast<double>(params->bec_droplets_z + 1);

                double cx = -x_tot / 2 + spacing_x * (i + 1);
                double cy = -y_tot / 2 + spacing_y * (j + 1);
                double cz = -z_tot / 2 + spacing_z * (k + 1);

                centers[idx] = {cx, cy, cz};
            }
        }
    }

    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int k = 1; k < nz - 1; k++) {
                double x = p_sctx->get_x(i);
                double y = p_sctx->get_y(j);
                double z = p_sctx->get_z(k);

                double sigma_x = (params->nx * params->dx) / 20.;
                double sigma_y = (params->ny * params->dy) / 20.;
                double sigma_z = (params->nz * params->dz) / 20.;

                double val = 0.0;
                for (int m = 0; m < n_maximas; m++) {
                    double xc = x - std::get<0>(centers[m]);
                    double yc = y - std::get<1>(centers[m]);
                    double zc = z - std::get<2>(centers[m]);

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

            int center_y_idx = n_maximas / 2;
            double y_offset  = (idx - center_y_idx) * (params->ny * params->dy) / (n_maximas + 1.);

            centers_y[idx] = y_offset;
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

                double sigma_x = (params->nx * params->dx) / 20.;
                double sigma_y = (params->ny * params->dy) / 20.;
                double sigma_z = (params->nz * params->dz) / 20.;

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

void DataInitializer::init_pote() {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    _pote.resize(nx, ny, nz);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                _pote(i, j, k) = _pote_func(i, j, k);
            }
        }
    }
}

void DataInitializer::set_pote_mexican() {
    _pote_func = [this](int ix, int iy, int iz) -> double {
        double x = p_sctx->get_x(ix);
        double y = p_sctx->get_y(iy);
        double z = p_sctx->get_z(iz);

        double vx = -params->b * std::pow(x, 2) + params->aa * std::pow(x, 4);

        double vy = 0.5 * params->m * std::pow(y, 2) * std::pow(params->omega_y, 2);
        double vz = 0.5 * params->m * std::pow(z, 2) * std::pow(params->omega_z, 2);

        return vx + vy + vz;
    };
}

void DataInitializer::set_pote_regular() {
    _pote_func = [this](int ix, int iy, int iz) -> double {
        double x = p_sctx->get_x(ix);
        double y = p_sctx->get_y(iy);
        double z = p_sctx->get_z(iz);

        double vx = params->aa * std::pow(x, 4);
        double vy = 0.5 * params->m * std::pow(y, 2) * std::pow(params->omega_y, 2);
        double vz = 0.5 * params->m * std::pow(z, 2) * std::pow(params->omega_z, 2);

        return vx + vy + vz;
    };
}

void DataInitializer::set_pote_cradle() {
    _pote_func = [this](int ix, int iy, int iz) -> double {
        double x = p_sctx->get_x(ix);
        double y = p_sctx->get_y(iy);
        double z = p_sctx->get_z(iz);

        double bec_spacing = params->nx * params->dx / (params->bec_droplets_x + 1);
        double alpha       = params->m * std::pow(params->omega_x, 2) / (std::pow(bec_spacing, 2));
        double beta        = params->m * std::pow(params->omega_x, 2) / 2;

        double vx = alpha * std::pow(x, 2);

        for (int i = 0; i < params->bec_droplets_x; i++) {
            double x_pos = -params->nx * params->dx / 2 + bec_spacing * (i + 1);
            vx -= beta * std::pow(x - x_pos, 2);
        }

        double vy = 0.5 * params->m * std::pow(y, 2) * std::pow(params->omega_y, 2);
        double vz = 0.5 * params->m * std::pow(z, 2) * std::pow(params->omega_z, 2);

        return vx + vy + vz;
    };
}
