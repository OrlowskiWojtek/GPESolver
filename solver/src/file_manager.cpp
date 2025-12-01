#include "include/file_manager.hpp"
#include "units.hpp"
#include <fstream>

FileManager::FileManager()
    : params(PhysicalParameters::getInstance()) {
}

void FileManager::set_data_pointer(StdMat3D<std::complex<double>> *data) {
    cpsi_data = data;
}

void FileManager::save_params() {
    nlohmann::json j;

    j["n_atoms"]            = params->n_atoms;
    j["load_initial_state"] = params->load_initial_state;
    j["m"]                  = UnitConverter::mass_au_to_Da(params->m);
    j["dd"]                 = UnitConverter::len_au_to_nm(params->dd);
    j["dx"]                 = UnitConverter::len_au_to_nm(params->dx);
    j["dy"]                 = UnitConverter::len_au_to_nm(params->dy);
    j["dz"]                 = UnitConverter::len_au_to_nm(params->dz);
    j["nx"]                 = params->nx;
    j["ny"]                 = params->ny;
    j["nz"]                 = params->nz;
    j["edd"]                = params->edd;

    std::ofstream file(PARAMS_FILENAME);
    file << j.dump(4);
    file.close();
}

void FileManager::load_params() {
    std::ifstream file(PARAMS_FILENAME);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open params file.");
    }

    nlohmann::json j;
    file >> j;
    file.close();

    params->n_atoms            = j["n_atoms"];
    params->load_initial_state = j["load_initial_state"];
    params->m                  = UnitConverter::mass_Da_to_au(j["m"]);
    params->dd                 = UnitConverter::len_nm_to_au(j["dd"]);
    params->edd                = j["edd"];
    params->dx                 = UnitConverter::len_nm_to_au(j["dx"]);
    params->dy                 = UnitConverter::len_nm_to_au(j["dy"]);
    params->dz                 = UnitConverter::len_nm_to_au(j["dz"]);
    params->nx                 = j["nx"];
    params->ny                 = j["ny"];
    params->nz                 = j["nz"];

    check_params();
}

void FileManager::save_initial_state() {
    std::ofstream file("initial_state.dat");

    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    if (!cpsi_data) {
        throw std::runtime_error("Data pointer not set.");
    }

    StdMat3D<std::complex<double>> &cpsi = *cpsi_data;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                file << cpsi(i, j, k).real() << "\t" << cpsi(i, j, k).imag() << "\n";
            }
        }
    }

    file.close();
}

void FileManager::load_initial_state() {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    if (!cpsi_data) {
        throw std::runtime_error("Data pointer not set.");
    }

    std::ifstream file("initial_state.dat");

    StdMat3D<std::complex<double>> &cpsi = *cpsi_data;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                double re, im;
                file >> re >> im;
                cpsi(i, j, k) = std::complex<double>(re, im);
            }
        }
    }

    file.close();
}

void FileManager::check_params() {
    if (params->nx <= 0 || params->ny <= 0 || params->nz <= 0) {
        throw std::runtime_error("Grid dimensions must be positive.");
    }
    if (params->dx <= 0 || params->dy <= 0 || params->dz <= 0) {
        throw std::runtime_error("Grid spacings must be positive.");
    }
    if (params->n_atoms <= 0) {
        throw std::runtime_error("Number of atoms must be positive.");
    }
}
