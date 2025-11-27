#include "include/file_manager.hpp"
#include <fstream>

FileManager::FileManager()
    : params(PhysicalParameters::getInstance()) {}

void FileManager::set_data_pointer(StdMat3D<std::complex<double>>* data) {
    cpsi_data = data;
}

void FileManager::save_params() {
    nlohmann::json j;

    j["n_atoms"] = params->n_atoms;
    j["m"]       = params->m;
    j["dx"]      = params->dx;
    j["dy"]      = params->dy;
    j["dz"]      = params->dz;
    j["nx"]      = params->nx;
    j["ny"]      = params->ny;
    j["nz"]      = params->nz;

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

    params->n_atoms = j["n_atoms"];
    params->m       = j["m"];
    params->dx      = j["dx"];
    params->dy      = j["dy"];
    params->dz      = j["dz"];
    params->nx      = j["nx"];
    params->ny      = j["ny"];
    params->nz      = j["nz"];
}

void FileManager::save_initial_state() {
    std::ofstream file("initial_state.dat");

    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    if(!cpsi_data) {
        throw std::runtime_error("Data pointer not set.");
    }

    StdMat3D<std::complex<double>>& cpsi = *cpsi_data;
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

    if(!cpsi_data) {
        throw std::runtime_error("Data pointer not set.");
    }

    std::ifstream file("initial_state.dat");

    StdMat3D<std::complex<double>>& cpsi = *cpsi_data;
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
