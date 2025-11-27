#include "include/file_manager.hpp"
#include <fstream>

FileManager::FileManager()
    : params(PhysicalParameters::getInstance()) {}


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


