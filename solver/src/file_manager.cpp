#include "include/file_manager.hpp"
#include "include/output.hpp"
#include "units.hpp"
#include <fstream>
#include "nlohmann/json.hpp"

FileManager::FileManager()
    : params(PhysicalParameters::getInstance()) {

    ene_file.open(ENE_FILENAME, std::ios::out);

    init_filesystem();
}

FileManager::~FileManager() {
    if (ene_file.is_open()) {
        ene_file.close();
    }
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
    j["calc_strategy"]      = params->calc_strategy.to_string();

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
    params->calc_strategy.from_string(j["calc_strategy"]);

    check_params();
}

void FileManager::save_initial_state() {
    std::ofstream file(INITIAL_STATE_FILENAME);

    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    if (!cpsi_data) {
        OutputFormatter::printError("Data pointer not set.");
        return;
    }

    if (!file.is_open()) {
        OutputFormatter::printError("Can't open initial state file for writing.");
        return;
    }

    file << nx << "\n" << ny << "\n" << nz << "\n";

    StdMat3D<std::complex<double>> &cpsi = *cpsi_data;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                file << cpsi(i, j, k).real() << "\t" << cpsi(i, j, k).imag() << "\n";
            }
        }
    }

    OutputFormatter::printInfo("Saved initial state to " + std::string(INITIAL_STATE_FILENAME));
    file.close();
}

void FileManager::load_initial_state() {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    if (!cpsi_data) {
        throw std::runtime_error("Data pointer not set.");
    }

    std::ifstream file(INITIAL_STATE_FILENAME);

    if (!file.is_open()) {
        throw std::runtime_error("Could not open initial state file for reading.");
    }

    file >> nx >> ny >> nz;

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

    OutputFormatter::printInfo("Loaded initial state from " + std::string(INITIAL_STATE_FILENAME));
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

void FileManager::save_last_state() {
    std::ofstream file(LAST_STATE_FILENAME, std::ios::out | std::ios::binary);

    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    if (!cpsi_data) {
        OutputFormatter::printError("Data pointer not set.");
        return;
    }

    if (!file.is_open()) {
        OutputFormatter::printError("Could not open last state file for writing.");
        return;
    }

    file.write(reinterpret_cast<char *>(&nx), sizeof(int));
    file.write(reinterpret_cast<char *>(&ny), sizeof(int));
    file.write(reinterpret_cast<char *>(&nz), sizeof(int));

    for (size_t i = 0; i < cpsi_data->size(); i++) {
        auto val    = (*cpsi_data)(i);
        double real = val.real();
        double imag = val.imag();

        file.write(reinterpret_cast<char *>(&real), sizeof(double));
        file.write(reinterpret_cast<char *>(&imag), sizeof(double));
    }

    OutputFormatter::printInfo("Saved data to " + std::string(LAST_STATE_FILENAME));
    file.close();
}

void FileManager::load_last_state() {
    std::ifstream file(LAST_STATE_FILENAME, std::ios::in | std::ios::binary);

    if (!file.is_open()) {
        throw std::runtime_error("Could not open last state file for reading.");
    }

    int nx_file, ny_file, nz_file;
    file.read(reinterpret_cast<char *>(&nx_file), sizeof(int));
    file.read(reinterpret_cast<char *>(&ny_file), sizeof(int));
    file.read(reinterpret_cast<char *>(&nz_file), sizeof(int));

    if (nx_file != params->nx || ny_file != params->ny || nz_file != params->nz) {
        throw std::runtime_error("Grid dimensions in the file do not match current parameters.");
    }

    if (!cpsi_data) {
        throw std::runtime_error("Data pointer not set.");
    }

    for (size_t i = 0; i < cpsi_data->size(); i++) {
        double real, imag;
        file.read(reinterpret_cast<char *>(&real), sizeof(double));
        file.read(reinterpret_cast<char *>(&imag), sizeof(double));

        (*cpsi_data)(i) = std::complex<double>(real, imag);
    }

    OutputFormatter::printInfo("Loaded data from " + std::string(LAST_STATE_FILENAME));
    file.close();
}

void FileManager::save_xy_to_file(int iter) {
    std::ofstream file(std::string(XY_CUT_FILENAME) + std::to_string(iter) + ".dat");

    int z_zero_idx = params->nz / 2 + 1;

    if (!cpsi_data) {
        OutputFormatter::printError("Data pointer not set.");
        return;
    }

    if (!file.is_open()) {
        OutputFormatter::printError("Could not open cut file for writing.");
        return;
    }

    file << "# Parameters:\n";
    file << "# Number of atoms: " << params->n_atoms << "\n";
    file << "# nx: " << params->nx << "\n";
    file << "# ny: " << params->ny << "\n";
    file << "# nz: " << params->nz << "\n";
    file << "# dx (nm): " << UnitConverter::len_au_to_nm(params->dx) << "\n";
    file << "# dy (nm): " << UnitConverter::len_au_to_nm(params->dy) << "\n";
    file << "# dz (nm): " << UnitConverter::len_au_to_nm(params->dz) << "\n";
    file << "# eps_dd: " << params->edd << "\n";
    file << "# X\tY\t|Psi|^2\n";

    StdMat3D<std::complex<double>> &cpsi = *cpsi_data;

    for (int i = 0; i < params->nx; i++) {
        for (int j = 0; j < params->ny; j++) {
            file << params->get_x(i) << "\t" << params->get_y(j) << "\t"
                 << std::norm(cpsi(i, j, z_zero_idx)) << "\n";
        }
    }

    OutputFormatter::printInfo("Saved XY cut to cut_xy_" + std::to_string(iter) + ".dat");
    file.close();
}

void FileManager::save_current_energies(int iter, BECEnergies &enes) {
    if (!ene_file.is_open()) {
        OutputFormatter::printError("Energy file not opened for writing.");
        return;
    }

    ene_file << iter << "\t" << enes.e_kin << "\t" << enes.e_pot << "\t" << enes.e_int << "\t"
             << enes.e_ext << "\t" << enes.e_bmf << "\t" << enes.e_total << std::endl;

    ene_file.close();
}

void FileManager::save_checkpoint(int iter) {
    std::ofstream file(std::string(CHECKPOUT_FILENAME) + std::to_string(iter) + ".bin",
                       std::ios::out | std::ios::binary);

    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    if (!cpsi_data) {
        OutputFormatter::printError("Data pointer not set.");
        return;
    }

    if (!file.is_open()) {
        OutputFormatter::printError("Could not open checkpoint file for writing.");
        return;
    }

    file.write(reinterpret_cast<char *>(&nx), sizeof(int));
    file.write(reinterpret_cast<char *>(&ny), sizeof(int));
    file.write(reinterpret_cast<char *>(&nz), sizeof(int));

    for (size_t i = 0; i < cpsi_data->size(); i++) {
        auto val    = (*cpsi_data)(i);
        double real = val.real();
        double imag = val.imag();

        file.write(reinterpret_cast<char *>(&real), sizeof(double));
        file.write(reinterpret_cast<char *>(&imag), sizeof(double));
    }

    OutputFormatter::printInfo("Saved checkpoint to " +
                               std::string(CHECKPOUT_FILENAME) + std::to_string(iter) + ".bin");
    file.close();
}

void FileManager::init_filesystem() {
    // Can't decide on a specific filesystem without external libraries,
    // now do nothing.
}
