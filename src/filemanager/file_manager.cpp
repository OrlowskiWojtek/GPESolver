#include "filemanager/file_manager.hpp"
#include "context/context.hpp"
#include "output.hpp"
#include "nlohmann/json.hpp"
#include "units.hpp"
#include <fstream>

const char FileManager::PARAMS_FILENAME[]       = "gpe_params.json";
const char FileManager::TEXT_FILE_EXTENSION[]   = ".gpe.dat";
const char FileManager::BINARY_FILE_EXTENSION[] = ".gpe.bin";
const char FileManager::ENERGIES_FILENAME[]     = "energy.gpe.dat";

FileManager::FileManager(AbstractSimulationMediator *mediator)
    : params(PhysicalParameters::getInstance())
    , mediator(mediator) {
}

FileManager::~FileManager() {
}

void FileManager::save_params() {
    OutputFormatter::printInfo("Saving simulation parameters to" + std::string(PARAMS_FILENAME));

    nlohmann::json j;
    j["n_atoms"]            = params->n_atoms;
    j["m"]                  = UnitConverter::mass_au_to_Da(params->m);
    j["dd"]                 = UnitConverter::len_au_to_nm(params->dd);
    j["dx"]                 = UnitConverter::len_au_to_nm(params->dx);
    j["dy"]                 = UnitConverter::len_au_to_nm(params->dy);
    j["dz"]                 = UnitConverter::len_au_to_nm(params->dz);
    j["nx"]                 = params->nx;
    j["ny"]                 = params->ny;
    j["nz"]                 = params->nz;
    j["edd"]                = params->edd;
    j["fftw_n_threads"]     = params->fftw_n_threads;
    j["calc_strategy"]      = params->calc_strategy.to_string();
    j["initial_maximas"]    = params->n_gauss_max;

    std::ofstream file(PARAMS_FILENAME);
    file << j.dump(4);
    file.close();
}

void FileManager::load_params() {
    OutputFormatter::printInfo("Loading simulation parameters from: " + std::string(PARAMS_FILENAME));

    std::ifstream file(PARAMS_FILENAME);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open params file.");
    }

    nlohmann::json j;
    file >> j;
    file.close();

    params->n_atoms            = j["n_atoms"];
    params->m                  = UnitConverter::mass_Da_to_au(j["m"]);
    params->dd                 = UnitConverter::len_nm_to_au(j["dd"]);
    params->edd                = j["edd"];
    params->dx                 = UnitConverter::len_nm_to_au(j["dx"]);
    params->dy                 = UnitConverter::len_nm_to_au(j["dy"]);
    params->dz                 = UnitConverter::len_nm_to_au(j["dz"]);
    params->nx                 = j["nx"];
    params->ny                 = j["ny"];
    params->nz                 = j["nz"];
    params->fftw_n_threads     = j["fftw_n_threads"];
    params->n_gauss_max        = j["initial_maximas"];
    params->calc_strategy.from_string(j["calc_strategy"]);

    mediator->on_params_loaded();
    check_params();
}

void FileManager::save_to_text_file(const wavefunction_t &psi,  std::string filename) {
    filename.append(TEXT_FILE_EXTENSION);
    OutputFormatter::printInfo("Saving to text file: " + std::string(filename));

    std::ofstream file(filename);

    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    if (!file.is_open()) {
        OutputFormatter::printError("Can't open initial state file for writing.");
        return;
    }

    file << nx << "\n" << ny << "\n" << nz << "\n";

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                file << psi(i, j, k).real() << "\t" << psi(i, j, k).imag() << "\n";
            }
        }
    }

    file.close();
}

void FileManager::load_from_text_file(std::string filename) {
    filename.append(TEXT_FILE_EXTENSION);
    OutputFormatter::printInfo("Loading from text file: " + filename);

    std::ifstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("Could not open initial state file for reading.");
    }

    int nx;
    int ny;
    int nz;
    file >> nx >> ny >> nz;

    if (nx != params->nx || ny != params->ny || nz != params->nz) {
        throw std::runtime_error("Grid dimensions in the file do not match current parameters.");
    }

    psi_loading_buffer.resize(nx, ny, nz);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                double re, im;
                file >> re >> im;
                psi_loading_buffer(i, j, k) = std::complex<double>(re, im);
            }
        }
    }

    mediator->on_data_loaded(psi_loading_buffer);
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

void FileManager::save_to_binary_file(const wavefunction_t &psi, std::string filename) {
    filename.append(BINARY_FILE_EXTENSION);
    OutputFormatter::printInfo("Saving to binary file: " + filename);

    std::ofstream file(filename, std::ios::out | std::ios::binary);

    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    if (!file.is_open()) {
        OutputFormatter::printError("Could not open last state file for writing.");
        return;
    }

    file.write(reinterpret_cast<char *>(&nx), sizeof(int));
    file.write(reinterpret_cast<char *>(&ny), sizeof(int));
    file.write(reinterpret_cast<char *>(&nz), sizeof(int));

    for (size_t i = 0; i < psi.size(); i++) {
        auto val    = psi(i);
        double real = val.real();
        double imag = val.imag();

        file.write(reinterpret_cast<char *>(&real), sizeof(double));
        file.write(reinterpret_cast<char *>(&imag), sizeof(double));
    }

    file.close();
}

void FileManager::load_from_binary_file(std::string filename) {
    filename.append(BINARY_FILE_EXTENSION);
    OutputFormatter::printInfo("Loading from binary file: " + filename);

    std::ifstream file(filename, std::ios::in | std::ios::binary);

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

    psi_loading_buffer.resize(nx_file, ny_file, nz_file);

    for (size_t i = 0; i < psi_loading_buffer.size(); i++) {
        double real, imag;
        file.read(reinterpret_cast<char *>(&real), sizeof(double));
        file.read(reinterpret_cast<char *>(&imag), sizeof(double));

        psi_loading_buffer(i) = std::complex<double>(real, imag);
    }

    mediator->on_data_loaded(psi_loading_buffer);
    file.close();
}

void FileManager::save_energies(const energies_container_t &energies) {
    std::string FILENAME = std::string(ENERGIES_FILENAME);
    OutputFormatter::printInfo("Saving energies to: " + FILENAME);

    std::ofstream file(FILENAME);

    int iter = 0;
    for (auto &enes : energies) {
        file << iter << "\t" << enes.e_kin << "\t" << enes.e_pot << "\t" << enes.e_int << "\t"
             << enes.e_ext << "\t" << enes.e_bmf << "\t" << enes.e_total << std::endl;
        iter++;
    }

    file.close();
}
