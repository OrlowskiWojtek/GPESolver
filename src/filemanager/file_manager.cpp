#include "filemanager/file_manager.hpp"
#include "context/context.hpp"
#include "nlohmann/json.hpp"
#include "output.hpp"
#include "units.hpp"
#include <fstream>

const char FileManager::PARAMS_FILENAME[]       = "gpe_params.json";
const char FileManager::TEXT_FILE_EXTENSION[]   = ".gpe.dat";
const char FileManager::BINARY_FILE_EXTENSION[] = ".gpe.bin";
const char FileManager::ENERGIES_FILENAME[]     = "energy.gpe.dat";

FileManager::FileManager(AbstractSimulationMediator *mediator)
    : mediator(mediator)
    , params(PhysicalParameters::getInstance()) {
}

FileManager::~FileManager() {
}

void FileManager::save_params() {
    OutputFormatter::printInfo("Saving simulation parameters to" + std::string(PARAMS_FILENAME));

    nlohmann::json j;

    j["n_atoms"]         = params->n_atoms;
    j["m"]               = UnitConverter::mass_au_to_Da(params->m);
    j["dd"]              = UnitConverter::len_au_to_nm(params->dd);
    j["dx"]              = UnitConverter::len_au_to_nm(params->dx);
    j["dy"]              = UnitConverter::len_au_to_nm(params->dy);
    j["dz"]              = UnitConverter::len_au_to_nm(params->dz);
    j["nx"]              = params->nx;
    j["ny"]              = params->ny;
    j["nz"]              = params->nz;
    j["edd"]             = params->edd;
    j["fftw_n_threads"]  = params->fftw_n_threads;
    j["calc_strategy"]   = params->calc_strategy.to_string();
    j["init_strategy"]   = params->init_strategy.to_string();
    j["load_filename"]   = params->load_filename;
    j["initial_maximas"] = params->n_gauss_max;
    j["iter_imag"]       = params->iter_imag;
    j["iter_real"]       = params->iter_real;
    j["omega_x"]         = UnitConverter::freq_au_to_Hz(params->omega_x);
    j["omega_y"]         = UnitConverter::freq_au_to_Hz(params->omega_y);
    j["omega_z"]         = UnitConverter::freq_au_to_Hz(params->omega_z);
    j["bec_droplets_x"]  = params->bec_droplets_x;
    j["bec_droplets_y"]  = params->bec_droplets_y;
    j["bec_droplets_z"]  = params->bec_droplets_z;

    std::ofstream file(PARAMS_FILENAME);
    file << j.dump(4);
    file.close();
}

void FileManager::load_params() {
    OutputFormatter::printInfo("Loading simulation parameters from: " +
                               std::string(PARAMS_FILENAME));

    std::ifstream file(PARAMS_FILENAME);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open params file.");
    }

    nlohmann::json j;
    file >> j;
    file.close();

    params->n_atoms = j["n_atoms"];
    params->m       = UnitConverter::mass_Da_to_au(j["m"]);

    params->dd = UnitConverter::len_nm_to_au(j["dd"]);
    params->dx = UnitConverter::len_nm_to_au(j["dx"]);
    params->dy = UnitConverter::len_nm_to_au(j["dy"]);
    params->dz = UnitConverter::len_nm_to_au(j["dz"]);
    params->nx = j["nx"];
    params->ny = j["ny"];
    params->nz = j["nz"];

    params->omega_x = UnitConverter::freq_Hz_to_au(j["omega_x"]);
    params->omega_y = UnitConverter::freq_Hz_to_au(j["omega_y"]);
    params->omega_z = UnitConverter::freq_Hz_to_au(j["omega_z"]);

    params->edd           = j["edd"];
    params->load_filename = j["load_filename"];
    params->iter_imag     = j["iter_imag"];
    params->iter_real     = j["iter_real"];

    params->fftw_n_threads = j["fftw_n_threads"];

    params->calc_strategy.from_string(j["calc_strategy"]);
    params->init_strategy.from_string(j["init_strategy"]);

    params->n_gauss_max    = j["initial_maximas"];
    params->bec_droplets_x = j["bec_droplets_x"];
    params->bec_droplets_y = j["bec_droplets_y"];
    params->bec_droplets_z = j["bec_droplets_z"];

    mediator->on_params_loaded();
    check_params();
}

void FileManager::save_to_text_file(const wavefunction_t &psi, std::string filename) {
    filename.append(TEXT_FILE_EXTENSION);
    OutputFormatter::printInfo("Saving to text file: " + std::string(filename));

    std::ofstream file(filename);

    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    int dx = UnitConverter::len_au_to_nm(params->dx);
    int dy = UnitConverter::len_au_to_nm(params->dy);
    int dz = UnitConverter::len_au_to_nm(params->dz);

    if (!file.is_open()) {
        OutputFormatter::printError("Can't open initial state file for writing.");
        return;
    }

    file << nx << "\n" << ny << "\n" << nz << "\n";
    file << dx << "\n" << dy << "\n" << dz << "\n";

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

    int dx;
    int dy;
    int dz;

    file >> dx >> dy >> dz;
    dx = UnitConverter::len_nm_to_au(dx);
    dy = UnitConverter::len_nm_to_au(dy);
    dz = UnitConverter::len_nm_to_au(dz);

    if (std::abs(dx - params->dx) > 1e-3 || std::abs(dy - params->dz) > 1e-3 ||
        std::abs(dy - params->dz) > 1e-3) {
        throw std::runtime_error("Grid size in the file do not match current parameters.");
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

    file.close();
    mediator->on_data_loaded(psi_loading_buffer);
}

void FileManager::load_from_text_file(std::string filename,
                                      std::function<void(wavefunction_t &)> callback) {
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

    int dx;
    int dy;
    int dz;

    file >> dx >> dy >> dz;
    dx = UnitConverter::len_nm_to_au(dx);
    dy = UnitConverter::len_nm_to_au(dy);
    dz = UnitConverter::len_nm_to_au(dz);

    if (std::abs(dx - params->dx) > 1e-3 || std::abs(dy - params->dz) > 1e-3 ||
        std::abs(dy - params->dz) > 1e-3) {
        throw std::runtime_error("Grid size in the file do not match current parameters.");
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

    file.close();
    if (callback) {
        callback(psi_loading_buffer);
    }
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
        if (params->n_atoms <= 0) {
            throw std::runtime_error("Number of atoms must be positive.");
        }
        if (params->m <= 0) {
            throw std::runtime_error("Mass must be positive.");
        }
        if (params->fftw_n_threads <= 0) {
            throw std::runtime_error("FFTW number of threads must be positive.");
        }
        if (params->n_gauss_max <= 0) {
            throw std::runtime_error("Number of Gaussian maxima must be positive.");
        }

        if (!is_fft_compatible(params->nx) || !is_fft_compatible(params->ny) ||
            !is_fft_compatible(params->nz)) {
            OutputFormatter::printWarning("Grid dimensions may be slow for FFTW. Consider using "
                                          "dimensions that factor into small primes (2, 3, 5, 7).");
        }
    }
}

bool FileManager::is_fft_compatible(int n) {
    // Check if n has only small prime factors (2, 3, 5, 7) which are efficient for FFTW
    int temp = n;
    while (temp % 2 == 0)
        temp /= 2;
    while (temp % 3 == 0)
        temp /= 3;
    while (temp % 5 == 0)
        temp /= 5;
    while (temp % 7 == 0)
        temp /= 7;
    return temp == 1;
}

void FileManager::save_to_binary_file(const wavefunction_t &psi, std::string filename) {
    filename.append(BINARY_FILE_EXTENSION);
    OutputFormatter::printInfo("Saving to binary file: " + filename);

    std::ofstream file(filename, std::ios::out | std::ios::binary);

    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    double dx = UnitConverter::len_au_to_nm(params->dx);
    double dy = UnitConverter::len_au_to_nm(params->dy);
    double dz = UnitConverter::len_au_to_nm(params->dz);

    if (!file.is_open()) {
        OutputFormatter::printError("Could not open last state file for writing.");
        return;
    }

    file.write(reinterpret_cast<char *>(&nx), sizeof(int));
    file.write(reinterpret_cast<char *>(&ny), sizeof(int));
    file.write(reinterpret_cast<char *>(&nz), sizeof(int));

    file.write(reinterpret_cast<char *>(&dx), sizeof(double));
    file.write(reinterpret_cast<char *>(&dy), sizeof(double));
    file.write(reinterpret_cast<char *>(&dz), sizeof(double));

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

    file.close();
    mediator->on_data_loaded(psi_loading_buffer);
}

void FileManager::load_from_binary_file(std::string filename,
                                        std::function<void(wavefunction_t &)> callback) {
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

    int dx_file, dy_file, dz_file;
    file.read(reinterpret_cast<char *>(&dx_file), sizeof(double));
    file.read(reinterpret_cast<char *>(&dy_file), sizeof(double));
    file.read(reinterpret_cast<char *>(&dz_file), sizeof(double));

    dx_file = UnitConverter::len_nm_to_au(dx_file);
    dy_file = UnitConverter::len_nm_to_au(dy_file);
    dz_file = UnitConverter::len_nm_to_au(dz_file);

    if (std::abs(dx_file - params->dx) > 1e-3 || std::abs(dy_file - params->dz) > 1e-3 ||
        std::abs(dy_file - params->dz) > 1e-3) {
        throw std::runtime_error("Grid size in the file do not match current parameters.");
    }

    psi_loading_buffer.resize(nx_file, ny_file, nz_file);

    for (size_t i = 0; i < psi_loading_buffer.size(); i++) {
        double real, imag;
        file.read(reinterpret_cast<char *>(&real), sizeof(double));
        file.read(reinterpret_cast<char *>(&imag), sizeof(double));

        psi_loading_buffer(i) = std::complex<double>(real, imag);
    }

    file.close();
    if (callback) {
        callback(psi_loading_buffer);
    }
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
