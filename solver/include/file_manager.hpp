#ifndef FILE_MANAGER_HPP
#define FILE_MANAGER_HPP

#include "include/context.hpp"
#include "include/params.hpp"
#include "mat3d/stdmat3d.hpp"
#include <complex>
#include <fstream>

class FileManager{
public:
    FileManager();
    ~FileManager();

    void set_data_pointer(StdMat3D<std::complex<double>>* data);

    void load_initial_state();
    void load_params();

    void save_params();
    void save_initial_state();

    void save_last_state();
    void load_last_state();

    void save_checkpoint(int iter);

    void save_xy_to_file(int iter);
    void save_current_energies(int iter, BECEnergies& enes);

private:
    static constexpr char PARAMS_FILENAME[] = "gpe_params.json";
    static constexpr char INITIAL_STATE_FILENAME[] = "initial_state.dat";  
    static constexpr char LAST_STATE_FILENAME[] = "last_state.bin";  
    static constexpr char ENE_FILENAME[] = "energy.dat";

    static constexpr char XY_CUT_FILENAME[] = "cut_xy_";
    static constexpr char CHECKPOUT_FILENAME[] = "checkpoint_";

    void init_filesystem();

    std::ofstream ene_file;

    StdMat3D<std::complex<double>>* cpsi_data;
    void check_params();

    PhysicalParameters* params;
};

#endif
