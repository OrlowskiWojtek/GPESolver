#ifndef FILE_MANAGER_HPP
#define FILE_MANAGER_HPP

#include "include/context.hpp"
#include "include/params.hpp"
#include "mat3d/stdmat3d.hpp"
#include <complex>
#include <fstream>

// Idea: implement into strategy pattern
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

    void load_from_different_mesh();

    void save_checkpoint(int iter);

    void save_xy_to_file(int iter);
    void save_current_energies(int iter, BECEnergies& enes);

private:
    static const char PARAMS_FILENAME[];
    static const char INITIAL_STATE_FILENAME[];  
    static const char LAST_STATE_FILENAME[];  
    static const char ENE_FILENAME[];

    static const char XY_CUT_FILENAME[];
    static const char CHECKPOUT_FILENAME[];
    static const char FORT_MESH_FILENAME[];

    void init_filesystem();

    std::ofstream ene_file;

    StdMat3D<std::complex<double>>* cpsi_data;
    void check_params();

    PhysicalParameters* params;
};

#endif
