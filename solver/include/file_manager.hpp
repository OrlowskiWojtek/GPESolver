#ifndef FILE_MANAGER_HPP
#define FILE_MANAGER_HPP

#include "include/params.hpp"
#include "mat3d/stdmat3d.hpp"
#include <complex>
#include <nlohmann/json.hpp>

class FileManager{
public:
    FileManager();

    void set_data_pointer(StdMat3D<std::complex<double>>* data);

    void load_initial_state();
    void load_params();

    void save_params();
    void save_initial_state();


private:
    static constexpr char PARAMS_FILENAME[] = "gpe_params.json";
    static constexpr char INITIAL_STATE_FILENAME[] = "initial_state.dat";  

    StdMat3D<std::complex<double>>* cpsi_data;
    void check_params();

    PhysicalParameters* params;
};

#endif
