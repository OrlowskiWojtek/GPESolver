#ifndef FILE_MANAGER_HPP
#define FILE_MANAGER_HPP

#include "include/params.hpp"
#include <nlohmann/json.hpp>

class FileManager{
public:
    FileManager();

    void load_initial_state();
    void load_params();

    void save_params();
    void save_initial_state();
private:

    static constexpr char PARAMS_FILENAME[] = "gpe_params.dat";
    static constexpr char INITIAL_STATE_FILENAME[] = "initial_state.dat";  

    PhysicalParameters* params;
};

#endif
