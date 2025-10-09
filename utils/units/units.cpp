#include "units.hpp"

double UnitConverter::ene_au_to_eV(double au){
    return au * 27.211407953;
}

double UnitConverter::ene_au_to_meV(double au){
    return au * 27211.407953;
}

double UnitConverter::ene_eV_to_au(double eV){
    return eV * 0.0367492929;
}

double UnitConverter::ene_meV_to_au(double meV){
    return meV * 0.0000367492929;
}


