#include "units.hpp"

double UnitConverter::ene_au_to_eV(double au) {
    return au * 27.211407953;
}

double UnitConverter::ene_au_to_meV(double au) {
    return au * 27211.407953;
}

double UnitConverter::ene_eV_to_au(double eV) {
    return eV * 0.0367492929;
}

double UnitConverter::ene_meV_to_au(double meV) {
    return meV * 0.0000367492929;
}

double UnitConverter::len_au_to_nm(double au) {
    return au * 0.05291772105;
}

double UnitConverter::len_nm_to_au(double nm) {
    return nm * 18.897261260649092;
}

double UnitConverter::mass_Da_to_au(double Da) {
    return Da / 0.000548579909;
}

double UnitConverter::mass_au_to_Da(double au) {
    return au * 0.000548579909;
}
