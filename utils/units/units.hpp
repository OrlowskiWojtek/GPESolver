#ifndef UNITS_HPP
#define UNITS_HPP

/*! class UnitConverter
 * \brief class containing static unit conversion instructions.
 *
 *  All instructions are meant to by static, without extra variables.
 */
class UnitConverter {
public:
    static double ene_eV_to_au(double eV);
    static double ene_meV_to_au(double meV);
    static double ene_au_to_eV(double au);
    static double ene_au_to_meV(double au);

    static double len_nm_to_au(double nm);
    static double len_au_to_nm(double au);

    //! Da - daltons (units aka 1/12 of C12)
    static double mass_Da_to_au(double Da);
    //! Da - daltons (units)
    static double mass_au_to_Da(double au);
};

#endif
