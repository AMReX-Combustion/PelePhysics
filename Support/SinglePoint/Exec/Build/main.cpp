#include <iostream>
#include <SPEval.H>

int
main (int   argc,
      char* argv[])
{
    std::string probin_file="probin";
    fort_data_init(probin_file);

    int nspecies = get_num_species();

    Real temperature = 300;
    Real density = 1.e-3;

    std::vector<Real> mass_fraction(nspecies, 0);
    mass_fraction[nspecies-1] = 1;

    std::vector<Real> mix_avg_diffusivity(nspecies);
    Real shear_viscosity, bulk_viscosity, conductivity;

    eval_single_point_transport(mass_fraction,temperature,density,
				mix_avg_diffusivity, shear_viscosity,
				bulk_viscosity, conductivity);

    std::cout << "Temperature: " << temperature << "(K)" << std::endl;
    std::cout << "Mass density: " << density << "(g/cm^3)" << std::endl;
    std::cout << "Mass Fractions (nspecies = " << nspecies << ")" << std::endl;
    Real sum = 0;
    for (int n=0; n<nspecies; ++n) {
	std::cout << "  Species " << n << ": " << mass_fraction[n] << std::endl;
	sum += mass_fraction[n];
    }
    std::cout << "-------------------------------" << std::endl;
    std::cout << "        Sum: " << sum << " (should equal 1)" << std::endl;
    std::cout << "===============================" << std::endl;
    std::cout << "      conductivity : " << conductivity << std::endl;
    std::cout << "   shear viscosity : " << shear_viscosity << std::endl;
    std::cout << "    bulk viscosity : " << bulk_viscosity << std::endl;
    for (int n=0; n<nspecies; ++n) {
	std::cout << "              D(" << n << ") : " << mix_avg_diffusivity[n] << std::endl;
    }

    return 0;
}




