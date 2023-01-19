#ifndef SPEVAL_H
#define SPEVAL_H

#include <vector>
#include <string>

// This normally comes in through BoxLibs REAL.H
// not sure what best to do here...
typedef double Real;

extern "C" {
    int get_num_species();

    void fort_data_init (const std::string& probin_file);

    void eval_single_point_transport(const std::vector<Real>& massFrac,
				     const Real&              temp,
				     const Real&              density,
				     std::vector<Real>&       mix_avg_diffusivity,
				     Real&                    shear_viscosity,
				     Real&                    bulk_viscosity,
				     Real&                    conductivity);
}

#endif /* SPEval_H */
