#include <SPEval.H>
#include <Transport_F.H>

#include <iostream>


extern "C" {
    void extern_init(const int* name, const int* namlen);
}

static bool fortran_data_initialized = false;

void
fort_data_init (const std::string& probin_file)
{
    if (!fortran_data_initialized)
    {
	// initialize the external runtime parameters -- these will
	// live in the probin

	std::cout << "reading extern runtime parameters ..." << std::endl;

	int probin_file_length = probin_file.length();
	std::vector<int> probin_file_name(probin_file_length);

	for (int i = 0; i < probin_file_length; i++)
	    probin_file_name[i] = probin_file[i];

	extern_init(&(probin_file_name[0]),&probin_file_length);

	fortran_data_initialized = true;
    }
}

void eval_single_point_transport(const std::vector<Real>& massFrac,
				 const Real&              temp,
				 const Real&              density,
				 std::vector<Real>&       mix_avg_diffusivity,
				 Real&                    shear_viscosity,
				 Real&                    bulk_viscosity,
				 Real&                    conductivity)
{
    int nspecies = get_num_species();

    mix_avg_diffusivity.resize(nspecies);

    // Build a trivial box, from ivt to ivt, to pass to transport coeff routine
    int ivt[3] = {0, 0, 0};

    get_transport_coeffs(ivt, ivt,
			 &(massFrac[0]),           ivt, ivt,
			 &temp,                    ivt, ivt,
			 &density,                 ivt, ivt,
			 &(mix_avg_diffusivity[0]),ivt, ivt,
			 &shear_viscosity,         ivt, ivt,
			 &bulk_viscosity,          ivt, ivt,
			 &conductivity,            ivt, ivt);
}
