.. highlight:: rst

.. _SprayInputs:

Spray Flags and Inputs
----------------------

* In the ``GNUmakefile``, specify ``USE_PARTICLES = TRUE`` and ``SPRAY_FUEL_NUM = N`` where ``N`` is the number of liquid species being used in the simulation

* Depending on the gas phase solver, spray solving functionality can be turned on in the input file using ``pelec.do_spray_particles = 1`` or ``pelelm.do_spray_particles = 1``

* There are many required ``particles.`` flags in the input file. For demonstration purposes, 2 liquid species of ``CH4`` and ``NC10H22`` will be used.

  * The liquid fuel species names are specified using ``particles.fuel_species = CH4 NC10H22``. The number of fuel species listed must match the number of liquid species specified in the ``GNUmakefile``. If the evaporated mass should contribute to a different gas phase species than what is modeled in the liquid phase, use ``particles.dep_fuel_species = SP3 SP3``, although this is not required or typical. All species specified must be present in the chemistry transport and thermodynamic data.

* The following table lists other inputs related to ``particles.``

.. table::
   :widths: 40 40 40 40

   +-----------------------+-------------------------------+-------------+--------------------------+
   |       Input           |        Description            | Per species |         Required         |
   |                       |                               |             |                          |
   |                       |                               |             |     (Default value)      |
   +=======================+===============================+=============+==========================+
   | ``fuel_ref_temp``     | Liquid reference temperature  |     No      |    Yes                   |
   +-----------------------+-------------------------------+-------------+--------------------------+
   | ``fuel_crit_temp``    | Critical temperature          |     Yes     |    Yes                   |
   +-----------------------+-------------------------------+-------------+--------------------------+
   | ``fuel_boil_temp``    | Boiling temperature at        |     Yes     |    Yes                   |
   |                       | atmospheric pressure          |             |                          |
   +-----------------------+-------------------------------+-------------+--------------------------+
   | ``fuel_cp``           | Liquid :math:`c_p` at         |     Yes     |    Yes                   |
   |                       | reference temperature         |             |                          |
   +-----------------------+-------------------------------+-------------+--------------------------+
   | ``fuel_latent``       | Latent heat at reference      |     Yes     |    Yes                   |
   |                       | temperature                   |             |                          |
   +-----------------------+-------------------------------+-------------+--------------------------+
   | ``mom_transfer``      | Couple momentum with gas      |     No      |    No (1)                |
   |                       | phase                         |             |                          |
   +-----------------------+-------------------------------+-------------+--------------------------+
   | ``mass_transfer``     | Evaporate mass and exchange   |     No      |    No (1)                |
   |                       | heat                          |             |                          |
   +-----------------------+-------------------------------+-------------+--------------------------+
   | ``fixed_parts``       | Fix particles in space        |     No      |    No (0)                |
   +-----------------------+-------------------------------+-------------+--------------------------+
   | ``parcel_size``       | Number of droplets per parcel |     No      |    No (1.)               |
   +-----------------------+-------------------------------+-------------+--------------------------+
   | ``write_ascii_files`` | Output ascii files of spray   |     No      |    No (0)                |
   |                       | data                          |             |                          |
   +-----------------------+-------------------------------+-------------+--------------------------+
   | ``init_function``     | Initialize with               |     No      |    No (1)                |
   |                       | ``InitSprayParticles()``      |             |                          |
   |                       |                               |             |                          |
   +-----------------------+-------------------------------+-------------+--------------------------+
   | ``init_file``         | Ascii file name to initialize |     No      |    No (0)                |
   |                       | sprays                        |             |                          |
   +-----------------------+-------------------------------+-------------+--------------------------+


* If an Antoine fit for saturation pressure is used, it must be specified for the individual species, ::

    particles.NC10H22_psat = 4.07857 1501.268 -78.67 1.E5

  where the numbers represent :math:`a`, :math:`b`, :math:`c`, and :math:`d`, respectively in:

  .. math::
     p_{\rm{sat}}(T) = d 10^{a - b / (T + c)}


* Templates to facilitate and simplify spray injection are available in `PeleMP`. To use them, changes must be made to the input and ``SprayParticlesInitInsert.cpp`` file. Inputs related to injection use the ``spray.`` parser name. To create a jet in the domain, modify the ``InitSprayParticles`` function in ``SprayParticleInitInsert.cpp``. Here is an example: ::

    void
    SprayParticleContainer::InitSprayParticles(
    const bool init_parts, ProbParm const& prob_parm)
    {
      amrex::ignore_unused(prob_parm);
      int num_jets = 1;
      m_sprayJets.resize(num_jets);
      std::string jet_name = "jet1";
      m_sprayJets[0] = std::make_unique<SprayJet>(jet_name, Geom(0));
      // Start without any particles
      m_injectVel = m_sprayJets[0]->jet_vel();
      return;
    }


  This creates a single jet that is named ``jet1``. This name will be used in the input file to reference this particular jet. For example, to set the location of the jet center for ``jet1``, the following should be included in the input file, ::

    spray.jet1.jet_cent = 0. 0. 0.

  If an injector is constructed using only a name and geometry, the injection parameters are read from the input file. Here is a list of injection related inputs:

.. table::
   :widths: 40 40 40

   +--------------------+--------------------------------+--------------------+
   |  Input             |  Description                   |   Required         |
   |                    |                                |                    |
   +====================+================================+====================+
   |  ``jet_cent``      | Jet center location            | Yes                |
   |                    |                                |                    |
   +--------------------+--------------------------------+--------------------+
   | ``jet_norm``       | Jet normal direction           |  Yes               |
   |                    |                                |                    |
   +--------------------+--------------------------------+--------------------+
   | ``jet_dia``        | Jet diameter                   |  Yes               |
   |                    |                                |                    |
   +--------------------+--------------------------------+--------------------+
   |    ``spread_angle``| Angle in degrees that droplets | Yes                |
   |                    | can vary from the normal       |                    |
   |                    | direction. This is the full    |                    |
   |                    | spread angle, meaning the      |                    |
   |                    | droplet will vary from         |                    |
   |                    | :math:`-\theta/2` to           |                    |
   |                    | :math:`\theta/2`.              |                    |
   +--------------------+--------------------------------+--------------------+
   |  ``T``             | Temperature of the injected    | Yes                |
   |                    | liquid                         |                    |
   +--------------------+--------------------------------+--------------------+
   | ``Y``              | Mass fractions of the injected | Yes, if            |
   |                    | liquid. Ordered based on       | ``SPRAY_FUEL_NUM`` |
   |                    | ``particles.fuel_species``     | > 1                |
   |                    |                                |                    |
   +--------------------+--------------------------------+--------------------+
   |   ``dist_type``    | Droplet diameter distribution  |  Yes               |
   |                    | type. Options are ``Uniform``, |                    |
   |                    | `` Normal``, ``LogNormal``,    |                    |
   |                    | ``Weibull``. Each distribution |                    |
   |                    | type has it's own required     |                    |
   |                    | inputs.                        |                    |
   |                    |                                |                    |
   +--------------------+--------------------------------+--------------------+




