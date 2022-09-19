.. highlight:: rst

.. _SprayInputs:

Spray Flags and Inputs
----------------------

* In the ``GNUmakefile``, specify ``USE_PARTICLES = TRUE`` and ``SPRAY_FUEL_NUM = N`` where ``N`` is the number of liquid species being used in the simulation

* Depending on the gas phase solver, spray solving functionality can be turned on in the input file using ``pelec.do_spray_particles = 1`` or ``pelelm.do_spray_particles = 1``

* There are many required ``particles.`` flags in the input file. For demonstration purposes, 2 liquid species of ``CH4`` and ``NC10H22`` will be used.

  * The liquid fuel species names are specified using ``particles.fuel_species = CH4 NC10H22``. The number of fuel species listed must match the number of liquid species specified in the ``GNUmakefile``. If the evaporated mass should contribute to a different gas phase species than what is modeled in the liquid phase, use ``particles.dep_fuel_species = SP3 SP3``, although this is not required or typical. All species specified must be present in the chemistry transport and thermodynamic data.

  * The following table lists other inputs related to ``particles.``

    +-----------------------+------------------------------+-------------+--------------------------+
    |       Input           |        Description           | Per species | Required                 |
    |                       |                              |             |                          |
    |                       |                              |             | (Default value)          |
    +=======================+==============================+=============+==========================+
    | ``fuel_ref_temp``     | Liquid reference temperature |     No      |    Yes                   |
    |                       |                              |             |                          |
    +-----------------------+------------------------------+-------------+--------------------------+
    | ``fuel_crit_temp``    | Critical temperature         |     Yes     |    Yes                   |
    +-----------------------+------------------------------+-------------+--------------------------+
    | ``fuel_boil_temp``    | Boiling temperature          |     Yes     |    Yes                   |
    |                       | at 1 atm                     |             |                          |
    +-----------------------+------------------------------+-------------+--------------------------+
    | ``fuel_cp``           | Liquid :math:`C_p` at        |     Yes     |    Yes                   |
    |                       | reference temperature        |             |                          |
    +-----------------------+------------------------------+-------------+--------------------------+
    | ``fuel_latent``       | Latent heat at               |     Yes     |    Yes                   |
    |                       | reference temperature        |             |                          |
    +-----------------------+------------------------------+-------------+--------------------------+
    | ``mom_transfer``      | Couple momentum with         |     No      |    No (1)                |
    |                       | gas phase                    |             |                          |
    +-----------------------+------------------------------+-------------+--------------------------+
    | ``mass_transfer``     | Evaporate mass and           |     No      |    No (1)                |
    |                       | exchange heat                |             |                          |
    +-----------------------+------------------------------+-------------+--------------------------+
    | ``fixed_parts``       | Fix particles in space       |     No      |    No (0)                |
    +-----------------------+------------------------------+-------------+--------------------------+
    | ``parcel_size``       | Number of droplets per       |     No      |    No (1.)               |
    |                       | parcel                       |             |                          |
    +-----------------------+------------------------------+-------------+--------------------------+
    | ``write_ascii_files`` | Output ascii                 |     No      |    No (0)                |
    |                       | files of spray data          |             |                          |
    +-----------------------+------------------------------+-------------+--------------------------+
    | ``init_function``     | Initialize with              |     No      |    No (1)                |
    |                       | ``InitSprayParticles()``     |             |                          |
    +-----------------------+------------------------------+-------------+--------------------------+
    | ``init_file``         | Ascii file name to initialize|     No      |    No (0)                |
    |                       | sprays                       |             |                          |
    +-----------------------+------------------------------+-------------+--------------------------+

  * If an Antoine fit for saturation pressure is used, it must be specified for the individual species, ::

      particles.NC10H22_psat = 4.07857 1501.268 -78.67 1.E5

    where the numbers represent :math:`a`, :math:`b`, :math:`c`, and :math:`d`, respectively in:

.. math::
   p_{\rm{sat}}(T) = d 10^{a - b / (T + c)}
