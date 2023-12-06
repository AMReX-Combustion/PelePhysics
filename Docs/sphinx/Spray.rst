.. highlight:: rst

.. _Spray:

*****
Spray
*****

Spray Equations
===============

This section outlines the modeling and mathematics for the spray routines.
Firstly, spray modeling relies on the following assumptions:

* Dilute spray: droplet volume inside an Eulerian cell is much smaller than the volume of the gas phase; the droplets can be modeled as Lagrangian point source terms relative to the Eulerian gas phase

* Infinite conductivity model: temperature within a droplet is temporally varying but spatially uniform

* One-third rule: the thermophysical properties in the film of an evaporating droplet can be approximated as a weighted average of the state at the droplet surface (weighted as 2/3) and the state of the surrounding gas (weighted as 1/3)

* Ideal equilibrium: the liquid and vapor state at the surface of the droplet are in equilibrium

* The radiation, Soret, and Dufour effects are neglected

The evaporation models follow the work by Abramzon and Sirignano [#abram]_ and the multicomponent evaporation is based on work by Tonini. [#ton]_ Details regarding the energy balance are provided in Ge et al. [#Ge]_

The subscript notation for this section is: :math:`d` relates to the liquid droplet, :math:`v` relates to the vapor state that is in equilibrium with the liquid and gas phase, :math:`L` relates to the liquid phase, and :math:`g` relates to the gas phase. The subscript :math:`r` relates to the reference state with which to approximate the thermophysical and transport properties. This reference state is assumed to be in the evaporating film that surrounds the droplet state and is approximated as

.. math::
   T_r &= T_d + A (T_g - T_d)

   Y_{r,n} &= Y_{v,n} + A (Y_{g,n} - Y_{v,n})

where :math:`A = 1/3` according the the one-third rule.
Additional nomenclature: :math:`M_n` is the molar mass of species :math:`n`, :math:`\overline{M}` is the average molar mass of a mixture, :math:`\mathcal{R}` is the universal gas constant, :math:`N_L` is the number of liquid species, and :math:`N_s` is the number of gas phase species. :math:`Y_n` and :math:`\chi_n` are the mass and molar fractions of species :math:`n`, respectively.
The user is required to provide a reference temperature for the liquid properties, :math:`T^*`, the critical temperature for each liquid species, :math:`T_{c,n}`, the boiling temperature for each liquid species at atmospheric pressure, :math:`T^*_{b,n}`, the latent heat and liquid specific heat at the reference temperature, :math:`h_{L,n}(T^*)` and :math:`c_{p,L,n}(T^*)`, respectively.
Note: this reference temperature is a constant value for all species and is not related to the reference state denoted by the subscript :math:`r`.

The equations of motion, mass, momentum, and energy for the Lagrangian spray droplet are:

.. math::
   \frac{d \mathbf{X}_d}{d t} &= \mathbf{u}_d,

   \frac{d m_d}{d t} &= \sum^{N_L}_{n=0} \dot{m}_n,

   m_d \frac{d Y_{d,n}}{d t} &= \dot{m}_n - Y_{d,n} \frac{d m_d}{d t},

   m_d \frac{d \mathbf{u}_d}{d t} &= \mathbf{F}_d + m_d \mathbf{g},

   m_d c_{p,L} \frac{d T_d}{d t} &= \sum^{N_L}_{n=0} \dot{m}_n h_{L,n}(T_d) + \mathcal{Q}_d.

where :math:`\mathbf{X}_d` is the spatial vector, :math:`\mathbf{u}_d` is the velocity vector, :math:`T_d` is the droplet temperature, :math:`m_d` is the mass of the droplet, :math:`\mathbf{g}` is an external body force (like gravity), :math:`\dot{m}` is evaporated mass, :math:`\mathcal{Q}_d` is the heat transfer between the droplet and the surrounding gas, and :math:`\mathbf{F}_d` is the momentum source term.
The density of the liquid mixture, :math:`\rho_d`, depends on the liquid mass fractions of the dropet, :math:`Y_{d,n}`,

.. math::
   \rho_d = \left( \sum^{N_L}_{n=0} \frac{Y_{d,n}}{\rho_{L,n}} \right)^{-1}

The droplets are assumed to be spherical with diameter :math:`d_d`. Therefore, the mass is computed as

.. math::
   m_d = \frac{\pi}{6} \rho_d d_d^3

The procedure is as follows for updating the spray droplet:

#. Interpolate the gas phase state to the droplet location using a trilinear interpolation scheme.
#. Compute the boiling temperature for species :math:`n` at the current gas phase pressure using the Clasius-Clapeyron relation

   .. math::
      T_{b,n} = \left(\log\left(\frac{p_{\rm{atm}}}{p_g}\right) \frac{\mathcal{R}}{M_n h_{L,n}(T^*_{b,n})} + \frac{1}{T^*_{b,n}}\right)

   The boiling temperature of the droplet is computed as

   .. math::
      T_{d,b} = \sum^{N_L}_{n=0} Y_{d,n} T_{b,n}

   Since we only have the latent heat at the reference condition temperature, we estimate the enthalpy at the boiling condition using Watson's law

   .. math::
      h_{L,n}(T^*_{b,n}) = h_{L,n}(T^*) \left(\frac{T_{c,n} - T^*}{T_{c,n} - T^*_{b,n}} \right)^{-0.38}

#. Compute the latent heat of the droplet using

   .. math::
      h_{L,n}(T_d) = h_{g,n}(T_d) - h_{g,n}(T^*) + h_{L,n}(T^*) - c_{p,L,n}(T^*) (T_d - T^*) \,.


   and the saturation pressure using either the Clasius-Clapeyron relation


   .. math::
      p_{{\rm{sat}}, n} = p_{\rm{atm}} \exp\left(\frac{h_{L,n}(T_d) M_n}{\mathcal{R}} \left(\frac{1}{T^*_{b,n}} - \frac{1}{T_d}\right)\right)

   or the Antoine curve fit

   .. math::
      p_{{\rm{sat}},n} = d 10^{a - b / (T_d + c)}

#. Estimate the mass fractions in the vapor state using Raoult's law

   .. math::
      Y_{v,n} &= \frac{\chi_{v,n} M_n}{\overline{M}_v + \overline{M}_g (1 - \chi_{v,{\rm{sum}}})} \; \forall n \in N_L

      \chi_{v,{\rm{sum}}} &= \sum^{N_L}_{n=0} \chi_{v,n}

      \chi_{v,n} &= \frac{\chi_{d,n} p_{{\rm{sat}},n}}{p_g}

      \chi_{d,n} &= \frac{Y_{d,n}}{M_n}\left(\sum^{N_L}_{k=0} \frac{Y_{d,k}}{M_k}\right)^{-1}

      \overline{M}_v &= \sum^{N_L}_{n=0} \chi_{v,n} M_n

   If :math:`\chi_{g,n} p_g > p_{{\rm{sat}},n}`, then :math:`\chi_{v,n} = Y_{v,n} = 0` for that particular species in the equations above, since that means the gas phase is saturated. The mass fractions in the reference state for the fuel are computed using the one-third rule and the remaining reference mass fractions are normalized gas phase mass fractions to ensure they sum to 1

   .. math::
      Y_{r,n} = \left\{\begin{array}{c l}
      \displaystyle Y_{v,n} + A (Y_{g,n} - Y_{v,n}) & {\text{If $Y_{v,n} > 0$}}, \\
      \displaystyle\frac{1 - \sum^{N_L}_{k=0} Y_{v,k}}{1 - \sum^{N_L}_{k=0} Y_{g,k}} Y_{g,n} & {\text{Otherwise}}.
      \end{array}\right. \; \forall n \in N_s.

#. The average molar mass, specific heat, and density of the reference state in the gas film are computed as

   .. math::
      \overline{M}_r &= \left(\sum^{N_s}_{n=0} \frac{Y_{r,n}}{M_n}\right)^{-1},

      c_{p,r} &= \sum^{N_s}_{n=0} Y_{s,n} c_{p,g,n}(T_r),

      \rho_r &= \frac{\overline{M}_r p_g}{\mathcal{R} T_r}.

#. Transport properties are computed using the reference state: dynamic viscosity, :math:`\mu_r`, thermal conductivity, :math:`\lambda_r`, and mass diffusion coefficient for species :math:`n`, :math:`D_{r,n}`.

#. It is important to note that `PelePhysics` provides mixture averaged mass diffusion coefficient :math:`\overline{(\rho D)}_{r,n}`, which is converted into the binary mass diffusion coefficient using

   .. math::
      (\rho D)_{r,n} = \overline{(\rho D)}_{r,n} \overline{M}_r / M_n.

   Mass diffusion coefficient is then normalized by the total fuel vapor molar fraction

   .. math::
      (\rho D)^*_{r,n} = \frac{\chi_{v,n} (\rho D)_{r,n}}{\chi_{v,{\rm{sum}}}} \; \forall n \in N_L

   and the total is

   .. math::
      (\rho D)_r = \sum_{n=0}^{N_L} (\rho D)_{r,n}^*

#. The momentum source is a function of the drag force

   .. math::
      \mathbf{F}_d = \frac{1}{2} \rho_r C_D A_d \left\|\Delta \mathbf{u}\right\| \Delta \mathbf{u}

   where :math:`\Delta \mathbf{u} = \mathbf{u}_g - \mathbf{u}_d`, :math:`A_d = \pi d_d^2/4` is the frontal area of the droplet, and :math:`C_D` is the drag coefficient for a sphere, which is estimated using the standard drag curve for an immersed sphere

   .. math::
      C_D = \frac{24}{{\rm{Re}}_d}\left\{\begin{array}{c l}
      1 & {\text{If Re$_d$ < 1}}, \\
      \displaystyle 1 + \frac{{\rm{Re}}^{2/3}_d}{6} & {\text{Otherwise}}.
      \end{array}\right.

   The droplet Reynolds number is defined as

   .. math::
      {\rm{Re}}_d = \frac{\rho_r d_d \left\|\Delta \mathbf{u}\right\|}{\mu_r}


#. The mass source term is modeled according to Abramzon and Sirignano (1989). The following non-dimensional numbers and factors are used:

   .. math::
      F(B) &= (1 + B)^{0.7}\frac{\log(1 + B)}{B}

      F_2 &= \max(1, \min(400, {\rm{Re}}_d)^{0.077})

      {\rm{Pr}}_r &= \frac{\mu_r c_{p,r}}{\lambda_r}

      {\rm{Sc}}_r &= \frac{\mu_r}{(\rho D)_r}

      {\rm{Sh}}_0 &= 1 + (1 + {\rm{Re}}_d {\rm{Sc}}_r)^{1/3} F_2

      {\rm{Nu}}_0 &= 1 + (1 + {\rm{Re}}_d {\rm{Pr}}_r)^{1/3} F_2

      {\rm{Sh}}^* &= 2 + \frac{{\rm{Sh}}_0 - 2}{F(B_M)}

      {\rm{Nu}}^* &= 2 + \frac{{\rm{Nu}}_0 - 2}{F(B_T)}

   * The Spalding numbers for mass transfer, :math:`B_M`, and heat transfer, :math:`B_T`, are computed using

     .. math::
        B_M &= \displaystyle\frac{\sum^{N_L}_{n=0} Y_{v,n} - \sum^{N_L}_{n=0} Y_{g,n}}{1 - \sum^{N_L}_{n=0} Y_{v,n}}

        B_T &= \left(1 + B_M\right)^{\phi} - 1

     where

     .. math::
        \phi = \frac{c_{p,r} (\rho D)_r {\rm{Sh}}^*}{\lambda_r {\rm{Nu}}^*}

     Note the dependence of :math:`{\rm{Nu}}^*` on :math:`B_T` means an iterative scheme is required to solve for both. The droplet vaporization rate and heat transfer become

     .. math::
        \dot{m}_n &= -\pi (\rho D)_{r,n}^* d_d {\rm{Sh}}^* \log(1 + B_M). \; \forall n \in N_L

        \mathcal{Q}_d &= \pi \lambda_r d_d (T_g - T_d) {\rm{Nu}}^* \frac{\log(1 + B_T)}{B_T}

   * If the gas phase is saturated for all liquid species, the equations for heat and mass transfer become

     .. math::
        \dot{m}_n &= 0

        \mathcal{Q}_d &= \pi \lambda_r d_d (T_g - T_d) {\rm{Nu}}_0

#. To alleviate conservation issues at AMR interfaces, each parcel only contributes to the gas phase source term of the cell containing it. The gas phase source terms for a single parcel to the cell are

    .. math::
       S_{\rho} &= \mathcal{C} \sum^{N_L}_{n=0} \dot{m}_n,

       S_{\rho Y_n} &= \mathcal{C} \dot{m}_n,

       \mathbf{S}_{\rho \mathbf{u}} &= \mathcal{C} \mathbf{F}_d,

       S_{\rho h} &= \mathcal{C}\left(\mathcal{Q}_d + \sum_{n=0}^{N_L} \dot{m}_n h_{g,n}(T_d)\right),

       S_{\rho E} &= S_{\rho h} + \frac{1}{2}\left\|\mathbf{u}_d\right\| S_{\rho} + \mathcal{C} \mathbf{F}_d \cdot \mathbf{u}_d

    where

    .. math::
       \mathcal{C} = -\frac{N_{d}}{V_{\rm{cell}}},

    :math:`N_{d}` is the number of droplets per computational parcel, and :math:`V_{\rm{cell}}` is the volume for the cell of interest. Note that the cell volume can vary depending on AMR level and if an EB is present.

Spray Flags and Inputs
======================

* In the ``GNUmakefile``, specify ``USE_PARTICLES = TRUE`` and ``SPRAY_FUEL_NUM = N`` where ``N`` is the number of liquid species being used in the simulation.

* Depending on the gas phase solver, spray solving functionality can be turned on in the input file using ``pelec.do_spray_particles = 1`` or ``peleLM.do_spray_particles = 1``.

* The units for `PeleLM` and `PeleLMeX` are MKS while the units for `PeleC` are CGS. This is the same for the spray inputs. E.g. when running a spray simulation coupled with `PeleC`, the units for ``particles.fuel_cp`` must be in erg/g.

* There are many required ``particles.`` flags in the input file. For demonstration purposes, 2 liquid species of ``NC7H16`` and ``NC10H22`` will be used.

  * The liquid fuel species names are specified using ``particles.fuel_species = NC7H16 NC10H22``. The number of fuel species listed must match ``SPRAY_FUEL_NUM``.

  * Many values must be specified on a per-species basis. Following the current example, one would have to specify ``particles.NC7H16_crit_temp = 540.`` and ``particles.NC10H22_crit_temp = 617.`` to set a critical temperature of 540 K for ``NC7H16`` and 617 K for ``NC10H22``.

  * Although this is not required or typical, if the evaporated mass should contribute to a different gas phase species than what is modeled in the liquid phase, use ``particles.dep_fuel_species``. For example, if we wanted the evaporated mass from both liquid species to contribute to a different species called ``SP3``, we would put ``particles.dep_fuel_species = SP3 SP3``. All species specified must be present in the chemistry transport and thermodynamic data.

* The following table lists other inputs related to ``particles.``, where ``SP`` will refer to a fuel species name

.. table::

   +-----------------------+-------------------------------+-------------+-------------------+
   |Input                  |Description                    |Required     |Default Value      |
   +=======================+===============================+=============+===================+
   |``fuel_species``       |Names of liquid species        |Yes          |None               |
   +-----------------------+-------------------------------+-------------+-------------------+
   |``dep_fuel_species``   |Name of gas phase species to   |Yes          |Inputs to          |
   |                       |contribute                     |             |``fuel_species``   |
   +-----------------------+-------------------------------+-------------+-------------------+
   |``fuel_ref_temp``      |Liquid reference temperature   |Yes          |None               |
   +-----------------------+-------------------------------+-------------+-------------------+
   |``SP_crit_temp``       |Critical temperature           |Yes          |None               |
   +-----------------------+-------------------------------+-------------+-------------------+
   |``SP_boil_temp``       |Boiling temperature at         |Yes          |None               |
   |                       |atmospheric pressure           |             |                   |
   +-----------------------+-------------------------------+-------------+-------------------+
   |``SP_cp``              |Liquid :math:`c_p` at reference|Yes          |None               |
   |                       |temperature                    |             |                   |
   +-----------------------+-------------------------------+-------------+-------------------+
   |``SP_latent``          |Latent heat at reference       |Yes          |None               |
   |                       |temperature                    |             |                   |
   +-----------------------+-------------------------------+-------------+-------------------+
   |``SP_rho``             |Liquid density                 |Yes          |None               |
   |                       |                               |             |                   |
   +-----------------------+-------------------------------+-------------+-------------------+
   |``SP_lambda``          |Liquid thermal conductivity    |No           |0.                 |
   |                       |(currently unused)             |             |                   |
   +-----------------------+-------------------------------+-------------+-------------------+
   |``SP_mu``              |Liquid dynamic viscosity       |No           |0.                 |
   |                       |(currently unused)             |             |                   |
   +-----------------------+-------------------------------+-------------+-------------------+
   |``mom_transfer``       |Couple momentum with gas phase |No           |``1``              |
   |                       |                               |             |                   |
   +-----------------------+-------------------------------+-------------+-------------------+
   |``mass_transfer``      |Evaporate mass and exchange    |No           |``1``              |
   |                       |heat with gas phase            |             |                   |
   +-----------------------+-------------------------------+-------------+-------------------+
   |``fixed_parts``        |Fix particles in space         |No           |``0``              |
   +-----------------------+-------------------------------+-------------+-------------------+
   |``parcel_size``        |:math:`N_{d}`; Number of       |No           |``1.``             |
   |                       |droplets per parcel            |             |                   |
   +-----------------------+-------------------------------+-------------+-------------------+
   |``write_ascii_files``  |Output ascii files of spray    |No           |``0``              |
   |                       |data                           |             |                   |
   +-----------------------+-------------------------------+-------------+-------------------+
   |``cfl``                |Particle CFL number for        |No           |``0.5``            |
   |                       |limiting time step             |             |                   |
   +-----------------------+-------------------------------+-------------+-------------------+
   |``init_file``          |Ascii file name to initialize  |No           |Empty              |
   |                       |droplets                       |             |                   |
   +-----------------------+-------------------------------+-------------+-------------------+


* If an Antoine fit for saturation pressure is used, it must be specified for individual species, ::

    particles.SP_psat = 4.07857 1501.268 -78.67 1.E5

  where the numbers represent :math:`a`, :math:`b`, :math:`c`, and :math:`d`, respectively in:

  .. math::
     p_{\rm{sat}}(T) = d 10^{a - b / (T + c)}

  * If no fit is provided, the saturation pressure is estimated using the Clasius-Clapeyron relation; see

* Temperature based fits for liquid density, thermal conductivity, and dynamic viscosity can be used; these can be specified as ::

    particles.SP_rho = 10.42 -5.222 1.152E-2 4.123E-7
    particles.SP_lambda = 7.243 1.223 4.223E-8 8.224E-9
    particles.SP_mu = 7.243 1.223 4.223E-8 8.224E-9

  where the numbers respresent :math:`a`, :math:`b`, :math:`c`, and :math:`d`, respectively in:

  .. math::
     \rho_L \,, \lambda_L = a + b T + c T^2 + d T^3

     \mu_L = a + b / T + c / T^2 + d / T^3

  If only a single value is provided, :math:`a` is assigned to that value and the other coefficients are set to zero, effectively using a constant value for the parameters.

Spray Injection
---------------

Templates to facilitate and simplify spray injection are available. To use them, changes must be made to the input and ``SprayParticlesInitInsert.cpp`` files. Inputs related to injection use the ``spray.`` parser name. To create a jet in the domain, modify the ``InitSprayParticles()`` function in ``SprayParticleInitInsert.cpp``. Here is an example: ::

  void
  SprayParticleContainer::InitSprayParticles(
  const bool init_parts, ProbParm const& prob_parm)
  {
    amrex::ignore_unused(prob_parm);
    int num_jets = 1;
    m_sprayJets.resize(num_jets);
    std::string jet_name = "jet1";
    m_sprayJets[0] = std::make_unique<SprayJet>(jet_name, Geom(0));
    return;
  }


This creates a single jet that is named ``jet1``. This name will be used in the input file to reference this particular jet. For example, to set the location of the jet center for ``jet1``, the following should be included in the input file, ::

  spray.jet1.jet_cent = 0. 0. 0.

No two jets may have the same name. If an injector is constructed using only a name and geometry, the injection parameters are read from the input file. Here is a list of injection related inputs:

.. table::
   :widths: 20 40 20

   +--------------------+--------------------------------+--------------------+
   |Input               |Description                     |Required            |
   |                    |                                |                    |
   +====================+================================+====================+
   |``jet_cent``        |Jet center location             |Yes                 |
   +--------------------+--------------------------------+--------------------+
   |``jet_norm``        |Jet normal direction            |Yes                 |
   +--------------------+--------------------------------+--------------------+
   |``jet_vel``         |Jet velocity magnitude          |Yes                 |
   +--------------------+--------------------------------+--------------------+
   |``jet_dia``         |Jet diameter                    |Yes                 |
   +--------------------+--------------------------------+--------------------+
   |``spread_angle``    |:math:`\theta_J`; Full spread   |Yes                 |
   |                    |angle in degrees from the jet   |                    |
   |                    |normal direction; droplets vary |                    |
   |                    |from                            |                    |
   |                    |:math:`[-\theta_J/2,\theta_J/2]`|                    |
   +--------------------+--------------------------------+--------------------+
   |``T``               |Temperature of the injected     |Yes                 |
   |                    |liquid                          |                    |
   +--------------------+--------------------------------+--------------------+
   |``Y``               |Mass fractions of the injected  |Yes, if             |
   |                    |liquid based on                 |``SPRAY_FUEL_NUM`` >|
   |                    |``particles.fuel_species``      |1                   |
   +--------------------+--------------------------------+--------------------+
   |``mass_flow_rate``  |:math:`\dot{m}_{\rm{inj}}`; Mass|Yes                 |
   |                    |flow rate of the jet            |                    |
   +--------------------+--------------------------------+--------------------+
   |``hollow_spray``    |Sets hollow cone injection with |No (Default: 0)     |
   |                    |angle :math:`\theta_J/2`        |                    |
   +--------------------+--------------------------------+--------------------+
   |``hollow_spread``   |:math:`\theta_h`; Adds spread to|No (Default: 0)     |
   |                    |hollow cone :math:`\theta_J/2\pm|                    |
   |                    |\theta_h`                       |                    |
   +--------------------+--------------------------------+--------------------+
   |``swirl_angle``     |:math:`\phi_S`; Adds a swirling |No (Default: 0)     |
   |                    |component along azimuthal       |                    |
   |                    |direction                       |                    |
   +--------------------+--------------------------------+--------------------+
   |``start_time`` and  |Beginning and end time for jet  |No                  |
   |``end_time``        |                                |                    |
   +--------------------+--------------------------------+--------------------+
   |``dist_type``       |Droplet diameter distribution   |Yes                 |
   |                    |type: ``Uniform``, ``Normal``,  |                    |
   |                    |``LogNormal``, ``Weibull``,     |                    |
   |                    |``ChiSquared``                  |                    |
   +--------------------+--------------------------------+--------------------+


.. figure:: /Visualization/inject_transform.png
   :align: center
   :figwidth: 60%

   Demonstration of injection angles. :math:`\phi_J` varies uniformly from :math:`[0, 2 \pi]`


Care must be taken to ensure the amount of mass injected during a time step matches the desired mass flow rate. For smaller time steps, the risk of over-injecting mass increases. To mitigate this issue, each jet accounts for three values: :math:`N_{P,\min}`, :math:`m_{\rm{acc}}`, and :math:`t_{\rm{acc}}` (labeled in the code as ``m_minParcel``, ``m_sumInjMass``, and ``m_sumInjTime``, respectively). :math:`N_{P,\min}` is the minimum number of parcels that must be injected over the course of an injection event; this must be greater than or equal to one. :math:`m_{\rm{acc}}` is the amount of uninjected mass accumulated over the time period :math:`t_{\rm{acc}}`. The injection routine steps are as follows:

#. The injected mass for the current time step is computed using the desired mass flow rate, :math:`\dot{m}_{\rm{inj}}` and the current time step

   .. math::
      m_{\rm{inj}} = \dot{m}_{\rm{inj}} \Delta t + m_{\rm{acc}}

#. The time period for the current injection event is computed using

   .. math::
      t_{\rm{inj}} = \Delta t + t_{\rm{acc}}

#. Using the average mass of an injected parcel, :math:`N_{d} m_{d,\rm{avg}}`, the estimated number of injected parcels is computed

   .. math::
      N_{P, \rm{inj}} = m_{\rm{inj}} / (N_{d} m_{d, \rm{avg}})

  * If :math:`N_{P, \rm{inj}} < N_{P, \min}`, the mass and time is accumulated as :math:`m_{\rm{acc}} = m_{\rm{inj}}` and :math:`t_{\rm{acc}} = t_{\rm{inj}}` and no injection occurs this time step.

  * Otherwise, :math:`m_{\rm{inj}}` mass is injected and convected over time :math:`t_{\rm{inj}}` and :math:`m_{\rm{acc}}` and :math:`t_{\rm{acc}}` are reset.

4. If injection occurs, the amount of mass injected, :math:`m_{\rm{actual}}`, is summed and compared with the desired mass flow rate. If :math:`m_{\rm{actual}} / t_{\rm{inj}} - \dot{m}_{\rm{inj}} > 0.05 \dot{m}_{\rm{inj}}`, then :math:`N_{P,\min}` is increased by one to reduce the liklihood of over-injecting in the future. A balance is necessary: the higher the minimum number of parcels, the less likely to over-inject mass but the number of time steps between injections can potentially grow as well.

Spray Validation
================

Single Droplet Tests
--------------------

Single droplet tests are performed and compared with computational or experimental results published in literature. These tests are setup in ``PeleProduction/PeleMPruns/single_drop_test``. To run a test case, simply open ``Validate.py`` and set the case name from the table below ::

  case = TestCaseName()

then do ``python Validate.py``.
The following table details the parameters of each test:

.. table::

   +---------------+-----------------+-----------------+-----------------+-----------------+-----------------------+-----------------+
   |Test Case Name | :math:`T_g` [K] |:math:`p_g` [bar]|:math:`T_d` [K]  |:math:`d_d` [um] | :math:`\Delta u` [m/s]|Ref              |
   |               |                 |                 |                 |                 |                       |                 |
   +===============+=================+=================+=================+=================+=======================+=================+
   |``Tonini_4_33``|1000             |1                |300              |200              |6.786                  |[#ton]_          |
   +---------------+-----------------+-----------------+-----------------+-----------------+-----------------------+-----------------+
   |``Abramzon``   |1500             |10               |300              |100              |15                     |[#abram]_        |
   +---------------+-----------------+-----------------+-----------------+-----------------+-----------------------+-----------------+
   |``Daif``       |348              |1                |294              |1334             |3.10                   |[#daif]_         |
   +---------------+-----------------+-----------------+-----------------+-----------------+-----------------------+-----------------+
   |``RungeHep``   |273              |1                |272              |500-570          |2.5                    |[#runge]_        |
   |``RungeDec``   |                 |                 |                 |                 |                       |                 |
   |``RungeMix``   |                 |                 |                 |                 |                       |                 |
   +---------------+-----------------+-----------------+-----------------+-----------------+-----------------------+-----------------+

.. figure:: /Visualization/ton_res.png
   :align: center
   :figwidth: 80%

   Droplet diameter, temperature, and n-octane mass fraction comparisons with Figure 4.33 in [#ton]_

.. figure:: /Visualization/abram_res.png
   :align: center
   :figwidth: 80%

   Droplet diameter and temperature comparisons with [#abram]_

.. figure:: /Visualization/daif_res.png
   :align: center
   :figwidth: 80%

   Droplet diameter and temperature comparisons with [#daif]_

.. [#ton] "Fuel spray modeling in direct-injection diesel and gasoline engines", S. Tonini, Dissertation, City University London (2006)

.. [#abram] "Droplet vaporization model for spray combustion calculations", B. Abramzon and W. A. Sirignano, Int. J. Heat Mass Transfer, Vol. 32, No. 9, pp. 1605-1618 (1989)

.. [#daif] "Comparison of multicomponent fuel droplet vaporization experiments in forced convection with the Sirignano model", A. Daı̈f and M. Bouaziz and X. Chesneau and A. Ali Chérif, Exp. Therm. Fluid Sci., Vol. 18, No. 4, pp. 282-290, Issn 0894-1777 (1998)

.. [#runge] "Low-temperature vaporization of JP-4 and JP-8 fuel droplets", T. Runge and M. Teske and C. E. Polymeropoulos, At. Sprays, Vol. 8, pp. 25-44 (1998)

.. [#Ge] "Development of a CPU/GPU portable software library for Lagrangian-Eulerian simulations of liquid sprays", W. Ge and R. Sankaran and J. H. Chen, Int. J. Multiph. Flow, Vol. 128 (2020)
