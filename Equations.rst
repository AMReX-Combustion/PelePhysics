.. highlight:: rst

.. _Equations:

Equations
=========

Spray Equations
---------------
This section outlines the modeling and mathematics for the spray routines in `PeleMP`.
Firstly, spray modeling in `PeleMP` relies on the following assumptions:

* Dilute spray: droplet volume inside an Eulerian cell is much smaller than the volume of the gas phase; the droplets can be modeled as Lagrangian point source terms relative to the Eulerian gas phase

* Infinite conductivity model: temperature is temporally varying but spatially uniform within a droplet

* 1/3 rule: the thermophysical properties in the film of an evaporating droplet can be approximated as an 2/3 vapor state and 1/3 gas phase state

* Ideal equilibrium: the liquid and vapor state at the surface of the droplet are in equilibrium

* The radiation, Soret, and Dufour effects are neglected

The evaporation models follow the work by Abramzon and Sirignano [#abram]_ and the multicomponent evaporation is based on work by Tonini. [#ton]_ Details regarding the energy balance are provided in Ge et al. [#Ge]_

The subscript notation for this section is: :math:`d` relates to the liquid droplet, :math:`v` relates to the vapor state that is in equilibrium with the liquid and gas phase, :math:`s` relates to the reference surface state, :math:`L` relates to the liquid phase, and :math:`g` relates to the gas phase.

The reference surface state is the state used to approximate the thermophysical and transport properties.
This state is approximated using the one-third rule.
Additional nomenclature: :math:`M_n` is the molar mass of species :math:`n`, :math:`\overline{M}` is the average molar mass of a mixture, :math:`\mathcal{R}` is the universal gas constant, :math:`N_L` is the number of liquid species, and :math:`N_s` is the number of gas phase species.
The user is required to provide a reference temperature for the liquid properties, :math:`T^*`, the critical temperature for each liquid species, :math:`T_{c,n}`, the boiling temperature for each liquid species at atmospheric pressure, :math:`T^*_{b,n}`, the latent heat and liquid specific heat at the reference temperature, :math:`h_{L,n}(T^*)` and :math:`c_{p,L,n}(T^*)`, respectively.
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

   Since we only have the latent heat at the reference temperature, we estimate the enthalpy at the boiling condition using Watson's law

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
      Y_{v,n} &= \frac{Y_{d,n} p_{{\rm{sat}}, n}}{\overline{M}_g(\phi_3 p_g - \phi_1) + \phi_2}

      \phi_1 &= \sum^{N_L}_{n=0} \frac{Y_{d,n} p_{{\rm{sat}},n}}{M_n}

      \phi_2 &= \sum^{N_L}_{n=0} Y_{d,n} p_{{\rm{sat}},n}

      \phi_3 &= \sum^{N_L}_{n=0} \frac{Y_{d,n}}{M_n}

   The reference surface mass fractions for the fuel are computed using the one-third rule and the remaining reference surface mass fractions are normalized gas phase mass fractions to ensure they sum to 1

   .. math::
      Y_{s,n} = \left\{\begin{array}{c l}
      \displaystyle\frac{2 Y_{v,n} + Y_{g,n}}{3} & {\text{If $Y_{v,n} > 0$}}, \\
      \displaystyle\frac{1 - \sum^{N_L}_{k=0} Y_{v,k}}{1 - \sum^{N_L}_{k=0} Y_{g,k}} Y_{g,n} & {\text{Otherwise}}.
      \end{array}\right. \; \forall n \in N_s.

#. The average molar mass, specific heat, and density of the reference surface state are computed as

   .. math::
      \overline{M}_s &= \left(\sum^{N_s}_{n=0} \frac{Y_{s,n}}{M_n}\right)^{-1},

      c_{p,s} &= \sum^{N_s}_{n=0} Y_{s,n} c_{p,g,n}(T_s),

      \rho_s &= \frac{\overline{M}_s p_g}{\mathcal{R} T_s}.

   where :math:`T_s = (2 T_d + T_g)/3`.

#. Transport properties are computed using the reference surface state: dynamic viscosity, :math:`\mu_s`, thermal conductivity, :math:`\lambda_s`, and mass diffusion coefficient for species :math:`n`, :math:`D_{s,n}`.

#. It is important to note that `PelePhysics` provides mixture averaged coefficient :math:`\rho_s \overline{D}_{s,n}`, which is converted into the binary coefficient with :math:`\rho_s D_{s,n} = \rho_s \overline{D}_{s,n} \overline{M}_s / M_n`. Additionally, the mass diffusion coefficient is normalized by the total fuel vapor molar fraction

   .. math::
      \rho_s D_{s,n}^* = \rho_s D_{s,n} \frac{Y_{d,n} p_{{\rm{sat}},n}/M_n}{\sum^{N_L}_{k=0}Y_{d,k} p_{{\rm{sat}},k} / M_k}

   and the total is

   .. math::
      \rho_s D_s = \sum_{n=0}^{N_L} \rho_s D_{s,n}^*

#. The momentum source is a function of the drag force

   .. math::
      \mathbf{F}_d = \frac{1}{2} \rho_s C_D A_d \left\|\Delta \mathbf{u}\right\| \Delta \mathbf{u}

   where :math:`\Delta \mathbf{u} = \mathbf{u}_g - \mathbf{u}_d`, :math:`A_d = \pi/4 d_d^2` is the frontal area of the droplet, and :math:`C_D` is the drag coefficient for a sphere, which is estimated using the standard drag curve for an immersed sphere

   .. math::
      C_D = \frac{24}{{\rm{Re}}_d}\left\{\begin{array}{c l}
      1 & {\text{If Re$_d$ < 1}}, \\
      \displaystyle 1 + \frac{{\rm{Re}}^{2/3}_d}{6} & {\text{Otherwise}}.
      \end{array}\right.

   The droplet Reynolds number is defined as

   .. math::
      {\rm{Re}}_d = \frac{\rho_s d_d \left\|\Delta \mathbf{u}\right\|}{\mu_s}


#. The mass source term is modeled according to Abramzon and Sirignano (1989). The following non-dimensional numbers and factors are used:

   .. math::
      F(B) &= (1 + B)^{0.7}\frac{\log(1 + B)}{B}

      F_2 &= \max(1, \min(400, {\rm{Re}}_d)^{0.077})

      {\rm{Pr}}_s &= \frac{\mu_s c_{p,s}}{\lambda_s}

      {\rm{Sc}}_s &= \frac{\mu_s}{\rho_s D_s}

      {\rm{Sh}}_0 &= 1 + (1 + {\rm{Re}}_d {\rm{Sc}}_s)^{1/3} F_2

      {\rm{Nu}}_0 &= 1 + (1 + {\rm{Re}}_d {\rm{Pr}}_s)^{1/3} F_2

      {\rm{Sh}}^* &= 2 + \frac{{\rm{Sh}}_0 - 2}{F(B_M)}

      {\rm{Nu}}^* &= 2 + \frac{{\rm{Nu}}_0 - 2}{F(B_T)}

   * The Spalding numbers for mass transfer, :math:`B_M`, and heat transfer, :math:`B_T`, are computed using

     .. math::
        B_M &= \displaystyle\frac{\sum^{N_L}_{n=0} Y_{v,n} - \sum^{N_L}_{n=0} Y_{g,n}}{1 - \sum^{N_L}_{n=0} Y_{v,n}}

        B_T &= \left(1 + B_M\right)^{\phi} - 1

     where

     .. math::
        \phi = \frac{c_{p,s} \rho_s D_s {\rm{Sh}}^*}{\lambda_s {\rm{Nu}}^*}

     Note the dependence of :math:`{\rm{Nu}}^*` on :math:`B_T` means an iterative scheme is required to solve for both. The droplet vaporization rate and heat transfer become

     .. math::
        \dot{m}_n &= -\pi \rho_s D_{s,n}^* d_d {\rm{Sh}}^* \log(1 + B_M). \; \forall n \in N_L

        \mathcal{Q}_d &= \pi \lambda_s d_d (T_g - T_d) {\rm{Nu}}^* \frac{\log(1 + B_T)}{B_T}

   * If :math:`X_{g,n} p_g > p_{{\rm{sat}},n}`, then the gas phase is saturated for species :math:`n`. In this case, :math:`Y_{v,n} = 0` and :math:`\dot{m}_n = 0`. If the gas phase is saturated for all liquid species, the equations for heat and mass transfer become

     .. math::
        \dot{m}_n &= 0

        \mathcal{Q}_d &= \pi \lambda_s d_d (T_g - T_d) {\rm{Nu}}_0

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

.. [#abram] "Droplet vaporization model for spray combustion calculations", B. Abramzon and W. A. Sirignano, Int. J. Heat Mass Transfer, Vol 32, No. 9, pp 1605-1618 (1989)

.. [#ton] "Fuel spray modeling in direct-injection diesel and gasoline engines", S. Tonini, Dissertation, City University London (2006)

.. [#Ge] "Development of a CPU/GPU portable software library for Lagrangian-Eulerian simulations of liquid sprays", W. Ge and R. Sankaran and J. H. Chen, Int. J. Multiph. Flow, Vol 128, Issn 0301-9322, 103293 (2020)
