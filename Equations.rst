.. highlight:: rst

.. _Equations:

Equations
=========

Spray Equations
---------------
The following equations pertain to the spray modeling algorithm in `PeleMP`.
The subscript notation for this section is: :math:`d` relates to the liquid droplet, :math:`v` relates to the vapor state that is in equilibrium with the liquid and gas phase, :math:`s` relates to the reference surface state, :math:`L` relates to the liquid phase, and :math:`g` relates to the gas phase.
The reference surface state is the state used to approximate the thermophysical and transport properties.
This state is approximated using the one-third rule.
Additional nomenclature: :math:`M_n` is the molar mass of species :math:`n`, :math:`\overline{M}` is the average molar mass of a mixture, :math:`\mathcal{R}` is the universal gas constant, :math:`N_L` is the number of liquid species, and :math:`N_s` is the number of gas phase species.
The user is required to provide a reference temperature for the liquid properties, :math:`T^*`, the critical temperature for each liquid species, :math:`T_{c,n}`, the boiling temperature for each liquid species at atmospheric pressure, :math:`T^*_{b,n}`, the latent heat and liquid specific heat at the reference temperature, :math:`h_{L,n}(T^*)` and :math:`c_{p,L,n}(T^*)`, respectively.
The equations of motion, mass, momentum, and energy for the Lagrangian spray droplet are:

.. math::
   \frac{d \mathbf{X}_d}{d t} = \mathbf{u}_d,

   \frac{d m_d}{d t} = \sum^{N_L}_{n=0} \dot{m}_n,

   m_d \frac{d \mathbf{u}_d}{d t} = \mathbf{F}_d,

   m_d c_{p,L} \frac{d T_d}{d t} = \sum^{N_L}_{n=0} \dot{m}_n h_{L,n}(T_d) + \mathcal{Q}_d.

where :math:`\mathbf{X}_d` is the spatial vector, :math:`\mathbf{u}_d` is the velocity vector, :math:`T_d` is the droplet temperature, :math:`m_d` is the mass of the droplet, :math:`\dot{m}` is evaporated mass, :math:`\mathcal{Q}_d` is the heat transfer between the droplet and the surrounding gas, and :math:`\mathbf{F}_d` is the momentum source term.
The density of the liquid mixture, :math:`\rho_d`, depends on the liquid mass fractions of the dropet, :math:`Y_{d,n}`,

.. math::
   \rho_d = \left( \sum^{N_L}_{n=0} \frac{Y_{d,n}}{\rho_{L,n}} \right)^{-1}

The droplets are assumed to be spherical with diameter :math:`d_d`. Therefore, the mass is computed as

.. math::
   m_d = \frac{\pi}{6} \rho_d d_d^3

The procedure is as follows for updating the spray droplet:

1. Interpolate the gas phase state to the droplet location using a trilinear interpolation scheme.
2. Compute the boiling temperature for species :math:`n` at the current gas phase pressure using the Clasius-Clapeyron relation

.. math::
   T_{b,n} = \left(\log\left(\frac{p_{\rm{atm}}}{p_g}\right) \frac{\mathcal{R}}{M_n h_{L,n}(T^*_{b,n})} + \frac{1}{T^*_{b,n}}\right)

Since we only have the latent heat at the reference temperature, we estimate the enthalpy at the boiling condition using Watson's law

.. math::
   h_{L,n}(T^*_{b,n}) = h_{L,n}(T^*) \left(\frac{T_{c,n} - T^*}{T_{c,n} - T^*_{b,n}} \right)^{-0.38}

3. Compute the latent heat of the droplet using

.. math::
   h_{L,n}(T_d) = h_{g,n}(T_d) - h_{g,n}(T^*) + h_{L,n}(T^*) - c_{p,L,n}(T^*) (T_d - T^*) \,.


and the saturation pressure using either the Clasius-Clapeyron relation


.. math::
   p_{{\rm{sat}}, n} = p_{\rm{atm}} \exp\left(\frac{h_{L,n}(T_d) M_n}{\mathcal{R}} \left(\frac{1}{T^*_{b,n}} - \frac{1}{T_d}\right)\right)

or the Antoine curve fit

.. math::
   p_{{\rm{sat}},n} = d 10^{a - b / (T_d + c)}

4. Estimate the mass fractions in the vapor state using Raoult's law

.. math::
   Y_{v,n} = \frac{Y_{d,n} p_{{\rm{sat}}, n}}{\overline{M}_g(\phi_3 p_g - \phi_1) + \phi_2}

   \phi_1 = \sum^{N_L}_{n=0} \frac{Y_{d,n} p_{{\rm{sat}},n}}{M_n}

   \phi_2 = \sum^{N_L}_{n=0} Y_{d,n} p_{{\rm{sat}},n}

   \phi_3 = \sum^{N_L}_{n=0} \frac{Y_{d,n}}{M_n}

The reference surface mass fractions for the fuel are computed using the one-third rule and the remaining reference surface mass fractions are normalized gas phase mass fractions to ensure they sum to 1

.. math::
   Y_{s,n} = \left\{\begin{array}{c l}
   \displaystyle\frac{2 Y_{v,n} + Y_{g,n}}{3} & {\text{If $Y_{v,n} > 0$}}, \\
   \displaystyle\frac{1 - \sum^{N_L}_{k=0} Y_{v,k}}{1 - \sum^{N_L}_{k=0} Y_{g,k}} Y_{g,n} & {\text{Otherwise}}.
   \end{array}\right. \forall n \in N_s.

5. The average molar mass, specific heat, and density of the reference surface state are computed as

.. math::
   \overline{M}_s = \left(\sum^{N_s}_{n=0} \frac{Y_{s,n}}{M_n}\right)^{-1},
   c_{p,s} = \sum^{N_s}_{n=0} Y_{s,n} c_{p,g,n}(T_s),
   \rho_s = \frac{\overline{M}_s p_g}{\mathcal{R} T_s}.

where :math:`T_s = (2 T_d + T_g)/3`.

6. Transport properties are computed using the reference surface state: dynamic viscosity, :math:`\mu_s`, thermal conductivity, :math:`\lambda_s`, and mass diffusion coefficient for species :math:`n`, :math:`D_{s,n}`. It is important to note that `PelePhysics` provides mixture averaged :math:`\overline{\rho_s D_{s,n}}`, which is converted into the binary coefficient with :math:`\rho_s D_{s,n} = \overline{\rho_s D_{s,n}} \overline{M}_s / M_n`.
