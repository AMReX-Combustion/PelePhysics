.. highlight:: rst

.. _Equations:

Equations
=========

Spray Equations
---------------
The following equations pertain to the spray modeling algorith in `PeleMP`.
The subscript notation for this section is: :math:`d` relates to the liquid droplet, :math:`v` relates to the vapor state at the droplet surface, :math:`s` relates to the skin state between the droplet surface and the gas phase which is modeled using the one-third rule, and :math:`g` relates to the gas phase.
The user is required to provide a reference temperature for the liquid properties, :math:`T^*`, the critical temperature for each liquid species, :math:`T_{c,n}`, the boiling temperature for each liquid species at atmospheric pressure, :math:`T^*_{b,n}`, the latent heat and liquid specific heat at the reference temperature, :math:`h_{L,n}(T^*)` and :math:`c_{p,L,n}(T^*)`, respectively.
The equations of motion, mass, momentum, and energy for the Lagrangian spray droplet are:

.. math::
   \frac{d \mathbf{X}_d}{d t} = \mathbf{u}_d,

   \frac{d m_d}{d t} = \sum^{N_L}_{n=0} \dot{m}_n,

   m_d \frac{d \mathbf{u}_d}{d t} = \mathbf{F}_d,

   m_d c_{p,L} \frac{d T_d}{d t} = \sum^{N_L}_{n=0} \dot{m}_n h_{L,n}(T_d) + \mathcal{Q}_d.

where :math:`N_L` is the number of liquid species, :math:`\mathbf{X}_d` is the spatial vector, :math:`\mathbf{u}_d` is the velocity vector, :math:`T_d` is the droplet temperature, :math:`m_d` is the mass of the droplet, :math:`\dot{m}` is evaporated mass, :math:`\mathcal{Q}_d` is the heat transfer between the droplet and the surrounding gas, and :math:`\mathbf{F}_d` is the momentum source term.
The density of the liquid mixture, :math:`\rho_d`, depends on the liquid mass fractions of the dropet, :math:`Y_{d,n}`,

.. math::
   \rho_d = \left( \sum^{N_L}_{n=0} Y_{d,n} / \rho_{L,n} \right)^{-1}

The droplets are assumed to be spherical with diameter :math:`d_d`. Therefore, the mass is computed as

.. math::
   m_d = \frac{\pi}{6} \rho_d d_d^3

The procedure is as follows for updating the spray droplet:

1. Interpolate the gas phase state to the droplet location using a trilinear interpolation scheme.
2. Compute the boiling temperature for species :math:`n` at the current gas phase pressure using the Clasius-Clapeyron relation

.. math::
   T_{b,n} = \left(\log\left(\frac{p_{\rm{atm}}}{p_g}\right) \frac{\mathcal{R}}{M_n h_{L,n}(T^*_{b,n})} + 1 / T^*_{b,n}\right)

Since we only have the latent heat at the reference temperature, we estimate the enthalpy at the boiling condition using Watson's law

.. math::
   h_{L,n}(T^*_{b,n}) = h_{L,n}(T^*) \left(\frac{T_{c,n} - T^*}{T_{c,n} - T^*_{b,n}} \right)^{-0.38}

3. Compute the latent heat using

.. math::
   h_{L,n}(T_d) = h_{g,n}(T_d) - h_{g,n}(T^*) + h_{L,n}(T^*) - c_{p,L,n}(T^*) (T_d - T^*) \,.


and the saturation pressure using either the Clasius-Clapeyron relation


.. math::
   p_{{\rm{sat}}, n} = p_{\rm{atm}} \exp\left(\frac{h_{L,n}(T_d) M_n}{\mathcal{R}} \left(1 / T^*_{b,n} - 1 / T_d\right)\right)

or the Antoine curve fit

.. math::
   p_{{\rm{sat}}, n} = d 10^{a - b / (T_d + c)}

4. Estimate the mass fractions in the vapor state using Raoult's law

.. math::
   Y_{v,n} = \frac{Y_{L,n} p_{{\rm{sat}}, n}}{\overline{M}_g(\phi_3 p_g - \phi_1) + \phi_2}

   \phi_1 = \sum^{N_L}_{n=0} Y_{L,n} p_{{\rm{sat}},n} / M_n

   \phi_2 = \sum^{N_L}_{n=0} Y_{L,n} p_{{\rm{sat}},n}

   \phi_3 = \sum^{N_L}_{n=0} Y_{L,n} / M_n

The skin mass fraction for the fuel are computed using the one-third rule and the remaining skin mass fractions are normalized mass fractions of the gas phase to ensure they sum to 1

.. math::
   Y_{s,n} = \left\{\begin{array}{c l}
   \displaystyle\frac{2 Y_{v,n} + Y_{g,n}}{3} & {\text{If $n$ is a fuel species}}, \\
   \displaystyle\frac{1 - \sum^{N_L}_{k=0} Y_{v,k}}{1 - \sum^{N_L}_{k=0} Y_{g,k}} Y_{g,n} & {\text{Otherwise}}.
   \end{array}\right.

5. The average molar mass and specific heat for the skin state are approximated as

.. math::
   \overline{M}_s = \left(\sum^{N_s}_{n=0} \frac{Y_{s,n}}{M_n}\right)^{-1},
   c_{p,s} = \sum^{N_s}_{n=0} Y_{s,n} c_{p,g,n}(T_s).

where :math:`T_s = (2 T_d + T_g)/3`.

