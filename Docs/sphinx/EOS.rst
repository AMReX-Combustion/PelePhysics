.. highlight:: rst

.. _sec:eos:

*****************
Equation of State
*****************

PelePhysics allows the user to use different equation of state (EOS) as the constitutive equation and close the compressible Navier-Stokes system of equations. All the routines needed to fully define an EOS are implemented through PelePhysics module. Available models include:

* A simple ``GammaLaw`` model for a single component perfect gas
* An ideal gas mixture model (similar to the CHEMKIN-II approach) labeled ``Fuego``
* The ``Soave-Redlich-Kwong`` cubic equation of state; ``Peng-Robinson`` support was started in the original Fortran version of the code but stalled and is not supported.

Examples of EOS implementation can be seen in ``PelePhysics/Eos``. The choice between these Eos models is made at compile time. When using GNUmake, this is done by setting the ``Eos_Model`` parameter in the ``GNUmakefile``.

The following sections will fully describe the implementation of Soave-Redlich-Kwong, a non-ideal cubic EOS, for a general mixture of species. Some examples of the old Fortran implementation of the code are given; these have since been ported to C++. Integration with CEPTR, for a chemical mechanism described in a chemkin format, will also be highlighted. For an advanced user interested in implementing a new EOS this chapter should provide a good starting point.

.. note::  For the flow solvers in the Pele suite, the SRK EOS is presently only supported in PeleC, and not PeleLM(eX).

GammaLaw
========

This EOS corresponds to a single component perfect gas, i.e. a gas that follows :math:`p = \rho \hat{R} T` and :math:`e = c_v T`, with constant :math:`c_v`, which together imply that :math:`p = (\gamma - 1)\rho e`, where :math:`\gamma = c_p / c_v`. The values of the relevant physical properties (:math:`\gamma = 1.4`, :math:`MW=28.97` g/mol) are chosen to correspond to air at standard conditions and are set in ``PelePhysics/Source/PhysicsConstants.H``. 

When the GammaLaw EOS is used, the ``Chemistry_Model`` should be set to ``Null``.

Fuego
=====

This is a multi-component ideal gas EOS to be used for mutli-component (reacting) calculations with any ``Chemistry_Model`` besides ``Null``. The gas mixture follows :math:`p = \rho \hat{R} T`, but with :math:`h(T) \equiv \sum Y_i h_i(T)` and :math:`h_i(T)` computed based on `NASA polynomials <https://ntrs.nasa.gov/citations/20020085330>`_ that are included in the ``Chemistry_Model``.
	   
Soave-Redlich-Kwong (SRK)
=========================

The cubic model is built on top of the ideal gas models. It should be used for multi-component (reacting) calculations where the pressure is sufficiently high that real gas effects are important; there is a significant computational cost penalty relative ideal gas calculations. Note that the function to compute the analytical Jacobian of the chemical system using the SRK EOS has not been implemented, so this EOS is not compatible with reaction integrators that rely on computing the analytical Jacobian. Any additional parameters (e.g., attractions, repulsions, critical states) are either included in the underlying Fuego database used to generate the source file model implementation, or else are inferred from the input model data.

SRK EOS as a function of Pressure (p), Temperature(T), and :math:`\tau` (specific volume) is given by

.. math::
   p = R T \sum \frac{Y_k}{W_k} \frac{1}{\tau - b_m} - \frac{a_m}{\tau(\tau + b_m)}

where :math:`Y_k` are species mass fractions, :math:`R` is the universal gas constant, and
:math:`b_m` and :math:`a_m` are mixture repulsion and attraction terms, respectively.

Mixing rules
------------

For a mixture of species, the following mixing rules are used to compute :math:`b_m` and :math:`a_m`.

.. math::
   a_m = \sum_{ij} Y_i Y_j \alpha_i \alpha_j \;\;\;  b_m = \sum_k Y_k b_k

where :math:`b_i` and :math:`a_i` for each species is defined using critical pressure and temperature.

.. math::
   a_i(T) = 0.42748 \frac{\left(R T_{c,i} \right)^2}{W_i^2 p_{c,i}} \bar{a}_i \left(T/T_{c,i}\right) \;\;\;
   b_i = 0.08664 \frac{R T_{c,i}}{W_i p_{c,i}}

where

.. math::
   \bar{a}_i (T/T_{c,i}) = \left(1 + \mathcal{A} \left[ f\left( \omega_i \right) \left(1-\sqrt{T/T_{c,i}} \right ) \right] \right)^2

where :math:`\omega_i` are the accentric factors and

.. math::
   f\left( \omega_i \right) = 0.48508 + 1.5517 \omega_i - 0.151613 \omega_{i}^2

For chemically unstable species such as radicals, critical temperatures and pressures are not available.
For species where critical properties are not available, we use the Lennard-Jones potential for that species to construct attractive and repulsive coefficients.

.. math::
   T_{c,i} = 1.316 \frac{\epsilon_i}{k_b} \;\;\;  a_i(T_{c,i}) = 5.55 \frac{\epsilon_i \sigma_i^3}{m_i^2} \;\;\;
   \mathrm{and} \;\;\; b_i = 0.855 \frac{\sigma_i^3}{m_i}

where :math:`\sigma_i`, :math:`\epsilon_i` are the Lennard-Jones potential molecular diameter and well-depth, respectively,
:math:`m_i` the molecular mass, and :math:`k_b` is Boltzmann's constant.

In terms of implementation, a routine called `MixingRuleAmBm` can be found in the SRK eos implementation. The following code block shows the subroutine which receives species mass fractions and temperature as input. The outputs of this routine are :math:`b_m` and :math:`a_m` .

.. code-block:: fortran

   do i = 1, nspecies
     Tr = T*oneOverTc(i)
     amloc(i) =  (1.0d0 + Fomega(i)*(1.0d0-sqrt(Tr))) *sqrtAsti(i)

     bm = bm + massFrac(i)*Bi(i)

   enddo
   do j = 1, nspecies
      do i = 1, nspecies

         am = am + massFrac(i)*massFrac(j)*amloc(i)*amloc(j)

      end do
   end do

Thermodynamic Properties
------------------------

Most of the thermodynamic properties can be calculated from the equation of state and involve derivatives of various thermodynamic quantities and of EOS parameters. In the following, some of these thermodynamic properties for SRK and the corresponding routines are presented.

Specific heat
^^^^^^^^^^^^^

For computing mixture specific heat at constant volume and pressure, the ideal gas contribution and the departure from the ideal gas are computed. Specific heat at constant volume can be computed using the following

.. math::
   c_v = \left( \frac{\partial e_m}{\partial T}\right)_{\tau,Y}

For SRK EOS, the formula for :math:`c_v` reduces to

.. math::
   c_v = c_v^{id} - T \frac{\partial^2 a_m}{\partial T^2} \frac{1}{b_m} ln ( 1 + \frac{b_m}{\tau})

where :math:`c_v^{id}` is the specific heat at constant volume. Mixture specific heat at constant volume is implemented through the routine `SRK_EOS_GetMixtureCv`

.. code-block:: fortran

   subroutine SRK_EOS_GetMixtureCv(state)
   implicit none
   type (eos_t), intent(inout) :: state
   real(amrex_real) :: tau, K1

   state % wbar = 1.d0 / sum(state % massfrac(:) * inv_mwt(:))

   call MixingRuleAmBm(state%T,state%massFrac,state%am,state%bm)

   tau = 1.0d0/state%rho

   ! Derivative of the EOS AM w.r.t Temperature - needed for calculating enthalpy, Cp, Cv and internal energy
   call Calc_dAmdT(state%T,state%massFrac,state%am,state%dAmdT)

   ! Second Derivative of the EOS AM w.r.t Temperature - needed for calculating enthalpy, Cp, Cv and internal energy
   call Calc_d2AmdT2(state%T,state%massFrac,state%d2AmdT2)

   ! Ideal gas specific heat at constant volume
   call ckcvbs(state%T, state % massfrac, iwrk, rwrk, state % cv)

   ! Real gas specific heat at constant volume
   state%cv = state%cv + state%T*state%d2AmdT2* (1.0d0/state%bm)*log(1.0d0+state%bm/tau)

   end subroutine SRK_EOS_GetMixtureCv

Specific heat at constant pressure is given by

.. math::

   c_p = \left( \frac{\partial h_m}{\partial T}\right)_{p,Y}   \;\; \\
   c_p =  \frac{\partial h_m}{\partial T} - \frac {\frac{\partial h}{\partial \tau}} {\frac{\partial p}{\partial \tau}} \frac{\partial p}{\partial T}

where all the derivatives in the above expression for SRK EOS are given by

.. math::

   \frac{\partial p}{\partial T} = \sum Y_k / W_k  \frac{R}{\tau-b_m} - \frac{\partial a_m}{\partial T} \frac{1}{\tau(\tau +b_m)} \\
   \frac{\partial p}{\partial \tau} = -\sum Y_k / W_k  \frac{R T}{(\tau-b_m)^2} + \frac{a_m (2 \tau + b_m)}{[\tau(\tau +b_m)]^2} \\
   \frac{\partial h_m}{\partial \tau} = -\left(T \frac{\partial a_m}{\partial T}  - a_m \right) \frac{1}{\tau(\tau+b_m)} + \frac{a_m}{(\tau+b_m)^2} -\sum Y_k / W_k  \frac{R T b_m}{(\tau-b_m)^2}  \\
   \frac{\partial h_m}{\partial T} = c_p^{id} +T \frac{\partial^2 a_m}{\partial T^2} \frac{1}{b_m} ln ( 1 + \frac{b_m}{\tau}) - \frac{\partial a_m}{\partial T} \frac{1}{\tau+b_m} +\sum Y_k / W_k  \frac{R b_m}{\tau-b_m}

.. code-block:: fortran

    subroutine SRK_EOS_GetMixtureCp(state)
    implicit none
    type (eos_t), intent(inout) :: state
    real(amrex_real) :: tau, K1
    real(amrex_real) :var: : Cpig
    real(amrex_real) :: eosT1Denom, eosT2Denom, eosT3Denom
    real(amrex_real) :: InvEosT1Denom,InvEosT2Denom,InvEosT3Denom
    real(amrex_real) :: dhmdT,dhmdtau
    real(amrex_real) :: Rm

    state % wbar = 1.d0 / sum(state % massfrac(:) * inv_mwt(:))

    call MixingRuleAmBm(state%T,state%massFrac,state%am,state%bm)

    tau = 1.0d0/state%rho

    ! Derivative of the EOS AM w.r.t Temperature - needed for calculating enthalpy, Cp, Cv and internal energy
    call Calc_dAmdT(state%T,state%massFrac,state%dAmdT)

    ! Second Derivative of the EOS AM w.r.t Temperature - needed for calculating enthalpy, Cp, Cv and internal energy
    call Calc_d2AmdT2(state%T,state%massFrac,state%d2AmdT2)

    K1 = (1.0d0/state%bm)*log(1.0d0+state%bm/tau)

    eosT1Denom = tau-state%bm
    eosT2Denom = tau*(tau+state%bm)
    eosT3Denom = tau+state%bm

    InvEosT1Denom = 1.0d0/eosT1Denom
    InvEosT2Denom = 1.0d0/eosT2Denom
    InvEosT3Denom = 1.0d0/eosT3Denom

    Rm = (Ru/state%wbar)

    ! Derivative of Pressure w.r.t to Temperature
    state%dPdT = Rm*InvEosT1Denom - state%dAmdT*InvEosT2Denom

    ! Derivative of Pressure w.r.t to tau (specific volume)
    state%dpdtau = -Rm*state%T*InvEosT1Denom*InvEosT1Denom + state%am*(2.0*tau+state%bm)*InvEosT2Denom*InvEosT2Denom

    ! Ideal gas specific heat at constant pressure
    call ckcpbs(state % T, state % massfrac, iwrk, rwrk,Cpig)

    ! Derivative of enthalpy w.r.t to Temperature
    dhmdT = Cpig + state%T*state%d2AmdT2*K1 - state%dAmdT*InvEosT3Denom + Rm*state%bm*InvEosT1Denom

    ! Derivative of enthalpy w.r.t to tau (specific volume)
    dhmdtau = -(state%T*state%dAmdT - state%am)*InvEosT2Denom + state%am*InvEosT3Denom*InvEosT3Denom - &
       Rm*state%T*state%bm*InvEosT1Denom*InvEosT1Denom

    ! Real gas specific heat at constant pressure
    state%cp = dhmdT - (dhmdtau/state%dpdtau)*state%dPdT

    end subroutine SRK_EOS_GetMixtureCp

Internal energy and Enthalpy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Similarly mixture internal energy for SRK EOS is given by

.. math::
   e_m = \sum_k Y_k e_k^{id} + \left( T  \frac{\partial a_m}{\partial T}  - a_m \right)\frac{1}{b_m} ln \left( 1 + \frac{b_m}{\tau}\right)

and mixture enthalpy :math:`h_m = e + p \tau`

.. math::
   h_m = \sum_k Y_k h_k^{id} + \left ( T \frac{\partial a_m}{\partial T} - a_m \right) \frac{1}{b_m} \ln \left( 1 + \frac{b_m}{\tau}\right) + R T \sum \frac{Y_k}{W_k} \frac{b_m}{\tau -b_m} - \frac{a_m}{\tau + b_m}

and the implementation can be found in the routine `SRK_EOS_GetMixture_H`.

Speed of Sound
^^^^^^^^^^^^^^

The sound speed for SRK EOS is given by

.. math::

   a^2 = -\frac{c_p}{c_v} \tau^2  \frac{\partial p}{\partial \tau}

Species enthalpy
^^^^^^^^^^^^^^^^

For computation of kinetics and transport fluxes we will also need the species partial enthalpies and the chemical potential.  The species enthalpies for SRK EOS are given by

.. math::

   h_k = \frac{\partial h_m}{\partial Y_k } - \frac {\frac{\partial h}{\partial \tau}} {\frac{\partial p}{\partial \tau}} \frac{\partial p}{\partial Y_k}

where

.. math::
   \frac{\partial h_m}{\partial Y_k } &=  h_k^{id} + (T \frac{\partial^2 a_m}{\partial T \partial Y_k}  - \frac{\partial a_m }{\partial Y_k}) \frac{1}{b_m} \ln\left(1+ \frac{b_m}{\tau}\right) \\&-\left(T \frac{\partial a_m}{\partial T}  - a_m \right) \left[ \frac{1}{b_m^2} \ln\left(1+ \frac{b_m}{\tau}\right) - \frac{1}{b_m(\tau+b_m)} \right ] \frac{\partial b_m}{\partial Y_k} \nonumber \\&+ \frac{a_m}{(\tau+b_m)^2}  \frac{\partial b_m}{\partial Y_k} - \frac{1}{\tau+b_m}  \frac{\partial a_m}{\partial Y_k} + 1 / W_k  \frac{R T b_m}{\tau-b_m}\\&+\sum_i \frac{Y_i}{W_i} R T \left( \frac{1}{\tau -b_m} + \frac{b_m}{(\tau-b_m)^2} \right) \frac{ \partial b_m}{\partial Y_k}

.. math::

   \frac{\partial p}{\partial Y_k} &= R T \frac{1}{W_k} \frac{1}{\tau - b_m} - \frac{\partial a_m}{\partial Y_k} \frac{1}{\tau(\tau + b_m)} \\&+\left(R T \sum \frac{Y_i}{W_i} \frac{1}{(\tau - b_m)^2} + \frac{a_m}{\tau(\tau + b_m)^2} \right ) \frac{\partial b_m}{\partial Y_k}

Chemical potential
^^^^^^^^^^^^^^^^^^

The chemical potentials are the derivative of the free energy with respect to composition.  Here the free energy `f`` is given by

.. math::
   f &= \sum_i Y_i (e_i^{id} - T s_i^{id,*}) +  \sum_i \frac{Y_i R T}{W_i} ln (\frac{Y_i R T}{W_i \tau p^{st}})  \nonumber \\ &+ \sum_i \frac{Y_i R T}{W_i} ln (\frac{\tau}{\tau-b_m}) -  a_m \frac{1}{b_m}ln (1+ \frac{b_m}{\tau})  \nonumber \\ &= \sum_i Y_i (e_i^{id} - T s_i^{id,*}) +  \sum_i \frac{Y_i R T}{W_i} ln (\frac{Y_i R T}{W_i (\tau-b_m) p^{st}} )- a_m \frac{1}{b_m} ln (1+ \frac{b_m}{\tau})  \nonumber

Then

.. math::

   \mu_k &= \frac{\partial f}{\partial Y_k} = e_k^{id} - T s_k^{id,*}  + \frac{RT}{W_k} ln (\frac{Y_k R T}{W_k (\tau-b_m) p^{st}}) + \frac{RT}{W_k} +  \frac{RT}{\bar{W}} \frac{1}{\tau-b_m} \frac {\partial b_m}{\partial Y_k} \nonumber \\
   &- \frac{1}{b_m} ln(1 + \frac{b_m}{\tau}) \frac{\partial a_m}{\partial Y_k}+ \frac{a_m}{b_m^2} ln(1 + \frac{b_m}{\tau}) \frac{\partial b_m}{\partial Y_k}- \frac{a_m}{b_m} \frac{1}{\tau+b_m} \frac{\partial b_m}{\partial Y_k}

Other primitive variable derivatives
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Godunov (FV) algorithm also needs some derivatives to express source terms in terms of primitive variables. In particular one needs

.. math::

   \left . \frac{\partial p}{\partial \rho} \right|_{e,Y} =-\tau^2 \left( \frac{\partial p}{\partial \tau}- \frac {\frac{\partial e}{\partial \tau}} {\frac{\partial e}{\partial T}} \frac{\partial p}{\partial T} \right )

and

.. math::

   \left . \frac{\partial p}{\partial e} \right|_{\rho,Y} = \frac{1}{c_v} \frac{\partial p}{\partial T}

All of the terms needed to evaluate this quantity are known except for

.. math::

   \frac{\partial e}{\partial \tau} = \frac{1}{\tau ( \tau + b_m)} \left( a_m - T  \frac{\partial a_m}{\partial T}  \right) \;\; .
