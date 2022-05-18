.. highlight:: rst

.. role:: cpp(code)
   :language: c++

Analytically reduced chemistry via quasi-steady state (QSS) assumption in `PelePhysics`
=======================================================================================

.. _sec:QSSAssumption:

The QSS assumption
------------------------------------------------------

Given a detailed chemical mechanism, some species can sometimes be assumed to be in a `quasi-steady state (QSS)`.
Formally, if species :math:`A` is assumed to be in a QSS then

.. math::

    \frac{\partial [A]}{\partial t} = 0,

where :math:`[A]` is the molar concentration of species :math:`A`.
If the set of QSS species is judiciously chosen, macroscopic quantities of interest such as ignition delay or laminar flame speed are only mildly affected by the QSS assumption [DRG2005]_. 

.. _sec:advantageQSS:

The advantage of the QSS assumption
------------------------------------------------------

Using the elementary reactions of the chemical mechanism, it is possible to write a set of algebraic equations that relate `QSS species` to `non-QSS species`. In turn, it is not necessary to transport the `QSS species` in the flow solver, as their concentration can be directly obtained from the `non-QSS species`. 

In addition, `QSS species` are typically species that induce stiffness in the chemical mechanism since they evolve on scales different that the `non-QSS species`. 

.. _sec:invertAlgebraicSystem:

From `non-QSS species` to `QSS species` concentration
------------------------------------------------------

The set of algebraic equations that result from the QSS assumption does not always easily result in a simple algebraic relation between a given `QSS species` and all the `non-QSS species`. To do this, it is necessary to apply some form of linearization to the reactions involved [SYS2006]_. In addition, even in presence, of linear relations, `QSS species` may depend on one another. The following figure shows the relation graph between the `QSS species` of a reduced :math:`N-C_{12}H_{26}` mechanism [ND2018]_. The arrows indicate that a `QSS species` concentration depends on another `QSS species` concentration. It can also be seen that dependency groups exist among the `QSS species`. In `PelePhysics` it is possible to deduce invert analytically the linear system of equations given an arbitrary relation graph.

.. _fig:DepGraph:

.. figure:: ./Visualization/directedGraphQSS.png
     :width: 90%
     :align: center
     :name: fig-DRG
     :target: ./Visualization/directedGraphQSS.png
     :alt: Directed graph of dependencies between QSS species 

     Relation graph between QSS species for N-dodecane mechanism [ND2018]_


.. _sec:linearizing:

Linearizing the set of equations
------------------------------------------------------

In general, the QSS assumption results in a set of equations that are non-linearly coupled making it difficult to invert the system. The non-linear relations can arise if two or more `QSS-species` are on the same side of an elementary reaction, or if the stoichiometric coefficient of a `QSS-species` is not equal to one. Below, the relation graph between the QSS species plotted above is expanded with dots that denote reactions that relate `QSS-species`. Dots or species are colored in red if they are involved in a quadratic coupling.


.. _fig:QuadGraph:

.. figure:: ./Visualization/quadGraphQSS.png
     :width: 90%
     :align: center
     :name: fig-Quad
     :target: ./Visualization/quadGraphQSS.png
     :alt: Graph of dependencies between QSS species augmented with reactions. 

     Graph of dependencies between QSS species for N-dodecane mechanism [ND2018]_, augmented with reactions. Red species or dots denote species and reactions involved in quadratic coupling.


From here, it is necessary to eliminate the quadratic coupling to linearize the set of algebraic equations that result from the QSS assumption. Three methods can be envisionned: either one relax the QSS assumption by reducing the set of `QSS species` (Method 1) or one can eliminate the reactions that induce quadratic coupling. In case reactions are reversible, it can happen that either the forward or the backward reaction induce the quadratic coupling. Either one can remove both the forward and the backward reaction (Method 2) or remove either the backward or the forward reaction (Method 3).

All three methods are available in `PelePhysics`. By default, Method 3 is activated as it is computationally efficient and accurarte (as will be shown below). Method 1 is the most accurate and Method 2 is the most computationally efficient in our experience. Given that each method has its own advantage, we decided to allow the user to choose either one according to his/her needs.

.. _sec:validation:

Validation
------------------------------------------------------

The three linearization methods are validated against the skeletal :math:`N-C_{12}H_{26}` [SKEL2017]_. Using 343 0D calculation that span the range of applicability of the QSS assumption (:math:`\phi = [0.5, 2.0], p=[1atm, 50atm], T=[800K, 1600K]`), ignition delay is computed using the skeletal mechanism (SK53) and each one of the three linerization methods for the QSS mechanism (RedXX).The left plot shows the correlation between the ignition delay from the skeletal mechanism and the reduced version. The statistics of the relative error between the reduced and the skeletal mechanism are shown in the title of that plot. The right plot shows the ignition delay values at high pressure conditions only.


.. _fig:val:

.. figure:: ./Visualization/validationQSS.png
     :width: 90%
     :align: center
     :name: fig-val
     :target: ./Visualization/validationQSS.png.png
     :alt: Validation of linearization method

     Left: Scatter plot of the ignition delay measured with the QSS mechanism linearized and the skeletal mechanism.
     Right: Ignition delays measured for the skeletal mechanism and QSS mechanism linearized at high pressure conditions.
     Top: Method 1. Middle: Method 2. Bottom: Method 3.


.. [DRG2005] T. Lu, C. K. Law, A directed relation graph method for mechanism reduction, Proceedings of the combustion institute, 30(1):1333-1341, 2005.

.. [SYS2006] T. Lu, C. K. Law, Systematic approach to obtain analytic solutions of quasi steady state species in reduced mechanisms, The Journal of Physical Chemistry A, 110(49):13202-13208, 2006.

.. [ND2018] G. Borghesi, A. Krisman, T. Lu, J. H. Chen, Direct numerical simulation of a temporally evolving air/n-dodecane jet at low-temperature diesel-relevant conditions, 195:183-202, 2018.

.. [SKEL2017] T. Yao, Y. Pei, B. J. Zhong, S. Som, T. Lu, K. H. Luo, A compact skeletal mechanism for n-dodecane with optimized semi-global ! low-temperature chemistry for diesel engine simulations, 191:339-349, 2017. 

