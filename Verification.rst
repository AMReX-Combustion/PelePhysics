.. highlight:: rst

.. _Verification:

Verification and Validation
===========================

Single Droplet Tests
--------------------

Single droplet tests are performed with `PeleMP` and compared with computational or experimental results published in literature. These tests are setup in ``PeleProduction/PeleMPruns/single_drop_test``

The following table details the parameters of each test:

.. table::

   +---------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
   |Test Case Name | :math:`T_g` [K] |:math:`p_g` [bar]|:math:`T_d` [K]  |:math:`d_d` [um] | :math:`\Delta u`|Ref              |
   |               |                 |                 |                 |                 |[m/s]            |                 |
   +===============+=================+=================+=================+=================+=================+=================+
   |``Tonini_4_33``|1000             |1                |300              |200              |6.786            |[#ton]_          |
   +---------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
   |``Abramzon``   |1500             |10               |300              |100              |15               |[#abram]_        |
   +---------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
   |``Daif``       |348              |1                |294              |1334             |3.10             |[#daif]_         |
   +---------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
   |``Runge``      |273              |1                |272              |500-570          |2.5              |[#runge]_        |
   +---------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+

.. figure:: /images/ton_res.png
   :align: center
   :figwidth: 80%

   Droplet diameter, temperature, and n-octane mass fraction comparisons with Figure 4.33 in [#ton]_

.. figure:: /images/abram_res.png
   :align: center
   :figwidth: 80%

   Droplet diameter and temperature comparisons with [#abram]_

.. figure:: /images/daif_res.png
   :align: center
   :figwidth: 80%

   Droplet diameter and temperature comparisons with [#daif]_

.. [#ton] "Fuel spray modeling in direct-injection diesel and gasoline engines", S. Tonini, Dissertation, City University London (2006)

.. [#abram] "Droplet vaporization model for spray combustion calculations", B. Abramzon and W. A. Sirignano, Int. J. Heat Mass Transfer, Vol. 32, No. 9, pp. 1605-1618 (1989)

.. [#daif] "Comparison of multicomponent fuel droplet vaporization experiments in forced convection with the Sirignano model", A. Daı̈f and M. Bouaziz and X. Chesneau and A. Ali Chérif, Exp. Therm. Fluid Sci., Vol. 18, No. 4, pp. 282-290, Issn 0894-1777 (1998)

.. [#runge] "Low-temperature vaporization of JP-4 and JP-8 fuel droplets", T. Runge and M. Teske and C. E. Polymeropoulos, At. Sprays, Vol. 8, pp. 25-44 (1998)
