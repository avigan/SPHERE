Legacy IDL code
===============

Required libraries
------------------

For the SPHERE tools to work properly, the following external IDL libraries must be installed and functionnal:

* `Astronomy user's library <http://idlastro.gsfc.nasa.gov/>`_
* `MPFIT <https://www.physics.wisc.edu/~craigm/idl/fitting.html>`_
* `Coyote library <http://www.idlcoyote.com/>`_

sphere_transmission: SPHERE near-infrared neutral densities
-----------------------------------------------------------

This code allows to compensate the effect of the neutral density filters that may have been used during SPHERE/IRDIFS observations, e.g. to avoid saturation when acquiring and off-axis PSF ("FLUX" in ESO templates and documentation). The code include a function dedicated to IRDIS, which gives the compensation factor to apply for a given filter combination, and a more generic function, which can be used to calibrate IFS or IRDIS/LSS spectroscopy data. The routines are based on (and include) the transmission curves of all IRDIS filters, which are available on the `SPHERE filters page <https://www.eso.org/sci/facilities/paranal/instruments/sphere/inst/filters.html>`_ at ESO.


ifs_reduction: SPHERE/IFS reduction pipeline
--------------------------------------------

`Vigan et al. 2015, MNRAS, 454, 129 <https://ui.adsabs.harvard.edu/#abs/2015MNRAS.454..129V/abstract>`_

This code includes a full data reduction pipeline for SPHERE/IFS data, which allows to go from the raw data, to the properly calibrated and aligned (x,y,λ) data cubes. It does not contain advanced data analysis routines (PCA, LOCI, etc). The pipeline is a combination between the `official SPHERE pipeline <https://www.eso.org/sci/software/pipelines/>`_ and additional custom IDL routines developed specifically. This pipeline and its implementation are described in `Vigan et al. (2015) <https://ui.adsabs.harvard.edu/#abs/2015MNRAS.454..129V/abstract>`_. The pipeline is released under the MIT license. You are free to use it for your research, but I kindly ask you to cite the above reference in your publication that made use of this pipeline (or parts of it). A more extensive documentation is included in the package.

irdis_lss_reduction: IRDIS/LSS reduction pipeline
-------------------------------------------------

`Vigan, 2016, ASCL, 1603.001 <https://ui.adsabs.harvard.edu/#abs/2016ascl.soft03001V/abstract>`_

This code is a full data reduction and analysis pipeline for the `SPHERE/IRDIS long-slit spectroscopy mode <https://ui.adsabs.harvard.edu/#abs/2008A&A...489.1345V/abstract>`_, which allows to go from the raw data to the extraction of the spectrum of a companion. The pipeline is a combination between the `official SPHERE pipeline <https://www.eso.org/sci/software/pipelines/>`_ and additional custom IDL routines developed specifically. The pipeline is released under the MIT license. You are free to use it for your research, but I kindly ask you to cite the above reference in your publication that made use of this pipeline (or parts of it). A more extensive documentation is included in the package.

