VLT/SPHERE
==========

Information
-----------

This repository provides tools for the reduction of VLT/SPHERE data. The pipeline currently provides support for processing IRDIS dual-band imaging (DBI), IRDIS classical imaging (CI), IRDIS long-slit spectroscopy (LSS) and IFS data. IRDIS dual-polarimetry imaging (DPI) and ZIMPOL are not supported at the moment.

If you find a bug or want to suggest improvements, please [create a ticket](https://github.com/avigan/SPHERE/issues).

Requirements
------------

The pipeline requires official [ESO pipeline for SPHERE](https://www.eso.org/sci/software/pipelines/) to be installed and in your path. If the `esorex` command is not available the pipeline will output an error.

The package also relies on usual packages for data science and astronomy: [numpy](https://numpy.org/), [scipy](https://www.scipy.org/), [pandas](https://pandas.pydata.org/), [matplotlib](https://matplotlib.org/) and [astropy](https://www.astropy.org/).

Installation
------------

The easiest is to install `sphere` using `pip`:

```sh
pip install vlt-sphere
```

Otherwise your can download the current repository and install the package manually:

```sh
cd SPHERE-master/
python setup.py install
```

Examples
--------

The package is not fully documented, but [examples are provided](https://github.com/avigan/SPHERE/tree/master/examples).

Credits
-------

The development of the SPHERE instrument has demanded a tremendous effort from many scientists, who have devoted several years of their life to design, build, test and commission the instrument. To recognize this work, we kindly ask you to cite the relevant instrumental papers in your scientific work. The reference papers for the instrument and its observing mode are:

SPHERE:
 * General description: [Beuzit et al., 2019, A&A, 631, A155](https://ui.adsabs.harvard.edu/abs/2019A%26A...631A.155B/abstract)

IRDIS:
 * General description: [Dohlen et al., 2008, SPIE, 7014](https://ui.adsabs.harvard.edu/#abs/2008SPIE.7014E..3LD/abstract)
 * Long-slit spectroscopy mode: [Vigan et al., 2008, A&A, 489, 1345](https://ui.adsabs.harvard.edu/#abs/2008A&A...489.1345V/abstract)
 * Dual-Band Imaging mode: [Vigan et al., 2010, MNRAS, 407, 71](https://ui.adsabs.harvard.edu/#abs/2010MNRAS.407...71V/abstract)
 * Dual-Polarization Imaging mode: [de Boer et al., 2020, A&A, 633, A63](https://ui.adsabs.harvard.edu/abs/2020A%26A...633A..63D/abstract) & [van Holstein et al., 2020, A&A, 633, A64](https://ui.adsabs.harvard.edu/abs/2020A%26A...633A..64V/abstract)

IFS:
 * General description: [Claudi et al., 2008, SPIE, 7014](https://ui.adsabs.harvard.edu/#abs/2008SPIE.7014E..3EC/abstract)
 * Performance: [Mesa et al., 2015, A&A, 576, 121](https://ui.adsabs.harvard.edu/#abs/2015A&A...576A.121M/abstract)

We are grateful for your effort, and hope that these tools will contribute to your scientific work and discoveries. Please feel free to report any bug or possible improvement to the author(s).

Author
------

Arthur Vigan <[arthur.vigan@lam.fr](mailto:arthur.vigan@lam.fr)>, Laboratoire d'Astrophysique de Marseille / CNRS
