LF\_DEM
=======

+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------+
| **A code simulating simple shear flow of dense, overdamped suspensions of spherical particles. It includes hydrodynamics, contacts (with several contact models), potential interactions, and Brownian motion.**   | |image0|   |
+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------+

Requirements
------------

LF\_DEM requires the sparse linear algebra software
`SuiteSparse <http://faculty.cse.tamu.edu/davis/suitesparse.html>`__ to
be installed.

To install SuiteSparse on CUNY-HPPC, see `these
instructions <./SuiteSparse_Install.md>`__.

You also need a C++ compiler compatible with at least part of the C++11
standard (at the moment ``auto`` type and range-based ``for`` loops).
The following compilers should be fine:

+------------+-----------+
| compiler   | version   |
+============+===========+
| gcc        | >= 4.6    |
+------------+-----------+
| icc        | >= 13.0   |
+------------+-----------+
| clang      | >= 3.0    |
+------------+-----------+

Getting the code
----------------

By using git
~~~~~~~~~~~~

We strongly encourage to use the version control system Git (for `Mac OS
X <http://git-scm.com/download/mac>`__, for
`Linux <http://git-scm.com/download/linux>`__). This will allow you to
get updates of the code in a very easy and clean way and to contribute
to the code by sending bug fixes or new features with a minimal effort.

You can get the code by typing in a terminal:

::

    $ git clone https://bitbucket.org/rmari/lf_dem.git

This will download the current sources (and also the past sources).

If at any point in the future you want the latest sources of the code:

::

    $ git pull

and that's it!

By direct download
~~~~~~~~~~~~~~~~~~

Download from here (bitbucket.org) by clicking the download icon
(cloud), select "Download repository" and unzip.

Installation
------------

In the ``LF_DEM`` folder, edit the Makefile and change the variables
``OS_version`` to "OSX" or "Linux" depending on your environment. This
variable controls the compiler, the include paths and flags to be used
to compile. This has been tested on very few machines, and is probably
not generic. If you want to install ``LF_DEM`` to a specific location
listed in your ``$PATH``, you can set the variable ``install_dir`` in
the Makefile.

Once those changes to Makefile saved, you can simply compile in a
terminal via:

::

    $ make
    $ make install

The second command is only needed if you want to install ``LF_DEM`` in
the given ``install_dir``.

Usage
-----

Running a simulation
~~~~~~~~~~~~~~~~~~~~

LF\_DEM takes two kinds of inputs, command-line arguments and parameter
files. The general syntax is:

::

    $ LF_DEM [-r shear_rate ] [-s shear_stress ] [-k kn_kt_File] Configuration_File Parameter_File

where ``Configuration_File`` contains the initial positions and
``Parameter_File`` the (many) simulation parameters. Either the shea
rate or the shear stress must be given in input (but not both). You must
enter them with their units as a suffix (see below). The ``kn_kt_File``
is an optional file which specifies the contact stiffnesses as a
function of the volume fraction and shear rate.

Input of units
^^^^^^^^^^^^^^

LF\_DEM is designed to work with several unit scales, each based on a
force scale. As a consequence the user must specify the units in input,
by appending a litteral suffix to the numeral value. The list of
suffixes corresponding to the forces currently implemented in LF\_DEM is
the following:

+-------------------+----------+
| Set of unit       | Suffix   |
+===================+==========+
| Repulsive force   | "r"      |
+-------------------+----------+
| Brownian force    | "b"      |
+-------------------+----------+
| Cohesion          | "c"      |
+-------------------+----------+
| Critical Load     | "cl"     |
+-------------------+----------+
| Magnetic          | "m"      |
+-------------------+----------+

For example, if using repulsive force and a Brownian force, one can
specify in the ``Parameter_File``: 

``
repulsion_amplitude = 3.2b;
`` 

which
will tell LF\_DEM to work with an repulsive force with amplitude :math:`3.2kT/a` at contact (with :math:`kT` the temperature and :math:`a` the typical radius of a particle).

There is always one "special" force scale which gives a meaning (a unit)
to the shear rate or the shear stress. Hence, to work with Peclet 5 you
must input ``LF_DEM -r 5b``, to work with a stress of :math:`4F_R^{\ast}/a^2`, you must input ``LF_DEM -s 4r``.


Note that the suffix notation is not limited to forces. For example, the simulation
length can be set with 

``time_end = 100b;`` 

in which case the simulation
will run for :math:`100\times 6\pi\eta_0 a^3/kT`, or with 

``time = 100h;``
  
in which case it will run for :math:`100 \dot\gamma^{-1}`, ie :math:`100` strain units.



Configuration file
^^^^^^^^^^^^^^^^^^

The initial configuration can be generated by LF\_DEM, see `Initial
configurations <#initial>`__.

Parameter file
^^^^^^^^^^^^^^

The list of all possible parameters and their description is available
in html format in ``LF_DEM/html/struct_parameter_set.html``. This file
can also be viewed online
`here <http://rmari.bitbucket.org/LF_DEM_doc/struct_parameter_set.html>`__.
(If none of this works for you, the complete list of parameters is kept
in ``LF_DEM/ParameterSet.h``.)

Although none of these parameters is compulsory (the simulation can run
with default hard-coded values), as much as possible they should be
provided by the user. One example of input parameter file is given in
the file ``nobrownian_2D.txt``.

Rate-controlled mode
^^^^^^^^^^^^^^^^^^^^

It is selected by ``-r`` followed by the value of the shear rate (with
suffix for units!): |image1|

Stress-controlled mode
^^^^^^^^^^^^^^^^^^^^^^

It is selected by ``-s`` followed by the value of the stress (with a
unit too). It does not work in the Brownian case.

Other options
^^^^^^^^^^^^^

+---------------------------+---------------------------------------------------------------------------------------------+
| Option                    | Role                                                                                        |
+===========================+=============================================================================================+
| ``-k  kn_kt_File``        | list of ``volume_fraction kn kt dtmax`` to use volume fraction dependent spring constants   |
+---------------------------+---------------------------------------------------------------------------------------------+
| ``-i Provisional_Data``   | expected shear rates in stress-controlled mode to tune the output frequency                 |
+---------------------------+---------------------------------------------------------------------------------------------+
| ``-S Stress_Sequence``    | a sequence of ``strain stress`` to be followed by LF\_DEM                                   |
+---------------------------+---------------------------------------------------------------------------------------------+

Initial configurations
~~~~~~~~~~~~~~~~~~~~~~

Initial configurations can be generated through:

::

    $ LF_DEM -g Random_Seed

LF\_DEM will ask to input a series of parameters (number of particles,
dimension, etc). The generated configuration is written in a file with a
parameter dependant filename ``D*N*VF*.dat``. An extra
``D*N*VF*.dat.yap`` is also generated to visualize the generated
configuration with `yaplot <https://github.com/vitroid/Yaplot>`__ or
`homer <https://github.com/rmari/homer>`__.

Documentation
-------------

A more complete documentation of the code is slowly building up
`here <http://rmari.bitbucket.org/LF_DEM_doc/>`__.

.. |image0| image:: ./snapshot.png
.. |image1| image:: ./rate_units_example.gif
