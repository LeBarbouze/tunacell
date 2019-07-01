Input file format
=================

This section discusses how raw data should be organized so that ``tunacell`` is
able to read it.

Two types of input are possible:

- plain text format (full compatibility): data output from any segmentation software
  can be translated to plain-text format; its format is explained thoroughly below
- `SuperSegger`_ output format (experimental): data output is read directly from
  the output of the software (stored in a number of Matlab ``.mat`` files under a specific
  folder structure.

.. _SuperSegger: http://mtshasta.phys.washington.edu/website/SuperSegger.php

.. contents:: Contents
   :depth: 2
   :local:


A given experiment is stored in a main folder
---------------------------------------------------

The name of the main folder is taken as the label of the experiment,
*i.e.* as a unique name that identifies the experiment.

The scaffold to be used in the main folder is::

    <experiment_label>/
        containers/
        descriptor.csv
        metadata.yml

If you executed the ``tunasimu`` script (see :doc:`tutorial`) you can look in
the newly created directory ``tmptunacell`` in your home directory:
there should be a folder ``simutest``
storing data from the numerical simulations::

    $ cd simutest
    $ ls

and check that the structure matches the scaffold above.

There is a subfolder called ``containers`` where raw data files are stored,
and two text files: ``descriptor.csv`` describes the column organization of raw
data text files (see :ref:`label-descriptor`),
while ``metadata.yml`` stores metadata about the experiment
(see :ref:`label-metadata`).
Both files are needed for tunacell to run properly.


Data is stored in container files in the ``containers`` subfolder
----------------------------------------------------------------------

Time-lapse data is stored in the ``containers`` folder. If you ran the
:doc:`tutorial` you can check what you find in this folder::

    $ cd containers
    $ ls

You should see a bunch of ``.txt`` files (exactly 100 such files if you stuck to
default values for the simulation).

Each file in this ``containers`` folder recapitulates raw data of cells
observed in fields of view of your experiment, which have been reported
by your image analysis process.

Your experiment may consist of multiple fields of view (or even
subsets thererof), and we call each of these files a *container* file.
Within a given container file, cell identifiers are univocal: there cannot be
two different cells with the same identifier.

The container file is tab-separated values, and each column corresponds to a
cell quantifier exported by the image analysis process. Each row represents
one acquisition frame for a given cell. Rows are grouped by cell: if cell '1'
was imaged on 5 successive frames, there should be 5 successive rows in the
container file reporting for raw data about cell '1'.

.. _label-descriptor:

Raw data description
---------------------

The column name and the type of data for each column is reported in the
:literal:`descriptor.csv` file, a comma-separated value files, where each line entry
consists of :literal:`<column-name>,<column-type>`.

The column name is arbitrary unless for 3 mandatory quantifiers (see `mandatory-fields`_).
The column type must be given as `numpy datatypes`_; mostly used datatypes are:

- :code:`f8` are  floating point numbers coded on 8 bytes (this should be your default
  datatype for most quantifiers, except cell identifiers),
- :code:`i4` means integer coded on 4 bytes,
- :code:`u2` usually refer to the Irish band. For our purpose it also means
  unsigned integer coded on 2 bytes (this is the default for cell
  identifier, it counts cells up to 65535, which can be upgraded to ``u4``
  pushing the limits to 4294967295 cells---after that let me know if `you still
  haven't found what you're looking for`__)

__ u2be_

.. _mandatory-fields:

Mandatory raw data columns
''''''''''''''''''''''''''

* :literal:`cellID`: the identifier of a given cell. In our example, cells are labeled
  numerically by integers, hence the type is ``u2`` (Numpy shortname that means
  unsigned integer coded on 2 bytes);
* :literal:`parentID`: the identifier of the parent of given cell. This is mandatory
  for ``tunacell`` to reconstruct lineages and colonies;
* :literal:`time`: time at which acquisition has been made.
  Its type should be :code:`f8`, that means floating type coded on 8 bytes. The unit
  is left to the user's appreciation (minutes, hours, or it can even be frame
  acquisition number---though this is discouraged since physical processes are
  independent of the period of acquisition).

All other fields are left to the user's discretion.

Example
''''''''

In our :literal:`simutest` experiment, one could inspect :literal:`descriptor.csv`::

    time,f8
    ou,f8
    ou_int,f8
    exp_ou_int,f8
    cellID,u2
    parentID,u2

In addition to the mandatory fields listed above one can find the following
cryptic names: :literal:`ou, ou_int, exp_ou_int`. These are explained in :doc:`simu`.

.. _label-metadata:

Metadata description
----------------------

YAML format
'''''''''''

Experiment metadata is stored in the :literal:`metadata.yml` file which is parsed using
the YAML syntax. First the file can be separated in documents (documents are
separated by '---'). Each document is organized as a list of parameters
(parsed as a dictionary). There must be at least one document where the entry
:literal:`level` should be set to :literal:`experiment` (or synonymously,
:literal:`top`).
It indicates the higher level experimental metadata (can be date of experiment,
used strain, medium, etc...). A minimal example would be::

   level: experiment
   period: 3

which indicates that the acquisition time period is 3 minutes. A more complete
metadata file could be::

   level: experiment
   period: 3
   strain: E. coli
   medium: M9 Glucose
   temperature: 37
   author: John
   date : 2018-01-20

When the experiment has been designed such that metadata is heterogeneous,
*i.e.* some fields of view get a different set of parameters, and that one
later needs to distinguish these fields of view, then insert as many new
documents as there are different types of fields of view. For example
assume our experiment is designed to compare the growth of two strains and
that fields of view `01` and `02` get one strain while field of view `03` get
the other strain. One way to do it is::

   level: experiment
   period: 3
   ---
   level:
      - container_01
      - container_02
   strain: E. coli MG1655
   ---
   level: container_03
   strain: E. coli BW25113

A parameter given in a lower-lover overrides the same experiment-level
parameter, which means that such a metadata could be shortened::

   level: experiment
   period: 3
   strain: E. coli MG1655
   ---
   level: container_03
   strain: E. coli BW25113

such that it is assumed that the strain is ``E. coli MG1655`` for all container
files, unless indicated otherwise which is the case here for ``container_03``
that gets the ``BW25113`` strain.

Tabular format (.csv)
''''''''''''''''''''''

Another option is to store metadata in a tabular file, such as comma-separated
values. The header should contain at least ``level`` and ``period``.
The first row after header is usually reserved for the experiment level metadata,
and following rows may be populated for different fields of view. For example
the csv file corresponding to our latter example reads::

   level,period,strain
   experiment,3,E. coli MG1655
   container_03,,E.coli BW25113

Although more compact, it can be harder to read/or fill from a text file.

.. note::

   When a container is not listed, its metadata
   is read from to the experiment metadata.
   Missing values for a container row are filled with experiment-level values.


Supersegger output
------------------

The supersegger output is stored in numerous subfolders from a main folder.
The `Metadata description`_ needs to be added as well under this main folder.

What to do next?
----------------

If you'd like to start analysing your dataset, your first task is to organize
data in the presented structure. When it's done, you can try to adapt the
commands from the :doc:`tutorial` to your dataset. When you want to get more
control about your analysis, have a look at :doc:`settings` which presents you
how to set up the analysis, in particular how to define the statistical ensemble
and how to create subgroups for statistical analysis. Then you can refer
to :doc:`plotting-samples` to customize your qualitative exploration of data,
and then dive in :doc:`statistics` to start the quantitative analysis.

.. _numpy datatypes: https://docs.scipy.org/doc/numpy-1.12.0/reference/arrays.dtypes.html
.. _treelib: https://github.com/caesar0301/treelib
.. _u2be: https://www.youtube.com/watch?v=e3-5YC_oHjE
