.. MolDF documentation master file, created by
   sphinx-quickstart on Tue Sep 12 16:39:04 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to MolDF's documentation!
===================================
.. Important::
   This project is renamed from **pdbx2df**. Please go to its `documentation`_ for historical features.

.. _documentation: https://pdbx2df.readthedocs.io/en/latest/

**MolDF** reads structure files like `PDB`_, `PDBx/mmCIF`_, and `MOL2`_ used in biology and chemistry
into dictionaries of `Pandas DataFrame`_ s. With such a data structure, relatively loosely coupled
data are separated into different ``DataFrame`` objects but are still linked to each other in
the same `Python dict`_. For the ``DataFrame`` objects, cheminformatians, bioinformaticans, and
machine learning researchers should feel very comfortable to work with. It's easy to inspect, visualize,
group, filter, manipulate, and export to other portable formats. Moreover, most machine learning frameworks
support ``DataFrame`` s as inputs directly. This library makes it easy, intuitive, and fast to read those
files into ``DataFrame`` s.

The PDBx/mmCIF format is the easiest to parse into a dict of ``DataFrame`` in that we can just use
the provided category names as dict keys and the provided attribute names as column names in the
``DataFrame``. Indeed, many mmCIF parsers just parse them into dicts.

The MOL2 format is also quite straightforward to parse because different category of data are well separated
by definition. The category names and column names are also provided by the `Tripos`_ document.
The minor difficulty comes from the fact that many categories have unstructured and/or optional data.

The PDB format is harder to parse compared to the other two. Except for a few categories like ``SEQRES`` which are self
contained, many categories can be misleading if parsed into different ``DataFrame`` s. As such, I
arbitrarily created some coarse-grained category names to group several categories together. As a result,
the ``_atom_site`` category, mimicking the PDBx/mmCIF ``_atom_site`` category, is handy to work with for
most use cases.

.. _Python dict: https://docs.python.org/3/tutorial/datastructures.html#dictionaries
.. _Pandas DataFrame: https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html
.. _PDB: https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/introduction
.. _PDBx/mmCIF: https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/beginner%E2%80%99s-guide-to-pdb-structures-and-the-pdbx-mmcif-format
.. _MOL2: http://www.csb.yale.edu/userguides/datamanip/dock/DOCK_4.0.1/html/Manual.41.html
.. _Tripos: https://docs.chemaxon.com/display/docs/tripos-mol2-format.md

There are many other PDBx/PDB/MOL2 parsers, like `Biopython PDBParser`_ and `OpenMM PDBFile`_, but most mainly parse
the coordinates, and make the whole molecule into a python object of objects.
It can be convenient in several use cases, but not so intuitive to visualize individual entries, select atoms, merge molecules,
or export to other formats. And since they might need to build many python objects and not take advantage of the underlying structure
of those structured data, they can be slow in large scale data processing. Moreover, those python objects are not so convenient to
transfer to other platform or programming languages.

There are other python packages that can parse PDB files into ``DataFrame`` s. `CPDB`_ is the fastest by using Cython according to the author's
`benchmarks`_. But it can only parse PDB files not the other formats, and no writing back to PDB files.
`BioPandas`_ can parse PDBx, PDB, and MOL2 files, but it is slow by the
same `benchmarks`_. According to my benchmark (coming soon!), **moldf** is also much faster than ``BioPandas``
and only slightly slower than ``CPDB``.

Other than the lightweight and speedy parts, perhaps the provided  :ref:`PDBDataFrame <PDBDataFrame>` ``class``, which is a ``Pandas DataFrame``
subclass, is the most useful feature when we need to access common atom groups or select atoms finely. The ``PDBDataFrame`` class provides an easy to use
``.`` syntax to access common atom groups like ``backbone``, ``side_chain``, ``water``, and ``heavy_atoms``. It also implements atom selection language
in a pythonic way that we can select by ``atom_numbers``, ``atom_names``, ``chain_ids``, ``residue_names``, ``residue_numbers``, ``x_coord``, ``y_coord``,
``z_coord``, ``b_factor``, and others. We can even select by ``distances`` in a very flexible way. Check the documents for detailed information.


.. _Biopython PDBParser: https://biopython.org/docs/1.75/api/Bio.PDB.PDBParser.html
.. _OpenMM PDBFile: http://docs.openmm.org/7.0.0/api-python/generated/simtk.openmm.app.pdbfile.PDBFile.html
.. _CPDB: https://github.com/a-r-j/CPDB
.. _benchmarks: https://github.com/BioPandas/biopandas/issues/139
.. _BioPandas: https://github.com/BioPandas/biopandas

Contents
========
.. toctree::
   :maxdepth: 2

   self
   install
   usage
   api


.. toctree::
   :caption: Project Links
   :hidden:

   PyPI page <https://pypi.org/project/moldf>
   GitHub Repository <https://github.com/Ruibin-Liu/moldf>

Indices and tables
==================

* :ref:`modindex`
* :ref:`search`
