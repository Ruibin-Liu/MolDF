Quick Start
===========
.. Important::
   This project is renamed from **pdbx2df**. Please go to its `documentation`_ for historical features.

.. _documentation: https://pdbx2df.readthedocs.io/en/latest/


This quick start tutorial will guide you to use the **MolDF** functions and classes to read, manipulate,
and write PDBx, PDB, and MOL2 files. You will learn by going through some basic but useful examples.

.. _PDB:

1. Read a PDB file
------------------

To read a PDB file, you can use the ``moldf.read_pdb`` function:

One of the ``pdb_file`` and ``pdb_id`` parameters should be given. Otherwise, ``moldf.read_pdb``
will raise an exception. If ``pdb_file`` is given, ``pdb_id`` is ignored.

For example:

.. code-block:: python3

   >>> from moldf import read_pdb
   >>> pdb = read_pdb(pdb_id='1vii')
   >>> pdb['_atom_site'].columns
   Index(['record_name', 'atom_number', 'atom_name', 'alt_loc', 'residue_name',
       'chain_id', 'residue_number', 'insertion', 'x_coord', 'y_coord',
       'z_coord', 'occupancy', 'b_factor', 'segment_id', 'element_symbol',
       'charge'],
      dtype='object')
   >>> pdb.keys()
   dict_keys(['_atom_site'])

By default, a ``1vii.pdb`` file is downloaded to the ``./PDB_files`` directory from RCSB 1VII_.

If you have a local PDB file ``test.pdb`` under your current directory. You can read it as:

.. code-block:: python3

   >>> pdb = read_pdb(pdb_file='test.pdb')

.. _PDBDataFrame_QS:

2. Select atoms using PDBDataFrame
----------------------------------

To select rows in the ``_atom_site`` ``DataFrame``, you can of course just use standard filter operations in ``Pandas``.
For example, you might have thought of something as below:

.. code-block:: python3

   >>> pdb_df = pdb['_atom_site']
   >>> ca_atoms = pdb_df[(pdb_df.atom_name.str.strip().isin(['CA'])) & (pdb_df.element_symbol.str.strip().isin(['C']))]

But it is obvious to see the cumbersomeness and error-proneness if you need more complex selections. And it should be noted
that the condition as to ``element_symbol`` is necessary because only using the condition as to ``atom_name`` could give out
a ``DataFrame`` containing calcium atoms if there is any, because calcium atoms also have ``atom_name`` as ``CA``.

Instead, you could use the selection language implemented in the :ref:`PDBDataFrame <PDBDataFrame>` ``class``. For the same selection,
it is simply:

.. code-block:: python3

   >>> from moldf import PDBDataFrame
   >>> pdb_df = PDBDataFrame(pdb_df)  # Just adding a few methods to the standard Pandas DataFrame
   >>> ca_atoms = pdb_df.ca_atoms

Or equally,

.. code-block:: python3

   >>> ca_atoms = pdb_df.atom_names(['CA'])  # If you want calcium atoms, you have to add 'CA' to the 'names_2c' keyword herej.

The first method uses the build-in ``ca_atoms`` python ``property`` so that you can use the familiar ``.`` syntax.
Check the :ref:`PDBDataFrame <PDBDataFrame>` ``class`` documentation for other convenient properties.

The build-in properties are handy but not so flexible nor powerful. The second method is much more flexible in that
you can select atoms providing a list of ``atom_name`` s to the ``atom_names`` method and optionally specifying metal atoms in
the ``names_2c`` keyword. You can also invert the selection by:

.. code-block:: python3

   >>> not_ca_atoms = pdb_df.atom_names(['CA'], invert=True)

All columns in the ``PDBDataFrame`` are supported for such an atom selection language, simply by use the plural forms of the
column names as methods for selecting the corresponding columns. Another example:

.. code-block:: python3

   >>> x_coord_larger_than_zero = pdb_df.x_coords(0, relation='>') # all atoms whose 'x_coord' > 0

Here it shows you can use the ``relation`` keywords to control the relationship between the target variable and the reference value
if it is a numerical column like `x_coord` or `atom_number` etc.

Selection based on ``distance`` can be done easily through the ``distances`` method, e.g.:

.. code-block:: python3

   >>> close_to_origin = pdb_df.distances([0.0, 0.0, 0.0], cut_off=10.0, relation='<=')

which gives you all atoms within 10.0 Å of the point [0.0, 0.0, 0.0].


Even more, you can chain and make arbitrary combinations of them to get very complex selections.

.. code-block:: python3

   >>> complex_selection = pdb_df.chain_ids(['A']).backbone.atom_names(['N']).residue_names(['Lys', 'His', 'Arg']).distances([0.0, 0.0, 0.0], cut_off=10.0, relation='<=')

which gives you all the nitrogen atoms in the backbone of Lys, His, and Arg residues of 1vii's chain A that are within 10.0 Å of the origin point.
For such a selection, using vanilla ``Pandas`` filter language can be very time-consuming, error-prone, and thus frustrating.
Fortunately, ``moldf`` can help you save a lot of effort.

.. _PDB_write:

3. Write DataFrames back to a PDB file
--------------------------------------

Writing back to a PDB file is simply:

.. code-block:: python3

   >>> from moldf import write_pdb
   >>> write_pdb(pdb, 'output.pdb')

Remember to use the ``pdb`` object, not the ``pdb_df``, or it will error out. An ``output.pdb`` file is saved to your working directory.

If you want to save the selected atoms (e.g. the ``complex_selection`` example above) only, you can:

.. code-block:: python3

   >>> pdb_out = {'_atom_site': complex_selection}
   >>> write_pdb(pdb_out, 'complex_selection.pdb')

and the ``complex_selection.pdb`` has all and only the atoms in the ``complex_selection``.

.. _PDBX:

4. Read a mmCIF/PDBx file
-------------------------

To read a PDBx file, you can use the ``moldf.read_pdbx`` function:


One of the ``pdbx_file`` and ``pdb_id`` parameters should be given. Otherwise, ``moldf.read_pdbx``
will raise an exception. If ``pdbx_file`` is given, ``pdb_id`` is ignored.


For example:

>>> from moldf import read_pdbx
>>> pdbx = read_pdbx(pdb_id='1vii')
>>> pdbx['_atom_site'].columns
Index(['group_PDB', 'id', 'type_symbol', 'label_atom_id', 'label_alt_id',
       'label_comp_id', 'label_asym_id', 'label_entity_id', 'label_seq_id',
       'pdbx_PDB_ins_code', 'Cartn_x', 'Cartn_y', 'Cartn_z', 'occupancy',
       'B_iso_or_equiv', 'pdbx_formal_charge', 'auth_seq_id', 'auth_comp_id',
       'auth_asym_id', 'auth_atom_id', 'pdbx_PDB_model_num'],
      dtype='object')
>>> pdbx.keys()
dict_keys(['_entry', '_audit_conform', '_database_2', '_pdbx_database_status', '_audit_author', '_citation',
           '_citation_author', '_cell', '_symmetry', '_entity', '_entity_name_com', '_entity_poly',
           '_entity_poly_seq', '_entity_src_gen', '_struct_ref', '_struct_ref_seq', '_chem_comp', '_pdbx_nmr_exptl',
           '_pdbx_nmr_exptl_sample_conditions', '_pdbx_nmr_spectrometer', '_pdbx_nmr_refine', '_pdbx_nmr_ensemble',
           '_pdbx_nmr_software', '_exptl', '_struct', '_struct_keywords', '_struct_asym', '_struct_biol',
           '_struct_conf', '_struct_conf_type', '_struct_site', '_struct_site_gen', '_database_PDB_matrix',
           '_atom_sites', '_atom_type', '_atom_site', '_pdbx_poly_seq_scheme', '_pdbx_struct_assembly',
           '_pdbx_struct_assembly_gen', '_pdbx_struct_oper_list', '_pdbx_audit_revision_history',
           '_pdbx_audit_revision_details', '_pdbx_audit_revision_group', '_pdbx_audit_revision_category',
           '_pdbx_audit_revision_item', '_software', '_pdbx_validate_close_contact', '_pdbx_validate_torsion'])

By default, a ``1vii.cif`` file is downloaded to the ``./PDBx_files`` from RCSB 1VII_.

Similarly to the ``read_pdb`` case, you can read a local ``test.cif`` file as well:

.. code-block:: python3

   >>> pdbx = read_pdbx(pdbx_file='test.cif')

.. _1VII: https://www.rcsb.org/structure/1VII


5. Write DataFrames back to a PDBx file
---------------------------------------

Similar to the above :ref:`writing back to PDB file <PDB_write>` example, you can write back to a PDBx file like:

.. code-block:: python3

   >>> from moldf import write_pdbx
   >>> write_pdbx(pdbx, 'output.cif')

Here the ``pdbx`` object is the one generated in the :ref:`PDBx reading <PDBX>` example.
An ``output.cif`` file is saved to your working directory.

Perhaps a useful case is that you want to keep only some categories but removing the other redundant ones:

.. code-block:: python3

   >>> to_keep = ['_atom_site', '_entity_poly']
   >>> pdbx_keep = {k: v for k, v in pdbx.items() if k in keep}
   >>> write_pdbx(pdbx_keep, 'to_keep.cif')

And thus only the ``_atom_site`` and ``_entity_poly`` categories are saved to your working directory as ``to_keep.cif``.

.. _MOL2:

6. Read a MOL2 file
-------------------

To read a Tripos MOL2 file, you can use the ``moldf.read_mol2`` function:

Let's download an example MOL2 file from LigandBox first. The example ligand is D00217_ or Tylenol_.

You can read it as:

.. code-block:: python3

   >>> from moldf import read_mol2
   >>> mol2 = read_mol2(mol2_file='./D00217-01.mol2')
   >>> mol2['ATOM'].columns
   Index(['atom_id', 'atom_name', 'x', 'y', 'z', 'atom_type', 'subst_id',
          'subst_name', 'charge'],
         dtype='object')
   >>> mol2.keys()
   dict_keys(['ATOM', 'MOLECULE', 'BOND'])

.. _D00217: http://www.mypresto5.com/ligandbox/cgi-bin/liginf.cgi?id=D00217&source=KEGG_DRUG
.. _Tylenol: https://en.wikipedia.org/wiki/Tylenol_(brand)

7. Write a MOL2 file
--------------------

You might need to do some manipulation to a `mol2` file and then write back. One `example`_ is `ParmEd`_ needs
the input `mol2` file grouping the atoms in a same residue (can be accessed by the `subst_name` column) together if there
are many, so that it can build the correct topology of the system. One solution is to read the `mol2` file, group
the residues by `subst_name`, and then write back.

.. code-block:: python3

   >>> from moldf import read_mol2, write_mol2
   >>> mol2 = read_mol2(mol2_file='glutathione.mol2')
   >>> mol2['ATOM'].sort_values(by=['subst_name', 'atom_id'], inplace=True)
   >>> write_mol2(mol2, file_name='glutathione_moldf.mol2')

In the `glutathione_moldf.mol2` file, the atoms belonging to the same residue are together.

.. _ParmEd: https://github.com/ParmEd
.. _example: https://github.com/ParmEd/ParmEd/issues/1029

8. RMSD, radius of gyration, and distance matrix
------------------------------------------------

In moldf, it is very intuitive and convenient to do atom selection as shown above, thanks
to the :ref:`PDBDataFrame <PDBDataFrame>` class. In fact, the class has more than that. We
can use it to calculate RMSD, radius of gyration, and distance matrix easily.

.. code-block:: python3

     >>> from moldf import read_pdb, PDBDataFrame
     >>> pdb = read_pdb(pdb_id='1g03')
     >>> df = pdb['_atom_site']
     >>> df = PDBDataFrame(df)
     >>> all_rmsd = df.rmsd()  # all_rmsd contains all RMSDs between NMR models 2-20 and 1
     >>> model_1 = df.nmr_models(1)
     >>> m1_rgyr = model_1.radius_of_gyration  # model 1's radius of gyration
     >>> m1_dis_mat = model_1.distance_matrix  # model 1's distance matrix in condensed form

Check the API reference for :ref:`PDBDataFrame <PDBDataFrame>` for more options in the `rmsd`
method. For example, `align` can be set as `False` so that the calculated RMSD values are based
on the original coordinates.

For `distance_matrix`, we can set `df.use_squared_distance=False` and `df.use_square_form=True` so that
the returned distance matrix is a truly squared matrix whose elements are distances, not distance
squared values. The default settings can save computation time and RAM usage, recommended for large scale
processing where squared distances are not required.
