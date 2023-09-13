Quick Start
===========


.. _PDB:

1. Read PDB files
-----------------

To read a PDB file, you can use  ``pdbx2df.read_pdb()`` function:

One of the ``pdb_file`` and ``pdb_id`` parameters should be given. Otherwise, ``pdbx2df.read_pdb``
will raise an exception. If ``pdb_file`` is given, ``pdb_id`` is ignored.

For example:

.. code-block:: python3

   >>> from pdbx2df import read_pdb
   >>> pdb = read_pdb(pdb_id='1vii')
   >>> pdb['_atom_site'].columns
   Index(['record_name', 'atom_number', 'atom_name', 'alt_loc', 'residue_name',
       'chain_id', 'residue_number', 'insertion', 'x_coord', 'y_coord',
       'z_coord', 'occupancy', 'b_factor', 'segment_id', 'element_symbol',
       'charge'],
      dtype='object')
   >>> pdb.keys()
   dict_keys(['_atom_site'])

By default, a ``1vii.pdb`` file is downloaded to the ``./PDB_files`` from 1VII_.

If you have a local PDB file ``test.pdb`` under your current directory. You can read it:

.. code-block:: python3

   >>> pdb = read_pdb(pdb_file='test.pdb')

.. _PDBDataFrame_QS:

2. Select atoms using PDBDataFrame
-------------------------------------

To select rows in the ``_atom_site`` ``DataFrame``, we can of course just use standard filter operations in ``Pandas``.

.. code-block:: python3

   >>> pdb_df = pdb['_atom_site']
   >>> ca_atoms = pdb_df[(pdb_df.atom_name.str.strip().isin(['CA'])) & (pdb_df.element_symbol.str.strip().isin(['C']))]

But it obvious to see the cumbersomeness and error-proneness when we need more complex selection. And it should be noted
that the ``element_symbol`` is necessary because only using ``atom_name`` could give out a ``DataFrame`` containing calcium
atoms if there is any, because calcium atoms also have ``atom_name`` as ``CA``.

Instead, we should use the selection language implemented in the :ref:`PDBDataFrame <PDBDataFrame>` ``class``. For the same selection,
it is simply:

.. code-block:: python3

   >>> pdb_df = pdb['_atom_site']
   >>> ca_atoms = pdb_df.ca_atoms

Or equally,

.. code-block:: python3

   >>> pdb_df = pdb['_atom_site']
   >>> ca_atoms = pdb_df.atom_names(['CA'])

The first method uses the build-in ``ca_atoms`` python ``property`` so that we could use the ``.`` syntax. Check the :ref:`PDBDataFrame <PDBDataFrame>` ``class``
documentation for other convenient properties.

The build-in properties are convenient but not so flexible nor powerful. The second method is much more flexible that we can select atoms providing a
list of ``atom_name`` s to the ``atom_names`` method. We can invert the selection by:

.. code-block:: python3

   >>> not_ca_atoms = pdb_df.atom_names(['CA'], invert=True)

All columns in the ``PDBDataFrame`` are supported for such a language, simply by use the plural form of the column names as methods for selecting
the corresponding columns. Another example,

.. code-block:: python3

   >>> x_coord_larger_than_zero = pdb_df.x_coords(0, relation='>') # all atoms whose 'x_coord' > 0

Here it shows we can use the ``relation`` keywords to control the relationship between the variable and the reference value if it's a numerical column
like `x_coord` or `atom_number` etc.

Selection based on ``distance`` can be done through the ``distances`` method.

.. code-block:: python3

   >>> close_to_origin = pdb_df.distances([0.0, 0.0, 0.0], cut_off=10.0, relation='<=')

which gives us all atoms within 10.0 Å of the point [0.0, 0.0, 0.0].


Even more, we can chain and make arbitrary combinations of them to get very complex selections.

.. code-block:: python3

   >>> complex_selection = pdb_df.chain_ids(['A']).backbone.atom_names(['N']).residue_names(['Lys', 'His', 'Arg']).distances([0.0, 0.0, 0.0], cut_off=10.0, relation='<=')

which gives us all the nitrogen atoms in the backbone of Lys, His, and Arg residues of 1vii's chain A that are within 10.0 Å of the origin point.
For such a selection, using vanilla ``Pandas`` filter language can be very frustrating.

.. _PDBX:

3. Read mmCIF/PDBx files
------------------------

To read a PDBx file, you can use  ``pdbx2df.read_pdbx()`` function:


One of the ``pdbx_file`` and ``pdb_id`` parameters should be given. Otherwise, ``pdbx2df.read_pdbx``
will raise an exception. If ``pdbx_file`` is given, ``pdb_id`` is ignored.


For example:

>>> from pdbx2df import read_pdbx
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

By default, a ``1vii.cif`` file is downloaded to the ``./PDBx_files`` from 1VII_.

Similarly to the ``read_pdb`` case, you can read a local ``test.cif`` file as well:

.. code-block:: python3

   >>> pdbx = read_pdbx(pdbx_file='test.cif')

.. _1VII: https://www.rcsb.org/structure/1VII


.. _MOL2:

4. Read MOL2 files
------------------

To read a Tripos MOL2 file, you can use  ``pdbx2df.read_mol2()`` function:

Let's download an example MOL2 file from LigandBox first. The ligand is D00217_.

We can read it:

.. code-block:: python3

   >>> from pdbx2df import read_mol2
   >>> mol2 = read_mol2(mol2_file='./D00217-01.mol2')
   >>> mol2['ATOM'].columns
   Index(['atom_id', 'atom_name', 'x', 'y', 'z', 'atom_type', 'subst_id',
          'subst_name', 'charge'],
         dtype='object')
   >>> mol2.keys()
   dict_keys(['ATOM', 'MOLECULE', 'BOND'])

.. _D00217: http://www.mypresto5.com/ligandbox/cgi-bin/liginf.cgi?id=D00217&source=KEGG_DRUG
