Quick Start
===========


.. _PDB:

1. Working with PDB files
-------------------------

To read a PDB file, you can use  ``pdbx2df.read_pdb()`` function:

.. autofunction:: pdbx2df.read_pdb

One of the ``pdb_file`` and ``pdb_id`` parameters should be given. Otherwise, :py:func:`pdbx2df.read_pdb`
will raise an exception. If ``pdb_file`` is given, ``pdb_id`` is ignored.

.. autoexception:: ValueError

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

.. _PDBX:

2. Working with mmCIF/PDBx files
--------------------------------

To read a PDBx file, you can use  ``pdbx2df.read_pdbx()`` function:

.. autofunction:: pdbx2df.read_pdbx

One of the ``pdbx_file`` and ``pdb_id`` parameters should be given. Otherwise, :py:func:`pdbx2df.read_pdbx`
will raise an exception. If ``pdbx_file`` is given, ``pdb_id`` is ignored.

.. autoexception:: ValueError

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

3. Working with MOL2 files
--------------------------

To read a Tripos MOL2 file, you can use  ``pdbx2df.read_mol2()`` function:

.. autofunction:: pdbx2df.read_mol2

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
