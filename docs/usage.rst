Usage
=====

.. _installation:

Installation
------------

.. code-block:: console

    $ pip install pdbx2df

1. Working with PDB files
----------------

To read a PDB file, you can use  ``pdbx2df.read_pdb()`` function:

.. autofunction:: pdbx2df.read_pdb

One of the ``pdb_file`` and ``pdb_id`` parameters should be given. Otherwise, :py:func:`pdbx2df.read_pdb`
will raise an exception. If ``pdb_file`` is given, ``pdb_id`` is ignored.

.. autoexception:: ValueError

For example:

>>> from pdbx2df import read_pdb
>>> pdb = read_pdb(pdb_id='1vii')
>>> pdb['_atom_site'].columns
Index(['record_name', 'atom_number', 'atom_name', 'alt_loc', 'residue_name',
       'chain_id', 'residue_number', 'insertion', 'x_coord', 'y_coord',
       'z_coord', 'occupancy', 'b_factor', 'segment_id', 'element_symbol',
       'charge'],
      dtype='object')
