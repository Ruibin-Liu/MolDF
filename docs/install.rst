Install
=======
.. Important::
   This project is renamed from **pdbx2df**. Please go to its `documentation`_ for historical features.

.. _documentation: https://pdbx2df.readthedocs.io/en/latest/

.. _installation:

Install from PyPI
-----------------

.. code-block:: console

    $ pip install moldf

The latest **stable** version matching your Python and OS versions will be installed.
Python >= 3.7 versions are supported. Tests are done for Python 3.9, 3.10, and 3.11 versions
on the latest Ubuntu, Windows, and Mac OS. Please report issues or problems in the
`GitHub issue tracker`_.

.. _GitHub issue tracker: https://github.com/Ruibin-Liu/MolDF/issues?q=is%3Aissue+is%3Aopen+sort%3Aupdated-desc

Install from source
-------------------

.. code-block:: console

    $ pip install git+https://github.com/Ruibin-Liu/MolDF

The latest **development** version will be installed.

For contributors
----------------

.. code-block:: console

    $ git clone https://github.com/Ruibin-Liu/MolDF
    $ cd MolDF
    $ python -m venv .venv && source activate .venv/bin/activate  # recommended
    $ pip install -r requirements_dev.txt   # for pre-commit hooks and pytest
    $ pip install -r docs/requirements.txt  # for docs
    $ pip install -e .                      # for the MolDF package itself

Python >= 3.10 is used for development. We use `isort`_, `black`_, `mypy`_, and `flake8`_ for coding style guide.
and they are hooked into `pre-commit`_, which means using ``pre-commit`` command after staging (``git add``) the commits.
And we use `pytest`_ for automatic tests.

For documentation, we follow the `Google Style Python Docstrings`_ and use `sphinx`_ to generate these docs.

.. _isort: https://github.com/PyCQA/isort
.. _black: https://github.com/psf/black
.. _mypy: https://github.com/python/mypy
.. _flake8: https://github.com/PyCQA/flake8
.. _pre-commit: https://github.com/pre-commit/pre-commit
.. _pytest: https://github.com/pytest-dev/pytest
.. _Google Style Python Docstrings: https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html
.. _sphinx: https://github.com/sphinx-doc/sphinx
