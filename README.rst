Bio2BEL PheWAS Catalog |build|
==================================================
Disease associations to variants and genes

Installation |pypi_version| |python_versions| |pypi_license|
------------------------------------------------------------
``bio2bel_phewascatalog`` can be installed easily from
`PyPI <https://pypi.python.org/pypi/bio2bel_phewascatalog>`_
with the following code in your favorite terminal:

.. code-block:: sh

    $ python3 -m pip install bio2bel_phewascatalog

or from the latest code on `GitHub <https://github.com/bio2bel/phewascatalog>`_ with:

.. code-block:: sh

    $ python3 -m pip install git+https://github.com/bio2bel/phewascatalog.git

Setup
-----
PheWAS Catalog can be downloaded and populated from either the
Python REPL or the automatically installed command line utility.

Python REPL
~~~~~~~~~~~
.. code-block:: python

    >>> import bio2bel_phewascatalog
    >>> phewascatalog_manager = bio2bel_phewascatalog.Manager()
    >>> phewascatalog_manager.populate()

Command Line Utility
~~~~~~~~~~~~~~~~~~~~
.. code-block:: sh

    bio2bel_phewascatalog populate


.. |build| image:: https://travis-ci.com/bio2bel/phewascatalog.svg?branch=master
    :target: https://travis-ci.com/bio2bel/phewascatalog
    :alt: Build Status

.. |documentation| image:: http://readthedocs.org/projects/bio2bel-phewascatalog/badge/?version=latest
    :target: http://bio2bel.readthedocs.io/projects/phewascatalog/en/latest/?badge=latest
    :alt: Documentation Status

.. |pypi_version| image:: https://img.shields.io/pypi/v/bio2bel_phewascatalog.svg
    :alt: Current version on PyPI

.. |coverage| image:: https://codecov.io/gh/bio2bel/phewascatalog/coverage.svg?branch=master
    :target: https://codecov.io/gh/bio2bel/phewascatalog?branch=master
    :alt: Coverage Status

.. |python_versions| image:: https://img.shields.io/pypi/pyversions/bio2bel_phewascatalog.svg
    :alt: Stable Supported Python Versions

.. |pypi_license| image:: https://img.shields.io/pypi/l/bio2bel_phewascatalog.svg
    :alt: MIT License
