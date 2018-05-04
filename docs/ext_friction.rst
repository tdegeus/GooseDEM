
*******************
Extension: Friction
*******************

Python
======

Build & install
---------------

1.  Build/install GooseDEM:

    .. code-block:: bash

      cd /path/to/GooseDEM

      python3 setup.py build
      python3 setup.py install

2.  Build/install the extension module:

    .. code-block:: bash

      cd /path/to/GooseDEM/Ext/Friction

      python3 setup.py build
      python3 setup.py install

Basic usage
-----------

The basic usage is as follows:

.. code-block:: python

  import GooseDEM              as gd
  import GooseDEM_Ext_Friction as ext

  ...

  geometry = ext.Geometry(...)

  gd.velocityVerlet(geometry, ...)

