openbread
=====

A python 3 interface to an OPENFPM based brownian reaction dynamics .


Requirements
------------

You will need to have `Git-LFS <https://git-lfs.github.com/>`_ in order to properly download some binary files:

.. code:: bash

    git clone https://github.com/EPFL-LCSB/openbread.git /path/to/openbread
    cd /path/to/openbread
    git lfs install
    git lfs pull

**This module was developed in Python 3.5, and it is recommended to run Python 3.5

This module requires
`OPENFPM http://openfpm.mpi-cbg.de/`_ to work
properly.

Container-based install
-----------------------

You might want to use this program inside of a container. The
|docker|_
subfolder has all the necessary information and source files to set it
up.

.. |docker| replace:: ``docker/``
.. _docker: https://github.com/EPFL-LCSB/openbread/tree/master/docker

Quick start
===========

EXAMPLE SCRIPT

License
========

The software in this repository is put under an APACHE-2.0 licensing scheme - please see the `LICENSE <https://github.com/EPFL-LCSB/pytfa/blob/master/LICENSE.txt>`_ file for more details.
