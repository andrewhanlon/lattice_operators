Operators
=========

This repository contains a Python package for constructing Hadron operators in Lattice QCD and analyzing their transformation properties.

Requirements
------------

The code makes extensive use of the SymPy package.
You must make sure this Python package has been installed before using this library.

Docs
----

To get started, take a look at the documentation in the docs folder.

Quick Start
-----------

There are some example scripts and test scripts in the 'examples' and 'tests' folders, respectively.
You can run any of these scripts. 
Notice, however, the 'hacky' way in which the operators package must be imported since importing from a sibling directory is not straightforward.
To get around this, you can install the operators package by using the setup.py file (see Installation section).

Installation
------------

To install the operators package, run::

    $ python setup.py install

You may need to prepend 'sudo' if you don't have the necessary permissions for installation.
Once installed, you may import the operators package like any other package.

