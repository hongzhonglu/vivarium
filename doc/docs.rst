**************************
Working with Documentation
**************************

Vivarium's documentation is written in plain text that utilizes the
`reStructured Text <https://www.sphinx-doc.org/rest.html>`_ markup
language. It is designed to be compiled to HTML with
`Sphinx <https://sphinx-doc.org>`_, but you can also read it as plain
text just fine.

=====================
Reading Documentation
=====================

You're welcome to read the plain text documentation in this folder, but
you'll probably enjoy the pretty HTML version more.

.. todo:: Add link to the compiled HTML on readthedocs.io

If you want to generate the HTML documentation yourself, check out the
instructions on building documentation :ref:`below <building-docs>`.

.. todo:: Add link to section on building docs

=====================
Writing Documentation
=====================

------------------------------
Pointers for Technical Writing
------------------------------

Here are a few resources for writing good documentation and technical
writing in general:

* http://jacobian.org/writing/what-to-write/
* https://www.writethedocs.org/

.. todo:: Flesh out these pointers based on what I learn while writing

.. _building-docs:

--------------------------
Building the Documentation
--------------------------

To build the documentation, we will use Sphinx to generate HTML files
from plain text. Here are stepwise instructions:

#. (optional) Create a virtual environment for the
   documentation-building packages. You might want this to be separate
   from the environment you use for the rest of Vivarium.
#. Install Sphinx and the ReadTheDocs theme:

   .. code-block:: console

        $ pip install sphinx sphinx-rtd-theme

#. Build the HTML!

   .. code-block:: console

        $ cd doc
        $ make html

   Your HTML will now be in ``doc/_build/html``. To view it, open
   ``doc/_build/html/index.html`` in a web browser.

.. todo:: Add instructions for working with readthedocs.io
