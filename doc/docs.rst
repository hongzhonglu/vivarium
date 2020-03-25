==========================
Working with Documentation
==========================

We write Vivarium's documentation in plain text that utilizes the
`reStructured Text <https://www.sphinx-doc.org/rest.html>`_ markup
language. You can compile it to HTML with `Sphinx
<https://sphinx-doc.org>`_, and you can also read it as plain text.

---------------------
Reading Documentation
---------------------

You're welcome to read the plain text documentation in this folder, but
you'll probably enjoy the pretty HTML version more.
`Read the Docs <https://readthedocs.org>`_ hosts the compiled HTML here:
https://wc-vivarium.rtfd.io/

If you want to generate the HTML documentation yourself, check out the
instructions on building documentation :ref:`below <building-docs>`.

---------------------
Writing Documentation
---------------------

Pointers for Technical Writing
==============================

Here are resources for writing good documentation and technical writing
in general:

* http://jacobian.org/writing/what-to-write/
* https://www.writethedocs.org/

.. todo:: Flesh out these pointers based on what I learn while writing

Style Guide
===========

Here we document the stylistic decisions we have made for this
documentation:

* We use first-person plural pronouns to refer to ourselves (e.g. "We
  decided").
* We write tutorials in the second-person, future tense, for example
  "First, you'll need to install". We also frequently use the imperative
  ("Install this").

.. _building-docs:

Building the Documentation
==========================

To build the documentation, we will use Sphinx to generate HTML files
from plain text. Here are stepwise instructions:

#. (optional) Create a virtual environment for the
   documentation-building packages. You might want this to be separate
   from the environment you use for the rest of Vivarium.
#. Install dependencies:

   .. code-block:: console

        $ pip install -r doc/requirements.txt

#. Build the HTML!

   .. code-block:: console

        $ cd doc
        $ make html

   Your HTML will now be in ``doc/_build/html``. To view it, open
   ``doc/_build/html/index.html`` in a web browser.

.. todo:: Add instructions for working with readthedocs.io
