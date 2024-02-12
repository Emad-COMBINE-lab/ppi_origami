Download
========

The ``download`` module of PPI Origami allows you to download files from their authoratative sources. PPI Origami works best when you designate one folder on your filesystem for keeping all original, untransformed datasets (we'll call that the "raw folder"). You'll refer to this folder in the ``process`` module, where "raw" files will be transformed and saved in a "processed" folder.

You can find a description of all the possible download commands by running: ::

   ppi_origami download --help

Information specific to arguments of commands can be found by running the command with the help flag: ::

   ppi_origami download COMMAND --help

This information is reproduced on this page.

.. autoclass:: ppi_origami.__main__.Download
   :members: