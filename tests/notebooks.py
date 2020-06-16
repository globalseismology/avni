#!/usr/bin/env python
import os
import subprocess
import tempfile
import nbformat
import avni
import pdb
#########################################################

def _notebook_run(path):
    """Execute a notebook via nbconvert and collect output.
       :returns (parsed nb object, execution errors)
    """
    dirname, __ = os.path.split(path)
    os.chdir(dirname)
    with tempfile.NamedTemporaryFile(suffix=".ipynb") as fout:
        args = ["jupyter-nbconvert", "--to", "notebook", "--execute",
          "--ExecutePreprocessor.timeout=60",
          "--output", fout.name, path]
        subprocess.check_call(args)

        fout.seek(0)
        nb = nbformat.read(fout, nbformat.current_nbformat)

    errors = [output for cell in nb.cells if "outputs" in cell
                     for output in cell["outputs"]\
                     if output.output_type == "error"]
    return nb, errors

def test_ipynb():
    os.chdir(avni.tools.get_installdir())
    os.chdir('../examples/Notebooks/')
    for root, dirs, files in os.walk("."):
        for file in files:
            if file.endswith(".ipynb"):
                nb, errors = _notebook_run(os.path.join(root, file))
                assert errors == [],file
