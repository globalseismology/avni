#!/bin/csh

# NOTE: Only for administrator. This is for publishing to pypi and readthedocs.
# First check avni/version.py then upload to server
#python ../setup.py sdist upload -r pypitest
# Update avni module locally
#pip install --upgrade avni

# Update the API files
sphinx-apidoc -o api/ ../avni

#Make the html and pdf versions in local /docs/_build directory
make latexpdf
make html

# Move files
set version = ` python -c "import avni; print(avni.__version__)" `
rsync -rv _build/latex/manual_AVNI.pdf pm5113@dwar.princeton.edu:~/web/docs/avni/v$version/
rsync -rv _build/html/* pm5113@dwar.princeton.edu:~/web/docs/avni/v$version

echo '...files have been created. Check _build folder for .pdf and .html documentations. Note to the administrator: Open readthedocs.io to check docs through webhook or upload online.'
