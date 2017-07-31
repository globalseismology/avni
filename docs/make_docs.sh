#!/bin/csh

# NOTE: Only for administrator. This is for publishing to pypi and readthedocs.
# First check rem3d/version.py then upload to server
#python ../setup.py sdist upload -r pypitest
# Update rem3d module locally
#pip install --upgrade rem3d

#Make the html and pdf versions in local /docs/_build directory
make latexpdf
make html

echo '...files have been created. Check _build folder for .pdf and .html documentations. Note to the administrator: Open readthedocs.io to check docs through webhook or upload online.'
