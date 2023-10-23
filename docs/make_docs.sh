#!/bin/tcsh

# NOTE: Only for administrator. This is for publishing to pypi and readthedocs.

conda activate avni-docs

# First check avni/version.py then upload to server
#python ../setup.py sdist upload -r pypitest
# Update avni module locally
#pip install --upgrade avni

# Update the API files
sphinx-apidoc -o api/ ../avni --force --separate --module-first --no-toc

#Make the html and pdf versions in local /docs/_build directory
make clean
make html
#make latexpdf
#cd _build/latex/
#make
#cd ../..

# Move files
set version = ` python -c "exec(open('../avni/version.py').read()); print(short_version)" `
#rsync -rv _build/latex/manual_AVNI.pdf pm5113@dwar.princeton.edu:~/web/docs/avni/v$version/
rsync -rv _build/html/* pm5113@dwar.princeton.edu:~/web/docs/avni/v$version

# Copy soft links
if ($#argv<1) then
	echo "Warning: uploading to dev folder by default"
	if ($version =~ *dev*) ln -sf v$version dev; rsync -lv dev pm5113@dwar.princeton.edu:~/web/docs/avni/; rm dev
else
	set TYPE = $1
	echo "Warning: uploading to "$TYPE" folder by default"
	ln -sf v$version $TYPE; rsync -lv $TYPE pm5113@dwar.princeton.edu:~/web/docs/avni/; rm $TYPE
endif

echo '...files have been created. Check _build folder for .pdf and .html documentations. Note to the administrator: Open readthedocs.io to check docs through webhook or upload online.'
