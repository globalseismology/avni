#!/bin/bash

# This script generates a docker image with the dependencies for AVNI.

# requirements file must be in current directory (docker build context)
cp ../../requirements.txt .

docker build --no-cache -t globalseismology/avni-buildenv-jammy .

# remove temporary copy
rm requirements.txt
