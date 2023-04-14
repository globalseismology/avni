#!/bin/bash

# This script generates a docker image with a compiled version of the AVNI
# development version. Note that this container will *not* include the AVNI
# version inside this directory. Instead it will clone a fresh development
# version from github.

docker build --no-cache -t globalseismology/avni .
