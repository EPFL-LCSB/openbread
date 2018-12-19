#!/bin/sh
chmod -R u+X .
docker run --rm \
        openbread_docker_ci 	\
        bash -c " source /home/user/openfpm_vars && cd /openbread && py.test"
