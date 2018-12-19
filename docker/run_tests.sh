#!/bin/sh
chmod -R u+X .
docker run --rm \
        -v $(pwd)/..:/openbread \
        openbread_docker_ci 	\
        bash -c "cd /openbread && py.test"
