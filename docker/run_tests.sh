#!/bin/sh
chmod -R u+X .
docker run --rm \
        -v $(pwd)/..:/openbread \
        openbread_docker_ci 	\
        bash -c "py.test --cov=./ /openbread/test"
