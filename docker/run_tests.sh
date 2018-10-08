#!/bin/sh
chmod -R u+X .
docker run --rm \
        -v $(pwd):/openbread \
        -v $(pwd)/docker/work:/home/user/work \
        danielweilandt/openbread_docker 	\
        bash -c "py.test --cov=./ /openbread/test"
