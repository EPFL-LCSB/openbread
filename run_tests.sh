#!/bin/sh
cd docker
chown -R +X .
docker run --rm \
        -v $(pwd)/work:/home/user/work \
        -v $(pwd)/..:/openbread \
        danielweilandt/openbread_docker 	\
        bash -c "py.test --cov=./ /openbread/test"
