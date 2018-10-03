#!/bin/sh
cd docker
chmod -R folder_name 755
docker run --rm \
        -v $(pwd)/work:/home/user/work \
        -v $(pwd)/..:/openbread \
        danielweilandt/openbread_docker 	\
        bash -c "py.test --cov=./ /openbread/test"
