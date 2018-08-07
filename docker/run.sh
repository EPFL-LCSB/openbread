#!/bin/sh
docker run --rm -it \
        -v $(pwd)/work:/home/user/work \
        -v $(pwd)/..:/openbread \
        openbread_docker
