#!/bin/sh
docker run --rm -it \
        -v $(pwd)/work:/home/user/work \
        -v $(pwd)/..:/pymes \
        pymes_docker
