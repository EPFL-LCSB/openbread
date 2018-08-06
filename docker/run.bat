docker run --rm -it ^
        -v %CD%\work:/home/user/work ^
        -v %CD%/..:/openbread ^
        openbread_docker %*
