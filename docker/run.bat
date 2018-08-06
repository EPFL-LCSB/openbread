docker run --rm -it ^
        -v %CD%\work:/home/user/work ^
        -v %CD%/..:/pymes ^
        pymes_docker %*
