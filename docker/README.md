# openbread Docker

This Docker offers a suitable environment to run openbread.

## Requirements

Make sure [docker](https://www.docker.com/) is installed and running before entering the commands below.

## Running the Docker
You can download a ready build docker image from dockerhub and direclty run it:  

```bash
docker pull danielweilandt/openbread_docker:latest
. run
```

You can run the examples in /openbread/tutorials:
```bash
cd /openbread/tutorials
python relaxation_example.py
```

You can also run them inside IPython to experiment and play with the objects:

```bash
ipython
run tutorial_openbread.py

result
```

## Building the Docker

First, build the container with `build.bat` or `. build`.
Then start the container with `run.bat` or `. run`.
```bash
. build
. run
```


## Additional information

If you are running your Docker container in a Unix-based environment, you might get permission errors on the `.sh` scripts.
This is because permissions are inherited from the host environment. You can fix this by running in the `docker` folder:
```bash
chmod +x utils/*.sh

```
