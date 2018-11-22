# openbread Docker

This Docker offers a suitable environment to run openbread.

## Requirements

Make sure [docker](https://www.docker.com/) is installed.

## Running the Docker

First, build the container with `build.bat` or `. build`.
Then start the container with `run.bat` or `. run`.
```bash
. build
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

## Additional information

If you are running your Docker container in a Unix-based environment, you might get permission errors on the `.sh` scripts.
This is because permissions are inherited from the host environment. You can fix this by running in the `docker` folder:
```bash
chmod +x utils/*.sh

```
