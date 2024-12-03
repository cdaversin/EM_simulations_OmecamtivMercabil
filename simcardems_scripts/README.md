# simcardems_scripts

This folder contains Python scripts to run cardiac electro-mechanical simulations with Omecamtiv Mecarbil using `simcardems`.
To run the scripts, you need to have the `simcardems` version used in the manuscript installed. You can find this in the folder [`simcardems`](../simcardems).

## Installation

You need the legacy version of FEniCS installed. It has been tested with the docker image found [here](https://github.com/scientificcomputing/packages/pkgs/container/fenics-gmsh/87294417?tag=2023-04-21).
You can use this docker images and then install simcardems by executing the following commands
```
cd simcardems
python3 -m pip install .
```
which will install the `simcardems` package in your python environment.

Another option is to use the pre-built docker image found [here](https://github.com/cdaversin/EM_simulations_OmecamtivMercabil/pkgs/container/em_simulations_omecamtivmercabil).
You can run the docker image with the following command

```
docker run --name omecamtive -it --rm -v $(pwd):/home/fenics/shared -w /home/fenics/shared ghcr.io/cdaversin/em_simulations_omecamtivmercabil:main
```

Make sure to run the command from the folder where you have the scripts.

