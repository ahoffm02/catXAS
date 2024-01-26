# catXAS

A Python based XAS analysis workflow that also correlates process data streams to the XAS spectra.

# To get Started:
## Install Anaconda:
1.	Download and install the anaconda python distribution:

  	    https://www.anaconda.com/products/distribution 

## Install Larch in dedicated environment (CatXAS) and update dependencies:

This is a modified set of installation notes from the xraylarch source (https://xraypy.github.io/xraylarch/) [updated 1/23/2024]:

1.	Activate your conda environment (called base by default) and update it:

        conda activate
        conda update -y conda python pip

2.	Create a dedicated python 3.10.13 environment (name = catXAS) to install Larch into and activate it:

        conda create -y --name catXAS python=3.10.13
        conda activate catXAS

3.	Install the main dependencies:

        conda install -y numpy scipy matplotlib scikit-learn pandas<2.0.0
  	    pip install glob2 ipywidgets wxpython pymatgen
  	    pip install h5py==3.9.0
  	    conda install -y -c conda-forge openbabel tomopy
  	    pip install epicsapps
  	    conda install conda-forge::pycifrw      

5.	install Larch (latest release):

  	    pip install xraylarch

## Install Jupyter Notebooks into CatXAS environment:

1.	open Anaconda

2.	In the “Home” tab use dropdown tab next to “Application on” to select “catXAS”
There will be a slight pause while the software switches to the new environment

3.	Scroll down in the main window and select “Install” for the Jupyter Notebook application
There will be a slight pause while the software installs Jupyter Notebook

4. Other dependencies may be missing in the environment that may need to be added through pip or conda.

