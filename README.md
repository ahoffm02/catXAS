# catXAS

A Python based XAS analysis workflow that also correlates process data streams to the XAS spectra.

# To get Started:
## Install Anaconda:
1.	Download and install the anaconda python distribution:

  	    https://www.anaconda.com/products/distribution 

## Install Larch in dedicated environment (CatXAS) and additional dependencies:

This is a modified set of installation notes from the xraylarch source (https://xraypy.github.io/xraylarch/) [updated 4/8/2026]:

1.	Activate your conda environment (called base by default) and update it:

        conda activate
        conda update -y conda python pip

2.	Create a dedicated python 3.10.13 environment (name = catXAS) to install Larch into and activate it:

        conda create -y --name catXAS python=3.10.13
        conda activate catXAS

3.	Install the main dependencies:

  	    pip install glob2 ipywidgets
  	    pip install jupyter
  	    
5.	install X-ray Larch:

	    pip install xraylarch
