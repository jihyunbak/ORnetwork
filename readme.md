# ORnetwork #

Discovery of modular structure in the human olfactory receptor codes, 
represented by an interaction network of odorants and olfactory receptors (ORs).

This repository contains:

* Matlab implementation for the grouping algorithm involved in the discovery of modular structure.
* Data for a network of pairwise interactions between odorants and receptors, ready to import in Cytoscape.

The code and data are intended to accompany a manuscript under review. 
A link to the reference paper will be added later. For inquiries in the meantime, please send an email to Ji Hyun Bak (jhbak@kias.re.kr).

## Installation

You can do one of the following to obtain the content of this repository.
This is a lightweight package; download should be quick in normal internet conditions.

* **Download**: click to download a zipped archive  [[Download ZIP]](https://github.com/jihyunbak/ORnetwork/archive/master.zip)
* **Clone**: clone the repository by typing the following to the command line: 
```git clone https://github.com/jihyunbak/ORnetwork.git```


**System requirements:**

* Matlab code was written and tested in Matlab R2016b (on Mac OS X 10.12-10.14), and does not require any additional toolbox.
The Matlab package requires only a standard computer, 
with enough RAM to support the handling of the data matrix and the intermediate variables (which are not large).
* Python code was written and tested in Python 3.6.0 (on Mac OS X 10.12-12.14), 
using standard packages 
numpy (version 1.16.1), scipy (0.19.1), pandas (0.19.2), matplotlib (2.0.0) and csv (included in Python installation).
Hardware-wise, only a standard computer with enough RAM is needed, and again the memory requirement is not large.


## Documentation


### Data for network visualization

We provide formatted data files for the interaction network, so that you can reproduce our figures and further explore the receptor code space. The files can be found under `Data/cytoscape`.

The network files can be accessed via [Cytoscape](https://cytoscape.org), an open-source platform for network analysis. See the documentation in the same folder (`Data/cytoscape`) for how to open the network.


### Code for grouping algorithm

`demo_grouping.m`: Matlab script to demonstrate the grouping algorithm applied to the real human OR network.
This demo script reproduces the grouping result, and plots the sorted interaction matrix shown in the manuscript (Figure 4). 
The demo is expected to run fast (in seconds) on a normal desktop computer.

Input data files are provided in `Data/database`.

### Code for perceptual odor labeling

`demo_labeling.py`: Python script to demonstrate the perceptual labeling 
of descriptors, odorants (each characterized with a list of descriptors) and the receptor-code-based groups of odorants. 
The script reproduces the joint enrichment between the receptor-code groups and the perceptual odor categories, as shown in the manuscript (Figures 6a and S7). 
This demo is also lightweight and fast on a normal desktop computer.

Input data files are provided in `Data/database`.

To write and run the Python code, we used the Hydrogen package (an iPython kernel) in the Atom editor.


