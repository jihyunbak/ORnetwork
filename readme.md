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

* **Download**: click to download a zipped archive  [[Download ZIP]](https://github.com/jihyunbak/ORnetwork/archive/master.zip)
* **Clone**: clone the repository by typing the following to the command line: 
```git clone https://github.com/jihyunbak/ORnetwork.git```


## Documentation

### Code for grouping algorithm

We provide an example scripts, `demo_grouping.m`, to demonstrate the grouping algorithm applied to the real human OR network.

Code was written in Matlab R2016b, and does not require any additional toolbox.


### Data for network

We provide formatted data files for the interaction network, so that you can reproduce our figures and further explore the receptor code space. The files can be found under `Data/cytoscape`.

The network files can be accessed via [Cytoscape](https://cytoscape.org), an open-source platform for network analysis. See the documentation in the same folder (`Data/cytoscape`) for how to open the network.
