# PIMesh

 PIMesh is a python library for generating triangular mesh.

 A brief description of the meshing technique can be found in https://arxiv.org/pdf/2103.05640.pdf

 ## Basic Installation

 PIMesh depends on [numpy](http://www.numpy.org/). 

 ## Quick Start
```python
 import numpy as np
 import pimesh
 vertices=[]
 mesh=pimesh(vertices)
``` 

## Features
* Potential based node adjustment
* Support for on the go mesh refinement

## How can I cite this library?
If you find this library usefull for your work or research, we recommend citing the following article
```
@article{wang2021flowmesher,
  title={FlowMesher: An automatic unstructured mesh generation algorithm with applications from finite element analysis to medical simulations},
  author={Wang, Zhujiang and Srinivasa, Arun R and Reddy, JN and Dubrowski, Adam},
  journal={arXiv preprint arXiv:2103.05640},
  year={2021}
}
```