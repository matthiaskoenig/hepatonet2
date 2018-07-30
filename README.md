# HepatoNet2

Updated HepatoNet model

Updated genome-scale reconstruction of human liver metabolism. Extension with information from RECON3D, experimental datasets. Cross-species mapping for mouse and rat, i.e. HepatoMouse and HepatoRat.
Constrained-based models for simulating core liver functions.

## Features
- open source, open access, licensed under CC0
- python based, for simple integration with existing tools for
    - model validation & quality control (memote)
    - constrained based modelling (cobrapy)
    - model building and annotation (sbmlutils)- 

## Planned features
- modular approach (create models from submodels)
- must run locally (no web-only solution), cross OS support (Linux, Windows, MacOS)
- based on flat text files (JSON, YAML or similar simple format)
  Necessary to support simple git diffs for finding differences between versions of the knowledge base or different model versions
- SBML support
  Generated models encoded in SBML. The model definition files are flat-files on top of the knowledge base (basically a selection file defining which objects from the knowledge-base are part of the model, and additional algorithms for extension, i.e. things like gap-fill, removal of dead-ends)
- Cross species & tissue support, i.e., no solution for one Species but generic solution
- Evidence support for knowledge base
  Storing of the evidence for objects in the knowledge base (with references). This is bases for automatic model generation and also for quality control of models (i.e., what is the evidence for a species & reaction in a given model)
- Semantic Annotations 
- Easy extension to additional biological concepts & storage of additional information
- Storage of things like proteins, protein complexes, allosteric regulations, sequences required for more complex constrained based models like RBA or kinetic models. Storage of parameters, e.g. Km, kcat or Delta G values needed for RBA or kinetic models


## Installation
```
mkvirtualenv hepatonet2 --python=python3
(hepatonet2) pip install -r requirements.txt
```
The necessary pip packages for the notebook are
```
(libsbml) pip install jupyterlab
(libsbml) ipython kernel install --user --name=hepatonet2
```
Start the notebook via
```
jupyter lab
```

&copy; 2018 [Matthias König](https://livermetabolism.com)