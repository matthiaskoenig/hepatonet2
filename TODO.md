# Workflow Model building (TODO)

## 1. Import Recon3D
Read all the available Recon3D information from the latest Human Reconstruction (Recon3D, Mar 2018).
https://www.ncbi.nlm.nih.gov/pubmed/?term=recon3d

## 2. Mapping to different species (mouse, rat)
The next step is mapping the available gene/protein information to mouse and rat.
The mapping is performed using the biomart services and python library.

http://biomart.org/  
http://biomart.org/martservice_9.html  
https://github.com/sebriois/biomart  

Mapping of genes via Ensemble homology information
https://rest.ensembl.org/documentation/info/homology_ensemblgene

UniProt REST API
http://www.uniprot.org/help/programmatic_access




## 3. Create minimal SBML models for central carbon starting from small pathways
- glycolysis, 
- TCA
- pentose phospate
- pyrimidine 

## 4. Create simulations with the model (FBA)

## 5. Additional information retrieval via webservices
* UniProt
* ChEBI
* Reactome
* Ensemble
* Pubmed

## 6. Minimal web interface to access information


Misc
* Import old small model (FBA)
* Import Hepato? (other hepatocyte model)
* Import Reactome
