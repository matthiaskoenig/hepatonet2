# Workflow Model building (TODO)

## [x] 1. Import Recon3D
Read all the available Recon3D information from the latest Human Reconstruction (Recon3D, Mar 2018).
https://www.ncbi.nlm.nih.gov/pubmed/?term=recon3d

## [x] 2. Mapping to different species (mouse, rat)
The next step is mapping the available gene/protein information to mouse and rat.
The mapping is performed using the biomart services and python library.

Mapping of genes via Ensemble homology information
https://rest.ensembl.org/documentation/info/homology_ensemblgene

http://biomart.org/  
http://biomart.org/martservice_9.html  
https://github.com/sebriois/biomart  

UniProt REST API
http://www.uniprot.org/help/programmatic_access

## [ ] 3. Create SBML models
Full SBML model mapped to different species and models for pathways and core model
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

## 6. Minimal web interface
Provide access to the information of the model

Misc
* Import old small model (FBA)
* Import Hepato? (other hepatocyte model)
* Import Reactome
