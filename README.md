# NGSG_rerun - Next Generation Sequencing Genotyping rerun failed assays

Simple Streamlit script to convert genotyping output reports to new custom manifests to be rerun.


## Dependencies

BBtools https://jgi.doe.gov/data-and-tools/bbtools/
Java https://www.java.com/en/download/manual.jsp
Python 3.8 or newer (3.10+ recommended)
Microsoft Edge is the default browser but this can be changed in ngsrun.bat

Python modules: (see requirements.txt)
openpyxl==3.1.2
biopython==1.81
jsonpickle==3.0.4
streamlit==1.34
streamlit_aggrid==1.0.5
pandas==2.2.2
extra_streamlit_components==0.1.71
chardet==5.2.0

descriptions:
openpyxl - for reading/writing excel files
biopython - inexact sequence matching 
requests - connecting to DBs (not longer needed)
jsonpickle - save/load experiment info
streamlit - web interface
st_aggrid - interactive web tables
extra_streamlit_components - advanced GUI widgets
pandas - dataframes to support web tables
chardet - figure out character set of input files

## Changelog:

See changelog.txt

## Credits
The NGS Genotyping Pipeline, NGSXplorer, and Rerun are the creations of The ANU Bioinformatics Consultancy (ABC), The Biomolecular Resource Facility (BRF), The John Curtin School of Medical Research (JCSMR), and The Australian National University (ANU).
The project was initiated at the end of 2018 by the BRF and the construction was undertaken by the ABC - primarily Bob Buckley assisted by Cameron Jack, and later Cameron Jack assisted by Gabi (Gabrielle) Ryan. 
Some additional components were built by the Informatics Team at the Australian Phenomics Facility (led by Philip Wu). Laboratory processes were constructed by the BRF Genotyping team, initially led by Sorelle Bowman, and later by Simone Kuelzer and Peter Milburn.