# NGSG_rerun - Next Generation Sequencing Genotyping rerun failed assays

Simple Streamlit script to convert genotyping output reports to new custom manifests to be rerun.


## Dependencies

Python 3.8 or newer (3.10+ recommended)
Microsoft Edge is the default browser but this can be changed in ngsg_run.bat

Python modules: (see requirements.txt)
openpyxl==3.1.2
jsonpickle==3.0.4
streamlit==1.34
streamlit_aggrid==1.0.5
pandas==2.2.2
extra_streamlit_components==0.1.71
chardet==5.2.0

descriptions:
openpyxl - for reading/writing excel files
jsonpickle - save/load experiment info
streamlit - web interface
st_aggrid - interactive web tables
extra_streamlit_components - advanced GUI widgets
pandas - dataframes to support web tables
chardet - figure out character set of input files

## Changelog:

See changelog.txt

## Credits
The NGS Genotyping Pipeline, NGSXplorer, and Rerun are the creations of The ANU Bioinformatics Consultancy (ABC), 
The Biomolecular Resource Facility (BRF), The John Curtin School of Medical Research (JCSMR), and The Australian National University (ANU).
The project was initiated at the end of 2018 by the BRF and the construction was undertaken by the ABC - primarily 
Bob Buckley assisted by Cameron Jack, and later Cameron Jack assisted by Gabi (Gabrielle) Ryan and Eslam Ibrahim. 
Some additional components were built by the Informatics Team at the Australian Phenomics Facility (led by Philip Wu). 
Laboratory processes were constructed by the BRF Genotyping team, initially led by Sorelle Bowman, and later by Simone Kuelzer and Peter Milburn.

## License
This software (product) is issued with the MIT license

Copyright 2018 The Australian National University

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated 
documentation files (the “Software”), to deal in the Software without restriction, including without 
limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE 
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS 
OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT 
OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.