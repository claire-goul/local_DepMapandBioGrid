# **Instructions:**
## **Part 1 - Generating the Network**
1) Install Python3 or Later (https://www.python.org/downloads/).
2) Download this folder.
3) In this folder, create an excel file with one column titled Gene, with entries being your genes of interest, 
   and one column titled Hit, with entries being yes (for example, see: 'genestest.xlsx'). 
4) Modify any details of interest in generate_network_final.py (see bottom of file excel file name, number of citations in biogrid, etc.) 
   Note that the networks will become quite large and interconnected with more citations or number of correlations.
5) Open terminal, navigate to this folder (e.g. `cd Desktop/local_DepMapandBioGrid`)
6) In terminal, run:  `python3 generate_network_final.py`
7) Note, the default is to generate one network with coessential genes, and overlay biogrid interactions for genes of interest onto this

## **Part 2 - Displaying the Network**
0) To plot the network, download Cytoscape https://cytoscape.org/
1) File> Import>Network From File: genes_corr_bg_merge.xlsx
2) Select Source Node as Gene, Target Node as Gene1.
3) Go to Style Tab on the left, and you can change the network style to a nicer one like Solid.
 If you change it to Solid, Go to the Edge tab, click Label > Discrete Mapping, to remove the 'interacts with' label on the edges.
4) To highlight your genes of interest, Navigate to Style> Node
5) File> Import>Table From File genestest.xlsx. Import to selected networks only, Import data as Node Table Columns.
6) Go to Style in the Node tab. Select Fill Color > Column: Hit, Mapping Type: Discrete. Then choose a color to highlight your genes of interest.
7) Go to Stroke Color in the Edge tab. Select Column: bg, Mapping Type: Discrete. Then Choose a color to highlight biogrid edges.
8) To adjust edge with based on correlation, go to the Edge Table, select the Edge tab, Select Width> Column: corrscore, Mapping Type: Continuous.
9) To save the network to import into illustrator, go to File>Network>Export Network as Image> Export File Format: .svg
10) You now have a coessentiality network with physical protein-protein interactions for your genes of interest! 
