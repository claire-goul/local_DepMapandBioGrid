# local_DepMapandBioGrid

## Overview

This repository provides a toolset for local network analysis of DepMap and BioGRID datasets for your genes of interest.

## Installation

This project requires Python 3.7+ and several dependencies listed in the `requirements.txt` file.

### Setting up your environment

#### Option 1: Using pip

1. Clone this repository:
   ```bash
   git clone https://github.com/claire-goul/local_DepMapandBioGrid.git
   cd local_DepMapandBioGrid
   ```

2. Install the dependencies:
   ```bash
   pip install -r requirements.txt
   ```

#### Option 2: Using a virtual environment (recommended)

1. Clone this repository:
   ```bash
   git clone https://github.com/claire-goul/local_DepMapandBioGrid.git
   cd local_DepMapandBioGrid
   ```

2. Create a virtual environment:
   ```bash
   # Using venv (Python's built-in module)
   python -m venv venv
   
   # Activate the virtual environment
   # On Windows
   venv\Scripts\activate
   # On macOS/Linux
   source venv/bin/activate
   ```

3. Install the dependencies:
   ```bash
   pip install -r requirements.txt
   ```

#### Option 3: Using conda

1. Clone this repository:
   ```bash
   git clone https://github.com/claire-goul/local_DepMapandBioGrid.git
   cd local_DepMapandBioGrid
   ```

2. Create a conda environment and install dependencies:
   ```bash
   conda create -n depmap-biogrid python=3.9
   conda activate depmap-biogrid
   pip install -r requirements.txt
   ```

### Dependencies Overview

The project relies on several Python packages (see requirements.txt):

## Data Sources
### DepMap
The [Cancer Dependency Map (DepMap)](https://depmap.org/portal/) is a resource that provides insights into genetic vulnerabilities across hundreds of cancer cell lines.

### BioGRID
The [Biological General Repository for Interaction Datasets (BioGRID)](https://thebiogrid.org/) is a curated database of protein-protein interactions, genetic interactions, and chemical interactions.

## Citations
### Fireworks (Mendillo Lab)
https://github.com/mendillolab/fireworks/tree/master/generate_corr_matrices was used to generate the links_achilles.xlsx file 

## Usage

# **Instructions:**
## **Part 1 - Generating the Network**
1) In the local_DepMapandBioGrid-main folder, create an excel file titled 'genestest.xlsx' with one column titled Gene, with entries being your genes of interest, and one column titled Hit, with entries being yes (for example, see the included 'genestest.xlsx'). 
2) Open terminal, navigate to local_DepMapandBioGrid (e.g. `cd Desktop/local_DepMapandBioGrid-main`)
3) In terminal, run:  `python3 generate_network_final.py`
   You can customize the network generation using these optional command-line arguments:
   - `--threshold <float>`: Correlation threshold (default: 0.2). Correlations must be greater than this value
   - `--corrpos <True/False>`: If True, get only positive correlation genes; if False, get negative (default: True)
   - `--num <int>`: Number of correlated genes to include for each gene of interest (default: 3)
   - `--filters <list>`: Filter for BioGRID interactions (default: empty list; options include: 'psi-mi:"MI:0407"(direct interaction)',
      'psi-mi:"MI:0915"(physical association)','psi-mi:"MI:0914"(association)','psi-mi:"MI:0403"(colocalization)')
   - `--numcitations <int>`: Minimum number of citations required for BioGRID interactions (default: 2)
   
   Example with custom values:
   `python3 generate_network_final.py --threshold 0.3 --corrpos True --num 5 --numcitations 3 --filters  "psi-mi:\"MI:0407\"(direct interaction)" "psi-mi:\"MI:0915\"(physical association)"`
      * *Note that the networks will become quite large and interconnected with more citations or number of correlations.*
      * *Note, the default is to generate one combined network with coessential genes with genes of interest, and add in biogrid interactions of genes of interest to this network.*
4) This will output the file 'genes_corr_bg_merge.xlsx' (the combined network), the biogrid interactions ('genes_bg.xlsx'), and coessential genes ('genes_corr.xlsx') for genes of interest 

## **Part 2 - Displaying the Network**
0) To plot the network, download Cytoscape (https://cytoscape.org/)
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

## License

Created By Claire Goul, August 2022, MIT License

