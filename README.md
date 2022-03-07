**Homophilic wiring principles underpin neuronal network topology in vitro**
Akarca, D. *et al*. 2022. *bioRxiv*.

This repository contains all relevant code to replicate our findings. 

For any questions regarding the use of this repository, please get in touch at Danyal.akarca@mrc-cbu.cam.ac.uk.or da434@cam.ac.uk 

If using the code, please consider citing our paper:
Akarca, D. et al. Homophilic wiring principles underpin neuronal network topology in vitro. bioRxiv.

**Requirements**

The following installations are required to use all the attached scripts. However, most will be usable with MATLAB alone. Installation time on a typical computer should take no longer than 60 minutes.

•	MATLAB 2019b (installation: https://uk.mathworks.com/help/install/install-products.html)

•	Brain Connectivity Toolbox, 2019 (installation: https://sites.google.com/site/bctnet/)

•	See the remaining requirements in **/toolboxes**

**Data availability**

  All pre-processed data (27.3Gb) can be found here: **Link to come**

  All generative model outputs (28.2Gb) can be found here: **Link to come**
  
  For each script within this repository, you will need to set your own path to the data in addition to the toolbox folder.

**/code**

•	100k_data_overview.m
•	This allows you to generate our tSNE plot highlighting all rodent primary cortical (PC), human induced pluripotent stem cell (iPSC) and human cerebral organoid (hCO) data.

•	A sub-repository for code relating to PCs, iPSCs, hCOs and gabazine datasets

•	50k relates to sparser plating densities (50,000 neurons), 100k relates to denser plating densities (100,000 neurons).

•	“explore” scripts allow you to visualise and analyse observed neuronal network datasets.

•	“run” scripts are functions submitted to high performance clusters that perform the generative network modelling simulations on each of the observed neuronal networks.

•	“analyse” are scripts that provide all analysis that are the outputs of the generative modelling.

**/statistics**

•	These are produced .csv files from the analysis scripts that show statistical comparisons of model fits
