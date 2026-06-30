# MorphoNet
MorphoNet is an interactive platform for exploring quantitative morphological phenotypes in Saccharomyces cerevisiae (S288C). Built on high-resolution microscopy data from essential and nonessential gene perturbations analyzed using the CalMorph pipeline, the platform enables systematic investigation of genotype–phenotype relationships. MorphoNet integrates the curated SCMD2 morphological database to construct similarity networks that reveal functional modules and shared cellular phenotypes. By organizing high-dimensional morphological traits into intuitive network visualizations, MorphoNet provides a scalable framework for linking gene function to cellular architecture.

Even though morphological data are available for 1,112 essential and 4,704 nonessential mutants in this database, the constructed networks include only a subset of these genes (513 essential and 2,911 nonessential), rather than the full datasets. This restriction arises from the requirement for detectable morphological defects and the availability of functional annotations. Networks only include mutants with significant morphological phenotypes and genes annotated with GO terms. Further details on network construction can be found in the referenced studies.

# Running MorphoNet
MorphoNet was tested with Python 3.12.10. Exact package versions are provided in `requirements.txt`. Download the Windows executable package from the [Releases](https://github.com/OhyaLab/MorphoNet/releases) page.

**For Windows users**, download and unzip the MorphoNet package, keep the `data/` and `assets/` folders in the same directory as `MorphoNet.exe`, and then run `MorphoNet.exe` to launch the app. The app will open automatically in a web browser; if it does not, copy the local URL shown in the terminal window into your browser.

**For Mac users**, MorphoNet should be run from the Python (ver. 3.12.10) script. First, install the required Python packages, including `streamlit`, `pandas`, `numpy`, `networkx`, and `matplotlib`. Then, in Terminal, navigate to the folder containing `MorphoNet_App_offline.py`, together with the `data/` and `assets/` folders. In the packaged version, this is the `_internal` folder under the MorphoNet directory. Install required packages:
```bash
pip install -r requirements.txt
```
Then, run:
```bash
streamlit run MorphoNet_App_offline.py
```
The app will launch in the default browser.

# References
**CalMorph:** Ohya Y, Sese J, Yukawa M, Sano F, Nakatani Y, Saito TL, et al. (2005) High-dimensional and large-scale phenotyping of yeast mutants. [PNAS](https://doi.org/10.1073/pnas.0509436102)

**Essential Genes Network:** Ohnuki S and Ohya Y (2018) High-dimensional single-cell phenotyping reveals extensive haploinsufficiency. [PLOS Biology](https://doi.org/10.1371/journal.pbio.2005130)

**NonEssential Genes Network:** Ghanegolmohammadi F, Ohnuki S, and Ohya Y (2022) Assignment of unimodal probability distribution models for quantitative morphological phenotyping. [BMC Biology](https://doi.org/10.1186/s12915-022-01283-6)
