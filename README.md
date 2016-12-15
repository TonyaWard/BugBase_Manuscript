## Analysis for the BugBase Manuscript

### Requirements:
* BugBase installed and able to run test dataset
* This git repository 
<br><br>

### Data Files:
#### HMP Body Sites
HMP OTU table was generated using QIIME v.1.8.0 (closed reference) with the GreenGenes 97% representative set v.13.5.8 and was filtered to keep only samples with at least 500 counts and samples from the tongue dorsum, stool and plaque
##### Reference: 
Human Microbiome Project Consortium. Structure, function and diversity of the healthy human microbiome. Nature 486, 207–214 (2012).
##### Files:
* HMP_Map.txt
* HMP_OTUS.biom
<br><br>

#### Western vs. Non-western Stool Samples
The OTU table for the western, non-western stool analysis was pulled from the Qitta websiteand was filtered to keep samples from only non-infants (>3 years of age) and samples with a minimum of 500 sequence counts
##### Reference:
Yatsunenko, T. et al. Human gut microbiome viewed across age and geography. Nature 486, 222–227 (2012).
##### Files:
* Yats_Adult_Map.txt
* Yats_Adult_OTUs.biom
<br><br>

#### Yellowstone Hotsprings Data
The OTU table for the Yellowstone hotspring analysis was pulled from the Qitta website and was filtered to keep only water samples and samples with a minimum of 500 sequence counts
##### Reference:
This data is not published elsewhere (yet!)
##### Files:
* Yellow_Stone_Map.txt
* Yellow_Stone_OTUs.biom
<br><br>

#### Soil pH
The OTU table for the soil analysis was pulled from the Qitta website and was filtered to keep samples with a minimum of 500 sequence counts
##### Reference:
Bartram, A. K. et al. Exploring links between pH and bacterial community composition in soils from the Craibstone Experimental Farm. FEMS Microbiol. Ecol. 87, 403–415 (2014).
##### Files:
* Soil_Map.txt
* Soil_OTUs.biom
<br><br>

#### Vaginal Samples
The vagina OTU table was generated using NINJA-OPS (closed reference) with the GreenGenes 97% representative set v13.5.8 all samples had a minimum of 300 sequences
##### Reference:
Ravel, J. et al. Vaginal microbiome of reproductive-age women. Proc. Natl. Acad. Sci. 108, 4680–4687 (2011).
##### Files:
* Vagina_Map.txt
* Vagina_OTUs.biom
<br><br>

### Running the Analysis:

To complete the BugBase analysis for all the work presented in the publication, you can run the commands listed in the `BugBase_Workflow.txt` file.

These commands should be run from within the BugBase_Manuscript directory.

