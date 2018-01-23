## README ##
This repo contains scripts used in the study "(Probabilistic enrichment of phosphopeptides by their mass defect)[https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2547851/]" by Bruce et al. (2006) in Anal Chem 78(13):4374-82. The study shows how the mass of tryptic peptides can be used to estimate probabilistically whether the peptide is phosphorylated or not. 

The Perl directory contains the following perl scripts used to generate the in-silico generated tryptic peptides and their characterization:

* digest.pl: Trypsin digestion human proteins from IPI, and calculate their masses.
    Generates all possible phosphorylations (up to 3 sites per peptide), and up to 2 missed 
    trypsin cuts per peptide. Peptides with methioinine are calculated with all possible variants 
    that are not/partially/fully oxidized. Cysteines are assumed to be carbamidomethylated. 
    In addition to the exact mass for each peptide, 4 corrupted mass values are generated, 
    corresponding to 1, 5 ,20, 100 ppm measurement error. 
* masscorrupt.pl: Script to determine how many peptides are incorrectly placed
    above or below the p=0.9 line, according to their phosphorylation status, if their masses 
    have a gaussian error of 1, 5, 20 or 100 ppm.
* find_neighbors.pl: For each peptide, calculates number of peptides having the same mass within 
    a given mass error (1, 5, 20, 100 ppm)

The Web directory contains the cgi-bin scripts, html and image file for the Web site mentioned in the original publication. The in-silico data used in this study can be queried with an interface that allows the user to enter monoisotopic, MH+, masses and measurement errors for a group of peptides and obtain output that lists for each peptide the number of phosphorylated and unphosphorylated peptides within the given mass measurement error, and its probability of being phosphorylated. The user can drill down deeper into the output to determine the sequences of candidate peptides whose mass falls within the mass measurement error, along with annotation data regarding the source. 

