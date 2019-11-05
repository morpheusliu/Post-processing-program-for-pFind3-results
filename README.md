# Post-processing-program-for-pFind3-results
For providing the site-centric quantification data, the output reports from pFind® 3 need to be further processed (i) to remove peptides matched to decoy protein hits and common contaminates; (ii) to remove peptides bearing more than one probe-derived modifications; (iii) to set a cut off value for sigma score to remove those quantification results interfered by coeluting isobaric signals during LC-MS/MS analysis; (iv) to assign the position number of cysteine site on protein sequences; (iv) to group the quantification ratios of those peptide matches (e.g., different charge states, normal and missed cleavage, with and without methionine oxidation) that can assign the same modified cysteine site, with the median value reported as the final ratio. To facilitate these processes, here we described an open-access algorithm based on pythonTM for automatic post-processing of the pFind® 3 searching results.
1. System requirements:
This algorithm was developed by python3 and tested on python 3.7.3. Windows® is the uniqie operating system for this algorithm.
2. Installation guide:
Python 3.7.3 or any version higher than 3.7.3 should be installed (https://www.python.org/). Please make sure that python has been added to the path when installed （please check the picture 'python installation' in folder 'Instructions'）. There is no need of installation of this algorithm. You can download this algotithm (two .py programs) in a new folder and run them follow to instructions.
3. Instructions: 
(1) Download and install python3;
(2) Download the two .py oprograms in a new folder;
(3) Creat two new subfolders, 'data' and 'data_protein', in chich folder .py programs in.
(4) Creat a new .txt file and nominate it 'database';
(5) Open the original database file in a .fasta format with a text editor, copy and paste all the contents into the text new created file ‘database’. Double click the program ‘construct_library.py’ to reconstruct the database. 
(6) Copy the file ‘pQuant.spectra.list’ from the ‘results’ folder in the folder of pfind result you have assigned and paste it into the subfolder ‘data’. Rename the file as ‘IPM[C]_XXX’.
(7) Copy the file ‘pFind.protein’ from the folder ‘results’ in the location setup in Step 75 and paste it into the subfolder ‘data_protein’. Rename the file as ‘IPM[C]_XXX_P’.
(8) Double click the program ‘pfind_post_processing_site_level.py’ and use the following command: 
-	‘Please enter the calculation style (lh for light/heavy & hl for heavy/light):’ Type ‘hl’ or ‘lh’ to proceed
-	‘Please enter the cut off value of interference score:’ Type ‘0.5’ to proceed. 
-	‘Please enter the number of modification in each spectrum (a for single modification & b for multi-modification):’. Type ‘a’ to proceed. 
 Press ‘Enter key’ and run the algorithm.
(9) Find the output file with the prefix ‘site_’ in the subfolder ‘results’ and import it into Excel.

There are some pictures for illustration in the folder 'instructions'.
