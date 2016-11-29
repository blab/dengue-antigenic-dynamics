inputRedirect = {};

/* specify input alignment (full path) */
inputRedirect["01"]="/home/sbell23/gard_denv2/denv2_gard.fasta";

/* specify nucleotide model (perhaps run NucModelCompare, or use GTR=012345 */
inputRedirect["02"]="012343";

/* site variation?  General Discrete is recommended, but other options are None (fastest) or Beta-Gamma */
inputRedirect["03"]="General Discrete";

/* How many distribution bins [2-32]? (unless None was chosen for site variation, in which case do not use this input) */
inputRedirect["04"]="3";

/* specify what output files will be called - there will be four outputs, with this stem in the name (give full path, otherwise it'll end up in /home/jayoung/malik_lab_shared/linux_gizmo/lib/hyphy/TemplateBatchFiles) */
inputRedirect["05"]="/home/sbell23/gard_denv2/denv2_gard.fasta.GARD";

/* do the analysis */
ExecuteAFile (HYPHY_BASE_DIRECTORY + DIRECTORY_SEPARATOR+"GARD.bf", inputRedirect);
