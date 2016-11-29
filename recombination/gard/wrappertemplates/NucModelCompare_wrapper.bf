inputRedirect = {};

/* Model options:  [Global] Model parameters are shared by all branches, branch lengths are estimated independently */
/*     other choices: Local - Global w/variation - Global w/variation+HM */
inputRedirect["01"]="Global";

/* Estimate Branch Lengths: [Once] Branch lengths obtained from the general reversible model are reused for subsequent models. */
/*     other choice: [Every Time] Branch lengths are reestimated for every model */
inputRedirect["02"]="Once";

/* specify input alignment (full path) */
inputRedirect["03"]="FULL_PATH_TO_ALIGNMENT";

/* specify input tree (full path) */
inputRedirect["04"]="FULL_PATH_TO_TREE_WITHOUT_BRANCH_LENGTHS";

/* model rejection level (0.05 is suggested) */
inputRedirect["05"]="0.05";

/* save each of the 203 fits to a separate file? */
inputRedirect["06"]="No";

/* specify what output file will be called (give full path, otherwise it'll end up in /home/jayoung/malik_lab_shared/linux_gizmo/lib/hyphy/TemplateBatchFiles) */
inputRedirect["07"]="NUCMODEL_OUTPUT_FILENAME_FULL_PATH";

/* do the analysis */
ExecuteAFile (HYPHY_BASE_DIRECTORY + DIRECTORY_SEPARATOR + "NucModelCompare.bf", inputRedirect);
