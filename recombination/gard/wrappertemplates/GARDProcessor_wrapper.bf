inputRedirect = {};

/* specify input alignment (full path) */
inputRedirect["01"]="FULL_PATH_TO_ALIGNMENT";

/* specify GARD splits file (name ends _splits, full path) */
inputRedirect["02"]="FULL_PATH_TO_GARD_SPLITS_OUTPUT";

/* note - output file will end up in the same directory as the input alignment and will be called alnFileName_multi.fit - might want to move it after running GARDProcessor*/

/* do the analysis */
ExecuteAFile (HYPHY_BASE_DIRECTORY + DIRECTORY_SEPARATOR + "GARDProcessor.bf", inputRedirect);
