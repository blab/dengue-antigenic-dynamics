#!/usr/bin/perl

use warnings;
use strict;
use Cwd;

####  given an input seq file and tree file, preps wrapper scripts for NucModelCompare.bf, GARD.bf, GARProcessor.bf

## give paths relative to current dir - if treeFile is "none" then we don't try to run NucModelCompare.bf

my $alnFile = "MxA_align.fasta";
#my $treeFile = "MxA_tree.nex";

my $treeFile = "none";

my $outgroupNames = "kangaroo_rat_MxA";
#my $outgroupNames = "zebra_finch";

my $numThreads = 12;

## specify templates
my $templateDir = "/home/sbell23/gard/wrappertemplates";
my $NucModelCompareWrapperTemplate = "NucModelCompare_wrapper.bf";
my $GARDwrapperTemplate = "GARD_wrapper.bf";
my $GARDProcessorWrapperTemplate = "GARDProcessor_wrapper.bf";
my $RscriptTemplate = "processGARDoutput.R";

#######################

## everything will go in the current working directory - double-check it's not the template dir
my $cwd = cwd();
if ($cwd =~ m/template_wrapper_scripts_JY/) {
    die "\n\nterminating - cannot run this script in the main templates directory\n\n";
}

## check template files exist
if (!-e "$templateDir/$NucModelCompareWrapperTemplate") {
    die "\n\nterminating - cannot find NucModelCompareWrapperTemplate $templateDir/$NucModelCompareWrapperTemplate\n\n";
}
if (!-e "$templateDir/$GARDwrapperTemplate") {
    die "\n\nterminating - cannot find GARDwrapperTemplate $templateDir/$GARDwrapperTemplate\n\n";
}
if (!-e "$templateDir/$GARDProcessorWrapperTemplate") {
    die "\n\nterminating - cannot find GARDProcessorWrapperTemplate $templateDir/$GARDProcessorWrapperTemplate\n\n";
}

## get full paths
my $fullAlnFile = "$cwd/$alnFile";
my $fullTreeFile = "$cwd/$treeFile";


## check those exist
if (!-e $fullAlnFile) {
    die "\n\nterminating - cannot find alignment file $fullAlnFile\n\n";
}
if ($treeFile ne "none") {
    if (!-e $fullTreeFile) {
        die "\n\nterminating - cannot find tree file $fullTreeFile\n\n";
    }
}

## I will also print suggested commands to run the wrapper scripts we create here:
open (COMMANDS, "> suggestedCommands.txt");

## open templates, print output wrapper scripts
# NucModelCompareWrapperTemplate
if ($treeFile ne "none") {
    my $nucModelOutFile = "$cwd/$alnFile.NucModelCompareOut";
    my $nucModelOutFile2 = "$cwd/$alnFile.NucModelCompareOut.BFout";
    open (OUT, "> this$NucModelCompareWrapperTemplate");
    open (IN, "< $templateDir/$NucModelCompareWrapperTemplate");
    while (<IN>) {
        my $line = $_;
        $line =~ s/FULL_PATH_TO_ALIGNMENT/$fullAlnFile/;
        $line =~ s/FULL_PATH_TO_TREE_WITHOUT_BRANCH_LENGTHS/$fullTreeFile/;
        $line =~ s/NUCMODEL_OUTPUT_FILENAME_FULL_PATH/$nucModelOutFile/;
        print OUT "$line";
    }
    close IN;
    close OUT;
    print COMMANDS "\n#####  NucModelCompare command ######\n\n";
    print COMMANDS "sbatch --exclude=/home/sbell23/gizmod.txt --time=20-0 --ntasks 12 --cpus-per-task 1 --wrap=\"mpirun -np $numThreads /home/sbell23/build/bin/hyphy BASEPATH=/home/sbell23/build/lib/hyphy/TemplateBatchFiles  $cwd/this$NucModelCompareWrapperTemplate > $nucModelOutFile\"2\n\n";
    print COMMANDS "    Then look in $nucModelOutFile to see what the best model is (AIC based winner). Edit the GARD.bf control file to use that model\n";
}

# GARDwrapperTemplate
open (OUT, "> this$GARDwrapperTemplate");
open (IN, "< $templateDir/$GARDwrapperTemplate");
my $gardOutputFileStem = "$cwd/$alnFile.GARD";
while (<IN>) {
    my $line = $_;
    $line =~ s/FULL_PATH_TO_ALIGNMENT/$fullAlnFile/;
    $line =~ s/OUTPUT_FILESTEM_FULL_PATH/$gardOutputFileStem/;
    print OUT "$line";
}
close IN;
close OUT;
print COMMANDS "\n#####  GARD command ######\n\n";
print COMMANDS "sbatch --exclude=/home/sbell23/gizmod.txt --time=20-0 --ntasks 12 --cpus-per-task 1 --wrap=\"mpirun -np $numThreads /home/sbell23/build/bin/hyphy BASEPATH=/home/sbell23/build/lib/hyphy/TemplateBatchFiles  $cwd/this$GARDwrapperTemplate > $gardOutputFileStem.BFout\"\n\n";
print COMMANDS "mv $alnFile.GARD $alnFile.GARD.html\n\n";

# GARDProcessorWrapperTemplate
open (OUT, "> this$GARDProcessorWrapperTemplate");
open (IN, "< $templateDir/$GARDProcessorWrapperTemplate");
my $gardOutputSplitFile = "$gardOutputFileStem"."_splits";
while (<IN>) {
    my $line = $_;
    $line =~ s/FULL_PATH_TO_ALIGNMENT/$fullAlnFile/;
    $line =~ s/FULL_PATH_TO_GARD_SPLITS_OUTPUT/$gardOutputSplitFile/;
    print OUT "$line";
}
close IN;
close OUT;
print COMMANDS "\n#####  GARDProcessor command ######\n\n";
print COMMANDS "sbatch --exclude=/home/sbell23/gizmod.txt --time=20-0 --ntasks 12 --cpus-per-task 1 --wrap=\"mpirun -np $numThreads /home/sbell23/build/bin/hyphy BASEPATH=/home/sbell23/build/lib/hyphy/TemplateBatchFiles  $cwd/this$GARDProcessorWrapperTemplate > $gardOutputFileStem"."ProcessorBFout\"\n\n";

# the R script
my $gardOutputSplitFileShort = $gardOutputSplitFile;
if ($gardOutputSplitFileShort =~ m/\//) {
    $gardOutputSplitFileShort = (split /\//, $gardOutputSplitFileShort)[-1];
}
my $gardProcOutFile = $gardOutputSplitFileShort;
$gardProcOutFile =~ s/_splits//;
$gardProcOutFile .= "ProcessorBFout";
open (OUT, "> this$RscriptTemplate");
open (IN, "< $templateDir/$RscriptTemplate");
while (<IN>) {
    my $line = $_;
    $line =~ s/GARD_SPLITS_FILE/$gardOutputSplitFileShort/;
    $line =~ s/OUTGROUP_NAME_OR_NAMES/$outgroupNames/;
    $line =~ s/GARD_PROCESSOR_OUTPUT_FILE/$gardProcOutFile/;
    print OUT "$line";
}
close IN;
close OUT;
print COMMANDS "\n#####  Rscript ######\n\n";
print COMMANDS "Now mess with this$RscriptTemplate to make some plots";

close COMMANDS;

print "\n\nDone - made three wrapper files (thisGARDProcessor_wrapper.bf, thisGARD_wrapper.bf and thisNucModelCompare_wrapper.bf) and an R script (thisprocessGARDoutput.R).\n\nCheck the file suggestedCommands.txt for how to run these wrappers\n\n";
