library(ape)
source("/home/jayoung/Rscripts/tree_rotate_all_nodes.R")

GARDsplitsFile <- "GARD_SPLITS_FILE"
GARDprocessedFile <- "GARD_PROCESSOR_OUTPUT_FILE"

titleToShow <- "THIS_TITLE"

## will try to reroot the tree using these as outgroups:
outgroups <- c("OUTGROUP_NAME_OR_NAMES")

## if "outgroups" is >1 species, and is not monophyletic, I will want a single outgroup
backupOutgroups <- "SINGLE_OUTGROUP"

scaleBarLength <- 0.1

################

#### look at the GARDProcessor output to get some summary stuff
GARDtextOut <- scan(GARDprocessedFile, what="character", sep="\n")

## get the delta AIC
deltaAIC <- GARDtextOut[grep("GARD vs the single tree/multiple partition model",GARDtextOut)+1]
deltaAIC <- round(as.numeric(gsub("Delta AIC = ","",deltaAIC)),1)

### get evidence ratio (I was previously calculating it from deltaAIC, but now it is given in the HyPhy output). evidenceRatio=exp((deltaAIC/2))  see http://en.wikipedia.org/wiki/Akaike_information_criterion https://github.com/veg/hyphy/issues/289#issuecomment-103138567 ~/malik_lab_shared/help/hyphy/GARDevidenceRatioAdviceFromSergei.pdf 
## evidenceRatio <- signif(exp((deltaAIC/2)),2)
evidenceRatio <- GARDtextOut[grep("GARD vs the single tree/multiple partition model",GARDtextOut)+2]
evidenceRatio <- as.numeric(gsub("Evidence ratio in favor of the GARD model = ","",evidenceRatio))
if (evidenceRatio>100) { evidenceRatio <- ">100" }

## get the breakpoint summary
topOfBreakpointTablePosition <- grep ("Breakpoint", GARDtextOut)[1]
significantBreakpointsPositions <-  grep ("At p = ", GARDtextOut)

### test for situations where there was no breakpoint table
breakpointWarnings <- "OK"
if ( (significantBreakpointsPositions[1]-1)==topOfBreakpointTablePosition) {
    breakpointWarnings <- "noBreakpoints"
} else {
    breakpointTable <- GARDtextOut[ (topOfBreakpointTablePosition+1): (significantBreakpointsPositions[1]-1)]
    breakpointTable <- strsplit( breakpointTable, "[ |]+")
    breakpointTable <- lapply(breakpointTable, as.numeric)
    breakpointTable <-  as.data.frame(t(as.data.frame(breakpointTable)))
    rownames(breakpointTable) <- NULL
    breakpointTable <- breakpointTable[,2:dim(breakpointTable)[2] ]
    colnames(breakpointTable) <- strsplit( GARDtextOut[topOfBreakpointTablePosition], " \\| ")[[1]]
    breakpointTable[,"signif"] <- ""
    breakpointTable[which(breakpointTable[,"LHS adjusted p"]<= 0.1 & 
                          breakpointTable[,"RHS adjusted p "]<=0.1 ),"signif"] <- "*"
    breakpointTable[which(breakpointTable[,"LHS adjusted p"]<= 0.05 & 
                          breakpointTable[,"RHS adjusted p "]<=0.05 ),"signif"] <- "**"
    breakpointTable[which(breakpointTable[,"LHS adjusted p"]<= 0.01 & 
                          breakpointTable[,"RHS adjusted p "]<=0.01 ),"signif"] <- "***"
}
### read in trees, convert to ape format
dat <- scan(GARDsplitsFile, what="character")
numTrees <- length(dat)/2
treeNames <- paste("seg",dat[ (2*1:numTrees-1) ], sep="")
trees <- dat[ (2*1:numTrees) ]
names(trees) <- treeNames
trees <- lapply( trees, function(x) {  paste(x, ";",sep="") } )
trees <- lapply( trees, function(x) {  read.tree(text=x) } )

### root using outgroup(s)
trees <- lapply( trees, function(x) {  
    outgroupsToUse <- outgroups
    if (! is.monophyletic(x, outgroupsToUse) ) {
        outgroupsToUse <- backupOutgroups
    }
    
    ## check whether the tree is ALREADY rooted with those species as outgroups
    commonAnc <- "undefined"
    if (length(outgroupsToUse)==2) {
        commonAnc <- mrca(x)[outgroupsToUse[1], outgroupsToUse[2]]
    }
    if (length(outgroupsToUse)==1) {
        ## get numeric ID of that species
        mySpecies <- which(x$tip.label==outgroupsToUse)
        commonAnc <- x$edge[ which(x$edge[,2]==mySpecies), 1]
    }
    needToReRoot <- "yes"
    if (commonAnc != "undefined") {
        mynumNodes <- sum( x$edge[,1] == commonAnc)
        if (mynumNodes == 3) { 
            needToReRoot <- "no" 
        }
    }
    if (needToReRoot == "yes") {
        outtree <- root(x, outgroup=outgroupsToUse, resolve.root=TRUE) 
    } else {
        outtree <- x
    }
    outtree
} )

### rotate, so they will be displayed upside down
treesRotated <- lapply( trees, function(x) {  rotate_all_nodes(x) } )

### I'm using treeDepths to figure out where to put the scale bar (at the right of each tree).  I don't know why this function works, but it seems like it does. https://stat.ethz.ch/pipermail/r-sig-phylo/2012-September/002267.html
treeDepths <- unlist(lapply( trees, function(x) { max(diag(vcv.phylo(x))) } ))

X11(height=6,width=10)

layout(matrix( c(rep(1,numTrees),(1+1:numTrees)), 2, numTrees, byrow=TRUE))
par(mar=c(5.1, 2.1, 2.1, 1.1), xpd=NA) ## xpd=NA is so that the scale bar is not hidden
#layout.show( (numTrees+1) )

### plot a diagram of the segments
bottomHeight <- 0
topHeight <- 0.2
segmentStarts <- as.numeric(sapply(strsplit(gsub("seg","",treeNames), "-"), "[[", 1))
segmentEnds <- as.numeric(sapply(strsplit(gsub("seg","",treeNames), "-"), "[[", 2))
plot( c(0,max(segmentEnds)), c(0,1), "n", xaxt="n", yaxt="n", bty="n", xlab="", ylab="")
for (i in 1:length(segmentStarts)) {
    polygon( c( segmentStarts[i], segmentEnds[i], segmentEnds[i], segmentStarts[i], segmentStarts[i]), c(bottomHeight, bottomHeight, topHeight, topHeight, bottomHeight), col=i, border=i )
}
## add breakpoint significance:
if (breakpointWarnings == "OK") {
    for(i in 1:dim(breakpointTable)[1]) {
        if( breakpointTable[i,"signif"] == "") {next}
        text( x=breakpointTable[i,"Breakpoint"],y=0.25, breakpointTable[i,"signif"], col=i, cex=2)
    }
}
## add title and deltaAIC
title(main=titleToShow, line=-8, cex.main=2)
title(paste("Delta AIC vs single tree/multiple partition model: ",deltaAIC, " (evidence ratio ",evidenceRatio,")", sep=""), line=-10, cex.main=1.25)

### plot the trees
for (i in 1:numTrees) { 
    plot( treesRotated[[i]], main=names(treesRotated)[i], font=1, cex=1, xpd=NA, edge.color=i, tip.color=i, col.main=i )
    add.scale.bar(x=0, y=0, length=scaleBarLength, lcol=i, col=i)
}

dev.print(pdf, file=paste(GARDsplitsFile,".plots.pdf", sep=""))
