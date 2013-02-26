###############################################
# It's a script to analyze splicing change in count data from high-throughput sequencing based on a rank product test.
# It is used to identify differential splicing events in exon/intron level, In particular, it will indentify up and down
# regulated exon/intron under one condition against another condition. eg case and control , patients and normal , two
# different tissue or female and male ...
# ZHANG bin -2013 - 02 - 18
###############################################

## creat percatage of splicing in object...
# countData1 should be the counts of splicing in reads while
# countData2 should be the counts of splicing out reads 
# conditions indicate the case or control of each sample
# rownames of the countData1 and countData should be the same (exon name or intron name)
newPSIData <- function(countData1,countData2,conditions){
    if(sum(dim(countData1) == dim(countData2)) != 2){
        stop("Dimensions of two countData are not equal.\n")
    }
    countData1 = as.matrix(countData1)
    countData2 = as.matrix(countData2)
    if((any(round(countData1) != countData1)) || (any(round(countData2) != countData2))){
        stop("The countData is not integer.\n")
    }
    sio = countData1 + countData2
    psi = countData1/sio
    psi[is.na(psi)] = 1
    sio = log2(sio + 1)
    conditions = factor(conditions)
    stopifnot(length(conditions) == ncol(countData1))
    PSIData = list(sio,psi,conditions)
    names(PSIData) = c('percentage_of_splicing_in','sum_splicing_in_and_out_reads_counts','conditions')
    PSIData
}



## calculate the percetage of splicing in between differnt conditions
## equalize the difference of PSI value zscore according to sum of splicing in and out reads
# it will calculate difference between condA and condB (condA - condB)

equalization = function(x,y,binNums){
    start = min(x)
    end = max(x)
    start = min(x)
    end = max(x)
    binSeq = seq(from = start,to = end,by = (end - start)/binNums)
    z = y[x == start]/(sd(y[x == start]))
    for(i in 1:binNums){
        binX = x[(x > binSeq[i]) & (x <= binSeq[i + 1])]
        binY = y[(x > binSeq[i]) & (x <= binSeq[i + 1])]
        zscore = binY/(sd(binY))
        names(zscore) = names(binY)
        #zscore[is.na(zscore)] = 0
        z = c(z,zscore)
    }
    z = z[order(names(z))]
    z[is.na(z)] = 0
    z
}

equalizeToZscore <- function(psiDiff,sioSums,binNums = 50){
    #cat(binNums,"\n")
    if(is.null(dim(psiDiff))){
        cat("There is only one comparison,it won't be examined splicing change using rank product test\n")
        x = sioSums
        y = psiDiff
        z = equalization(x,y,binNums)
    }else{
        allZ = vector()
        for(i in 1:ncol(psiDiff)){
            z = equalization(sioSums[,i],psiDiff[,i],binNums)
            allZ = cbind(allZ,z)
        }
    }
    colnames(allZ) = colnames(psiDiff)
    allZ = allZ[order(rownames(allZ)), ]
    allZ
}


equalizePSIChange <- function(PSIData,condA,condB) {
    conditions = PSIData$conditions
    psi = PSIData[[2]]
    sio = PSIData[[1]]
    stopifnot(condA %in% levels(conditions))
    stopifnot(condB %in% levels(conditions))
    colA <- conditions == condA
    colB <- conditions == condB
    psiA = psi[,colA]
    psiB = psi[,colB]
    sioA = sio[,colA]
    sioB = sio[,colB]
    sioSums = vector()
    psiDiff = vector()
    if(sum(colA) == 1){
        cat("There's only one ",condA,"sample\n")
        psiDiff = psiA - psiB
        sioSums = sioA + sioB
        colnames(psiDiff) = paste(condA,"-",condB,1:ncol(psiB),sep = "")
        colnames(sioSums) = paste(condA,"-",condB,1:ncol(psiB),sep = "")
    }else{
        for(i in 1:ncol(psiA)){
            psiA[,i] - psiB -> diff
            sioA[,i] + sioB -> sums
            colnames(diff) = paste(condA,i,"-",condB,1:ncol(psiB),sep = "")
            colnames(sums) = paste(condA,i,"-",condB,1:ncol(psiB),sep = "")
            psiDiff = cbind(psiDiff,diff)
            sioSums = cbind(sioSums,sums)
        }
    }
    psiDiff = psiDiff[order(rownames(psiDiff)),]
    sioSums = sioSums[order(rownames(sioSums)),]
    #equalizedPSIChange = equalizeToZscore(psiDiff,sioSums)
    equalizedPSIChange = list(sioSums,psiDiff,equalizeToZscore(psiDiff,sioSums))
    names(equalizedPSIChange) = c("Sum_Counts","PSI_Change","equalized_PSI_Change")
    equalizedPSIChange
}

### rank product test
# implent rank product test to examine splicing change
# the splicing change value should be equalized before using rank test
# the R package "RankProd" is required

# the last 2 columns of the matrix should be the rank product test "pfp" and "pvlaue"
getSig = function(mat,cutoff,method){
    n = ncol(mat)
    if(method == "pfp"){
        sig = mat[mat[,n - 1] < cutoff,]
    }else if(method == "pvalue"){
        sig = mat[mat[,n] < cutoff,]
    }
    sig
}

rankProdTest <- function(equalizedPSIChange,times = 100,cutoff = 0.1,method = "pfp"){
#rankProdTest <- function(equalizedPSIChange,times = 100){
    require(RankProd)
    cl = rep( 1 , ncol(equalizedPSIChange[[3]]))
    marker = apply( equalizedPSIChange[[3]], 1 , sum )
    upGroup = equalizedPSIChange[[3]][ marker >= 0 ,]
    downGroup = equalizedPSIChange[[3]][ marker < 0 ,]
    cat("rank product test for up regulate ones\n")
    rp.out = RP( upGroup ,cl, num.perm = times )
    upGroup = cbind(upGroup,rp.out$pfp[,2],rp.out$pval[,2])
    #top.gene.up = topGene(rp.out,cutoff = cutoff,method = "pfp")
    #if(method == "pfp"){
    #   top.gene.up = upGroup[rp.out$pfp[,2] < cutoff,]
    #}else if(method = "pvalue"){
    #   top.gene.up = upGroup[rp.out$pfp[,1] < cutoff,]
    #}
    #reverseDown = -downGroup
    cat("rank product test for down regulate ones\n")
    rp.out = RP(-downGroup,cl,num.perm = times)
    downGroup = cbind(downGroup,rp.out$pfp[,2],rp.out$pval[,1])
    #top.gene.down = topGene(rp.out,cutoff = cutoff,method = "pfp")
    #if(method = "pfp"){
    #    top.gene.down = downGroup[rp.out$pfp[,2] < cutoff,]
    #}else if(method == "pvalue"){
    #   top.gene.down = downGroup[rp.out$pfp[,1] < cutoff,]
    #}
    union = rbind(upGroup,downGroup)
    union = union[order(rownames(union)),]
    equalizedPSIChange[[1]] = cbind(equalizedPSIChange[[1]],union[,(ncol(union) - 1) : (ncol(union))])
    equalizedPSIChange[[2]] = cbind(equalizedPSIChange[[2]],union[,(ncol(union) - 1) : (ncol(union))])
    rankProdTestResult = list(equalizedPSIChange[[1]],equalizedPSIChange[[2]],union)
    #rankProdTestResult = lapply(rankProdTestResult,function(x) getSig(x,cutoff,method))
    names(rankProdTestResult) = c("Sum_Counts","PSI_Change","equalized_PSI_Change_pfp_pvalue")
    sig = lapply(rankProdTestResult,function(x) getSig(x,cutoff,method))
    rankProdTestResult = list(rankProdTestResult,sig)
    #rankProdTestResult = list(union,top.gene.up,top.gene.down)
    #names(rankProdTestResult) = c("flase_potive_percatage_and_pvalue","up_regulate","down_regulate")
    rankProdTestResult
}



### graphic

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    #require(ggplot2)
    require(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
        ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
        print(plots[[1]])
        
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
            layout.pos.col = matchidx$col))
        }
    }
}


# MA plot
plotMA = function(rankProdTestResult,title = "Case-Control",baseline = TRUE,sig = TRUE,point.col = "black",combine = FALSE,sig.col = "purple",ggplot = TRUE){
    Sum_Counts = rankProdTestResult[[1]]
    PSI_Change = rankProdTestResult[[2]]
    sig = "00"
    #sigUp = rankProdTestResult[[2]]
    #sigDown = rankProdTestResult[[3]]
    #Samples = factor(rep(colnames(Sum_Counts) , each = nrow(Sum_Counts)),levels = colnames(Sum_Counts))
    if(combine){
        Sum_Counts = apply(Sum_Counts,1,median)
        PSI_Change = apply(PSI_Change,1,median)
    }
    if(ggplot){
        cat("it will implement ggplot in plots\n")
        require(ggplot2)
        if(is.null(dim(Sum_Counts))){
            allPoints = data.frame(PSI_Change = PSI_Change, Sum_Counts = Sum_Counts)
            ggplot(data = allPoints,aes(x = Sum_Counts,y = PSI_Change)) + geom_point(point.col,size = .5) + labs(title = title) -> plots
            if(sig){
                sigPoints = data.frame(PSI_Change)
                geom_point()
            }
        }else{
            plotL = list()
            for(i in 1:ncol(Sum_Counts)){
                allPoints = data.frame(PSI_Change = PSI_Change[,i],Sum_Counts = Sum_Counts[,i])
                ggplot(data = allPoints,aes(x = Sum_Counts,y = PSI_Change)) + geom_point(point.col,size = .5) + labs(title = colnames(Sum_Counts)[i]) -> plotL[[i]]
            }
            multiplot(plotlist = plotL,cols = round(sqrt(ncol(Sum_Counts)))) -> plots
        }
    }else{
        
    }
}


