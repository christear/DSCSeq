#' @rdname DSCSeqData
#' @export

setClass("DSCSeqData",representation(splicingin="data.frame",splicingout="data.frame",condition="vector",expression="data.frame",psipir = "data.frame",me = "vector",sc="vector",zscore="vector",pval="vector",qval="vector"))

DSCSeqData=function(splicingin,splicingout,condition){
    psuodeC = 1
    if(sum(dim(splicingin) == dim(splicingout)) != 2){
        stop("Dimensions of splicing in and out matrix are not equal...\n")
    }
    
    if((any(round(splicingin) != splicingin)) || (any(round(splicingout) != splicingout))){
        stop("The count matrix is not integer...\n")
    }
    
    dsccd = new("DSCSeqData",splicingin = round(splicingin),splicingout = round(splicingout),condition = condition)
    cat("calculating the expression and PSI/PIR value...\n")
    dsccd@expression = log2(splicingin + splicingout + psuodeC)
    dsccd@expression = dsccd@expression[order(rownames(dsccd@expression)),]
    dsccd@psipir = splicingin/(splicingin + splicingout)
    dsccd@psipir = dsccd@psipir[order(rownames(dsccd@psipir)),]
    dsccd@condition = as.factor(condition)
    dsccd@splicingin = splicingin[order(rownames(splicingin)),]
    dsccd@splicingout = splicingout[order(rownames(splicingout)),]
    stopifnot(length(condition) == ncol(splicingin))
    return(dsccd)
}


### fill the NA with zero ....
fillNA = function(x){
    x[is.na(x)] = 0
    x
}

### calculate SD with NA value
calSd = function(x){
    psd = sd(x[!(is.na(x))])
    psd
}

## calculate the percetage of splicing in between differnt conditions
## equalize the difference of PSI value zscore according to sum of splicing in and out reads
# it will calculate difference between condA and condB (condA - condB)

equalization = function(x,y,binNums,method = "step_e"){
    psuodeC = 0.01
    start = min(x)
    end = max(x) + psuodeC
    if(method == "step_e"){
        binSeq = seq(from = start,to = end,by = (end - start)/binNums)
        z = y[x == start]/(calSd(y[x == start]))
        for(i in 1:binNums){
            binX = x[(x > binSeq[i]) & (x <= binSeq[i + 1])]
            binY = y[(x > binSeq[i]) & (x <= binSeq[i + 1])]
            zscore = binY/(calSd(binY))
            names(zscore) = names(binY)
            #zscore[is.na(zscore)] = 0
            z = c(z,zscore)
        }
        z = z[order(names(z))]
        z[is.na(z)] = 0
    }else if(method == "step_n"){
        sorted_x = x[order(x)]
        sorted_y = y[order(x)]
        binSeq = round(seq(1,length(x),length(x)/binNums))
        z = vector()
        for(i in 1:(length(binSeq) - 1)){
            binX = sorted_x[binSeq[i] : (binSeq[i + 1] - 1)]
            binY = sorted_y[binSeq[i] : (binSeq[i + 1] - 1)]
            zscore = binY/(sd(binY))
            names(zscore) = names(binY)
            z = c(z,zscore)
        }
        lastX = sorted_x[binSeq[length(binSeq)] : length(x)]
        lastY = sorted_y[binSeq[length(binSeq)] : length(y)]
        lastZ = lastY/(sd(lastY))
        z = c(z,lastZ)
        z = z[order(names(z))]
        z[is.na(z)] = 0
    }else{
        cat("equlize PSI method: ",method," has not been defined.\n")
    }
    z
}

equalizeToZscore <- function(psiDiff,sioSums,binNums = 50,method = "step_e"){
    #cat(binNums,"\n")
    if(is.null(dim(psiDiff))){
        cat("There is only one comparison,it won't be examined splicing change using rank product test...\n")
        x = sioSums
        y = psiDiff
        z = equalization(x,y,binNums,method)
    }else{
        allZ = vector()
        for(i in 1:ncol(psiDiff)){
            x = sioSums[,i]
            y = psiDiff[,i]
            names(x) = names(y) = rownames(sioSums)
            z = equalization(x,y,binNums,method)
            allZ = cbind(allZ,z)
        }
    }
    colnames(allZ) = colnames(psiDiff)
    allZ = allZ[order(rownames(allZ)), ]
    allZ
}

### convert the PSI/PIR value to z score
setGeneric("equalizePSI",function(object,method = NULL,condA = NULL,condB = NULL,paired = FALSE,binNums = 50) standardGeneric("equalizePSI"))
#' export
setMethod("equalizePSI","DSCSeqData",function(object,method = "step_e",condA = NULL,condB=NULL,paired = FALSE,binNums = 50){
    cat("converting the PSI/PIR value to z score...\n")
    psi = object@psipir
    sio = object@expression
    conditions = object@condition
    stopifnot(condA %in% levels(conditions))
    stopifnot(condB %in% levels(conditions))
    colA <- conditions == condA
    colB <- conditions == condB
    psiA = psi[,colA]
    psiB = psi[,colB]
    sioA = sio[,colA]
    sioB = sio[,colB]
    
    if(paired){
        cat("paired is true......\n")
        stopifnot(ncol(psiA) == ncol(psiB))
        object@sc = psiA - psiB
        object@me = sioA + sioB
        colnames(object@sc) = colnames(object@me) = paste(condA,"-",condB,1:ncol(psiA),sep = "")
    }else{
        if(sum(colA) == 1 | sum(colB) == 1){
            #cat("There's only one ",condA," or one",condB,"sample...\n")
            object@sc = psiA - psiB
            object@me = sioA + sioB
            if(sum(colB) == 1){
                cat("only one ",condB," sample...\n")
            }else{
                colnames(object@sc) = colnames(object@me) = paste(condA,"-",condB,1:ncol(psiB),sep = "")
            }
            if(sum(colA) == 1){
                cat("only one ",condA," sample...\n")
            }else{
                colnames(object@sc) = colnames(object@me) = paste(condA,1:ncol(psiA),"-",condB,sep = "")
            }
        }else{
            
            cat("merge",condB," as one group...\n")
            diff = psiA - apply(psiB,1,mean)
            sums = sioA + apply(sioB,1,mean)
            colnames(diff) = colnames(sums) = paste(paste(condA,1:ncol(psiA),sep = ""),"-",condB,"_mean",sep = "")
            object@sc = diff
            object@me = sums
        }
    }
    ### exclude the splicing events without supporting reads ... in any one of condition ....
    #NaNmarker = apply(object@sc,1,mean)
    #object@sc = object@sc[!(is.na(NaNmarker)),]
    #object@me = object@me[!(is.na(NaNmarker)),]
    object@zscore = equalizeToZscore(object@sc,object@me,binNums = binNums,method = method)
    ###
    return(object)
}
)

### applying rank product test ...
setGeneric("rankProdTest",function(object,permutime = NULL,method = NULL) standardGeneric("rankProdTest"))
#' export
setMethod("rankProdTest","DSCSeqData",function(object,permutime = 100,method = "rankproduct"){
    require(RankProd)
    cl = rep(1,ncol(object@zscore))
    #marker = apply(object@zscore,1,sum)
    #upGroup = object@zscore[marker > 0,]
    ###upGroup = object@zscore[marker >= 0 ,]
    #dnGroup = object@zscore[marker < 0,]
    #zeroGroup = object@zscore[marker == 0,]
    
    upmarker = apply(object@zscore,1,min)
    dnmarker = apply(object@zscore,1,max)
    
    #ambiGroup = object@zscore[upmarker < 0 & dnmarker > 0,]
    ambiGroup = object@zscore[upmarker * dnmarker < 0 | apply(object@zscore,1,median) == 0,]
    
    pval = vector()
    qval = vector()
    if(method == "rankproduct"){
        #rp.out = lapply(list(upGroup,-dnGroup),function(x) RP(x,cl,num.perm = permutime))
        rp.out = RP(object@zscore,cl,num.perm = permutime)
    }else if (method == "ranksum"){
        #rp.out = lapply(list(upGroup,-dnGroup),function(x) RSadvance(x,cl,1:ncol(zeroGroup),num.perm = permutime))
        rp.out = RSadvance(object@zscore,cl,1:ncol(object@zscore),num.perm = permutime)
    }else{
        cat("the rankprod test method: ", method, " has not been defined.\n")
    }

#uppval = rp.out[[1]]$pval[,2]
#    upqval = rp.out[[1]]$pfp[,2]
#names(uppval) = names(upqval) = rownames(upGroup)
    
    #    dnpval = rp.out[[2]]$pval[,2]
    #dnqval = rp.out[[2]]$pfp[,2]
    #   names(dnpval) = names(dnqval) = rownames(dnGroup)
    
    #zeropval = zeroqval = rep(1,nrow(zeroGroup))
    #names(zeropval) = names(zeroqval) = rownames(zeroGroup)

#pval = c(uppval[names(uppval) %in% setdiff(names(uppval),rownames(zeroGroup))],dnpval[names(dnpval) %in% setdiff(names(dnpval),rownames(zeroGroup))],zeropval)
#   qval = c(upqval[names(upqval) %in% setdiff(names(upqval),rownames(zeroGroup))],dnqval[names(dnqval) %in% setdiff(names(dnqval),rownames(zeroGroup))],zeroqval)

 pval = apply(rp.out$pval,1,min)
 qval = apply(rp.out$pfp,1,min)
    ### set the splicing change in different direction events p and q value to 1
    names(pval) = names(qval) = rownames(object@zscore)
    pval[names(pval) %in% rownames(ambiGroup)] = 1
    qval[names(qval) %in% rownames(ambiGroup)] = 1


    ###
    #  pval = pval[order(names(pval))]
    #qval = qval[order(names(qval))]
    qval[qval>1] = 1
    object@pval = pval
    object@qval = qval
    return(object)
}
)



### output the results
setGeneric("rpresults",function(object,permutime = NULL,method = NULL) standardGeneric("rpresults"))
setMethod("rpresults","DSCSeqData",function(object,method = "mean"){
    if(method == "mean"){
        rpr = cbind(apply(object@me,1,function(x) mean(x[!(is.na(x))])),apply(object@sc,1,function(x) mean(x[!(is.na(x))])),apply(object@zscore,1,function(x) mean(x[!(is.na(x))])),object@qval,object@pval)
    }else if(method == "median"){
        rpr = cbind(apply(object@me,1,function(x) median(x[!(is.na(x))])),apply(object@sc,1,function(x) median(x[!(is.na(x))])),apply(object@zscore,1,function(x) median(x[!(is.na(x))])),object@qval,object@pval)
    }
    colnames(rpr) = c("meanExpression","splicingchange","zscore","qval","pval")
    rpr
}
)

### plot the results
setGeneric("scplot",function(object,mecut = NULL,qcut = NULL,pcut = NULL,sccut = NULL,type = NULL,main = NULL,xlab = NULL,ylab = NULL,xlim = NULL,ylim = NULL,pch = NULL,cex = NULL,method = NULL) standardGeneric("scplot"))

setMethod("scplot","DSCSeqData",function(object,mecut = log2(400),qcut = 0.01,pcut = 0.001,sccut = 0,type = "density",main = NULL,xlab = "log2 Mean expression",ylab = "Delta PSI/PIR",xlim = NULL,ylim = NULL,pch = 19,cex = .8,method = "mean"){
    rpr = rpresults(object,method = method)
    subrpr = rpr[rpr[,1] >= mecut & (!(is.na(rpr[,2]))),]
    sigrpr = subrpr[subrpr[,4] < qcut & subrpr[,5] < pcut,]
    if(is.null(xlim)){
        xlim = c(min(subrpr[,1]),max(subrpr[,1] + 1))
    }
    if(is.null(ylim)){
        ylim = c(-max(abs(subrpr[,2])) - 0.1,max(abs(subrpr[,2])) + 0.1)
    }
    if(type == "density"){
        smoothScatter(subrpr[,1],subrpr[,2],main = main,xlab = xlab,ylab = ylab,xlim = xlim,ylim = ylim,col = "grey",pch = pch)
        #smoothScatter(sigrpr[,1],sigrpr[,2],col = "red",add = T)
    }else if(type == "point"){
        plot(subrpr[,1],subrpr[,2],main = main,xlab = xlab,ylab = ylab,xlim = xlim,ylim = ylim,col = "grey",pch = pch,cex = cex)
    }else{
        cat("the plot method: ",method," has not been defined...\n")
    }
    if(!(is.null(dim(sigrpr)))){
        points(sigrpr[,1],sigrpr[,2],col = "red",pch = pch,cex = cex)
    }
    abline(h = 0,lty = 3)
}
)

