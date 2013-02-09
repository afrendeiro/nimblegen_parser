


            


## parse the layout information from a pairs file
## assumed to be buried in the comment line
get.layout.info<-function(file) {
    con<-file(file, "r", blocking = FALSE)
    h<-readLines(con, n=1)
    gregexpr("designname=", h)[[1]]->s
    substr(h, s+11, 1000)->l
    gregexpr("\t", l)[[1]]->s
    substr(l, 1, s-1)->l
    close(con)
    l
}



## read information from pos and ndf file
## remove all non-experimental information
read.layout <- function(layout.name, verbose=T, prefix="") {
    file <- paste(prefix, layout.name, ".ndf", sep="")
    if (verbose) {
        cat(paste("reading", file,"\n"))
    }
    li <- .read.nimble(file, c("PROBE_ID","CONTAINER","PROBE_CLASS","PROBE_SEQUENCE","X","Y"),
        c("character","character", "character", "character","integer","integer"), skip=1)
    nd<-data.frame("probe.id"=li[[1]], "opt"=li[[2]], "class"=li[[3]],
        "sequence"=li[[4]], "x"=li[[5]], "y"=li[[6]])
    
    file <- paste(prefix, layout.name, ".pos", sep="")
    if (verbose) {
        cat(paste("reading", file,"\n"))
    }
    li <- .read.nimble(file, c("PROBE_ID", "CHROMOSOME", "POSITION", "LENGTH"),
        c("character", "character","integer","integer"), skip=1)
    pos<-data.frame("probe.id"=li[[1]], "chr"=sub("chr", "", li[[2]],
        ignore.case=T),
        "pos"=li[[3]], "length"=li[[4]])
    if (verbose) {
        cat("merging.. \n")
    }
    l<-merge(nd, pos)
    o<-order(l$probe.id)
    l<-l[o,]
    l<-l[l$opt!="RANDOM",c(1,4:9)]
    if (verbose) {
        cat(paste(dim(l)[1], "layout rows kept\n"))
    }
    l
}

## write continuous track information into an sgr file
## y - the value to be written for each probe defined in layout
## layout - the layout table with chromosome (chr) and position (pos) information
write.sgr <- function(y, layout, filename="target.sgr") {
    options(scipen=100)
    ex=data.frame(paste("chr", layout$chr,sep=""), 
    round((layout$pos+(layout$pos + layout$length))/2,0),
    round(as.vector(y), 2))
    flt <- !is.na(ex[,2]) & !is.na(ex[,3])
    con <- file(filename, open="wt")
    writeLines(paste("##", date(), sep=""), con, sep="\n")
    write.table(ex[flt,], row.names=F, col.names=F, quote=F, sep="\t", con)
    close(con)
}

## write continuous track information into a wiggle file
## y - the value to be written for each probe defined in layout
## layout - the layout table with chromosome (chr), position (pos) and length information
write.wiggle <- function(y, layout, filename="target.wig", descr="target") {
    options(scipen=100)
    ex=data.frame(paste("chr", layout$chr,sep=""), 
    layout$pos, layout$pos+layout$length, round(as.vector(y), 2))
    flt <- !is.na(ex[,2]) & !is.na(ex[,4])
    con <- file(filename, open="wt")
    writeLines(paste("##", date(), sep=""), con, sep="\n")
    head <- paste("track type=wiggle_0 name=\"",descr,"\" description=\"",descr,
        "\" visibility=full color=200,100,0 altColor=0,100,200 priority=20", sep="")
    writeLines(head, sep="\n", con)
    write.table(ex[flt,], row.names=F, col.names=F, quote=F, sep="\t", con)
    close(con)
}

## write continuous track information into an gff file
## y - the value to be written for each probe defined in layout
## layout - the layout table with chromosome (chr) and position (pos) information
write.gff.signals<-function(y, layout, filename="target.gff", descr="target") {
    gff=data.frame(layout$chr, "nimble", descr, layout$pos, 
        (layout$pos + layout$length), round(as.vector(y), 2), "+", "",
        paste("Name=ng",descr, sep=""))
    con <- file(filename, open="wt")
    writeLines("##gff-version 3", con, sep="\n")
    writeLines(paste("##", date(), sep=""), con, sep="\n")
    write.table(gff, row.names=F, col.names=F, quote=F, sep="\t", con)
    close(con)
}

## write region information into an gff file
## regions - a region table with chr (chromosome), start, and end columns
write.gff.regions<-function(regions, filename="regions.gff", descr="region") {
    gff=data.frame(regions$chr, "nimble", descr, regions$start, 
        regions$end, round(regions$score, 2), "+", " ",
        paste("Name=region_",descr, rownames(regions),sep=""))
    con <- file(filename, open="wt")
    writeLines("##gff-version 3", con, sep="\n")
    writeLines(paste("##", date(), sep=""), con, sep="\n")
    write.table(gff, row.names=F, col.names=F, quote=F, sep="\t", con)
    close(con)
}

###############################################################################
## QUALITY CONTROL 
###############################################################################

## quality plots for pre-processed signals
## operates on a NimbleBatch object
quality.control <- function(signals, samples, sample.name) {
    require(affy)
    require(vsn)
    inp <- signals[,which(samples$sample.name==sample.name &
        samples$sample.type=="i")]
    enr <- signals[,which(samples$sample.name==sample.name &
        samples$sample.type=="e")]
    par(mfrow=c(3,2), oma = 0.1+c(0,0,4,0)) 
    plot(density(log(enr)), main="IP")
    plot(density(log(inp)), main="Input")
    qqnorm(log2(enr), main="IP", pch=".")
    qqline(log2(enr), lty=2, col=2)
    qqnorm(log2(inp), main="Input", pch=".")
    qqline(log2(inp), lty=2, col=2)
    A = (log2(enr) + log2(inp))/2
    M = (log2(enr) - log2(inp))
    ma.plot(A, M, plot.method="smoothScatter", cex=0.8)
    meanSdPlot(cbind(log2(enr),log2(inp)), ranks=TRUE)
    mtext(sample.name, line = 1, outer = TRUE, cex=1.5, font=2)
}

## reconstruct array image
nimble.image <- function(signal, layout) {
    require(lattice)
    require(grid)
    print(levelplot(log(signals)~layout$y*layout$x, 
        cuts=16, col.regions=rainbow(17), aspect="iso",
        xlab="Y position", ylab="X position", scales=list(cex=1),
        panel=function(...) {
            grid.rect(gp=gpar(col=NA, fill="black"))
            panel.levelplot(...)
        }))
} 

###############################################################################
## PAIRWISE CORRELATIONS
###############################################################################

# pair wise correlations of all columns in signals
# scatterplots and spearman correlation coefficients
correlate.batch <- function(signals) {
    require(geneplotter)
    panel.cor <- function(x, y, digits=2, prefix="", cex.cor) {
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(0,1,0,1))
        r <- (cor(x, y,method="spearman"))
        txt <- format(c(r, 0.123456789), digits=digits)[1]
        txt <- paste(prefix, txt, sep="") 
        if(missing(cex.cor))
        text(0.5, 0.5, txt) 
    } 
    pairs(signals , 
        upper.panel=function(...) {par(new=TRUE);smoothScatter(..., nrpoints=0)}, 
        lower.panel=panel.cor)
}

###############################################################################
## PROBE LEVEL STATS
###############################################################################

# calculate the probe mean ratios of all replicates
probe.mean <- function(signals, sample.type) {
    require(st)
    X = t(signals)
    L <- ifelse(sample.type=="i", 2, 1)
    diffmean.stat(X, L)    
}

# calculate probe summary statiscs. default is sam
probe.stat <- function(signals, sample.type, algorithm="sam") {
    if (algorithm=="sam") {
        L <- ifelse(sample.type=="i", 1, 2)
    } else {
        L <- ifelse(sample.type=="i", 2, 1)
    }
    X = t(signals) 
    do.call(paste(algorithm, ".stat", sep=""), list(X, L))
}


###############################################################################
## HMM
###############################################################################

# computes hmm based states for each probe
# stat = test statistic for each probe
# layout = layout description of the probes (chr, pos and length columns required)
# frag.size = mean fragment size of chromatin
# returns a list with regions and probe.states 
stateHMM <- function(stat, layout, frag.size=700) { 
    require(tileHMM)
    hmm.init<-hmm.setup(stat, state=c("non-enriched","enriched"),
        pos.state=2, probe.region=median(layout$length), frag.size=frag.size)
    max.gap <- 400 
    gap <- diff(layout$pos) > max.gap 
    gap.idx <- which(gap) 
    start <- c(1, gap.idx + 1) 
    end <- c(gap.idx, length(stat)) 
    lst <- mapply(function(s, e, data) data[s:e], start, end, MoreArgs = list(stat)) 
    lengths <- unlist(lapply(lst, length))
    flt <- lengths > 100
    hmm.opt <- viterbiEM(hmm.init, lst[flt], df = 9, verbose = 2) 
    post <- lapply(lst, posterior, hmm.opt)
    state.seq <- lapply(post, apply, 2, which.max)
    state.seq <- states(hmm.opt)[c(state.seq, recursive = TRUE)] 
    probe.state <- ifelse(state.seq=="enriched", 1, 0)
    regions.idx <- region.position(state.seq, region = "enriched")
    regions.pos <- matrix(layout[regions.idx, 6], nrow = 2, ncol = dim(regions.idx)[2])
    post.enriched <- lapply(post,"[",2,)
    post.enriched <- exp(c(post.enriched, recursive = TRUE))
    region.score <- apply(regions.idx, 2, function(reg, post) mean(post[reg[1]:reg[2]]), 
        post.enriched)
    regions.clean <- remove.short(regions.idx,post.enriched, 
        data.frame(chromosome=layout$chr, position=layout$pos), min.length=400, min.score=0.8)
    gff <- reg2gff(regions.clean, post.enriched, data.frame(chromosome=layout$chr, 
        position=layout$pos))
    regions <- data.frame(chr=gff$chr,start=gff$start, end=gff$end, score=gff$score)
    l <- list(regions=regions, probe.state=probe.state)
    l
}


