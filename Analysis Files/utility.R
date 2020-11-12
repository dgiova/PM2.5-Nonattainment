

## puts the results of a pair match in a nice form
summarize.match <- function(dat, ms, ps.name="prop", keep.mset=FALSE) {
    adat <- dat
    adat$mset <- ms
    adat <- adat[!is.na(adat$mset),]
    adat.treat <- adat[adat$z==1, ]
    adat.ctrl <- adat[adat$z==0, ]

    adat.m <- merge(adat.treat, adat.ctrl, by="mset", suffixes=c(".1", ".0"))

    if(!keep.mset) {
        adat.m <- adat.m[, -which(names(adat.m) %in% c("z.1", "z.0", "mset"))]
    } else {
        adat.m <- adat.m[, -which(names(adat.m) %in% c("z.1", "z.0"))]        
    }
    adat.m <- adat.m[, sort(names(adat.m), index.return=TRUE)$ix]
    
    p0.name <- paste0(ps.name,".", 0)
    p1.name <- paste0(ps.name,".",1)

    adat.m.tmp.1 <- adat.m[, -which(names(adat.m) %in% c(p0.name, p1.name))]
    adat.m.tmp.2 <- adat.m[, c(p0.name, p1.name)]

    adat.m <- cbind(adat.m.tmp.1, adat.m.tmp.2)
    
    return(adat.m)
}

## preprocesses the results of pair matching for an analysis
## using `senm'.
cast.senm <- function(dat, ms.arg, two.outcomes=FALSE) {
    ms <- as.vector(ms.arg)

    y <- dat$pmbase[!is.na(ms)]
    mset <- ms[!is.na(ms)]
    z <- dat$a[!is.na(ms)]
    
    dico.names <- unique(mset)
    dico <- seq(length(dico.names))
    names(dico) <- dico.names
    mset <- as.integer(dico[mset])

    if(two.outcomes==FALSE) {
        return(list(y=y, mset=mset, z=z))
    } else {
        y2 <- dat$y2[!is.na(ms)]
        return(list(y=y, y2=y2, mset=mset, z=z))
    }
}
