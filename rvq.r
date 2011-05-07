rvq <- function(x, centers, iter.max = 10, nstart = 1)
{
    dyn.load("rvq.so")
    do_one <- function() {
        Z <- .C("kmeans_Rvq", as.double(x), as.integer(m),
                               as.integer(ncol(x)),
                               centers = as.double(centers), as.integer(k),
                               c1 = integer(m), iter = as.integer(iter.max),
                               nc = integer(k), wss = double(k))
                       if(Z$iter > iter.max)
                           warning("did not converge in ",
                                  iter.max, " iterations", call.=FALSE)
                       if(any(Z$nc == 0))
                           warning("empty cluster: try a better set of initial centers", call.=FALSE)
                       Z
    }
    x <- as.matrix(x)
    m <- nrow(x)
    if(missing(centers))
	stop("'centers' must be a number or a matrix")
    if(length(centers) == 1L) {
	k <- centers
        ## we need to avoid duplicates here
        if(nstart == 1)
            centers <- x[sample.int(m, k), , drop = FALSE]
        if(nstart >= 2 || any(duplicated(centers))) {
            cn <- unique(x)
            mm <- nrow(cn)
            if(mm < k)
                stop("more cluster centers than distinct data points.")
            centers <- cn[sample.int(mm, k), , drop=FALSE]
        }
    } else {
	centers <- as.matrix(centers)
        if(any(duplicated(centers)))
            stop("initial centers are not distinct")
        cn <- NULL
	k <- nrow(centers)
        if(m < k)
            stop("more cluster centers than data points")
    }
    if(iter.max < 1) stop("'iter.max' must be positive")
    if(ncol(x) != ncol(centers))
	stop("must have same number of columns in 'x' and 'centers'")
    Z <- do_one()
    best <- sum(Z$wss)
    if(nstart >= 2 && !is.null(cn))
	for(i in 2:nstart) {
	    centers <- cn[sample.int(mm, k), , drop=FALSE]
	    ZZ <- do_one()
	    if((z <- sum(ZZ$wss)) < best) {
		Z <- ZZ
		best <- z
	    }
	}
    centers <- matrix(Z$centers, k)
    dimnames(centers) <- list(1L:k, dimnames(x)[[2L]])
    cluster <- Z$c1
    if(!is.null(rn <- rownames(x)))
        names(cluster) <- rn
    totss <- sum(scale(x, scale = FALSE)^2)
    structure(list(cluster = cluster, centers = centers, totss = totss,
                   withinss = Z$wss, tot.withinss = best,
                   betweenss = totss - best, size = Z$nc),
	      class = "rvq")
}

## modelled on print methods in the cluster package
print.rvq <- function(x, ...)
{
    cat("RVQ clustering with ", length(x$size), " clusters of sizes ",
        paste(x$size, collapse=", "), "\n", sep="")
    cat("\nCluster means:\n")
    print(x$centers, ...)
    cat("\nClustering vector:\n")
    print(x$cluster, ...)
    cat("\nWithin cluster sum of squares by cluster:\n")
    print(x$withinss, ...)
    cat(sprintf(" (between_SS / total_SS = %5.1f %%)\n",
		100 * x$betweenss/x$totss),
	"Available components:\n", sep="\n")
    print(names(x))
    invisible(x)
}
