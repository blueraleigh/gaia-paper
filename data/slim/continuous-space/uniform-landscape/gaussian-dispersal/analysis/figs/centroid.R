
x = read.csv("../ancestor-estimates-v2.csv", header=FALSE)

png(file="ancestor-centroid-error.png", width=6.5, height=10, units="in", res=300)
SIGMA = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0)
par(oma=c(4,4,0,0), mar=c(1,1,1,1))
layout(matrix(1:10, 5, 2))
for (i in 1:10) {
    idx = x[,1] == SIGMA[i]
    plot(x[idx,7], x[idx,5], log='', las=1, bty="l",
        xlab="Estimated distance from sample centroid",
        ylab="Ancestor location error",
        type="n", xaxt="n", yaxt="n", ylim=c(0,0.8))
    axis(1)
    #if (i == 5 || i == 10)
    #    axis(1, at=c(0,100,200,300,400,500,600))
    #else
    #    axis(1, at=c(0,100,200,300,400,500,600), labels=c('','','','','','',''))
    if (i %in% 1:5)
        axis(2, las=1)
    else
        axis(2, seq(0,0.8,.2),labels=rep('',5))
    if (i %in% 1:5)
        mtext("Ancestor location error", 2, line=2.5, cex=0.8)
    if (i == 5 || i == 10)
        mtext("Estimated distance from sample centroid", 1, line=2.5, cex=0.8)
    myhsp(x[idx,7], x[idx,5], colpal="heat", log='',pch=19)
    legend("topright", legend=bquote(paste(sigma, " = ", .(SIGMA[i]))), bty="n")
}
dev.off()



png(file="ancestor-centroid-estimates.png", width=6.5, height=10, units="in", res=300)
SIGMA = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0)
par(oma=c(4,4,0,0), mar=c(1,1,1,1))
layout(matrix(1:10, 5, 2))
for (i in 1:10) {
    idx = x[,1] == SIGMA[i]
    plot(x[idx,6], x[idx,7], log='', las=1, bty="l",
        xlab="True distance from sample centroid",
        ylab="Estimated distance from sample centroid",
        type="n", xaxt="n", yaxt="n")
    axis(1)
    #if (i == 5 || i == 10)
    #    axis(1, at=c(0,100,200,300,400,500,600))
    #else
    #    axis(1, at=c(0,100,200,300,400,500,600), labels=c('','','','','','',''))
    
    axis(2, las=1)
    
    if (i == 3)
        mtext("Estimated distance from sample centroid", 2, line=2.5, cex=0.8)
    if (i == 5 || i == 10)
        mtext("True distance from sample centroid", 1, line=2.5, cex=0.8)
    myhsp(x[idx,6], x[idx,7], colpal="heat", log='',pch=19)
    abline(0,1)
    legend("topleft", legend=bquote(paste(sigma, " = ", .(SIGMA[i]))), bty="n")
}
dev.off()

png(file="ancestor-centroid-age-true.png", width=6.5, height=10, units="in", res=300)
SIGMA = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0)
par(oma=c(4,4,0,0), mar=c(1,1,1,1))
layout(matrix(1:10, 5, 2))
for (i in 1:10) {
    idx = x[,1] == SIGMA[i]
    plot(x[idx,4], x[idx,6], log='x', las=1, bty="l",
        xlab="Ancestor age",
        ylab="True distance from sample centroid",
        type="n", xaxt="n", yaxt="n")
    axis(1)
    #if (i == 5 || i == 10)
    #    axis(1, at=c(0,100,200,300,400,500,600))
    #else
    #    axis(1, at=c(0,100,200,300,400,500,600), labels=c('','','','','','',''))
    
    axis(2, las=1)
    
    if (i == 3)
        mtext("True distance from sample centroid", 2, line=2.5, cex=0.8)
    if (i == 5 || i == 10)
        mtext("Ancestor age", 1, line=2.5, cex=0.8)
    myhsp(x[idx,4], x[idx,6], colpal="heat", log='x',pch=19)
    legend("topright", legend=bquote(paste(sigma, " = ", .(SIGMA[i]))), bty="n")
}
dev.off()

library(LSD)
# modified from LSD::heatscatterpoints
myhsp = function(x, y, pch = 19, cexplot = 0.5, nrcol = 30, grid = 100, 
    colpal = "heat", simulate = FALSE, daltonize = FALSE, cvd = "p", 
    alpha = NULL, rev = FALSE, xlim = NULL, ylim = NULL, only = "none", 
    add.contour = FALSE, nlevels = 10, color.contour = "black", 
    greyscale = FALSE, log = "", ...) 
{
    if (!is.vector(x) | !is.vector(y)) 
        stop("First two argument must be numeric vectors!")
    if (length(x) != length(y)) 
        stop("Data vectors must be of the same length!")
    sound = which((!(is.na(x) | is.nan(x) | (x == Inf) | (x == 
        -Inf))) & (!(is.na(y) | is.nan(y) | (y == Inf) | (y == 
        -Inf))))
    if (length(sound) == 0) 
        stop("There are no valid point pairs to plot!")
    x = x[sound]
    y = y[sound]
    if (!is.null(xlim)) {
        cut = x >= xlim[1] & x <= xlim[2]
        x = x[cut]
        y = y[cut]
    }
    if (!is.null(ylim)) {
        cut = y >= ylim[1] & y <= ylim[2]
        y = y[cut]
        x = x[cut]
    }
    colpal = colorpalette(colpal, nrcol, simulate = simulate, 
        daltonize = daltonize, cvd = cvd, alpha = alpha, rev = rev)
    if (greyscale) {
        colpal = convertgrey(colpal)
    }
    todiscrete = function(t, tmin, tmax, bins) {
        erg = round((t - tmin)/(tmax - tmin) * bins + 0.5)
        erg = pmin(pmax(erg, 1), bins)
        return(erg)
    }
    kde2d.adj = function(x, y, h, n = 25, lims = c(range(x), 
        range(y)), only = "none") {
        nx = length(x)
        gx = seq.int(lims[1], lims[2], length.out = n)
        gy = seq.int(lims[3], lims[4], length.out = n)
        bandwidth.nrd.adj = function(x) {
            r = quantile(x, c(0.25, 0.75))
            h = (r[2] - r[1])/1.34
            return(4 * 1.06 * min(sqrt(var(x)), h) * length(x)^(-1/5))
        }
        if (missing(h)) {
            bx = bandwidth.nrd.adj(x)
            by = bandwidth.nrd.adj(y)
            if (all(c(bx, by) == 0)) {
                h = rep(0.01, 2)
            }
            else if (any(c(bx, by) == 0)) {
                h = rep(max(bx, by), 2)
            }
            else {
                h = c(bx, by)
            }
        }
        else h = rep(h, length.out = 2)
        h = h/4
        ax = outer(gx, x, "-")/h[1]
        ay = outer(gy, y, "-")/h[2]
        norm.ax = dnorm(ax)
        norm.ay = dnorm(ay)
        if (only == "x") {
            norm.ay = rep(1, length(ay))
        }
        if (only == "y") {
            norm.ax = rep(1, length(ax))
        }
        z = tcrossprod(matrix(norm.ax, , nx), matrix(norm.ay, 
            , nx))/(nx * h[1] * h[2])
        list(x = gx, y = gy, z = z)
    }
    if (log == "") {
        xlog = x
        ylog = y
    }
    else if (log == "x") {
        xlog = log(x, 10)
        ylog = y
    }
    else if (log == "y") {
        xlog = x
        ylog = log(y, 10)
    }
    else if (log %in% c("xy", "yx")) {
        xlog = log(x, 10)
        ylog = log(y, 10)
    }
    d = kde2d.adj(xlog, ylog, n = grid, only = only)
    xdiscrete = todiscrete(xlog, min(xlog), max(xlog), bins = grid)
    ydiscrete = todiscrete(ylog, min(ylog), max(ylog), bins = grid)
    getfrommat = function(a) {
        d$z[a[1], a[2]]
    }
    heatvec = unlist(apply(cbind(xdiscrete, ydiscrete), 1, getfrommat))
    coldiscrete = todiscrete(heatvec, min(d$z), max(d$z), bins = nrcol)
    ord = order(coldiscrete)
    points(x[ord], y[ord], col = colpal[coldiscrete[ord]], pch = pch, cex = cexplot, 
        ...)
    if (add.contour) {
        contour(d, add = TRUE, nlevels = nlevels, col = color.contour)
    }
}
