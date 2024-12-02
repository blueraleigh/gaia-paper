foo = function(x) {
    node_age = x[, 1]
    node_err = x[, 2]
    ages = c(seq(0,10,1),seq(20,100,10),seq(200,1000,100),seq(2000,10000,1000))
    q025 = tapply(
        node_err, 
        findInterval(node_age, ages, rightmost.closed=TRUE, left.open=TRUE),
        quantile, prob=0.025
    )
    q975 = tapply(
        node_err, 
        findInterval(node_age, ages, rightmost.closed=TRUE, left.open=TRUE),
        quantile, prob=0.975
    )
    cbind(ages[-1], q025, q975)
}

x = read.csv("../ancestor-estimates-linear.csv", header=FALSE)
y = read.csv("../ancestor-estimates-linear-pareto.csv", header=FALSE)

png(file="ancestor-estimates-linear.png", width=6.5, height=3.5, units="in", res=300)
par(oma=c(4,4,0,0), mar=c(1,1,1,1))
layout(matrix(1:2, 1, 2))
plot(x[,3], x[,4], log='x', las=1, bty="l",
    xlab="Ancestor age",
    ylab="Ancestor location error",
    type="n", xaxt="n", yaxt="n", ylim=c(0,0.8))
axis(1, at=c(1,10,100,1000,10000,''))
axis(2, las=1)
mtext("Ancestor location error", 2, line=2.5, cex=0.8)
mtext("Ancestor age (generations)", 1, line=2.5, cex=0.8)
myhsp(x[,3], x[,4], colpal="heat", log='x',pch=19)
p = foo(x[, 3:4])
polygon(c(p[,1], rev(p[,1])), c(p[,3],rev(p[,2])), lwd=1)
legend("topleft", legend="Gaussian kernel", bty="n")
plot(y[,3], y[,4], log='x', las=1, bty="l",
    xlab="Ancestor age",
    ylab="Ancestor location error",
    type="n", xaxt="n", yaxt="n", ylim=c(0,0.8))
axis(1, at=c(1,10,100,1000,10000,''))
axis(2, las=1, at=seq(0,0.8,0.2), labels=rep('', 5))
mtext("Ancestor age (generations)", 1, line=2.5, cex=0.8)
myhsp(y[,3], y[,4], colpal="heat", log='x',pch=19)
p = foo(y[, 3:4])
polygon(c(p[,1], rev(p[,1])), c(p[,3],rev(p[,2])), lwd=1)
legend("topleft", legend="Pareto kernel", bty="n")
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
