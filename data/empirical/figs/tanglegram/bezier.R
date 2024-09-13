# Construct a Bezier curve
#
# x       The set of control points ordered by increasing time
# times   Times associated with each control point.
bezier = function(x, times)
{
    n = nrow(x)
    stopifnot(n >= 2)
    tt = (times - min(times)) / (max(times) - min(times))
    if (n == 2)
    {
        function(t) {
            tmp = sapply(t, function(tt) {
                (1 - tt)*x[1,,drop=FALSE] + tt*x[2,,drop=FALSE]
            })
            ti = apply(tmp, 2, function(p) which.min(sqrt(rowSums((sweep(x, 2, p))^2))))
            list(x=tmp[1,], y=tmp[2,], time=times[ti])
        }
    }
    else
    {
        d = n-1
        function(t) {
            tmp = sapply(t, function(tt) {
                p = (1 - tt)^d * x[1,]
                for (i in 2:(n-1)) {
                    p = p + choose(d, i-1) * (1-tt)^(d-i+1) * tt^(i-1) * x[i,]
                }
                p = p + tt^d * x[n,]
                p
            })
            ti = apply(tmp, 2, function(p) which.min(sqrt(rowSums((sweep(x, 2, p))^2))))
            list(x=tmp[1,], y=tmp[2,], time=times[ti])
        }
    }
}
