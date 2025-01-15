library(sf)
library(stars)

stem.cost = function(node, x, y, parent, time, node_cost, stem_cost) {
    u = node + 1
    if (parent[u] == -1) return (stem_cost)
    b = time[parent[u]+1] - time[u]
    num_children = length(which(parent == node))
    if (num_children > 0)
    {
        p1 = node_cost[u, 1] / (b*node_cost[u, 1] + 1)
        p2 = node_cost[u, 2] / (b*node_cost[u, 2] + 1)
        p3 = node_cost[u, 3] / (b*node_cost[u, 1] + 1)
        p4 = node_cost[u, 4] / (b*node_cost[u, 2] + 1)
        p5 = node_cost[u, 5] - 
            ((node_cost[u, 3])^2 / (4*(node_cost[u, 1]+ 1/b))) -
            ((node_cost[u, 4])^2 / (4*(node_cost[u, 2]+ 1/b)))
    }
    else
    {
        p1 = 1 / b
        p2 = 1 / b
        p3 = -(2*x[u]) / b
        p4 = -(2*y[u]) / b
        p5 = (x[u]^2 + y[u]^2) / b
    }
    stem_cost[u,] = c(p1,p2,p3,p4,p5)
    return (stem_cost)
}

node.cost = function(node, x, y, parent, time, node_cost, stem_cost) {
    u = node + 1
    children = which(parent == node)
    num_children = length(children)
    if (num_children == 0) return (node_cost)
    node_cost[u,] = colSums(stem_cost[children,,drop=FALSE])
    return (node_cost)
}

final.cost = function(node, x, y, parent, time, node_cost, stem_cost,
    final_cost)
{
    u = node + 1
    if (parent[u] == -1) {
        final_cost[u,] = node_cost[u,]
        return (final_cost)
    }
    v = u
    u = parent[v] + 1
    b = time[u] - time[v]
    p = final_cost[u, ] - stem_cost[v, ]
    final_cost[v,1] = p[1] / (b*p[1]+1)
    final_cost[v,2] = p[2] / (b*p[2]+1)
    final_cost[v,3] = p[3] / (b*p[1]+1)
    final_cost[v,4] = p[4] / (b*p[2]+1)
    final_cost[v,5] = p[5] - 
            ((p[3])^2 / (4*(p[1]+ 1/b))) -
            ((p[4])^2 / (4*(p[2]+ 1/b)))
    final_cost[v,] = final_cost[v,] + node_cost[v,]
    return (final_cost)
}

mpr = function(x, y, parent, time, postorder, num_nodes) {
    node_cost = matrix(0, num_nodes, 5)
    stem_cost = matrix(0, num_nodes, 5)
    final_cost = matrix(0, num_nodes, 5)
    for (node in postorder)
    {
        node_cost = node.cost(node,x,y,parent,time,node_cost,stem_cost)
        stem_cost = stem.cost(node,x,y,parent,time,node_cost,stem_cost)
    }
    for (node in rev(postorder))
    {
        final_cost = final.cost(
            node,x,y,parent,time,node_cost,stem_cost,final_cost)
    }
    return (structure(final_cost, node_cost=node_cost, stem_cost=stem_cost))
}


plot.tree = function(parent, time, postorder, num_nodes, 
    xlim=c(1,3), ylim=c(0,1), ...)
{
    coords = matrix(NA, num_nodes, 2)
    plot.new()
    plot.window(xlim=xlim, ylim=ylim)
    t = 1
    for (node in postorder) {
        u = node + 1
        children = which(parent == node)
        if (length(children)) {
            x0 = mean(coords[children, 1])
            y0 = time[u]
            segments(coords[children[1], 1], y0,
                coords[children[2], 1], y0, ...)
            if (parent[u] != -1) {
                x1 = x0
                y1 = time[parent[u]+1]
                segments(x0,y0,x1,y1, ...)
            }
            coords[u,] = c(x0, y0)
        } else {
            x0 = t
            y0 = 0
            x1 = t
            y1 = time[parent[u]+1]
            segments(x0,y0,x1,y1, ...)
            coords[u,] = c(x0, y0)
            t = t + 1
        }
    }
    invisible(coords)
}

node.labels = function(nodes, coords, ...) {
    u = nodes + 1
    points(coords[u, 1], coords[u, 2], ...)
}

node.text = function(nodes, labels, coords, ...) {
    u = nodes + 1
    text(coords[u, 1], coords[u, 2], labels, ...)
}

# Example tree sequence:
#
#      [0, 20)                     [20, 80)                  [80, 100)
#
#    +----6----+
#    |         |                                           +----5----+
#    |         |                                           |         |
#    |     +---4---+              +---4---+                |     +---4---+
#    |     |       |              |       |                |     |       |
#    |     |       |           +--3--+    |                |     |       |  
#    |     |       |           |     |    |                |     |       |
#    0     1       2           0     2    1                0     1       2
#
#
# node_id    time
#  0          0
#  1          0
#  2          0
#  3          0.15
#  4          0.60
#  5          0.80
#  6          1.00
#


t1 = c(6,4,4,-1,6,-1,-1)
t2 = c(3,4,3,4,-1,-1,-1)
t3 = c(5,4,4,-1,5,-1,-1)
time = c(0,0,0,0.15,0.6,0.8,1.0)
x = c(0.5, 1.25, 2)
y = c(2.7, .41, 1.5)

s1 = mpr(x, y, t1, time, c(1,2,4,0,6), 7)
s2 = mpr(x, y, t2, time, c(0,2,3,1,4), 7)
s3 = mpr(x, y, t3, time, c(1,2,4,0,5), 7)
s4 = .2*s1[5,] + .6*s2[5,] + .2*s3[5,]


scoreq = function(x, y, fn) {
    fn[1]*x*x + fn[2]*y*y + fn[3]*x + fn[4]*y + fn[5]
}

minq = function(fn) {
    c(-fn[3]/(2*fn[1]), -fn[4]/(2*fn[2]))
}


vv = seq(0,3,.01)
v = (vv[-301] + vv[-1]) / 2


layout.mat = matrix(
c(1,2,3,4,4,
  8,8,8,4,4,
  5,6,7,4,4), ncol=5, byrow=TRUE
)

pdf(file="concept-fig.pdf", width=7.25, height=4.25, colormodel="cmyk")
#layout.show(layout(layout.mat, heights=c(1,.1,.6)))
layout(layout.mat, heights=c(1,.1,.6))

# local trees

par(mar=c(1,1,1,1), oma=c(1,1,1,1), xpd=NA, ps=9)
tc = plot.tree(t1, time, c(1,2,4,0,6), 7)
node.labels(c(4, 6), tc, pch=19, cex=1)
node.text(c(4, 6), c(4,6), tc, cex=1, adj=c(-.5,-.5))
node.labels(0:2, tc, pch=22, bg="white", cex=1)
node.text(0:2, 0:2, tc, cex=1, pos=1)

tc = plot.tree(t2, time, c(0,2,3,1,4), 7)
node.labels(c(3, 4), tc, pch=19, cex=1)
node.text(c(3, 4), c(3,4), tc, cex=1, adj=c(-.5,-.5))
node.labels(0:2, tc, pch=22, bg="white", cex=1)
node.text(0:2, 0:2, tc, cex=1, pos=1)

tc = plot.tree(t3, time, c(1,2,4,0,5), 7)
node.labels(c(4, 5), tc, pch=19, cex=1)
node.text(c(4, 5), c(4,5), tc, cex=1, adj=c(-.5,-.5))
node.labels(0:2, tc, pch=22, bg="white", cex=1)
node.text(0:2, 0:2, tc, cex=1, pos=1)
axis(4,las=1,line=1, tcl=-0.2, mgp=c(3,.5,0),cex.axis=1)
text(3,1,"time", adj=c(-1,-1),cex=1)

# cost surfaces

colpal = LSD::colorpalette("heat")
colpal = colorRampPalette(colpal[c(1:4)])(18)

par(mar=c(0,0,0,0))
z = outer(v, v, scoreq, fn=s2[4,])
pmat = persp(
    x=seq(0,3,.01),
    y=seq(0,3,.01),
    matrix(0.15,301,301),zlim=c(0,1),
    col="white"
    ,border=NA,box=FALSE,ticktype="detailed",expand=1.75,phi=-3,theta=75,r=2*sqrt(3),d=1,
    xlab="X", ylab="Y",zlab="Time"
)
cl = st_contour(st_as_stars(z))
for (i in length(cl[[4]]):1)
{
    cr = st_coordinates(cl[[4]][[i]])
    crs = split(as.data.frame(cr), cr[,4])
    for (j in 1:length(crs))
    {
        tr = trans3d(crs[[j]]$X/100, crs[[j]]$Y/100, 0.15, pmat)
        polygon(tr$x, tr$y, col=colpal[i],
            border=colpal[i])
    }
}
z = outer(v, v, scoreq, fn=s4)
cl = st_contour(st_as_stars(z))
for (i in length(cl[[4]]):1)
{
    cr = st_coordinates(cl[[4]][[i]])
    crs = split(as.data.frame(cr), cr[,4])
    for (j in 1:length(crs))
    {
        tr = trans3d(crs[[j]]$X/100, crs[[j]]$Y/100, 0.6, pmat)
        polygon(tr$x, tr$y, col=colpal[i],
            border=colpal[i])
    }
}
z = outer(v, v, scoreq, fn=s3[6,])
cl = st_contour(st_as_stars(z))
for (i in length(cl[[4]]):1)
{
    cr = st_coordinates(cl[[4]][[i]])
    crs = split(as.data.frame(cr), cr[,4])
    for (j in 1:length(crs))
    {
        tr = trans3d(crs[[j]]$X/100, crs[[j]]$Y/100, 0.8, pmat)
        polygon(tr$x, tr$y, col=colpal[i],
            border=colpal[i])
    }
}
z = outer(v, v, scoreq, fn=s1[7,])
cl = st_contour(st_as_stars(z))
for (i in length(cl[[4]]):1)
{
    cr = st_coordinates(cl[[4]][[i]])
    crs = split(as.data.frame(cr), cr[,4])
    for (j in 1:length(crs))
    {
        tr = trans3d(crs[[j]]$X/100, crs[[j]]$Y/100, 1, pmat)
        polygon(tr$x, tr$y, col=colpal[i],
            border=colpal[i])
    }
}

# ARG nodes

points(trans3d(minq(s2[4,])[1], minq(s2[4,])[2], 0.15, pmat),pch=19,cex=1)
text(
    trans3d(minq(s2[4,])[1], minq(s2[4,])[2], 0.15, pmat)$x,
    trans3d(minq(s2[4,])[1], minq(s2[4,])[2], 0.15, pmat)$y, 3, pos=4)
points(trans3d(minq(s4)[1], minq(s4)[2], .6, pmat),pch=19,cex=1)
text(
    trans3d(minq(s4)[1], minq(s4)[2], .6, pmat)$x,
    trans3d(minq(s4)[1], minq(s4)[2], .6, pmat)$y, 4, pos=4)
points(trans3d(minq(s3[6,])[1], minq(s3[6,])[2], 0.8, pmat),pch=19,cex=1)
text(
    trans3d(minq(s3[6,])[1], minq(s3[6,])[2], 0.8, pmat)$x,
    trans3d(minq(s3[6,])[1], minq(s3[6,])[2], 0.8, pmat)$y, 5, pos=4)
points(trans3d(minq(s1[7,])[1], minq(s1[7,])[2], 1, pmat),pch=19,cex=1)
text(
    trans3d(minq(s1[7,])[1], minq(s1[7,])[2], 1, pmat)$x,
    trans3d(minq(s1[7,])[1], minq(s1[7,])[2], 1, pmat)$y, 6, pos=4)

# ARG edges

segments(
    trans3d(minq(s1[7,])[1], minq(s1[7,])[2], 1, pmat)$x,
    trans3d(minq(s1[7,])[1], minq(s1[7,])[2], 1, pmat)$y,
    trans3d(minq(s4)[1], minq(s4)[2], .6, pmat)$x,
    trans3d(minq(s4)[1], minq(s4)[2], .6, pmat)$y
)
segments(
    trans3d(minq(s3[6,])[1], minq(s3[6,])[2], .8, pmat)$x,
    trans3d(minq(s3[6,])[1], minq(s3[6,])[2], .8, pmat)$y,
    trans3d(minq(s4)[1], minq(s4)[2], .6, pmat)$x,
    trans3d(minq(s4)[1], minq(s4)[2], .6, pmat)$y
)
segments(
    trans3d(minq(s4)[1], minq(s4)[2], .6, pmat)$x,
    trans3d(minq(s4)[1], minq(s4)[2], .6, pmat)$y,
    trans3d(minq(s2[4,])[1], minq(s2[4,])[2], .15, pmat)$x,
    trans3d(minq(s2[4,])[1], minq(s2[4,])[2], .15, pmat)$y
)
segments(
    trans3d(minq(s4)[1], minq(s4)[2], .6, pmat)$x,
    trans3d(minq(s4)[1], minq(s4)[2], .6, pmat)$y,
    trans3d(x[2], y[2], 0, pmat)$x,
    trans3d(x[2], y[2], 0, pmat)$y
)
segments(
    trans3d(minq(s4)[1], minq(s4)[2], .6, pmat)$x,
    trans3d(minq(s4)[1], minq(s4)[2], .6, pmat)$y,
    trans3d(x[3], y[3], 0, pmat)$x,
    trans3d(x[3], y[3], 0, pmat)$y
)
segments(
    trans3d(minq(s2[4,])[1], minq(s2[4,])[2], .15, pmat)$x,
    trans3d(minq(s2[4,])[1], minq(s2[4,])[2], .15, pmat)$y,
    trans3d(x[3], y[3], 0, pmat)$x,
    trans3d(x[3], y[3], 0, pmat)$y
)
segments(
    trans3d(minq(s2[4,])[1], minq(s2[4,])[2], .15, pmat)$x,
    trans3d(minq(s2[4,])[1], minq(s2[4,])[2], .15, pmat)$y,
    trans3d(x[1], y[1], 0, pmat)$x,
    trans3d(x[1], y[1], 0, pmat)$y
)
segments(
    trans3d(minq(s1[7,])[1], minq(s1[7,])[2], 1, pmat)$x,
    trans3d(minq(s1[7,])[1], minq(s1[7,])[2], 1, pmat)$y,
    trans3d(x[1], y[1], 0, pmat)$x,
    trans3d(x[1], y[1], 0, pmat)$y
)
segments(
    trans3d(minq(s3[6,])[1], minq(s3[6,])[2], .8, pmat)$x,
    trans3d(minq(s3[6,])[1], minq(s3[6,])[2], .8, pmat)$y,
    trans3d(x[1], y[1], 0, pmat)$x,
    trans3d(x[1], y[1], 0, pmat)$y
)

# present time slice

polygon(
    trans3d(c(0,0,3,3), c(0,3,3,0), c(0,0), pmat)$x,
    trans3d(c(0,0,3,3), c(0,3,3,0), c(0,0), pmat)$y)


points(trans3d(x, y, 0, pmat),pch=22, bg="white",cex=1)
text(
    trans3d(x, y, 0, pmat)$x,
    trans3d(x, y, 0, pmat)$y,
    0:2, pos=1, cex=1
)

# axes

# run x axis in from (3,3) to (0,3)
arrows(
    trans3d(3,3,0,pmat)$x,
    trans3d(3,3,0,pmat)$y,
    trans3d(0,3,0,pmat)$x,
    trans3d(0,3,0,pmat)$y, length=0.05,lwd=1.1
)
text(trans3d(0,3,0,pmat)$x,trans3d(0,3,0,pmat)$y, "x", pos=4, cex=1)
# run y axis in from (3,3) to (3,0)
arrows(
    trans3d(3,3,0,pmat)$x,
    trans3d(3,3,0,pmat)$y,
    trans3d(3,0,0,pmat)$x,
    trans3d(3,0,0,pmat)$y, length=0.05,lwd=1.4
)
text(trans3d(3,0,0,pmat)$x,trans3d(3,0,0,pmat)$y, "y", pos=1, cex=1)
arrows(
    trans3d(3,3,0,pmat)$x,
    trans3d(3,3,0,pmat)$y,
    trans3d(3,3,1,pmat)$x,
    trans3d(3,3,1,pmat)$y, length=0.05, lwd=1.4
)
text(trans3d(3,3,1,pmat)$x,trans3d(3,3,1,pmat)$y, "time", pos=4, cex=1)


text(trans3d(3,3,-0.08,pmat)$x,
    trans3d(3,3,-0.08,pmat)$y, 
"Average minimum cost dispersal surfaces", cex=1,xpd=NA, adj=c(1,0))

par(mar=c(2,0,2,0))

plot.new()
plot.window(xlim=c(0,1),ylim=c(0,1))

text(0.5, 0.5,
bquote(f[4](x,y) == .(format(s1[5,1],digits=2))*x^2 + 
                       .(format(s1[5,2],digits=2))*y^2 -
                       .(format(abs(s1[5,3]),digits=2))*x - 
                       .(format(abs(s1[5,4]),digits=2))*y + .(format(s1[5,5],digits=2))),
cex=1)

text(0.5, 0.0,
bquote(f[6](x,y) == .(format(s1[7,1],digits=2))*x^2 + 
                       .(format(s1[7,2],digits=2))*y^2 -
                       .(format(abs(s1[7,3]),digits=2))*x - 
                       .(format(abs(s1[7,4]),digits=2))*y + .(format(s1[7,5],digits=2))),
cex=1)

text(0.5, 0.25, "--")
text(0.5, 0.75, "--")

segments(0,0.95,1,0.95)
segments(0,-0.2,1,-0.2)


plot.new()
plot.window(xlim=c(0,1),ylim=c(0,1))

text(0.5, 0.75,
bquote(f[3](x,y) == .(format(s2[4,1],digits=2))*x^2 + 
                       .(format(s2[4,2],digits=2))*y^2 -
                       .(format(abs(s2[4,3]),digits=2))*x - 
                       .(format(abs(s2[4,4]),digits=2))*y + .(format(s2[4,5],digits=2))),
cex=1)

text(0.5, 0.5,
bquote(f[4](x,y) == .(format(s2[5,1],digits=2))*x^2 + 
                       .(format(s2[5,2],digits=2))*y^2 -
                       .(format(abs(s2[5,3]),digits=2))*x - 
                       .(format(abs(s2[5,4]),digits=2))*y + .(format(s2[5,5],digits=2))),
cex=1)

text(0.5, 0.25, "--")
text(0.5, 0.0, "--")

segments(0,0.95,1,0.95)
segments(0,-0.2,1,-0.2)

text(0.5,0.95,"Minimum cost dispersal surfaces",cex=1,adj=c(0.5,-0.8))


plot.new()
plot.window(xlim=c(0,1),ylim=c(0,1))

text(0.52, 0.5,
bquote(f[4](x,y) == .(format(s3[5,1],digits=2))*x^2 + 
                       .(format(s3[5,2],digits=2))*y^2 -
                       .(format(abs(s3[5,3]),digits=2))*x - 
                       .(format(abs(s3[5,4]),digits=2))*y + .(format(s3[5,5],digits=2))),
cex=1)

text(0.52, 0.25,
bquote(f[5](x,y) == .(format(s3[6,1],digits=2))*x^2 + 
                       .(format(s3[6,2],digits=2))*y^2 -
                       .(format(abs(s3[6,3]),digits=2))*x - 
                       .(format(abs(s3[6,4]),digits=2))*y + .(format(s3[6,5],digits=2))),
cex=1)

text(0.5, 0.75, "--")
text(0.5, 0.0, "--")

segments(0,0.95,1,0.95)
segments(0,-0.2,1,-0.2)

par(mar=c(0,0,0,0))
plot.new()
plot.window(xlim=c(0,1), ylim=c(0,1))

segments(0.2,0,0.31,1,lty=2)
segments(0.31,1,0.31,3.5,lty=2,xpd=NA)

segments(0.8,0,0.69,1,lty=2)
segments(0.69,1,0.69,3.5,lty=2,xpd=NA)

axis(1, at=c(0,0.2,0.8,1), tcl=-0.2, mgp=c(3,.25,0), cex.axis=1)
text(1,0,"Genome position",adj=c(1,-0.2),cex=1)

dev.off()
