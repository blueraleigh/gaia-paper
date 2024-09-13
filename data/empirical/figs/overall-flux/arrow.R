# Add polygon arrows to a plot
#
# x0,y0     The starting point for the arrow base.
# x1,y1     The ending point for the arrow head.
# awd       The maximum width of the arrow shaft. The 
#           width of arrow head is twice this quantity.
# length    The length of the arrow head.
# col       The color of the arrow interior.
# border    The color of the arrow border.
# f         The width of the arrow base as a fraction of awd.
# ...       Furth arguments passed to polygon.
arrow = function(x0, y0, x1, y1, 
    awd=0.02, length=0.02, col=1, border=1, f=1, ...)
{
    x11 = x1
    y11 = y1

    x1 = 0
    y1 = 0
    x0 = x0 - x11 
    y0 = y0 - y11

    m = (y1 - y0) / (x1 - x0)
    b = y1 - m*x1

    d = sqrt((x1-x0)^2 + (y1-y0)^2)
    r = d - length

    x3 = (1 - r/d)*x0 + (r/d)*x1
    y3 = (1 - r/d)*y0 + (r/d)*y1

    pm = -1 / m
    
    pb1 = y3 - pm*x3 # arrow head perpendicular y-intercept
    pb2 = y0 - pm*x0 # arrow base perpendicular y-intercept

    root1 = -pb1 / pm # arrow head perpendicular x-intercept
    root2 = -pb2 / pm # arrow base perpendicular x-intercept

    # distance from (x3, y3) to perpendicular y-intercept
    d1 = sqrt((x3 - 0)^2 + (y3 - pb1)^2)

    # distance from (x0, y0) to perpendicular y-intercept
    d2 = sqrt((x0 - 0)^2 + (y0 - pb2)^2)

    # distance from (x3, y3) to perpendicular x-intercept
    d3 = sqrt((x3 - root1)^2 + (y3 - 0)^2)
    
    # distance from (x0, y0) to perpendicular x-intercept
    d4 = sqrt((x0 - root2)^2 + (y0 - 0)^2)

    xv = c(
        (1 - 0.5*awd*f/d2)*x0 + (0.5*awd*f/d2)*0
        , (1 - 0.5*awd/d1)*x3 + (0.5*awd/d1)*0
        , (1 - awd/d1)*x3 + (awd/d1)*0
        , x1
        , (1-awd/d3)*x3 + (awd/d3)*root1
        , (1-0.5*awd/d3)*x3 + (0.5*awd/d3)*root1
        , (1 - 0.5*awd*f/d4)*x0 + (0.5*awd*f/d4)*root2
    )
    yv = c(
        (1 - 0.5*awd*f/d2)*y0 + (0.5*awd*f/d2)*pb2
        , (1 - 0.5*awd/d1)*y3 + (0.5*awd/d1)*pb1
        , (1 - awd/d1)*y3 + (awd/d1)*pb1
        , y1
        , (1-awd/d3)*y3 + (awd/d3)*0
        , (1-0.5*awd/d3)*y3 + (0.5*awd/d3)*0
        , (1 - 0.5*awd*f/d4)*y0 + (0.5*awd*f/d4)*0
    )
    
    polygon(xv + x11, yv + y11, col=col, border=border, ...)
}
