
# https://community.rstudio.com/t/creating-contour-plots-in-r/111024/2

#generate data
a <- rnorm(100,10,10)
b <- rnorm(100,5,5)
c <- rnorm(100,1,1)
d <- data.frame(a,b,c)
d = as.matrix(d)

image(d)

d <- data.frame(a,b,c)
DENS <- MASS::kde2d(d$a,d$b)
contour(DENS)

filled.contour(DENS,plot.axes = {
   axis(1)
   axis(2)
   contour(DENS,add = TRUE)})


#load library
library(akima)

#create plot

fld <- with(DF, interp(x = X1, y = X2, z = total))

filled.contour(x = fld$x,
               y = fld$y,
               z = fld$z,
               color.palette = colorRampPalette(c("white", "blue")),
               xlab = "X1",
               ylab = "X2",
               main = "Contour Plot",
               key.title = title(main = "total ", cex.main = 1))


library(fields)
image.plot(volcano)
contour(volcano, add = TRUE)

# https://subscription.packtpub.com/book/data/9781783988785/9/ch09lvl1sec99/creating-filled-contour-plots
# 
# https://r-charts.com/correlation/contour-plot/
#    
# https://r-charts.com/correlation/hexbin-chart-ggplot2/
#    
# https://r-graph-gallery.com/hexbin-map.html










