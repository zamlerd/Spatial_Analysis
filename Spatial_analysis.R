install.packages("spatstat")
library(spatstat)
library(tidyverse)

dat.spat <- read.csv("~/Downloads/object_results.csv", )
mypattern <- ppp(dat.spat$XMax, dat.spat$YMax, xrange = c(3,1349),yrange = c(2.5, 1009))
plot(mypattern)
summary(mypattern)
plot(Kest(mypattern))
plot(density(mypattern))
marx <- nonadat.spat[,c(10,11,12,13,14,40,41,42,43,44,45)]
plot(Smooth(nonamypattern))
mypattern <- na.omit(mypattern)
is.marked(mypattern)
nonadat.spat <- na.omit(dat.spat)
nonamypattern <- ppp(nonadat.spat$XMax, nonadat.spat$YMax, xrange = c(3,1349),yrange = c(2.5, 1009), marks = nonadat.spat[,c(15,21,26)])
plot(nonamypattern)
poly <- cbind(d1.spat$XMax, d1.spat$YMax)
spat.levels <- levels(d1.spat$Image.Location)
d1.spat <- read.csv("~/Documents/IACS/MS-FIT/Object_Data/D1_Object_Data.csv")
cut.ppp(d1.spat, z=marks(d1.spat$Image.Location))
mypattern <- ppp(d1.spat$XMax, d1.spat$YMax, poly = poly)
plot(mypattern)
summary(mypattern)
plot(Kest(mypattern))
plot(density(mypattern))


###########
#WORKSPACE#
###########

myeloid.spat <- read.csv("~/Documents/IACS/MS-FIT/Myeloid.csv", )
mye.list = list.files(path = "~/Documents/IACS/MS-FIT/Mye/", pattern = "*.csv")
setwd("~/Documents/IACS/MS-FIT/Mye/")
mye = lapply(mye.list, read.delim)
tcell.spat <- read.csv("~/Documents/IACS/MS-FIT/T-cell.csv", )

mye.lpa3.c11.spat <- read.csv("~/Documents/IACS/MS-FIT/Mye/Mye_D9_C11.csv", )
mye.lpa3.c12.spat <- read.csv("~/Documents/IACS/MS-FIT/Mye/Mye_D9_C12.csv", )
mye.lpa3.c14.spat <- read.csv("~/Documents/IACS/MS-FIT/Mye/Mye_D9_C14.csv", )
mye.lpa3.c15.spat <- read.csv("~/Documents/IACS/MS-FIT/Mye/Mye_D9_C15.csv", )
mye.393p.A1.spat <- read.csv("~/Documents/IACS/MS-FIT/Mye/Mye_D10_A1.csv", )
mye.393p.A5.spat <- read.csv("~/Documents/IACS/MS-FIT/Mye/Mye_D10_A5.csv", )
mye.393p.A6.spat <- read.csv("~/Documents/IACS/MS-FIT/Mye/Mye_D10_A6.csv", )
mye.bp.A1.spat <- read.csv("~/Documents/IACS/MS-FIT/Mye/Mye_D11_A3.csv", )
mye.bp.A8.spat <- read.csv("~/Documents/IACS/MS-FIT/Mye/Mye_D11_A8.csv", )
mye.yum31.A8.spat <- read.csv("~/Documents/IACS/MS-FIT/Mye/Mye_D12_A6.csv", )

mye.393p.A1.pat <- ppp(mye.393p.A1.spat$XMax, mye.393p.A1.spat$YMax, range(mye.393p.A1.spat$XMax),range(mye.393p.A1.spat$YMax), marks = mye.393p.A1.mark)
mye.393p.A1.mark <- mye.393p.A1.spat[,c(6,7,8,9,11)]
plot(mye.393p.A1.pat.micro)

mye.393p.A1.pat.micro <- subset.ppp(mye.393p.A1.pat, mye.393p.A1.pat$marks$Microglia > 0)

myeloid.pattern <- ppp(myeloid.spat$XMax, myeloid.spat$YMax, range(myeloid.spat$XMax), range(myeloid.spat$YMax), marks = marx)
split.mye.pattern <- cut.ppp(myeloid.pattern, z = marks(myeloid.pattern), breaks = 5)
plot(myeloid.pattern)
summary(mypattern)
plot(Kest(mypattern))
plot(density(mypattern))
plot(Smooth(nonamypattern))
micro.pattern <- subset.ppp(myeloid.pattern, myeloid.pattern$marks$Microglia > 0)
cd163.pattern <- subset.ppp(myeloid.pattern, c(myeloid.pattern$marks$CD163..Microglia > 0, myeloid.pattern$marks$CD163..Macrophages > 0 ))
summary(myeloid.pattern)
tcell.pattern <- ppp(tcell.spat$XMax, tcell.spat$YMax)
plot(tcell.pattern)

#bring marks about cytoplasm and nuclei over
na.myeloid.spat <- na.omit(myeloid.spat)
marx <- na.myeloid.spat[,c(6, 7, 8, 9, 10)]



