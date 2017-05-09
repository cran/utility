# Plots of aggregation functions
# ==============================

if ( !require(utility) ) { install.packages("utility"); library(utility) }

dir.fig <- "./utility/man/figures"

pdf.width1  <- 4.5
pdf.width2  <- 9
pdf.height2 <- 9
pdf.height1 <- 4.5
pdf.mar     <- c(4.5,4,3.5,1) + 0.1  # c(bottom, left, top, right)

png.width1  <- 360
png.width2  <- 720
png.height2 <- 720
png.height1 <- 360
png.mar     <- c(4.5,4,3.5,1) + 0.1  # c(bottom, left, top, right)


obj1 <- utility.endnode.intpol1d.create(name.node   = "v1", 
                                        name.attrib = "v1", 
                                        range       = c(0,1), 
                                        x           = c(0,1), 
                                        u           = c(0,1), 
                                        utility     = FALSE)
obj2 <- utility.endnode.intpol1d.create(name.node   = "v2", 
                                        name.attrib = "v2", 
                                        range       = c(0,1), 
                                        x           = c(0,1), 
                                        u           = c(0,1), 
                                        utility     = FALSE)

# add:
# ----

AddAggregation1 <- utility.aggregation.create(
  name.node = "Additive Aggregation (w1=w2=0.5)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.add", 
  par       = c(1,1))
AddAggregation2 <- utility.aggregation.create(
  name.node = "Additive Aggregation (w1=0.25,w2=0.75)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.add", 
  par       = c(1,3))

pdf(paste(dir.fig,"aggregationadd.pdf",sep="/"),width=pdf.width2,height=pdf.height1)
par(mfrow=c(1,2),mar=pdf.mar)
plot(AddAggregation1,type="nodes",node="Additive Aggregation (w1=w2=0.5)")
plot(AddAggregation2,type="nodes",node="Additive Aggregation (w1=0.25,w2=0.75)")
dev.off()

png(paste(dir.fig,"aggregationadd.png",sep="/"),width=png.width2,height=png.height1)
par(mfrow=c(1,2),mar=png.mar)
plot(AddAggregation1,type="nodes",node="Additive Aggregation (w1=w2=0.5)")
plot(AddAggregation2,type="nodes",node="Additive Aggregation (w1=0.25,w2=0.75)")
dev.off()


# min:
# ----

MinAggregation <- utility.aggregation.create(
  name.node = "Minimum Aggregation", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.min",
  par       = numeric(0))

pdf(paste(dir.fig,"aggregationmin.pdf",sep="/"),width=pdf.width1,height=pdf.height1)
par(mfrow=c(1,1),mar=pdf.mar)
plot(MinAggregation,type="nodes",node="Minimum Aggregation")
dev.off()

png(paste(dir.fig,"aggregationmin.png",sep="/"),width=png.width1,height=png.height1)
par(mfrow=c(1,1),mar=png.mar)
plot(MinAggregation,type="nodes",node="Minimum Aggregation")
dev.off()


# max:
# ----

MaxAggregation <- utility.aggregation.create(
  name.node = "Maximum Aggregation", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.max",
  par       = numeric(0))

pdf(paste(dir.fig,"aggregationmax.pdf",sep="/"),width=pdf.width1,height=pdf.height1)
par(mfrow=c(1,1),mar=pdf.mar)
plot(MaxAggregation,type="nodes",node="Maximum Aggregation")
dev.off()

png(paste(dir.fig,"aggregationmax.png",sep="/"),width=png.width1,height=png.height1)
par(mfrow=c(1,1),mar=png.mar)
plot(MaxAggregation,type="nodes",node="Maximum Aggregation")
dev.off()


# mult:
# -----

MultAggregation1 <- utility.aggregation.create(
  name.node = "Multiplicative Aggregation (w1=w2=0.5)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.mult", 
  par       = c(0.5,0.5))
MultAggregation2 <- utility.aggregation.create(
  name.node = "Multiplicative Aggregation (w1=w2=1)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.mult", 
  par       = c(1,1))
MultAggregation3 <- utility.aggregation.create(
  name.node = "Multiplicative Aggregation (w1=w2=0.3)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.mult", 
  par       = c(0.3,0.3))
MultAggregation4 <- utility.aggregation.create(
  name.node = "Multiplicative Aggregation (w1=0.3,w2=1)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.mult", 
  par       = c(0.3,1))

pdf(paste(dir.fig,"aggregationmult.pdf",sep="/"),width=pdf.width2,height=pdf.height2)
par(mfrow=c(2,2),mar=pdf.mar)
plot(MultAggregation1,type="nodes",node="Multiplicative Aggregation (w1=w2=0.5)")
plot(MultAggregation2,type="nodes",node="Multiplicative Aggregation (w1=w2=1)")
plot(MultAggregation3,type="nodes",node="Multiplicative Aggregation (w1=w2=0.3)")
plot(MultAggregation4,type="nodes",node="Multiplicative Aggregation (w1=0.3,w2=1)")
dev.off()

png(paste(dir.fig,"aggregationmult.png",sep="/"),width=png.width2,height=png.height2)
par(mfrow=c(2,2),mar=png.mar)
plot(MultAggregation1,type="nodes",node="Multiplicative Aggregation (w1=w2=0.5)")
plot(MultAggregation2,type="nodes",node="Multiplicative Aggregation (w1=w2=1)")
plot(MultAggregation3,type="nodes",node="Multiplicative Aggregation (w1=w2=0.3)")
plot(MultAggregation4,type="nodes",node="Multiplicative Aggregation (w1=0.3,w2=1)")
dev.off()


# geo/geooff:
# -----------

GeoAggregation1 <- utility.aggregation.create(
  name.node = "Geometric Aggregation (w1=w2=0.5)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.geo", 
  par       = c(1,1))
GeoAggregation2 <- utility.aggregation.create(
  name.node = "Geometric Aggregation (w1=0.25,w2=0.75)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.geo", 
  par       = c(1,3))
GeoAggregation3 <- utility.aggregation.create(
  name.node = "Geometric Offset Agg. (w1=w2=0.5,d=0.5)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.geooff", 
  par       = c(1,1,0.5))
GeoAggregation4 <- utility.aggregation.create(
  name.node = "Geometric Offset Agg. (w1=0.25,w2=0.75,d=0.5)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.geooff", 
  par       = c(1,3,0.5))

pdf(paste(dir.fig,"aggregationgeo.pdf",sep="/"),width=pdf.width2,height=pdf.height2)
par(mfrow=c(2,2),mar=pdf.mar)
plot(GeoAggregation1,type="nodes",node="Geometric Aggregation (w1=w2=0.5)")
plot(GeoAggregation2,type="nodes",node="Geometric Aggregation (w1=0.25,w2=0.75)")
plot(GeoAggregation3,type="nodes",node="Geometric Offset Agg. (w1=w2=0.5,d=0.5)")
plot(GeoAggregation4,type="nodes",node="Geometric Offset Agg. (w1=0.25,w2=0.75,d=0.5)")
dev.off()

png(paste(dir.fig,"aggregationgeo.png",sep="/"),width=png.width2,height=png.height2)
par(mfrow=c(2,2),mar=png.mar)
plot(GeoAggregation1,type="nodes",node="Geometric Aggregation (w1=w2=0.5)")
plot(GeoAggregation2,type="nodes",node="Geometric Aggregation (w1=0.25,w2=0.75)")
plot(GeoAggregation3,type="nodes",node="Geometric Offset Agg. (w1=w2=0.5,d=0.5)")
plot(GeoAggregation4,type="nodes",node="Geometric Offset Agg. (w1=0.25,w2=0.75,d=0.5)")
dev.off()


# revgeo/revgeooff:
# -----------------

RevGeoAggregation1 <- utility.aggregation.create(
  name.node = "Rev. Geometric Aggregation (w1=w2=0.5)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.revgeo", 
  par       = c(1,1))
RevGeoAggregation2 <- utility.aggregation.create(
  name.node = "Rev. Geometric Aggregation (w1=0.25,w2=0.75)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.revgeo", 
  par       = c(1,3))
RevGeoAggregation3 <- utility.aggregation.create(
  name.node = "Rev. Geometric Offset Agg. (w1=w2=0.5,d=0.5)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.revgeooff", 
  par       = c(1,1,0.5))
RevGeoAggregation4 <- utility.aggregation.create(
  name.node = "Rev. Geometric Offset Agg. (w1=0.25,w2=0.75,d=0.5)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.revgeooff", 
  par       = c(1,3,0.5))

pdf(paste(dir.fig,"aggregationrevgeo.pdf",sep="/"),width=pdf.width2,height=pdf.height2)
par(mfrow=c(2,2),mar=pdf.mar)
plot(RevGeoAggregation1,type="nodes",node="Rev. Geometric Aggregation (w1=w2=0.5)")
plot(RevGeoAggregation2,type="nodes",node="Rev. Geometric Aggregation (w1=0.25,w2=0.75)")
plot(RevGeoAggregation3,type="nodes",node="Rev. Geometric Offset Agg. (w1=w2=0.5,d=0.5)")
plot(RevGeoAggregation4,type="nodes",node="Rev. Geometric Offset Agg. (w1=0.25,w2=0.75,d=0.5)")
dev.off()

png(paste(dir.fig,"aggregationrevgeo.png",sep="/"),width=png.width2,height=png.height2)
par(mfrow=c(2,2),mar=png.mar)
plot(RevGeoAggregation1,type="nodes",node="Rev. Geometric Aggregation (w1=w2=0.5)")
plot(RevGeoAggregation2,type="nodes",node="Rev. Geometric Aggregation (w1=0.25,w2=0.75)")
plot(RevGeoAggregation3,type="nodes",node="Rev. Geometric Offset Agg. (w1=w2=0.5,d=0.5)")
plot(RevGeoAggregation4,type="nodes",node="Rev. Geometric Offset Agg. (w1=0.25,w2=0.75,d=0.5)")
dev.off()


# harmo/harmooff:
# ---------------

HarmoAggregation1 <- utility.aggregation.create(
  name.node = "Harmonic Aggregation (w1=w2=0.5)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.harmo", 
  par       = c(1,1))
HarmoAggregation2 <- utility.aggregation.create(
  name.node = "Harmonic Aggregation (w1=0.25,w2=0.75)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.harmo", 
  par       = c(1,3))
HarmoAggregation3 <- utility.aggregation.create(
  name.node = "Harmonic Offset Agg. (w1=w2=0.5,d=0.5)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.harmooff", 
  par       = c(1,1,0.5))
HarmoAggregation4 <- utility.aggregation.create(
  name.node = "Harmonic Offset Agg. (w1=0.25,w2=0.75,d=0.5)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.harmooff", 
  par       = c(1,3,0.5))

pdf(paste(dir.fig,"aggregationharmo.pdf",sep="/"),width=pdf.width2,height=pdf.height2)
par(mfrow=c(2,2),mar=pdf.mar)
plot(HarmoAggregation1,type="nodes",node="Harmonic Aggregation (w1=w2=0.5)")
plot(HarmoAggregation2,type="nodes",node="Harmonic Aggregation (w1=0.25,w2=0.75)")
plot(HarmoAggregation3,type="nodes",node="Harmonic Offset Agg. (w1=w2=0.5,d=0.5)")
plot(HarmoAggregation4,type="nodes",node="Harmonic Offset Agg. (w1=0.25,w2=0.75,d=0.5)")
dev.off()

png(paste(dir.fig,"aggregationharmo.png",sep="/"),width=png.width2,height=png.height2)
par(mfrow=c(2,2),mar=png.mar)
plot(HarmoAggregation1,type="nodes",node="Harmonic Aggregation (w1=w2=0.5)")
plot(HarmoAggregation2,type="nodes",node="Harmonic Aggregation (w1=0.25,w2=0.75)")
plot(HarmoAggregation3,type="nodes",node="Harmonic Offset Agg. (w1=w2=0.5,d=0.5)")
plot(HarmoAggregation4,type="nodes",node="Harmonic Offset Agg. (w1=0.25,w2=0.75,d=0.5)")
dev.off()


# revharmo/revharmooff:
# ---------------------

RevHarmoAggregation1 <- utility.aggregation.create(
  name.node = "Rev. Harmonic Aggregation (w1=w2=0.5)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.revharmo", 
  par       = c(1,1))
RevHarmoAggregation2 <- utility.aggregation.create(
  name.node = "Rev. Harmonic Aggregation (w1=0.25,w2=0.75)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.revharmo", 
  par       = c(1,3))
RevHarmoAggregation3 <- utility.aggregation.create(
  name.node = "Rev. Harmonic Offset Agg. (w1=w2=0.5,d=0.5)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.revharmooff", 
  par       = c(1,1,0.5))
RevHarmoAggregation4 <- utility.aggregation.create(
  name.node = "Rev. Harmonic Offset Agg. (w1=0.25,w2=0.75,d=0.5)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.revharmooff", 
  par       = c(1,3,0.5))

pdf(paste(dir.fig,"aggregationrevharmo.pdf",sep="/"),width=pdf.width2,height=pdf.height2)
par(mfrow=c(2,2),mar=pdf.mar)
plot(RevHarmoAggregation1,type="nodes",node="Rev. Harmonic Aggregation (w1=w2=0.5)")
plot(RevHarmoAggregation2,type="nodes",node="Rev. Harmonic Aggregation (w1=0.25,w2=0.75)")
plot(RevHarmoAggregation3,type="nodes",node="Rev. Harmonic Offset Agg. (w1=w2=0.5,d=0.5)")
plot(RevHarmoAggregation4,type="nodes",node="Rev. Harmonic Offset Agg. (w1=0.25,w2=0.75,d=0.5)")
dev.off()

png(paste(dir.fig,"aggregationrevharmo.png",sep="/"),width=png.width2,height=png.height2)
par(mfrow=c(2,2),mar=png.mar)
plot(RevHarmoAggregation1,type="nodes",node="Rev. Harmonic Aggregation (w1=w2=0.5)")
plot(RevHarmoAggregation2,type="nodes",node="Rev. Harmonic Aggregation (w1=0.25,w2=0.75)")
plot(RevHarmoAggregation3,type="nodes",node="Rev. Harmonic Offset Agg. (w1=w2=0.5,d=0.5)")
plot(RevHarmoAggregation4,type="nodes",node="Rev. Harmonic Offset Agg. (w1=0.25,w2=0.75,d=0.5)")
dev.off()


# mix:
# ----

MixAggregation1 <- utility.aggregation.create(
  name.node = "Mixture Agg. (w1=w2=0.5,wa=1,wm=0,wg=0)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.mix", 
  par       = c(1,1,1,0,0))
MixAggregation2 <- utility.aggregation.create(
  name.node = "Mix. Agg. (w1=w2=0.5,wa=0.33,wm=0.33,wg=0.33)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.mix", 
  par       = c(1,1,1,1,1))
MixAggregation3 <- utility.aggregation.create(
  name.node = "Mixture Agg. (w1=w2=0.5,wa=0.5,wm=0,wg=0.5)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.mix", 
  par       = c(1,1,1,0,1))
MixAggregation4 <- utility.aggregation.create(
  name.node = "Mix. (w1=0.25,w2=0.75,wa=0.33,wm=0.33,wg=0.33)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.mix", 
  par       = c(1,3,1,1,1))

pdf(paste(dir.fig,"aggregationmix.pdf",sep="/"),width=pdf.width2,height=pdf.height2)
par(mfrow=c(2,2),mar=pdf.mar)
plot(MixAggregation1,type="nodes",node="Mixture Agg. (w1=w2=0.5,wa=1,wm=0,wg=0)")
plot(MixAggregation2,type="nodes",node="Mix. Agg. (w1=w2=0.5,wa=0.33,wm=0.33,wg=0.33)")
plot(MixAggregation3,type="nodes",node="Mixture Agg. (w1=w2=0.5,wa=0.5,wm=0,wg=0.5)")
plot(MixAggregation4,type="nodes",node="Mix. (w1=0.25,w2=0.75,wa=0.33,wm=0.33,wg=0.33)")
dev.off()

png(paste(dir.fig,"aggregationmix.png",sep="/"),width=png.width2,height=png.height2)
par(mfrow=c(2,2),mar=png.mar)
plot(MixAggregation1,type="nodes",node="Mixture Agg. (w1=w2=0.5,wa=1,wm=0,wg=0)")
plot(MixAggregation2,type="nodes",node="Mix. Agg. (w1=w2=0.5,wa=0.33,wm=0.33,wg=0.33)")
plot(MixAggregation3,type="nodes",node="Mixture Agg. (w1=w2=0.5,wa=0.5,wm=0,wg=0.5)")
plot(MixAggregation4,type="nodes",node="Mix. (w1=0.25,w2=0.75,wa=0.33,wm=0.33,wg=0.33)")
dev.off()


# addmin:
# -------

AddMinAggregation1 <- utility.aggregation.create(
  name.node = "Additive-Minimum Aggregation (w1=w2=0.5,a=0.9)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.addmin", 
  par       = c(1,1,0.9))
AddMinAggregation2 <- utility.aggregation.create(
  name.node = "Additive-Minimum Aggregation (w1=w2=0.5,a=0.5)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.addmin", 
  par       = c(1,1,0.5))
AddMinAggregation3 <- utility.aggregation.create(
  name.node = "Additive-Minimum Aggregation (w1=w2=0.5,a=0.1)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.addmin", 
  par       = c(1,1,0.1))
AddMinAggregation4 <- utility.aggregation.create(
  name.node = "Additive-Minimum Agg. (w1=0.25,w2=0.75,a=0.5)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.addmin", 
  par       = c(1,3,0.5))

pdf(paste(dir.fig,"aggregationaddmin.pdf",sep="/"),width=pdf.width2,height=pdf.height2)
par(mfrow=c(2,2),mar=pdf.mar)
plot(AddMinAggregation1,type="nodes",node="Additive-Minimum Aggregation (w1=w2=0.5,a=0.9)")
plot(AddMinAggregation2,type="nodes",node="Additive-Minimum Aggregation (w1=w2=0.5,a=0.5)")
plot(AddMinAggregation3,type="nodes",node="Additive-Minimum Aggregation (w1=w2=0.5,a=0.1)")
plot(AddMinAggregation4,type="nodes",node="Additive-Minimum Agg. (w1=0.25,w2=0.75,a=0.5)")
dev.off()

png(paste(dir.fig,"aggregationaddmin.png",sep="/"),width=png.width2,height=png.height2)
par(mfrow=c(2,2),mar=png.mar)
plot(AddMinAggregation1,type="nodes",node="Additive-Minimum Aggregation (w1=w2=0.5,a=0.9)")
plot(AddMinAggregation2,type="nodes",node="Additive-Minimum Aggregation (w1=w2=0.5,a=0.5)")
plot(AddMinAggregation3,type="nodes",node="Additive-Minimum Aggregation (w1=w2=0.5,a=0.1)")
plot(AddMinAggregation4,type="nodes",node="Additive-Minimum Agg. (w1=0.25,w2=0.75,a=0.5)")
dev.off()


# addpower:
# ---------

AddPowerAggregation1 <- utility.aggregation.create(
  name.node = "Additive Power Aggregation (w1=w2=0.5,a=2)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.addpower", 
  par       = c(1,1,2))
AddPowerAggregation2 <- utility.aggregation.create(
  name.node = "Additive Power Aggregation (w1=w2=0.5,a=0.5)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.addpower", 
  par       = c(1,1,0.5))
AddPowerAggregation3 <- utility.aggregation.create(
  name.node = "Add. Power Aggregation (w1=0.25,w2=0.75,a=2)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.addpower", 
  par       = c(1,3,0.5))
AddPowerAggregation4 <- utility.aggregation.create(
  name.node = "Add. Power Aggregation (w1=0.25,w2=0.75,a=0.5)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.addpower", 
  par       = c(1,3,0.5))

pdf(paste(dir.fig,"aggregationaddpower.pdf",sep="/"),width=pdf.width2,height=pdf.height2)
par(mfrow=c(2,2),mar=pdf.mar)
plot(AddPowerAggregation1,type="nodes",node="Additive Power Aggregation (w1=w2=0.5,a=2)")
plot(AddPowerAggregation2,type="nodes",node="Additive Power Aggregation (w1=w2=0.5,a=0.5)")
plot(AddPowerAggregation3,type="nodes",node="Add. Power Aggregation (w1=0.25,w2=0.75,a=2)")
plot(AddPowerAggregation4,type="nodes",node="Add. Power Aggregation (w1=0.25,w2=0.75,a=0.5)")
dev.off()

png(paste(dir.fig,"aggregationaddpower.png",sep="/"),width=png.width2,height=png.height2)
par(mfrow=c(2,2),mar=png.mar)
plot(AddPowerAggregation1,type="nodes",node="Additive Power Aggregation (w1=w2=0.5,a=2)")
plot(AddPowerAggregation2,type="nodes",node="Additive Power Aggregation (w1=w2=0.5,a=0.5)")
plot(AddPowerAggregation3,type="nodes",node="Add. Power Aggregation (w1=0.25,w2=0.75,a=2)")
plot(AddPowerAggregation4,type="nodes",node="Add. Power Aggregation (w1=0.25,w2=0.75,a=0.5)")
dev.off()


# revaddpower:
# ------------

RevAddPowerAggregation1 <- utility.aggregation.create(
  name.node = "Rev. Add. Power Aggregation (w1=w2=0.5,a=2)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.revaddpower", 
  par       = c(1,1,2))
RevAddPowerAggregation2 <- utility.aggregation.create(
  name.node = "Rev. Add. Power Aggregation (w1=w2=0.5,a=0.5)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.revaddpower", 
  par       = c(1,1,0.5))
RevAddPowerAggregation3 <- utility.aggregation.create(
  name.node = "Rev. Add. Power Agg. (w1=0.25,w2=0.75,a=2)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.revaddpower", 
  par       = c(1,3,0.5))
RevAddPowerAggregation4 <- utility.aggregation.create(
  name.node = "Rev. Add. Power Agg. (w1=0.25,w2=0.75,a=0.5)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.revaddpower", 
  par       = c(1,3,0.5))

pdf(paste(dir.fig,"aggregationrevaddpower.pdf",sep="/"),width=pdf.width2,height=pdf.height2)
par(mfrow=c(2,2),mar=pdf.mar)
plot(RevAddPowerAggregation1,type="nodes",node="Rev. Add. Power Aggregation (w1=w2=0.5,a=2)")
plot(RevAddPowerAggregation2,type="nodes",node="Rev. Add. Power Aggregation (w1=w2=0.5,a=0.5)")
plot(RevAddPowerAggregation3,type="nodes",node="Rev. Add. Power Agg. (w1=0.25,w2=0.75,a=2)")
plot(RevAddPowerAggregation4,type="nodes",node="Rev. Add. Power Agg. (w1=0.25,w2=0.75,a=0.5)")
dev.off()

png(paste(dir.fig,"aggregationrevaddpower.png",sep="/"),width=png.width2,height=png.height2)
par(mfrow=c(2,2),mar=png.mar)
plot(RevAddPowerAggregation1,type="nodes",node="Rev. Add. Power Aggregation (w1=w2=0.5,a=2)")
plot(RevAddPowerAggregation2,type="nodes",node="Rev. Add. Power Aggregation (w1=w2=0.5,a=0.5)")
plot(RevAddPowerAggregation3,type="nodes",node="Rev. Add. Power Agg. (w1=0.25,w2=0.75,a=2)")
plot(RevAddPowerAggregation4,type="nodes",node="Rev. Add. Power Agg. (w1=0.25,w2=0.75,a=0.5)")
dev.off()


# addsplitpower:
# --------------

AddSplitPowerAggregation1 <- utility.aggregation.create(
  name.node = "Add. Split Power Aggregation (w1=w2=0.5,a=2)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.addsplitpower", 
  par       = c(1,1,2,0.5))
AddSplitPowerAggregation2 <- utility.aggregation.create(
  name.node = "Add. Split Power Aggregation (w1=w2=0.5,a=0.5)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.addsplitpower", 
  par       = c(1,1,0.5,0.5))
AddSplitPowerAggregation3 <- utility.aggregation.create(
  name.node = "Add. Split Power Agg. (w1=0.25,w2=0.75,a=2)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.addsplitpower", 
  par       = c(1,3,2,0.5))
AddSplitPowerAggregation4 <- utility.aggregation.create(
  name.node = "Add. Split Power Agg. (w1=0.25,w2=0.75,a=0.5)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.addsplitpower", 
  par       = c(1,3,0.5,0.5))

pdf(paste(dir.fig,"aggregationaddsplitpower.pdf",sep="/"),width=pdf.width2,height=pdf.height2)
par(mfrow=c(2,2),mar=pdf.mar)
plot(AddSplitPowerAggregation1,type="nodes",node="Add. Split Power Aggregation (w1=w2=0.5,a=2)")
plot(AddSplitPowerAggregation2,type="nodes",node="Add. Split Power Aggregation (w1=w2=0.5,a=0.5)")
plot(AddSplitPowerAggregation3,type="nodes",node="Add. Split Power Agg. (w1=0.25,w2=0.75,a=2)")
plot(AddSplitPowerAggregation4,type="nodes",node="Add. Split Power Agg. (w1=0.25,w2=0.75,a=0.5)")
dev.off()

png(paste(dir.fig,"aggregationaddsplitpower.png",sep="/"),width=png.width2,height=png.height2)
par(mfrow=c(2,2),mar=png.mar)
plot(AddSplitPowerAggregation1,type="nodes",node="Add. Split Power Aggregation (w1=w2=0.5,a=2)")
plot(AddSplitPowerAggregation2,type="nodes",node="Add. Split Power Aggregation (w1=w2=0.5,a=0.5)")
plot(AddSplitPowerAggregation3,type="nodes",node="Add. Split Power Agg. (w1=0.25,w2=0.75,a=2)")
plot(AddSplitPowerAggregation4,type="nodes",node="Add. Split Power Agg. (w1=0.25,w2=0.75,a=0.5)")
dev.off()


# revaddsplitpower:
# -----------------

RevAddSplitPowerAggregation1 <- utility.aggregation.create(
  name.node = "Rev. Add. Split Power Agg. (w1=w2=0.5,a=2)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.revaddsplitpower", 
  par       = c(1,1,2,0.5))
RevAddSplitPowerAggregation2 <- utility.aggregation.create(
  name.node = "Rev. Add. Split Power Agg. (w1=w2=0.5,a=0.5)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.revaddsplitpower", 
  par       = c(1,1,0.5,0.5))
RevAddSplitPowerAggregation3 <- utility.aggregation.create(
  name.node = "Rev. Add. Split Pow. Agg. (w1=0.25,w2=0.75,a=2)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.revaddsplitpower", 
  par       = c(1,3,2,0.5))
RevAddSplitPowerAggregation4 <- utility.aggregation.create(
  name.node = "Rev. Add. Split Pow. Agg. (w1=0.25,w2=0.75,a=0.5)", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.revaddsplitpower", 
  par       = c(1,3,0.5,0.5))

pdf(paste(dir.fig,"aggregationrevaddsplitpower.pdf",sep="/"),width=pdf.width2,height=pdf.height2)
par(mfrow=c(2,2),mar=pdf.mar)
plot(RevAddSplitPowerAggregation1,type="nodes",node="Rev. Add. Split Power Agg. (w1=w2=0.5,a=2)")
plot(RevAddSplitPowerAggregation2,type="nodes",node="Rev. Add. Split Power Agg. (w1=w2=0.5,a=0.5)")
plot(RevAddSplitPowerAggregation3,type="nodes",node="Rev. Add. Split Pow. Agg. (w1=0.25,w2=0.75,a=2)")
plot(RevAddSplitPowerAggregation4,type="nodes",node="Rev. Add. Split Pow. Agg. (w1=0.25,w2=0.75,a=0.5)")
dev.off()

png(paste(dir.fig,"aggregationrevaddsplitpower.png",sep="/"),width=png.width2,height=png.height2)
par(mfrow=c(2,2),mar=png.mar)
plot(RevAddSplitPowerAggregation1,type="nodes",node="Rev. Add. Split Power Agg. (w1=w2=0.5,a=2)")
plot(RevAddSplitPowerAggregation2,type="nodes",node="Rev. Add. Split Power Agg. (w1=w2=0.5,a=0.5)")
plot(RevAddSplitPowerAggregation3,type="nodes",node="Rev. Add. Split Pow. Agg. (w1=0.25,w2=0.75,a=2)")
plot(RevAddSplitPowerAggregation4,type="nodes",node="Rev. Add. Split Pow. Agg. (w1=0.25,w2=0.75,a=0.5)")
dev.off()


# bonusmalus:
# -----------

BonusMalusAggregation1 <- utility.aggregation.create(
  name.node = "Bonus-Malus Aggregation (par=c(1,NA,1))", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.bonusmalus", 
  par       = c(1,NA,1))
BonusMalusAggregation2 <- utility.aggregation.create(
  name.node = "Bonus-Malus Aggregation (par=c(1,1,NA))", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.bonusmalus", 
  par       = c(1,1,NA))
BonusMalusAggregation3 <- utility.aggregation.create(
  name.node = "Bonus-Malus Aggregation (par=c(1,NA,-1))", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.bonusmalus", 
  par       = c(1,NA,-1))
BonusMalusAggregation4 <- utility.aggregation.create(
  name.node = "Bonus-Malus Aggregation (par=c(1,NA,2))", 
  nodes     = list(obj1,obj2),
  name.fun  = "utility.aggregate.bonusmalus", 
  par       = c(1,NA,2))

pdf(paste(dir.fig,"aggregationbonusmalus.pdf",sep="/"),width=pdf.width2,height=pdf.height2)
par(mfrow=c(2,2),mar=pdf.mar)
plot(BonusMalusAggregation1,type="nodes",node="Bonus-Malus Aggregation (par=c(1,NA,1))")
plot(BonusMalusAggregation2,type="nodes",node="Bonus-Malus Aggregation (par=c(1,1,NA))")
plot(BonusMalusAggregation3,type="nodes",node="Bonus-Malus Aggregation (par=c(1,NA,-1))")
plot(BonusMalusAggregation4,type="nodes",node="Bonus-Malus Aggregation (par=c(1,NA,2))")
dev.off()

png(paste(dir.fig,"aggregationbonusmalus.png",sep="/"),width=png.width2,height=png.height2)
par(mfrow=c(2,2),mar=png.mar)
plot(BonusMalusAggregation1,type="nodes",node="Bonus-Malus Aggregation (par=c(1,NA,1))")
plot(BonusMalusAggregation2,type="nodes",node="Bonus-Malus Aggregation (par=c(1,1,NA))")
plot(BonusMalusAggregation3,type="nodes",node="Bonus-Malus Aggregation (par=c(1,NA,-1))")
plot(BonusMalusAggregation4,type="nodes",node="Bonus-Malus Aggregation (par=c(1,NA,2))")
dev.off()

