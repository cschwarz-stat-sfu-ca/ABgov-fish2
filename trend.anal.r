# analysis of  BLTR, RNTR, and BKTR for a trend analysis

library(lmerTest)
library(ggforce)
library(reshape2)

doParallel <- TRUE  # should I set up parallel processing of the emodels.
if(doParallel) {
  library(doMC)  # for parallel model fitting
  #library(foreach)

  # see http://viktoriawagner.weebly.com/blog/five-steps-to-parallel-computing-in-r
  detectCores() 
  cl <- makeCluster(4)
  # May need to export some libraries to the cluster
  # see http://stackoverflow.com/questions/18981932/logging-with-plyr-in-parallel-true
  # clusterEvalQ(cl, library(unmarked))
  registerDoMC(4) 
} 


source("read.data.R")
source("slr-power-stroup.r")

# BLTR, BKTR, and RNTR analysis
age.classes.csv <- textConnection(
  "Species, AgeClass, LBound, UBound
  BLTR, Immature, 0, 150
  BLTR, Mature, 150, Inf
  RNTR, Immature, 0, 142
  RNTR, Mature, 142, Inf
  BKTR, AllAges, 0,  Inf")   # BKTR not divided because of sparseness

age.classes <- read.csv(age.classes.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)

# select only species of interest above
xtabs(~Species, data=fish.data, exclude=NULL, na.action=na.pass)
fish.data <- fish.data[ fish.data$Species %in% age.classes$Species,]
xtabs(~Species, data=fish.data, exclude=NULL, na.action=na.pass)


# Divide fish by age class
fish.data <- plyr::ddply(fish.data, "Species", function(x, age.classes){
   # match the age class with species of interest
   age.classes <- age.classes[ age.classes$Species == x$Species[1],]
   age.class <- apply(outer(x$ForkLength, age.classes$LBound, ">=")  &
                      outer(x$ForkLength, age.classes$UBound, "<" ), 1, which)
   x$AgeClass <- age.classes$AgeClass[age.class]
   x
  
}, age.classes=age.classes)

xtabs(~Species+AgeClass, data=fish.data, exclude=NULL, na.action=na.pass)

# Check if actual number of fish differs from nominal number of fish.
species.check.n <- plyr::ddply(fish.data, 
                      c("Watershed","Stream","Year","Site","Species"), plyr::summarize,
                      nominal.fish= Captured_C[1],
                      actual.fish  = length(Species))
head(species.check.n)
species.check.n[ species.check.n$nominal.fish != species.check.n$actual.fish,]


# Count number of fish in each age class
species.n <- plyr::ddply(fish.data, 
                      c("Watershed","Stream","Year","Site","Species","AgeClass"), plyr::summarize,
                      n.fish=length(Species),
                      distance=min(Distance))
head(species.n)

# We need to add back Stream-Year-Sites where 0 fish were captured
# Note that RNTR only occurs in the Athbathsak/NSask watersheds and in Mackenzie and Moon streams
# so we need to remove all other watersheds/streams for this species.
dim(species.n)
xtabs(~Stream+Year, data=species.n, exclude=NULL, na.action=na.pass)

# expand the site table for each species of interest
dim(site.table)
expanded.site.table <- plyr::ddply(age.classes[,c("Species","AgeClass")], 
                                                c("Species","AgeClass"), function(SA, site.table){
   site.table$Species <- SA$Species
   site.table$AgeClass<- SA$AgeClass
   site.table
}, site.table=site.table)
dim(expanded.site.table)

expanded.site.table.complete <- expanded.site.table # has all streams

xtabs(~interaction(Watershed,Stream,drop=TRUE)+Species, data=expanded.site.table)
# Remove RNTR from all but MacKenzie and Moon streams
select <- (!expanded.site.table$Stream %in% c("Mackenzie","Moon")) & 
            expanded.site.table$Species %in% c("RNTR")
expanded.site.table[select,c("Watershed","Stream","Species")]
expanded.site.table <- expanded.site.table[!select,]
xtabs(~interaction(Watershed,Stream,drop=TRUE)+Species, data=expanded.site.table)




dim(species.n)
species.n <- merge(species.n, cbind(expanded.site.table[,c("Watershed","Stream","Year","Site","Species","AgeClass","siteDistance")]),
                   all.y=TRUE)
dim(species.n)
xtabs(~Stream+Year, data=species.n, exclude=NULL, na.action=na.pass)
xtabs(~interaction(Watershed,Stream, drop=TRUE)+Species, data=species.n)

# replace all Distance = NA by distance=0 
species.n$distance[ is.na(species.n$distance)] <- species.n$siteDistance[ is.na(species.n$distance)]
# replace all NA  for the fish counts by 0
species.n$n.fish[ is.na(species.n$n.fish)] <- 0



# compute and plot the CUE (fish/300 m)
species.n$CUE <- species.n$n.fish  /species.n$distance * 300

# Figure out the 1/2 of smallest non-zero CUE to use an offset
offset <- plyr::ddply(species.n, c("Species","AgeClass"), plyr::summarize,
                      offset=0.5*min(CUE[CUE>0]))
offset

species.n <- merge(species.n, offset)

# We create plotdata so that all streams are shown even if no data are present
# For example, for RNTR, fish only occur in Mackenzie and Moon, but we want
# panels for all streams
plotdata <- species.n
dim(plotdata)
plotdata <- merge(plotdata, expanded.site.table.complete, all.y=TRUE)
dim(plotdata)


plyr::d_ply(plotdata, c("Species","AgeClass"), function (x){
  plot1 <- ggplot(data=x, aes(x=Year, y=CUE))+
    ggtitle(paste("Raw CUE/300 m for ", x$Species[1]," ", x$AgeClass[1],sep=""))+
    geom_point(position=position_jitter(h=.05))+
    facet_wrap(~Stream, scales="free_y")+
    scale_x_continuous(breaks=seq(1900,2100,2))
  plot(plot1)

  ggsave(plot1, 
         file=file.path("Report",paste('raw-CUE-',x$Species[1],'-',x$AgeClass[1],'.png',sep="")),
         h=6, w=6, units="in", dpi=300)
})


plyr::d_ply(plotdata, c("Species","AgeClass"), function (x){
  plot1 <- ggplot(data=x, aes(x=Year, y=log(CUE+offset)))+
    ggtitle(paste("Raw log( CUE/300 m) for ", x$Species[1], " ", x$AgeClass[1],sep=""))+
    geom_point(position=position_jitter(h=.05))+
    facet_wrap(~Stream, scales="free_y", ncol=3)+
    scale_x_continuous(breaks=seq(1900,2100,2))

  plot(plot1)
  ggsave(plot1, 
         file=file.path("Report",paste('log-CUE-',x$Species[1],'-',x$AgeClass[1],'.png',sep="")),
        h=6, w=6, units="in", dpi=300)
})


# What is average CUE for each type of fish
# Many of these differ from the report because JR included ALL fish. 
mean.cue <- plyr::ddply(species.n, c("Watershed","Stream","Species","Year","AgeClass"), plyr::summarize,
            n.sites          = length(CUE),
            mean.CUE= mean(CUE),
            se.CUE = sd(CUE)/sqrt(n.sites))
mean.cue <- merge(mean.cue, offset)
head(mean.cue)
mean.cue$log.mean.cue    <- log(mean.cue$mean.CUE+mean.cue$offset)
mean.cue$log.mean.cue.se <- mean.cue$se.CUE / mean.cue$mean.CUE

mean.cue <- mean.cue[ order(mean.cue$AgeClass,mean.cue$Watershed, mean.cue$Stream),]
mean.cue

# again, we want all of the plots to look the same
plotdata <- mean.cue
dim(plotdata)
plotdata <- merge(plotdata, 
     unique(expanded.site.table.complete[,c("Species","AgeClass","Watershed","Stream","Year")]),
     all.y=TRUE)
dim(plotdata)


plyr::d_ply(plotdata, c("Species"), function (x){
  plot.mean.cue <- ggplot(data=x, aes(x=Year, y=log(mean.CUE), color=Stream))+
   ggtitle(paste("Mean CUE for ",x$Species[1],sep=""))+
   geom_point()+
   geom_line(aes(group=Stream))+
   ylab("log(Mean CUE (fish/300 m))")+
   facet_grid(Watershed~AgeClass, scales="free_y")+
   scale_x_continuous(breaks=seq(1900,2100,2))
 plot(plot.mean.cue)
 ggsave(plot.mean.cue,
       file=file.path("Report",paste("mean.cue.",x$Species[1],".png")),
       h=6, w=6, units="in", dpi=300)
})


# Estimate sources of variation on log(CUE) for use in Regression power analysis

# are there any sites measured in more than one year?
n.site.year <- plyr::ddply(unique(species.n[, c("Species","Stream","Site","Year")]),
                                  c("Species","Stream","Site"), plyr::summarize,
                                  n.year = length(Year)
                           )
n.site.year[n.site.year$n.year > 1,]


# Is every stream measured in more than one year
xtabs(~Stream+Year, data=species.n, exclude=NULL, na.action=na.pass)
xtabs(~interaction(Stream,Watershed, drop=TRUE)+Year+Species, data=species.n, exclude=NULL, na.action=na.pass)





fits <- plyr::dlply(species.n, c("Watershed","Species","AgeClass"), function (x){
  # Fit a mixed linear model to log(CUE + offset.)
  cat("\n\n\n*** Starting ", as.character(x$Watershed[1]), x$Species[1], as.character(x$AgeClass[1]), "\n")
  x$SiteF <- factor(x$Site)
  x$YearF <- factor(x$Year)
  x$StreamF<-factor(x$Stream)
  x$Year  <- x$Year - 2000
  #if(x$Watershed[1]=="Ram")browser()
  if(length(unique(x$Year)) >1  & length(unique(x$Stream)) >1 )fit <- lmerTest::lmer(log(CUE+offset) ~ Year + 
                                                                                       (1|YearF) + (1|StreamF)+(1|StreamF:YearF), data=x)
  if(length(unique(x$Year))==1  & length(unique(x$Stream)) >1 )fit <- lmerTest::lmer(log(CUE+offset) ~  (1|StreamF), data=x)
  if(length(unique(x$Year))==1  & length(unique(x$Stream))==1 )fit <- lm(log(CUE+offset) ~  1, data=x)
  if(length(unique(x$Year))==2  & length(unique(x$Stream))==1 )fit <- lmerTest::lmer(log(CUE+offset) ~ 1 + (1|YearF), data=x)
  list(Watershed=x$Watershed[1],
       Species  =x$Species  [1],
       AgeClass =x$AgeClass [1],
       fit      =fit)  
  
  
})

# Extract the variance components
vc <- plyr::ldply(fits, function(x){
  sd.Year   <- NA
  sd.Stream <- NA
  sd.Resid  <- NA
  sd.StreamYear<-NA
  if(class(x$fit)=="lmerModLmerTest"){
    vc <- as.data.frame(VarCorr(x$fit))
    if(any(vc$grp=="YearF"        ))sd.Year      <- vc[vc$grp=="YearF"   ,"sdcor"]
    if(any(vc$grp=="StreamF"      ))sd.Stream    <- vc[vc$grp=="StreamF" ,"sdcor"]
    if(any(vc$grp=="StreamF:YearF"))sd.StreamYear<- vc[vc$grp=="StreamF:YearF", "sdcor"]
    if(any(vc$grp=="Residual"     ))sd.Resid     <- vc[vc$grp=="Residual","sdcor"]
  }
  if(class(x$fit)=="lm"){
    sd.Resid  <- summary(x$fit)$sigma
  }
  #browser()
  
  data.frame(Watershed=x$Watershed,
             Species  =x$Species,
             AgeClass =x$AgeClass,
             sd.Year  =sd.Year,
             sd.Stream=sd.Stream,
             sd.StreamYear=sd.StreamYear,
             sd.Resid =sd.Resid, 
             stringsAsFactors=FALSE)
  
  
  
})

vc <- vc[ order(vc$Species, vc$AgeClass, vc$Watershed),]
vc

temp <- vc
temp[, 4:7] <- round(temp[,4:7],3) 
temp
write.csv(temp,
          file=file.path("Report","vc.csv"), row.names=FALSE)

  




# find the mean VC and use it to impute any missing values
vc2 <- plyr::ddply(vc, c("Species","AgeClass"), function(x){
  mean.sd.Year      <- mean(x$sd.Year,       na.rm=TRUE)
  mean.sd.Stream    <- mean(x$sd.Stream,     na.rm=TRUE)
  mean.sd.StreamYear<- mean(x$sd.StreamYear, na.rm=TRUE)
  mean.sd.Resid     <- mean(x$sd.Resid     , na.rm=TRUE)
  x$sd.Year      [ is.na(x$sd.Year)]       <- mean.sd.Year
  x$sd.Stream    [ is.na(x$sd.Stream) ]    <- mean.sd.Stream
  x$sd.StreamYear[ is.na(x$sd.StreamYear)] <- mean.sd.StreamYear
  x$sd.Resid     [ is.na(x$sd.Resid)]      <- mean.sd.Resid
  x
})
vc2


##############################################################################
##############################################################################
##############################################################################
##############################################################################
# Power analysis for trends

source("slr-power-stroup.r")

# The process error is generally small. I will use the average from the streams 
# where the process error is present

scenario <- expand.grid(
            trend= c(seq(0.02,.10,.02),.15, .25, .50, 1), 
            n.years=c(5,7,9,11,13,15,17,20),
            n.sites.per.year=c(5,10,15),
            alpha=c(0.05, 0.10, 0.15),
            stringsAsFactors=FALSE)
scenarios <- merge(vc2, scenario)  # get every combination


power <- plyr::adply(scenarios, 1, function (x){
  
  power <- slr.power.stroup( Trend=x$trend,
                             Xvalues=rep(1:x$n.years, x$n.sites.per.year),
                        Process.SD=x$sd.Year, 
                        Sampling.SD=x$sd.Resid,          alpha=x$alpha)
  x <- cbind(x, power)  
  x
  
}, .parallel=doParallel)


# graph the power curves
pdf(file.path("Report",paste("Figure07-power-plot-trend.pdf")), h=6, w=6)
plyr::d_ply(power, c("trend","alpha","Species"),function(x){
   cat("Plotting ", x$trend[1], x$alpha[1], x$Species[1], "\n")
   power.plot <- ggplot(data=x, aes(x=n.years, y=os.power1,  color=as.factor(n.sites.per.year)))+
      ggtitle(paste("TREND power to detect a ", 100*x$trend,'% increase in mean CUE/year',
                    "\n",x$Species[1],
                    "; alpha= ",formatC(x$alpha[1],format="f", digits=2),sep=""))+
      #geom_point()+
      geom_line()+
      facet_grid(Watershed~AgeClass)+
      xlab("Number of years after restoration")+ylab("Power to detect trend in mean CUE")+ylim(0,1)+
      geom_hline(yintercept=0.80)+
      scale_color_discrete(name="Sites\nper\nyear")+
      scale_x_continuous(breaks=seq(0,100,2))
   plot(power.plot) 
  # ggsave(power.plot, 
  #        file=file.path("Report",paste("power-plot-BACI-",
  #                                      exp(x$baci[1]),"x effect;",
  #                                      x$n.years.b[1],"-years-before-",
  #                                      x$Species[1],"-alpha-",formatC(100*x$alpha[1],width=2, flag="0"),".png",sep="")),
  #        h=6,w=6, units="in", dpi=300)
  
})
dev.off()

