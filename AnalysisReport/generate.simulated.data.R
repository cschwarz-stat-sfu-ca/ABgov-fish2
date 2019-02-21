# Generate simulated data for use in the report showing how to analyze the data.
# The simulated data sets will use the 
#    - key species BLTR, RNTR, and BKTR
#      two age classes - immature, mature, both
# and will also include imputed 0 for species not seen at a year-site-visit.
# The existing data will only be used to estimate the variance components
# and simulated data will be used to "replace" the existing data


# Get the existing data
#

source(file.path("..",'read.data.r'), chdir=TRUE)

library(arm)
library(lmerTest)


# This doesn't have the 0 records and we must impute them

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


# Check now that we have all of the data we need
dim(species.n)
names(species.n)
head(species.n)
xtabs(~interaction(Stream,Site,drop=TRUE)+Year, data=species.n, exclude=NULL, na.action=na.pass)

species.n[ species.n$Stream == "Fall",]

# All of the years in the data base are "Before"
species.n$BA <- "Before"

# Specify the Treated or Control Streams
stream.treatment.csv <- textConnection(
  "Stream,  CI
Elk, C
Cutoff, C
Rocky, T
Limestone, C
Moon, C,
Mackenzie, T
Thistle, C
Fall, T")

stream.treatment <- read.csv(stream.treatment.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)

# check that streams are all coded
setdiff(stream.treatment$Stream, species.n$Stream)
setdiff(species.n$Stream, stream.treatment$Stream)

# add CI to the data
species.n <- merge(species.n, stream.treatment)

# All streams to be mesured in 2018 onwards (after period) in every year
# Generate 5 years of data with similar numbers of sites measured/year
# and similar year effects etc.
# This is done by finding the variance components and using to generate the data 
# This is done for each species


set.seed(234234)
new.years=2018:2023

stream.year.sites <- unique(species.n[, c("Stream","Year","Site")])
stream.year.n.sites <- plyr::ddply(stream.year.sites, c("Stream","Year"), plyr::summarize, n.sites=length(Site))

# how many sites are currently sampled
xtabs(n.sites~Stream+Year, data=stream.year.n.sites, exclude=NULL, na.action=na.pass)


# How many sites will be sampled?
new.stream.year.n.sites <- expand.grid(Stream=unique(stream.year.sites$Stream), Year=new.years, stringsAsFactors=FALSE)
new.stream.year.n.sites$n.sites <- pmax(2,sample(stream.year.n.sites$n.sites, nrow(new.stream.year.n.sites), replace=TRUE))

new.stream.year.n.sites$n.sites.perm <- 4 # number of permanent sites sampled. Not all permanent sites will be sampled in a particular year

# How many simulated sites sampled?
xtabs(n.sites~Stream+Year, data=new.stream.year.n.sites, exclude=NULL, na.action=na.pass)

head(new.stream.year.n.sites)

species.n.sim <- plyr::ddply(species.n, c("Species","AgeClass"), function(x,new.years, stream.year.n.sites){
   x$WatershedF <- factor(x$Watershed)
   x$StreamF    <- factor(x$Stream)
   x$YearF      <- factor(x$Year)
   x$SiteF      <- factor(x$Site)
   offset <- 0.5 * min(x$CUE[x$CUE>0])
   #browser()
   fit <- lmerTest::lmer(log(CUE+offset) ~ 
                               (1|StreamF)   +
                               (1|YearF)+
                               (1|StreamF:YearF), data=x)
   rand.eff <- ranef(fit)
   new.data <- expand.grid(Stream=unique(x$Stream), Year=new.years, stringsAsFactors=FALSE)
   # Generate multiple sites at each stream-year. We sample from observed data.
   new.data <- merge(new.data, stream.year.n.sites)
   #browser()
   # expand for each site
   new.data <- merge(new.data, unique(x[,c("Stream","Watershed","CI")]), all.x=TRUE)
   new.data$BA <- 'After'
   # Some of these sites are permanent and some are temporary
   new.data <- plyr::ddply(new.data, c("Watershed","Stream","Year","CI","BA"), function(x){
       x <- x[rep(1,x$n.sites),]  # generate rows for each site
       #browser()
       site.names <- paste("SP",1:min(x$n.sites.perm[1],nrow(x)),sep="")
       if(nrow(x)>x$n.sites.perm[1]){site.names <- c(site.names, 
                     paste("ST",".",x$Year[1],".",(1+x$n.sites.perm[1]):nrow(x) ,sep=""))}
       x$Site <- site.names
       #browser()
       x
   })
   #browser()
   # take the existing data and strip off the current data values
   old.data <- x[,c("Stream","Year","Watershed","CI","BA","Site")]
   # add the old data to the new data
   new.data <- plyr::rbind.fill(old.data, new.data)

   # add back stream effect
   stream.eff <- data.frame(stream.eff=unlist(rand.eff$StreamF$"(Intercept)"),
                            Stream = row.names(rand.eff$StreamF))
   new.data <- merge(new.data, stream.eff, all.x=TRUE)

   # generate new year effects by sampling from existing values 
   new.years <- c(new.years,unique(x$Year)) # for the old and new years
   year.eff <- data.frame(Year=new.years, 
                          year.eff=sample(unlist(rand.eff$YearF$"(Intercept)"), length(new.years), replace=TRUE))
   new.data <- merge(new.data, year.eff, all.x=TRUE)

   # generate stream x year interactions
   vc <- as.data.frame(VarCorr(fit))
   sd.stream.year <- max(.1, vc[ vc$grp=="StreamF:YearF","sdcor"])
   stream.year.eff <- unique(new.data[,c("Stream","Year")])
   stream.year.eff$stream.year.eff <- rnorm(nrow(stream.year.eff), sd=sd.stream.year)
   new.data <- merge(new.data, stream.year.eff, all.x=TRUE)


   # Add everything together, generate the log(CUE+offset) and unlog
   new.data$logCUE <- fixef(fit) + new.data$stream.eff + new.data$year.eff + new.data$stream.year.eff +
                    rnorm(nrow(new.data), sd=vc[ vc$grp=="Residual","sdcor"])
   new.data$CUE <- pmax(0, round(exp(new.data$logCUE)-offset,2))
   
   new.data[,c("Watershed","Stream","Year","CI","BA","Site","CUE")]
   
}, new.years=new.years, stream.year.n.sites=new.stream.year.n.sites)

head(species.n.sim)

min(species.n.sim$CUE)
unique(species.n.sim$Year)
unique(species.n$Year)

#species.n <- plyr::rbind.fill(species.n, species.n.sim) # sp
species.n <- species.n.sim  # we now simulate the past data as well

xtabs(~BA+Year, data=species.n, exclude=NULL, na.action=na.pass)
xtabs(~Stream+CI, data=species.n, exclude=NULL,na.action=na.pass)
xtabs(~Stream+Year, data=species.n, exclude=NULL,na.action=na.pass)

# Thistle Creek doesn't have two years before so we make 2018 also a before year
species.n$BA[ species.n$Year==2018] <- "Before"
xtabs(~Stream+paste(Year,'.',substr(BA,1,1),sep=""), data=species.n, exclude=NULL,na.action=na.pass)


# find the grand mean for each species x stream
species.mean <- plyr::ddply(species.n, c("Stream","Species"), plyr::summarize,
          species.mean=mean(log(CUE+.02)))
species.mean

species.n <- merge(species.n, species.mean, all.x=TRUE)
head(species.n)

# Detrend each dataset to have a slope of 0 both before and after prior to adding in simulated effects below
species.n <- plyr::ddply(species.n, c("Stream","Species","BA"), function(x){
    # fit the trend using ordinary least squares
    cat(x$Stream[1], x$Species[1], x$BA[1], "\n")
    #browser()
    fit <- lm(log(CUE+.02) ~ Year, data=x)
    pred.values <- fitted(fit)
    x$logCUE <- log(x$CUE+.02)-pred.values+x$species.mean
    x$CUE <- exp(x$logCUE)
    x
})

plyr::ddply(species.n, c("Stream","Species"), plyr::summarize,
          slope=coef(lm(log(CUE+.02)~Year))[2])

sum(is.na(species.n$CUE))


# compute the year of restoration
year.restore <- plyr::ddply(species.n, c("Stream"), plyr::summarize, 
                            year.restore=min(Year[BA=="After"]))
year.restore

species.n <- merge(species.n, year.restore, all.x=TRUE)

# The above generated data has NO effect of restoration either in the mean or in the slope.
# We will add shift in the meanfor the three streams in the watershet for BKTR
# These are on the log-scale
effects.size.mean.csv <- textConnection(
  "Stream, BKTR, BLTR, RNTR
Fall,      2, 0, 0
Mackenzie, 1.5, 0, 0
Rocky,     2.5, 0, 0")
effects.size.mean <- read.csv(effects.size.mean.csv, header=TRUE, strip.white=TRUE, as.is=TRUE)
write.csv(effects.size.mean, file="effects.size.mean.csv", row.names=FALSE)

effects.size.mean.long <- reshape2::melt(effects.size.mean,
                                         id.var="Stream",
                                         variable.name="Species",
                                         value.name="Rest.add.effect")
effects.size.mean.long$BA <- "After"

species.n <- merge(species.n, effects.size.mean.long, all.x=TRUE)
species.n$Rest.add.effect[ is.na(species.n$Rest.add.effect)] <- 0  # for all other cases the restoration effect (log scale)

species.n$CUE <- species.n$CUE*exp(species.n$Rest.add.effect)


# We will add shift in the slope for the three streams in the watershet for BLTR
# These are on the log-scale
effects.size.slope.csv <- textConnection(
  "Stream, BKTR, BLTR, RNTR
Fall,      0, .4, 0
Mackenzie, 0, .45, 0
Rocky,     0, .35, 0")
effects.size.slope <- read.csv(effects.size.slope.csv, header=TRUE, strip.white=TRUE, as.is=TRUE)
write.csv(effects.size.slope, file="effects.size.slope.csv", row.names=FALSE)

effects.size.slope.long <- reshape2::melt(effects.size.slope,
                                         id.var="Stream",
                                         variable.name="Species",
                                         value.name="Rest.slope.effect")
effects.size.slope.long$BA <- "After"

species.n <- merge(species.n, effects.size.slope.long, all.x=TRUE)
species.n$Rest.slope.effect[ is.na(species.n$Rest.slope.effect)] <- 0  # for all other cases the restoration effect (log scale)

head(species.n)
head(species.n[ species.n$Species=="BLTR" &species.n$BA=="After" & species.n$Stream=="Fall",])

species.n$CUE <- species.n$CUE*exp(species.n$Rest.slope.effect*(pmax(0,species.n$Year-species.n$year.restore)) )

head(species.n)
head(species.n[ species.n$Species=="BLTR" &species.n$BA=="After" & species.n$Stream=="Fall",])



write.csv(species.n[, c("Species","AgeClass","Watershed","Stream","Year",
                        "CI","BA","Site","CUE")],
      file="simulated.data.csv", row.names=FALSE)

