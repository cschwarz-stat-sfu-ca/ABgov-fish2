# Read in the data from the FWMIS system

library(ggplot2)
library(plyr)
library(readxl)

# Notice that Thistle Creek has a different format


file.list.csv <- textConnection(
  "file, sheet
  Cutoff.xlsx,    raw data
  Elk.xlsx,       raw data
  Fall.xlsx,      Sheet1
  FallCreek2008.xlsx, Sheet1
  Limestone.xlsx, raw data
  Mackenzie.xlsx, raw data
  Moon.xlsx,      raw data
  Moon2008.xlsx,  Sheet1
  Rocky.xlsx,     raw data")

file.list <- read.csv(file.list.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
file.list

all.data <- plyr::alply(file.list,1, function(x){
     # read in the data
     fish.data <- readxl::read_excel(file.path("RawData",x$file), sheet=x$sheet)
     fish.data$Stream <- gsub(".xlsx","",x$file)
     list(stream   =x$file,
          fish.data=fish.data
          )
})

# The Thistle site is in a different format

fish.data.T <- readxl::read_excel(file.path("RawData","Thistle2.xlsx"), sheet="Sheet1")
names(fish.data.T) <- make.names(names(fish.data.T))
names(fish.data.T)

fish.data.T$Stream <- "Thistle"
fish.data.T$Distance <- 300
fish.data.T$ForkLength <- fish.data.T$Fork.Length..mm.
fish.data.T$Fork.Length..mm. <- NULL
fish.data.T$Location.Code       <- NULL
fish.data.T$Species.Code     <- NULL

length(all.data)
all.data <- append(all.data, list(list(stream="Thistle", fish.data=fish.data.T)))
length(all.data)



# check that all of the names are each file
field.names <- plyr::ldply(all.data, function(x){
    Stream = x$stream
    var.names=names(x$fish.data)
    data.frame(Stream=Stream,
               var.names=var.names, stringsAsFactor=FALSE)
})
xtabs(~var.names+Stream, data=field.names, exclude=NULL, na.action=na.pass)

# put all of the data together
fish.data <- plyr::rbind.fill(plyr::llply(all.data, function(x){x$fish.data}))

# convert UTM to Lat/Long
# https://stackoverflow.com/questions/30018098/how-to-convert-utm-coordinates-to-lat-and-long-in-r?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa

library(rgdal)
select <- !is.na(fish.data$UTM.Easting)
sputm <- SpatialPoints(fish.data[select,c("UTM.Easting","UTM.Northing")], proj4string=CRS("+proj=utm +zone=11 +datum=WGS84"))  
spgeo <- spTransform(sputm, CRS("+proj=longlat +datum=WGS84"))

fish.data$Long[select ] <- spgeo@coords[,1]
fish.data$Lat [select ] <- spgeo@coords[,2]


# One site at moon needs to be discarded because no site information available
select <- fish.data$Stream=="Moon" & fish.data$Distance==0
fish.data[select,]
dim(fish.data)
fish.data <- fish.data[ !select,]
dim(fish.data)

# some basic data checking
unique(fish.data$Date)
str(fish.data$Date)
fish.data$Date <- as.Date(fish.data$Date)
str(fish.data$Date)

fish.data$Year <- as.numeric(format(fish.data$Date, "%Y"))
fish.data$Year.Month <- fish.data$Year + as.numeric(format(fish.data$Date, "%m"))/100
xtabs(~Year+Stream, data=fish.data, exclude=NULL, na.action=na.pass)

# stream vs official names
fish.data$OfficialNa <- toupper(fish.data$OfficialNa)
fish.data$Stream[ fish.data$Stream=="FallCreek2008"] <- "Fall"
fish.data$Stream[ fish.data$Stream=="Moon2008"]      <- "Moon"

xtabs(~Stream+OfficialNa, data=fish.data, exclude=NULL, na.action=na.pass)



# Assign watershed
unique(fish.data$Stream)

watershed.csv <- textConnection(
  "Stream, Watershed
  Cutoff, Clearwater
  Elk,    Clearwater
  Fall,   Ram
  Limestone, Clearwater
  Mackenzie, Ath/NSask
  Moon,      Ath/NSask
  Rocky,     Clearwater
  Thistle,   Ath/NSask")

watershed <- read.csv(watershed.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)

dim(fish.data)
fish.data <- merge(fish.data, watershed)
dim(fish.data)

xtabs(~Stream+Watershed, data=fish.data, exclude=NULL, na.action=na.pass)
xtabs(~interaction(Stream,Watershed, drop=TRUE)+Species, data=fish.data, exclude=NULL, na.action=na.pass)



# define a Site as a unique combination of lat long
site.locations <- ggplot2::ggplot(data=fish.data, aes(x=Long, y=Lat, color=as.factor(Year)))+
  geom_point()+
  facet_wrap(~Stream, ncol=3, scale="free")+
  theme(legend.justification=c(1,0), legend.position=c(1,-.1))+
  scale_color_discrete(name="Year")+
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  guides(color=guide_legend(ncol=2))
site.locations
ggsave(site.locations,
       file=file.path("Report","site-locations.png"),
       h=6, w=6, units="in", dpi=300)

fish.data$Site <- paste(fish.data$Long, "--", fish.data$Lat, sep="") 

site.table <- plyr::ddply(fish.data, c("Watershed","Stream","Year","Year.Month","Site"), plyr::summarize,
                          n.fish=length(Year),
                          siteDistance=Distance[1])
xtabs(n.fish~Stream+Year, data=site.table)
head(site.table)
site.table[ is.na(site.table$siteDistance),]

temp <- xtabs(~interaction(Stream,Watershed, drop=TRUE)+Year.Month, data=site.table)
temp
temp[ temp==0] <- NA
temp

write.csv(temp, file=file.path("Report","n.sites.year.csv"), na="", row.names=TRUE)


# Number of sites
# This includes sites with no fish captured
xtabs(~Stream+Year, data=site.table) # number of sites
xtabs(~Stream+Year.Month, data=site.table)

# Exclude all sites not measured in June-> October
dim(fish.data)
fish.data <- fish.data[ as.numeric(format(fish.data$Date, "%m"))%in% 6:10,]
dim(fish.data)

site.table <- plyr::ddply(fish.data, c("Watershed","Stream","Year","Site"), plyr::summarize,
                          siteDistance=Distance[1])
xtabs(~Stream+Year, data=site.table) # number of sites


# any stream measured more than once?
n.site.year <- plyr::ddply(unique(site.table[, c("Stream","Site","Year")]),
                                  c("Stream","Site"), plyr::summarize,
                                  n.year = length(Year)
                           )
repeated.sites <- n.site.year[n.site.year$n.year > 1,]
site.table[ site.table$Site %in% repeated.sites$Site,]


# check species names
# The NA are NULL records (i.e. a visit to a site, but no fish captured)
xtabs(~Stream+Species, data=fish.data, exclude=NULL, na.action=na.pass)




# Look at distances
xtabs(~Distance+Stream, data=fish.data, exclude=NULL, na.action=na.pass)
# distances = 0 typically indicate problems and so are removed
dim(fish.data)
fish.data <- fish.data[ fish.data$Distance > 0,]
dim(fish.data)
xtabs(~Distance+Stream, data=fish.data, exclude=NULL, na.action=na.pass)

# check that distance doesn't vary within a Stream-Year-Site
distance.sd <- plyr::ddply(fish.data, c("Stream","Year","Site"), plyr::summarize,
                           distance.mean=mean(Distance),
                           distance.sd  =sd  (Distance))
distance.sd[ distance.sd$distance.sd >0 & !is.na(distance.sd$distance.sd),]



# check the length values to see if can divide in mature and immmature by length
xtabs(~ForkLength+Species, data=fish.data, exclude=NULL, na.action=na.pass)
sum(is.na(fish.data$ForkLength))
# which BLTR have length ==0
fish.data[ fish.data$Species=="BLTR" & !is.na(fish.data$Species=="BLTR") &
             fish.data$ForkLength==0,
           c("Stream","Year","Site","Species","ForkLength")]
# exclude all fish with a fork length of 0
dim(fish.data)
fish.data <- fish.data[ fish.data$ForkLength > 0,]
dim(fish.data)



# At this point we only have POSITIVE records, i.e .fish captured at a particular site-year combination
# Any further analysis needs to add back the sites with 0 counts.






