# Author: Andrew Zammit Mangion
# Date: April 2012
# Description: Compares the AWD to the ACLED, ANSO and GTD datasets (in terms of count frequency)
# 
# Requires: the below packages + the datasets acledData, gtdData and origData (enclosed)
#	    
# Please set the current working directory to the parent directory of where this file resides. Output pngs are saved in child directory.

# Packages
library(rgdal)		# for spTransform
library(zoo)		# for yearmon
library(spatstat)	
library(maptools)
library(ggplot2)
library(RecordLinkage)

Extractdata = F  # If you wish to extract the data from source you will need to obtain the raw data. 

if (Extractdata == T) { 

	# Sort AWD data
	origData<-read.csv(paste(getwd(), "/TextS5_Corroboration/DataFiles/afg.csv", sep=""), header=FALSE, as.is = TRUE)
	colnames(origData)<-c("reportkey", "date", "type", "category", "tracking_n", "title", "summary", "region", "attackon", "complex", "reportingu", "unitname", "typeofunit", "FriendlyWo", "FriendlyKi", "HostNatWou", "HostNatKil", "CivWounded", "CivKill", "EnemyWound", "EnemyKille", "EnemyDetai", "MGRS", "lat", "long", "Originator", "UpdatedByG", "CCIR", "sigact", "affiliatio", "dcolor", "classification")

	# Sort set by date
	origData$Date_ <- as.Date(origData$date, format="%m/%d/%Y")
	origData <- origData[order(origData$Date_),]

	# Add week number and month number to set
	year <- as.POSIXlt(origData$date)$year + 1900
	month <- as.POSIXlt(origData$date)$mon
	origData <- cbind(origData,year=year)
	origData <- cbind(origData,month=month)


	# Load shapefiles
	afgprov.outline.unfort <- readShapePoly(paste(getwd(),"/Shapefiles/admin2_poly_32.shp",sep=""))
	afg.outline.unfort <- readShapePoly(paste(getwd(),"/Shapefiles/admin1_poly_32.shp",sep=""))
	afg.outline <- fortify.SpatialPolygons(afg.outline.unfort)

	# Sort provinces alphabetically so that order matches factor levels (ascending alphabetically)
	afgprov.outline.unfort <- afgprov.outline.unfort[order(afgprov.outline.unfort$PRV_NAME),]

	# For setting outlines we need to set this as the polygon files we've got have some overlaps
	spatstat.options(checkpolygons = FALSE)

	# Separate data by provinces
	sepfun <- function(x) {
	  afgprov.outline <- fortify.SpatialPolygons(x)
	  provdata = subset(origData, point.in.polygon(origData$lon,origData$lat,afgprov.outline$long,afgprov.outline$lat)==1)
	  provdata <- cbind(provdata,Province=x$PRV_NAME)
	  provwin = as(x,"owin")
	  return(list(provdata=provdata,provwin=provwin,provoutline=afgprov.outline))
	}

	afgprov.win <- vector("list",32)
	afgprov.data <- vector("list",32)
	afgprov.outline <- vector("list",32)
	for (i in 1:32) {temp <- sepfun(afgprov.outline.unfort[i,])  
			 afgprov.win[[i]] = temp$provwin
			 afgprov.outline[[i]] = temp$provoutline
			 afgprov.data[[i]] = temp$provdata
			 }

	# Merge sets together
	origData <- do.call("rbind",  afgprov.data)
	# Sort again by date
	origData <- origData[order(origData$Date_),]
	save(origData, file="origData")
	
	
	
	# Sort ACLED data
	acledData<-read.csv(paste(getwd(), "/TextS5_Corroboration/DataFiles/acled.csv", sep=""), header=TRUE, as.is = TRUE)

	# Sort set by date
	acledData$Date_ <- as.Date(acledData$EVENT_DATE, format="%m/%d/%Y")
	acledData <- acledData[order(acledData$Date_),]
	
	# Add week number and month number to set
	year <- as.POSIXlt(acledData$Date_)$year + 3900  # 2000 added because year is in "00" format
	acledData <- cbind(acledData,year=year)
	
	# Separate data by provinces
	sepfun <- function(x) {
		  afgprov.outline <- fortify.SpatialPolygons(x)
		  provdata = subset(acledData, point.in.polygon(acledData$LONGITUDE,acledData$LATITUDE,afgprov.outline$long,afgprov.outline$lat)==1)
		  provdata <- cbind(provdata,Province=x$PRV_NAME)
		  provwin = as(x,"owin")
		  return(list(provdata=provdata,provwin=provwin,provoutline=afgprov.outline))
		}
	
	afgprov.win <- vector("list",32)
	afgprov.data <- vector("list",32)
	afgprov.outline <- vector("list",32)
	for (i in 1:32) {temp <- sepfun(afgprov.outline.unfort[i,])  
			 afgprov.win[[i]] = temp$provwin
			 afgprov.outline[[i]] = temp$provoutline
			 afgprov.data[[i]] = temp$provdata
			 }
	
	# Merge sets together
	acledData <- do.call("rbind",  afgprov.data)
	# Sort again by date
	acledData <- acledData[order(acledData$Date_),]
	
	save(acledData, file="acledData")
	
	
	
	# Sort GTD data
	gtdData1<-read.csv(paste(getwd(), "/TextS5_Corroboration/DataFiles/gtd_2004_2007.csv", sep=""), header=TRUE, as.is = TRUE)
	gtdData2<-read.csv(paste(getwd(), "/TextS5_Corroboration/DataFiles/gtd_2008_2009.csv", sep=""), header=TRUE, as.is = TRUE)
	gtdData<-rbind(gtdData1,gtdData2)

	# Sort set by date
	gtdData$Date_ <- as.Date(gtdData$DATE, format="%Y-%m-%d")
	gtdData <- gtdData[order(gtdData$Date_),]
	
	# Add week number and month number to set
	year <- as.POSIXlt(gtdData$Date_)$year + 1900  # 2000 added because year is in "00" format
	month <- as.POSIXlt(gtdData$Date_)$mon
	gtdData <- cbind(gtdData,month=month)
	gtdData <- cbind(gtdData,year=year)
	
	# Load city shapefiles
	afgcities <- readShapePoints(paste(getwd(),"/Shapefiles/07_03_settlements.shp",sep=""))
	
	ClosestMatch2 = function(string, stringVector){
	  distance = levenshteinSim(string, stringVector);
	  stringVector[distance == max(distance)]
	}
	
	CITY2 <- list()
   	length(CITY2) <- length(gtdData$DATE)
	# Make city list
	for (i in 1:length(gtdData$DATE)) {
		CITY2[i] <- as.character(ClosestMatch2(toupper(gtdData$CITY[i]),toupper(afgcities$NAME[afgcities$TYPE>0]))[1])
	}
	
	gtdData <- cbind(gtdData,CITY2=unlist(CITY2))
	
	Province <- list()
	length(Province) <- length(gtdData$DATE)
	for (i in 1:length(gtdData$DATE))
	{
	 Province[i] <-  as.character(afgcities$PROV_NAM[toupper(afgcities$NAME)==unlist(CITY2[i])])
	}
	gtdData <- cbind(gtdData,Province=unlist(Province))
	gtdData <- gtdData[order(gtdData$Date_),]
	
	save(gtdData, file="gtdData")
} else {
graphics.off()

# Load shapefiles
afgprov.outline.unfort <- readShapePoly(paste(getwd(),"/Shapefiles/admin2_poly_32.shp",sep=""))
# Sort provinces alphabetically so that order matches factor levels (ascending alphabetically)
afgprov.outline.unfort <- afgprov.outline.unfort[order(afgprov.outline.unfort$PRV_NAME),]

# For setting outlines we need to set this as the polygon files we've got have some overlaps
spatstat.options(checkpolygons = FALSE)

load("./TextS5_Corroboration/origData",.GlobalEnv) 
load("./TextS5_Corroboration/acledData",.GlobalEnv) 
load("./TextS5_Corroboration/gtdData",.GlobalEnv)

studyyear <- function(data1,data2,year) {
		data1_stat <- vector("list",32)
		data2_stat <- vector("list",32)
		for (i in 1:32) {
			data1_stat[[i]] <- length(subset(data1,(toupper(data1$Province) == toupper(afgprov.outline.unfort[i,]$PRV_NAME)))$Province)
			data2_stat[[i]] <- length(subset(data2,(toupper(data2$Province) == toupper(afgprov.outline.unfort[i,]$PRV_NAME)))$Province)
		}
		return(rbind(data1_stat,data2_stat))
	}

plot_with_stat <- function(my_stat,i,dataname)
{
	   #dev.new(width=5, height=4)
	   png(paste(getwd(),'/TextS5_Corroboration/',dataname,'_',i,'.png',sep=""))  
	   my_stat[1,is.infinite(unlist(my_stat[1,]))] <- NaN
	   my_stat[2,is.infinite(unlist(my_stat[2,]))] <- NaN
	   plot(my_stat[1,],my_stat[2,],main=paste('Year: ',i),xlab=dataname,ylab="AWD",xlim=c(0, max(unlist(my_stat[1,]))),ylim=c(0, max(unlist(my_stat[2,]))))
	   abline(lm(unlist(my_stat[2,]) ~ unlist(my_stat[1,])))
	   x <- cor.test(unlist(my_stat[1,]),unlist(my_stat[2,]), method="pearson",na.rm = T)
	   x_text <- paste('r=',format(x$estimate,digits =2),', p-value=',format(x$p.value,digits=4))
	   text(max(unlist(my_stat[1,]),na.rm=T)/40,max(unlist(my_stat[2,]),na.rm=T)*15/16,x_text,pos=4,font=8)
	   
	   if ((dataname == "ANSO") || (dataname == "GTD per province"))
	   {
	   x <- cor.test(unlist(my_stat[1,c(seq(1,9),seq(11,32))]),unlist(my_stat[2,c(seq(1,9),seq(11,32))]), method="pearson",na.rm = T)
	   x_text <- paste('(w/o Helmand) r=',format(x$estimate,digits =2),', p-value=',format(x$p.value,digits=4))
	   text(max(unlist(my_stat[1,]),na.rm=T)/40,max(unlist(my_stat[2,]),na.rm=T)*14/16,x_text,pos=4,font=8)
	   }
	   dev.off()
	 
}




for (i in seq(2008,2009))
	{  my_stat <- studyyear(subset(acledData,subset=(year==i)),subset(origData,subset=(year==i)),i)
	      plot_with_stat(my_stat,i,"ACLED")
	}
	

ANSO <- c(45,79,45,22,15,114,56,502,27,448,114,28,174,740,118,338,723,125,119,191,260,77,90,253,230,48,5,2,25,108,285,249)
ANSO <- rbind(ANSO,c(43,239,101,88,22,162,137,461,83,620,227,48,177,970,116,478,1318,292,149,188,295,135,85,379,180,52,10,23,65,196,414,259))
colnames(ANSO) <- afgprov.outline.unfort$PRV_NAME
ANSO <- as.data.frame(ANSO)

for (i in seq(2008,2009))
	{  my_stat <- studyyear(subset(acledData,subset=(year==i)),subset(origData,subset=(year==i)),i)
	my_stat[1,] <- ANSO[i-2007,]
	plot_with_stat(my_stat,i,"ANSO")
	}



data1 <- gtdData
data2 <- origData
data1_stat <- vector("list",6)
data2_stat <- vector("list",6)
 for (i in seq(2004,2009)) {
	data1_stat[[i-2003]] <- length(subset(data1,(data1$year == i))$year)
	data2_stat[[i-2003]] <- length(subset(data2,(data2$year == i))$year)
}
my_stat <- rbind(data1_stat,data2_stat)
png(paste(getwd(),'/TextS5_Corroboration/GTD.png',sep=""))  
plot_with_stat(my_stat,i,"GTD per year")
dev.off()


data1 <-gtdData
data2 <- origData
data1_stat <- vector("list",32)
data2_stat <- vector("list",32)
for (i in 1:32) {
			data1_stat[[i]] <- length(subset(data1,(toupper(data1$Province) == toupper(afgprov.outline.unfort[i,]$PRV_NAME)))$Province)
			data2_stat[[i]] <- length(subset(data2,(toupper(data2$Province) == toupper(afgprov.outline.unfort[i,]$PRV_NAME)))$Province)
		}
my_stat <- rbind(data1_stat,data2_stat)
plot_with_stat(my_stat,"All","GTD per province")





data1_stat <- vector("list",32)
data2_stat <- vector("list",32)
data1a <- subset(origData,subset=(year==2008))
data1b <- subset(origData,subset=(year==2009))
for (i in 1:32) {
		count1a <- length(subset(data1a,(data1a$Province == afgprov.outline.unfort[i,]$PRV_NAME))$Province)
		count1b <- length(subset(data1b,(data1b$Province == afgprov.outline.unfort[i,]$PRV_NAME))$Province)
		count2a <- ANSO[1,i]	
		count2b <- ANSO[2,i]
			
			
		data1_stat[[i]] <- count1b/count1a
		data2_stat[[i]] <- count2b/count2a
		}
		
my_stat <- rbind(data1_stat,data2_stat)
plot_with_stat(my_stat,"2008_2009","ANSO ratio")
}