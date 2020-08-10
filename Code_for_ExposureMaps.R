##------------------------------------------------------------------------##
## This script reproduces Figure 1 of Exposure maps to environmental      ## 
## variables used for analysis of CI using 1992-2017 data-set            ##
##------------------------------------------------------------------------##

library(rgdal)
library(rgeos)
library(raster)
library(ggplot2)
library(ggspatial)
library(maps)
library(mapdata)
library(ggmap)
library(grid)
library(ggpubr)
library(egg)

##----------------------------
#Load data-set used on Analyses 
##----------------------------
CIdata=read.csv("CIdata_full.csv")


##-----------------------------------------------------------------------------------------
#Load GBR-Map shapefiles
##Download spatial layers from: http://www.gbrmpa.gov.au/geoportal/catalog/main/home.page
##-----------------------------------------------------------------------------------------

land <- readOGR(dsn="GBRMPA_Data", layer="Great_Barrier_Reef_Features")	
land1 <- gSimplify(land, tol = 0.00001)##if error
crs(land1)
projection(land1)<-CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")

MapAll<-CIdata

coordinates(MapAll)=~LON+LAT
#crop map area using the whole latitudinal extent of the study
sw.b=extent(MapAll@bbox[1]-1,MapAll@bbox[3]+0.5,MapAll@bbox[2]-0.5, MapAll@bbox[4]+0.5)
land.sub <- crop(land1, sw.b)


## identifying reefs per region
##As an inset map
sitesmap2=ggplot(CIdata, aes(x=LON,y=LAT, fill=NewRegion))+
  geom_polygon(data =land.sub,aes(group = group, x = long, y = lat),colour="grey79", fill=NA,alpha=0.5)+
  geom_point(data=CIdata, aes(colour=NewRegion), size=1.5,show.legend = FALSE)+
  scale_color_manual(values = c("steelblue4","magenta4"))+ 
  theme_classic() +
  coord_equal()+
  theme_bw(base_size=7) + 
  theme(plot.title = element_text(size = 12, face = "bold"),
        panel.background = element_rect(fill="white"),  
  axis.title.y = element_blank(),
  axis.title.x = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.background = element_blank())+
  ggtitle("GBR-regions")+
  annotate("text", x = 148.5, y = -17, label = expression("Central"), color="steelblue4")+
  annotate("text", x = 148, y = -22.1, label = expression("Southern"), color="magenta4")


##make inset map of Australia and GBRMPA
aus<-map("worldHires", "Australia", fill=TRUE, xlim=c(110,160),
         ylim=c(-45,-5), mar=c(0,0,0,0))

aumap=ggplot(fortify(aus), aes(y=lat, x=long, group=group))+ 
  geom_rect(mapping=aes(xmin=143, xmax=152, ymin=-23, ymax=-12), color="red", fill=alpha("grey79",0.1))+ 
  geom_polygon()+ 
  theme_classic()+
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA))

##-----------------------------------------
##Plot Exposure to Environmental variables
##-----------------------------------------

####Add Towns to map
LON=c(145.7206,149.181966,151.250934) 
LAT=c(-16.9180,-21.144880,-23.843128)
town=c("Cairns","Mackay","Gladstone")
gbrtowns=data.frame(town,LON,LAT, stringsAsFactors=FALSE)


## Chlorophyll a
chl=ggplot(CIdata, aes(x=LON,y=LAT))+
  geom_polygon(data =land.sub,aes(group = group, x = long, y = lat),colour="grey", fill=NA,alpha=0.9)+
 geom_point(data=CIdata, aes(colour=Chla), size=2.5)+
  scale_color_gradient(low="chartreuse", high="darkgreen")+
  theme_classic()+ylab("Latitude")+xlab("Longitude")+
  theme_classic() +
  coord_equal()+
  theme(axis.text.x=element_text(margin=unit(c(t=-6,r=0,b=0,l=0),"mm")),
        axis.text.y=element_text(margin=unit(c(t=0,r=-1,b=0,l=0),"mm")),
        legend.position= c(0.35,0.8), 
        legend.background = element_rect(fill="transparent",size=8), 
        plot.title = element_text(hjust=0.5,size = 10, face = "bold"),legend.title=element_blank())+
  annotate("text", x = 146.5, y = -14, label = expression("Chl"~italic(a)), size=5)+ 
  annotate("text", x = 148, y = -17.1, label = expression("mg m"^-3))

chl=chl+geom_point(data=gbrtowns,aes(x=LON,y=LAT),colour="black",size=3)+
  geom_text(data=gbrtowns,aes(label=town), hjust = 1,  vjust = 0)


  
##First map with inset of AU
p1=chl+annotation_custom(ggplotGrob(aumap+theme_inset(background_image(fill="transparent"))), 
                         xmin=144.4, xmax=148.5, ymin=-24.1, ymax=-20.5)+
  annotation_north_arrow(height=unit(1,"cm"), width=unit(1,"cm"),
                         pad_x = unit(1, "cm"),pad_y = unit(1, "cm"),which_north="true", location="bl")


## Salinity
sal=ggplot(CIdata, aes(x=LON,y=LAT))+
  geom_polygon(data =land.sub,
               aes(group = group, x = long, y = lat),colour="grey", fill=NA,alpha=0.9)+
  geom_point(data=CIdata, aes(colour=log(Sal+1)), size=2.5)+
  scale_color_gradient(low="cyan", high="cyan4")+
  #ylab("Latitude")+xlab("Longitude")+
  theme_classic() +
  coord_equal()+
  theme(axis.text.x=element_text(margin=unit(c(t=-6,r=0,b=0,l=0),"mm")),
        axis.text.y=element_text(margin=unit(c(t=0,r=-1,b=0,l=0),"mm")),
        legend.position= c(0.37,0.8),
        legend.background = element_rect(fill="transparent",size=8), 
        plot.title = element_text(hjust=0.5,size = 10, face = "bold"),legend.title=element_blank())+
  annotate("text", x = 146.5, y = -14.7, label = expression("Salinity\nIndex"), size=5)+ 
  annotate("text", x = 149.2, y = -18.3, label = expression("log-Days\nbelow\noptima (30)"))

##Fine Seds
TSS=ggplot(CIdata, aes(x=LON,y=LAT))+
  geom_polygon(data =land.sub,
               aes(group = group, x = long, y = lat),colour="grey", fill=NA,alpha=0.9)+
  geom_point(data=CIdata, aes(colour=sqrt(Mud)), size=2.5)+
  scale_color_gradient(low="navajowhite3", high="lightsalmon4")+
  #ylab("Latitude")+xlab("Longitude")+
  theme_classic() +
  coord_equal()+
  theme(axis.text.x=element_text(margin=unit(c(t=-6,r=0,b=0,l=0),"mm")),
        axis.text.y=element_text(margin=unit(c(t=0,r=-1,b=0,l=0),"mm")),
        legend.position= c(0.46,0.8),
        legend.background = element_rect(fill="transparent",size=8), 
        plot.title = element_text(hjust=0.5,size = 10, face = "bold"),legend.title=element_blank())+
  annotate("text", x = 147, y = -14.7,label = expression("Fine suspended\n  sediments"), size=5)+ 
  annotate("text", x =149.1, y = -17.3, label = expression("sqrt- kg m"^-3))

##DINe
DINe=ggplot(CIdata, aes(x=LON,y=LAT))+
  geom_polygon(data =land.sub,
               aes(group = group, x = long, y = lat),colour="grey", fill=NA,alpha=0.9)+
  geom_point(data=CIdata, aes(colour=log(DINe+1)), size=2.5)+
  scale_color_gradient(low="khaki1", high="khaki4")+
  ylab("Latitude")+xlab("Longitude")+
  theme_classic() +
  coord_equal()+
  theme(axis.text.x=element_text(margin=unit(c(t=-6,r=0,b=0,l=0),"mm")),
        axis.text.y=element_text(margin=unit(c(t=0,r=-1,b=0,l=0),"mm")),
        legend.position= c(0.35,0.8), 
        legend.background = element_rect(fill="transparent",size=8), 
        plot.title = element_text(hjust=0.5,size = 10, face = "bold"),legend.title=element_blank())+
  annotate("text", x = 146, y = -14,label = expression("DIN"), size=5)+ 
  annotate("text", x =148, y = -17.3,label = expression("log mg m"^-3))


##PAR
PAR=ggplot(CIdata, aes(x=LON,y=LAT))+
  geom_polygon(data =land.sub,
               aes(group = group, x = long, y = lat),colour="grey", fill=NA,alpha=0.9)+
  geom_point(data=CIdata, aes(colour=PAR),size=2.5)+ 
  scale_color_gradient(low="navajowhite4", high="blue")+
  ylab("Latitude")+xlab("Longitude")+
  theme_classic()+
  coord_equal()+
  theme(axis.text.x=element_text(margin=unit(c(t=-6,r=0,b=0,l=0),"mm")),
        axis.text.y=element_text(margin=unit(c(t=0,r=-1,b=0,l=0),"mm")),
        legend.position= c(0.35,0.8), 
        legend.background = element_rect(fill="transparent",size=8), 
        plot.title = element_text(hjust=0.5,size = 10, face = "bold"),legend.title=element_blank())+
  annotate("text", x = 146, y = -14,label = expression("PAR"), size=5)+ 
  annotate("text", x = 148.2, y = -17.3, label = expression("mol photons m"^-2*d ^-1))


##River DIN
DINr=ggplot(CIdata, aes(x=LON,y=LAT))+
  geom_polygon(data =land.sub,
               aes(group = group, x = long, y = lat),colour="grey", fill=NA,alpha=0.9)+
  geom_point(data=CIdata, aes(colour=log(DINriv+1)),size=2.5)+
  scale_color_gradient(low="khaki1", high="khaki4")+
  #ylab("Latitude")+xlab("Longitude")+
  theme_classic() +
  coord_equal()+
  theme(axis.text.x=element_text(margin=unit(c(t=-6,r=0,b=0,l=0),"mm")),
        axis.text.y=element_text(margin=unit(c(t=0,r=-1,b=0,l=0),"mm")),
        legend.position= c(0.35,0.8), 
        legend.background = element_rect(fill="transparent",size=8), 
        plot.title = element_text(hjust=0.5,size = 10, face = "bold"),legend.title=element_blank())+
  annotate("text", x = 146, y = -14,label = expression("River DIN"), size=5)+ 
  annotate("text", x = 148.1, y = -17.5, label = expression("log mg m"^-3))

p6=DINr+annotation_custom(ggplotGrob(sitesmap2), xmin=149, xmax=153, ymin=-14, ymax=-19)


## Arrange all maps
fplot2=egg::ggarrange(p1+rremove("axis")+rremove("ticks")+rremove("xlab")+rremove("x.text"),sal+rremove("axis")+rremove("ticks")+rremove("xylab")+rremove("xy.text"),
                     TSS+rremove("axis")+rremove("ticks")+rremove("xylab")+rremove("xy.text"),DINe+rremove("axis")+rremove("ticks"),
                     PAR+rremove("axis")+rremove("ticks")+rremove("ylab")+rremove("y.text"),p6+rremove("axis")+rremove("ticks")+rremove("ylab")+rremove("y.text"), nrow=2,ncol=3)

ggsave("Figure1_maps.png", fplot2, width=15, height=15,dpi=300)
