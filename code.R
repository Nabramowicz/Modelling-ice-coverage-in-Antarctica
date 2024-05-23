library(dplyr)
library(sf)
library(terra)
library(tidyverse)
library(xts)
library(lubridate)


#zaimportowanie danych
ice_data<-read.csv("daily_ice_edge.csv")

#zmienna z datami
ice_date<-ice_data[, 1]

#zmianna z długościami geograficznymi i odpowiadającymi im szarokościami, dla poszczególnych lat
ice_nodate<-ice_data[, 2:(ncol(ice_data)-1)]


#0 - minimalny zasieg lodu

min_lat <- abs(apply(ice_nodate, 2, FUN=min))
lon <- matrix(0:359, ncol=ncol(ice_nodate))

min_lat_lon <- as.data.frame(t(rbind(min_lat, lon)))
colnames(min_lat_lon) <- c("min_lat", "lon")

Antarktyda_min<-ggplot(min_lat_lon, aes(x=lon, y=min_lat))+
  geom_path(colour="green")+coord_polar()+ylim(90, 50)+
  xlab("Longitude")+ylab("Latitude [°S]")+
  scale_x_continuous(breaks=c(0, 45, 90, 135, 180, 225, 270, 315), labels=c("0°", "45°E", "90°E", "135°E", "180°", "135°W", "90°W", "45°W"))+
  ggtitle("Minimalny zasięg lodu")


#zasieg matemtyczny

#1 - modelowanie zasiegu

mat_lat_lon<-data.frame("Date"=(0:11160))
predict_mat_lat_lon<-mat_lat_lon

ice_date=dmy(ice_date) #zmiana zapisu daty
daty=as.Date(ice_date) #przekonwertowanie stringa na date
liczba_dni=daty-daty[1] #obliczenie dni od pierwszej daty
ice_data_2<-ice_data
ice_data_2$Date=as.numeric(liczba_dni) #zamiana dni na liczby

lat_lon_days<-left_join(mat_lat_lon, ice_data_2, by="Date")
lat_lon_days<-na.approx(lat_lon_days)
lat_lon_days<-as.data.frame(lat_lon_days)

lat_lon_days_nodate<-lat_lon_days[, 2:(ncol(lat_lon_days)-1)]



#2 - zmiana zasiegu lodu

#stworzenie tymczasowego katalogu na wykresy
dir_out2 <- file.path(tempdir(), "wykresy")
dir.create(dir_out2, recursive = TRUE)

#tworzenie wykresow
dzien_poczatkowy<-ice_date[1]
colors <- c("Zasięg minimalny" = "green", "Zasięg rzeczywisty" = "blue", "Zasięg modelowany" = "red")
for(i in 1:11161){
  dzien<-lat_lon_days[i, 1]
  dzien_iteracji<-dzien_poczatkowy+dzien
  
  #rzeczywisty zasieg
  ice_i<-lat_lon_days_nodate[i, ]
  lat_i <- abs(apply(ice_i, 2, FUN=min))
  lat_lon_i<-as.data.frame(t(rbind(lat_i, lon)))
  colnames(lat_lon_i) <- c("lat", "lon")
  
  #matematyczny zasieg
  mat_lat<-sin(2*pi*lat_lon_i$lat/180)
  mat_lon<-cos(2*pi*lat_lon_i$lon/360)
  
  mat_model<-lm(lat_lon_i$lat~mat_lat+mat_lon)
  predict_mat<-predict(mat_model, as.data.frame(lat_lon_i$lon))
  
  mat_lat_lon_i<-as.data.frame(t(rbind(predict_mat, lon)))
  colnames(mat_lat_lon_i) <- c("lat", "lon")
  
  Antarktyda<-ggplot(min_lat_lon, aes(x=lon, y=min_lat, color="Zasięg minimalny"))+
    geom_path()+
    geom_line(data=lat_lon_i, aes(x=lon, y=lat, color="Zasięg rzeczywisty"))+
    geom_line(data=mat_lat_lon_i, aes(x=lon, y=lat, color="Zasięg modelowany"))+
    coord_polar()+ylim(90, 50)+
    scale_x_continuous(breaks=c(0, 45, 90, 135, 180, 225, 270, 315), labels=c("0°", "45°E", "90°E", "135°E", "180°", "135°W", "90°W", "45°W"))+
    ggtitle(paste(" Zasięg lodu na Antarktydzie "), dzien_iteracji)+
    labs(x="Longitude",
         y="Latitude [°S]",
         color="Legenda")+
    scale_color_manual(values = colors)+
    theme(plot.title = element_text(hjust = 1.2, vjust=4, size=20),
          legend.background = element_rect(fill = "lightblue", size=0.5, linetype="solid", color ="darkblue"))
  
  #zapis wykresow do katalogu
  fp <- file.path(dir_out2, paste0(i, ".png"))
  ggsave(plot = Antarktyda, filename = fp, device = "png", height=7, width=7,dpi=150)
}



#zasieg drugi

#1 - modelowanie zasiegu

mod_lat_lon<-data.frame("Date"=(0:11160))
predict_lat_lon<-mod_lat_lon

ice_date=dmy(ice_date) #zmiana zapisu daty
daty=as.Date(ice_date) #przekonwertowanie stringa na date
liczba_dni=daty-daty[1] #obliczenie dni od pierwszej daty
ice_data_2<-ice_data
ice_data_2$Date=as.numeric(liczba_dni) #zamiana dni na liczby

lat_lon_days<-left_join(mod_lat_lon, ice_data_2, by="Date")
lat_lon_days<-na.approx(lat_lon_days)
lat_lon_days<-as.data.frame(lat_lon_days)

lat_lon_days_nodate<-lat_lon_days[, 2:(ncol(lat_lon_days)-1)]

for(i in 2:361){
  ice_1<-lat_lon_days[, c(1, i)]
  mod_lat_lon<-left_join(mod_lat_lon, ice_1, by="Date")
  
  #obliczenia do dopasowania lini regresji
  spect<-spectrum(mod_lat_lon[, i]) #obliczenie widma dla danych zasięgu lodu
  per<-1/spect$freq[spect$spec==max(spect$spec)] #określenie okresu na podstawie widma
  
  rg<-diff(range(mod_lat_lon[, i])) #zakres danych zasięgu lodu
  plot(mod_lat_lon[, 1], type = "n", ylim = c(min(mod_lat_lon[, i]) - 0.1 * rg, max(mod_lat_lon[, i]) + 0.1 * rg))
  
  reslm2<-lm(mod_lat_lon[, i] ~ sin(2*pi/per*mod_lat_lon[, 1])+cos(2*pi/per*mod_lat_lon[, 1])+sin(4*pi/per*mod_lat_lon[, 1])+cos(4*pi/per*mod_lat_lon[, 1])) #dopasowanie modelu regresji
  summary(reslm2)
  lines(fitted(reslm2)~mod_lat_lon[, 1],col=3) #dodanie linii regresji do wykresu
  
  #przewidywanie wartosci na podstawie lini regresji
  predict_lat<-predict(reslm2, as.data.frame(mod_lat_lon[, 1]))
  predict_lat_lon_i<-as.data.frame(t(rbind(mod_lat_lon[, 1], predict_lat)))
  colnames(predict_lat_lon_i) <- c("Date", "mod_lat")
  
  #pozbywanie sie dni bez pomiarow
  predict_lat_lon<-merge(predict_lat_lon, predict_lat_lon_i, by="Date")
  colnames(predict_lat_lon)[i]<-paste(i-2)
}

#dociecie danych
mod_ice_nodate<-predict_lat_lon[, 2:(ncol(mod_lat_lon))]


#2 - zmiana zasiegu lodu

#stworzenie tymczasowego katalogu na wykresy
dir_out <- file.path(tempdir(), "wykresy")
dir.create(dir_out, recursive = TRUE)

#tworzenie wykresow
dzien_poczatkowy<-ice_date[1]
colors <- c("Zasięg minimalny" = "green", "Zasięg rzeczywisty" = "blue", "Zasięg modelowany" = "red")
for(i in 1:11161){
  dzien<-lat_lon_days[i, 1]
  dzien_iteracji<-dzien_poczatkowy+dzien
  
  #rzeczywisty zasieg
  ice_i<-lat_lon_days_nodate[i, ]
  lat_i <- abs(apply(ice_i, 2, FUN=min))
  lat_lon_i<-as.data.frame(t(rbind(lat_i, lon)))
  colnames(lat_lon_i) <- c("lat", "lon")
  
  #modelowany zasieg
  model_ice_i<-mod_ice_nodate[i, ]
  model_lat_i <- abs(apply(model_ice_i, 2, FUN=min))
  model_lat_lon_i<-as.data.frame(t(rbind(model_lat_i, lon)))
  colnames(model_lat_lon_i) <- c("lat", "lon")
  
  Antarktyda<-ggplot(min_lat_lon, aes(x=lon, y=min_lat, color="Zasięg minimalny"))+
    geom_path()+
    geom_line(data=lat_lon_i, aes(x=lon, y=lat_i, color="Zasięg rzeczywisty"))+
    geom_line(data=model_lat_lon_i, aes(x=lon, y=model_lat_i, color="Zasięg modelowany"))+
    coord_polar()+ylim(90, 50)+
    scale_x_continuous(breaks=c(0, 45, 90, 135, 180, 225, 270, 315), labels=c("0°", "45°E", "90°E", "135°E", "180°", "135°W", "90°W", "45°W"))+
    ggtitle(paste(" Zasięg lodu na Antarktydzie "), dzien_iteracji)+
    labs(x="Longitude",
         y="Latitude [°S]",
         color="Legenda")+
    scale_color_manual(values = colors)+
    theme(plot.title = element_text(hjust = 1.2, vjust=4, size=20),
          legend.background = element_rect(fill = "lightblue", size=0.5, linetype="solid", color ="darkblue"))
  
  #zapis wykresow do katalogu
  fp <- file.path(dir_out, paste0(i, ".png"))
  ggsave(plot = Antarktyda, filename = fp, device = "png", height=7, width=7,dpi=150)
}

