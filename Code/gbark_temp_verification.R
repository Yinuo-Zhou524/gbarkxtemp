setwd("F://LSU//research//gbark_temperature")
#gbark_temp_verif about heat loss through evaporation
temperature<-read.csv("gbark_temperature_Goff_1957.csv")
C55<-temperature[temperature$treatment=="55C",]
C45<-temperature[temperature$treatment=="45C",]

#T1=Ta-LE/hc
#55C L=42.69 KJ/mol
deltaC55<-42.69*C55$E/(29.3*0.135*sqrt(2/(C55$diam*10^(-3))))

#45C L=43.13
deltaC45<-43.13*C45$E/(29.3*0.135*sqrt(2/(C45$diam*10^(-3))))

