require(moments)
require(corrgram)
require(pear)
require(gridExtra)
require(knitr)
require(tcltk)
source('funcoes.R')

input_1 <- tk_choose.files()#Arquivo de serie historica
input_2 <- tk_choose.files()#Arquivo de vazao sintetica

entrada <- data.frame(read.csv(input_1, sep=";", dec=",")) #Leitura de dados historicos mensais
serie_sintetica <- data.frame(read.csv(input_2, header = F,sep=";", dec=",")) #Leitura de serie sintetica gerada para essa mesma serie

serie_historica <- div_mensais(entrada)

serie_desagregada <- desagregacao_parametrica(serie_sintetica,serie_historica)




