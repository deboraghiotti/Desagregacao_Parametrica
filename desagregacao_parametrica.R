require(moments)
require(corrgram)
require(pear)
require(gridExtra)
require(knitr)
require(tcltk)
source('funcoes.R')
#source('C:\Users\debor\Desktop\IC-2019-2020\Desagregacao\Desagregacao_Parametrica\Funcoes_abr2020.R')

input_1 <- tk_choose.files()#Arquivo de serie historica
input_2 <- tk_choose.files()#Arquivo de vazao sintetica
input_3 = scan(what = integer(), nmax = 1)

entrada <- data.frame(read.csv(input_1, sep=";", dec=",")) #Leitura de dados historicos mensais
serie_sintetica <- data.frame(read.csv(input_2, header = F,sep=";", dec=",")) #Leitura de serie sintetica gerada para essa mesma serie
serie_sintetica_normalizada <- apply(log(serie_sintetica), 2, function(x) (x-mean(x)))

SeriesDadosHist <- div_mensais(entrada)
SeriesDadosHist_normalizada <- apply(log(SeriesDadosHist),2, function(x) (x-mean(x)))

Parametro_Hist <- desag_param_info()
#ParametroP_Hist <- desag_param_info_padronizado()

#DesagregadoP = desag_param_mult(Parametro_Hist)
DesagregadoP <- desagregacao_parametrica_mult(Parametro_Hist)
for(i in 1:input_3){
  DesagregadoP[[i]] <- abs(DesagregadoP[[i]])
}
