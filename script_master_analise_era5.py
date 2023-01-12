#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 26 17:12:16 2022
@author: everson
"""

from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from netCDF4 import num2date
from siphon.catalog import TDSCatalog

def delay(dia_inicial,dia_final,int_forecast=17,fator_multi=1):
    #int_forecast = intervalo final das rodadas
    #fator_multi = fator multiplicativo do limite da rodada x 8
    delta_t = dia_inicial - dia_final
    count = 0
    
    if 
    soma_fator = 8*1
    int_forecast = int_forecast + soma_fator
    dados_gfs = TDSCatalog(f'https://rda.ucar.edu/thredds/catalog/files/g/ds084.1/2022/202205{dia_inicial}/catalog.xml')
    validade = dados_gfs.datasets[:int_forecast:2]
    delta_t = dia_inicial - dia_final
    
    for i in range()
    print(validade)


# dia = 15
# count = 0
# forecast = 41 #+8
# soma_for = +8
# # count_dia =

















# 'print(validade)
# # Contador de dias com delay
# while count <= 13:
#     # Abre o catalogo
    
    
#     for indice in range(len(validade)):
#         print (indice)
#         for rodada in validade[indice]:
#             print(rodada)
        

# # count = count + 1
# # forecast = forecast + 8    