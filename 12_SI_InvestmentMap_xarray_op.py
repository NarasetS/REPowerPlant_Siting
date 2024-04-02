import linopy
import pandas as pd
import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
from pathlib import Path
from shapely import wkt
import numpy as np


coarsenscale = 5
lccs_resolution = 300 * coarsenscale #m
areapergrid = (lccs_resolution/1000) ** 2 ## km2

############################################## Solar Wind Config #########################################################################################################################
mwperkm2_wind = 4.5 ## originally 9 MW/km2 but deduct by 50% of the technically available from IEA's Thailand CET
mwperkm2_solar = 15 ## originally 30 MW/km2 but deduct by 50% of the technically available from IEA's Thailand CET

mwpergrid_wind = np.round(areapergrid * mwperkm2_wind,2)
mwpergrid_solar =  np.round(areapergrid * mwperkm2_solar,2)

print('areapergrid = ',areapergrid)
print('mwpergrid_wind = ',mwpergrid_wind)
print('mwpergrid_solar = ',mwpergrid_solar)

maxcapacityfor_spp = 90
maxcapacityfor_vspp = 10

rollingwindow_spp_wind = int(np.ceil(np.sqrt(maxcapacityfor_spp/mwpergrid_wind)))
# if ( rollingwindow_spp_wind % 2 == 0 ) : rollingwindow_spp_wind += 1
rollingwindow_vspp_wind = int(np.ceil(np.sqrt(maxcapacityfor_vspp/mwpergrid_wind)))
# if ( rollingwindow_vspp_wind % 2 == 0 ) : rollingwindow_vspp_wind += 1

print('rollingwindow_spp_wind = ',rollingwindow_spp_wind,' * ',rollingwindow_spp_wind,' = ', mwpergrid_wind * (rollingwindow_spp_wind**2))
print('rollingwindow_vspp_wind = ',rollingwindow_vspp_wind,' * ',rollingwindow_vspp_wind,' = ', mwpergrid_wind * (rollingwindow_vspp_wind**2))

rollingwindow_spp_solar = int(np.ceil(np.sqrt(maxcapacityfor_spp/mwpergrid_solar)))
# if ( rollingwindow_spp_solar % 2 == 0 ) : rollingwindow_spp_solar += 1
rollingwindow_vspp_solar = int(np.ceil(np.sqrt(maxcapacityfor_vspp/mwpergrid_solar)))
# if ( rollingwindow_vspp_solar % 2 == 0 ) : rollingwindow_vspp_solar += 1

print('rollingwindow_spp_solar = ',rollingwindow_spp_solar,' * ',rollingwindow_spp_solar,' = ', mwpergrid_solar * (rollingwindow_spp_solar**2))
print('rollingwindow_vspp_solar = ',rollingwindow_vspp_solar,' * ',rollingwindow_vspp_solar,' = ', mwpergrid_solar * (rollingwindow_vspp_solar**2))
############################################## Solar Wind Config #########################################################################################################################

############################################## Fuel-based Plant Config #########################################################################################################################
suitablearea_biomass = 150 ## km2
suitablearea_bgec = 150 ## km2
suitablearea_bgww = 150 ## km2
suitablearea_msw = 150 ## km2
suitablearea_iew = 150 ## km2

rollingwindow_spp_biomass = int(np.ceil(np.sqrt(suitablearea_biomass * 9 /areapergrid)))
rollingwindow_vspp_biomass = int(np.ceil(np.sqrt(suitablearea_biomass/areapergrid)))

rollingwindow_spp_bgec = int(np.ceil(np.sqrt(suitablearea_bgec * 9 /areapergrid)))
rollingwindow_vspp_bgec = int(np.ceil(np.sqrt(suitablearea_bgec/areapergrid)))

rollingwindow_spp_bgww = int(np.ceil(np.sqrt(suitablearea_bgww * 9 /areapergrid)))
rollingwindow_vspp_bgww = int(np.ceil(np.sqrt(suitablearea_bgww/areapergrid)))

rollingwindow_spp_msw = int(np.ceil(np.sqrt(suitablearea_msw * 9 /areapergrid)))
rollingwindow_vspp_msw = int(np.ceil(np.sqrt(suitablearea_msw/areapergrid)))

rollingwindow_spp_iew = int(np.ceil(np.sqrt(suitablearea_iew * 9 /areapergrid)))
rollingwindow_vspp_iew = int(np.ceil(np.sqrt(suitablearea_iew/areapergrid)))

print('rollingwindow_spp_biomass = ',rollingwindow_spp_biomass,' * ',rollingwindow_spp_biomass)
print('rollingwindow_vspp_biomass = ',rollingwindow_vspp_biomass,' * ',rollingwindow_vspp_biomass)

print('rollingwindow_spp_bgec = ',rollingwindow_spp_bgec,' * ',rollingwindow_spp_bgec)
print('rollingwindow_vspp_bgec = ',rollingwindow_vspp_bgec,' * ',rollingwindow_vspp_bgec)

print('rollingwindow_spp_bgww = ',rollingwindow_spp_bgww,' * ',rollingwindow_spp_bgww)
print('rollingwindow_vspp_bgww = ',rollingwindow_vspp_bgww,' * ',rollingwindow_vspp_bgww)

print('rollingwindow_spp_msw = ',rollingwindow_spp_msw,' * ',rollingwindow_spp_msw)
print('rollingwindow_vspp_msw = ',rollingwindow_vspp_msw,' * ',rollingwindow_vspp_msw)

print('rollingwindow_spp_iew = ',rollingwindow_spp_iew,' * ',rollingwindow_spp_iew)
print('rollingwindow_vspp_iew = ',rollingwindow_vspp_iew,' * ',rollingwindow_vspp_iew)
############################################## Fuel-based Plant Config #########################################################################################################################

############################################## Coarsen #########################################################################################################################
xr_final_SI = xr.open_dataset('Output\\xr_final_SI_all.nc')
xr_final_SI = xr_final_SI.drop_vars('ADM1_EN')
xr_final_SI = xr_final_SI.coarsen(lat = coarsenscale, lon= coarsenscale, boundary='pad').sum()
xr_final_SI['SI_Wind'] = xr_final_SI['SI_Wind'] / (coarsenscale**2)
xr_final_SI['SI_Solar'] = xr_final_SI['SI_Solar'] / (coarsenscale**2)
xr_final_SI['SI_Biomass'] = xr_final_SI['SI_Biomass'] / (coarsenscale**2)
xr_final_SI['SI_BGEC'] = xr_final_SI['SI_BGEC'] / (coarsenscale**2)
xr_final_SI['SI_BGWW'] = xr_final_SI['SI_BGWW'] / (coarsenscale**2)
xr_final_SI['SI_MSW'] = xr_final_SI['SI_MSW'] / (coarsenscale**2)
xr_final_SI['SI_IEW'] = xr_final_SI['SI_IEW'] / (coarsenscale**2)
print("AVA Wind = ",xr_final_SI['AVA_Wind'].sum())
print("AVA Solar = ",xr_final_SI['AVA_Solar'].sum())
print("A_Biomass = ",xr_final_SI['A_Biomass'].sum())
print("A_BGEC = ",xr_final_SI['A_BGEC'].sum())
print("A_BGWW = ",xr_final_SI['A_BGWW'].sum())
print("A_MSW = ",xr_final_SI['A_MSW'].sum())
print("A_IEW = ",xr_final_SI['A_IEW'].sum())
print('Max SI Wind = ',xr_final_SI['SI_Wind'].max())
print('Max SI Solar = ',xr_final_SI['SI_Solar'].max())
print('Max SI SI_Biomass = ',xr_final_SI['SI_Biomass'].max())
print('Max SI SI_BGEC = ',xr_final_SI['SI_BGEC'].max())
print('Max SI SI_BGWW = ',xr_final_SI['SI_BGWW'].max())
print('Max SI SI_MSW = ',xr_final_SI['SI_MSW'].max())
print('Max SI SI_IEW = ',xr_final_SI['SI_IEW'].max())
############################################## Coarsen #########################################################################################################################

############################################## Assign Region #########################################################################################################################
region = pd.read_csv('Data\\Region.csv',index_col=False)
thailandmap = gpd.read_file('Data\\tha_admbnda_adm1_rtsd_20220121\\tha_admbnda_adm1_rtsd_20220121.shp')

thailandmap.crs = {'init': 'epsg:4326'}
list_region = []

count = 0
for i in thailandmap['ADM1_TH']:
    r = region['region'].loc[region['province'] == i]
    try : 
        # print(i,r.values[0])
        list_region.append(r.values[0])
    except :
        print(i,'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
        list_region.append('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')       
    
thailandmap['region'] = list_region
thailandmap['center'] = thailandmap['geometry'].centroid
thailandmap = thailandmap.set_geometry('center')

thailandmap = thailandmap.drop(columns=['Shape_Leng',
                                        'Shape_Area',
                                        'ADM1_PCODE',
                                        'ADM1_REF',
                                        'ADM1ALT1EN',
                                        'ADM1ALT2EN',
                                        'ADM1ALT1TH',
                                        'ADM1ALT2TH',
                                        'ADM0_EN',
                                        'ADM0_TH',
                                        'ADM0_PCODE',
                                        'date',
                                        'validOn',
                                        'validTo'
                                        ,'geometry'
                                        ])

df_final_SI = xr_final_SI.to_dataframe()
df_final_SI.reset_index(inplace=True)
df_final_SI = gpd.GeoDataFrame(df_final_SI, geometry =gpd.points_from_xy(df_final_SI['lon'],df_final_SI['lat']))
df_final_SI.crs = {'init': 'epsg:4326'}
df_final_SI.reset_index(inplace= True, drop = False)

df_final_SI_2 = gpd.sjoin_nearest(df_final_SI,thailandmap,how = 'left')
df_final_SI_2 = df_final_SI_2.drop(columns=['ADM1_TH','geometry','index_right'])

df_final_SI_2 = df_final_SI_2.drop_duplicates('index')
df_final_SI_2 = df_final_SI_2.drop(columns=['index'])

df_final_SI_2.reset_index(inplace= True, drop = True)
df_final_SI_2 = df_final_SI_2.set_index(['lat', 'lon'])

xr_final_SI = xr.Dataset.from_dataframe(df_final_SI_2)
############################################## Assign Region #########################################################################################################################

############################################## Summary #########################################################################################################################
print("AVA Wind = ",xr_final_SI['AVA_Wind'].sum())
print("AVA Solar = ",xr_final_SI['AVA_Solar'].sum())
print("A_Biomass = ",xr_final_SI['A_Biomass'].sum())
print("A_BGEC = ",xr_final_SI['A_BGEC'].sum())
print("A_BGWW = ",xr_final_SI['A_BGWW'].sum())
print("A_MSW = ",xr_final_SI['A_MSW'].sum())
print("A_IEW = ",xr_final_SI['A_IEW'].sum())

print('Max SI Wind = ',xr_final_SI['SI_Wind'].max())
print('Max SI Solar = ',xr_final_SI['SI_Solar'].max())
print('Max SI SI_Biomass = ',xr_final_SI['SI_Biomass'].max())
print('Max SI SI_BGEC = ',xr_final_SI['SI_BGEC'].max())
print('Max SI SI_BGWW = ',xr_final_SI['SI_BGWW'].max())
print('Max SI SI_MSW = ',xr_final_SI['SI_MSW'].max())
print('Max SI SI_IEW = ',xr_final_SI['SI_IEW'].max())
xr_ref = xr_final_SI
############################################## Summary #########################################################################################################################

# ######## PDP ############################################################
# SPP_quota_wind_total =  6335
# SPP_quota_wind_R0 = 0
# SPP_quota_wind_R1 = 0
# SPP_quota_wind_R2 = 6075
# SPP_quota_wind_R3 =  260
# SPP_quota_wind_R4 =  0

# VSPP_quota_wind_total =  0
# VSPP_quota_wind_R0 =  0
# VSPP_quota_wind_R1 =  0
# VSPP_quota_wind_R2 =  0
# VSPP_quota_wind_R3 =  0
# VSPP_quota_wind_R4 = 0

# SPP_quota_solar_total =  18810
# SPP_quota_solar_R0 =  0
# SPP_quota_solar_R1 =  3150
# SPP_quota_solar_R2 =  3060
# SPP_quota_solar_R3 = 8100
# SPP_quota_solar_R4 =  4500

# VSPP_quota_solar_total = 5372
# VSPP_quota_solar_R0 = 16
# VSPP_quota_solar_R1 =  983
# VSPP_quota_solar_R2 =  1929
# VSPP_quota_solar_R3 =   725
# VSPP_quota_solar_R4 =  1719

# SPP_quota_biomass_total =  0
# SPP_quota_biomass_R0 =  0
# SPP_quota_biomass_R1 =  0
# SPP_quota_biomass_R2 =  0
# SPP_quota_biomass_R3 = 0
# SPP_quota_biomass_R4 =  0

# VSPP_quota_biomass_total = 1010
# VSPP_quota_biomass_R0 = 0
# VSPP_quota_biomass_R1 =  213
# VSPP_quota_biomass_R2 =  361
# VSPP_quota_biomass_R3 =   90
# VSPP_quota_biomass_R4 =  346

# SPP_quota_bgec_total =  0
# SPP_quota_bgec_R0 =  0
# SPP_quota_bgec_R1 =  0
# SPP_quota_bgec_R2 =  0
# SPP_quota_bgec_R3 = 0
# SPP_quota_bgec_R4 =  0

# VSPP_quota_bgec_total = 1042
# VSPP_quota_bgec_R0 = 0
# VSPP_quota_bgec_R1 =  232
# VSPP_quota_bgec_R2 =  452
# VSPP_quota_bgec_R3 =   67
# VSPP_quota_bgec_R4 =  291

# SPP_quota_msw_total =  0
# SPP_quota_msw_R0 =  0
# SPP_quota_msw_R1 =  0
# SPP_quota_msw_R2 =  0
# SPP_quota_msw_R3 = 0
# SPP_quota_msw_R4 =  0

# VSPP_quota_msw_total = 342
# VSPP_quota_msw_R0 = 80
# VSPP_quota_msw_R1 =  77
# VSPP_quota_msw_R2 =  100
# VSPP_quota_msw_R3 =   35
# VSPP_quota_msw_R4 =  50
# ######## PDP ############################################################

######## PDP ############################################################
SPP_quota_wind_total =  0
SPP_quota_wind_R0 = 0
SPP_quota_wind_R1 = 0
SPP_quota_wind_R2 = 0
SPP_quota_wind_R3 =  0
SPP_quota_wind_R4 =  0

VSPP_quota_wind_total =  0
VSPP_quota_wind_R0 =  0
VSPP_quota_wind_R1 =  0
VSPP_quota_wind_R2 =  0
VSPP_quota_wind_R3 =  0
VSPP_quota_wind_R4 = 0

SPP_quota_solar_total =  0
SPP_quota_solar_R0 =  0
SPP_quota_solar_R1 =  0
SPP_quota_solar_R2 =  0
SPP_quota_solar_R3 = 0
SPP_quota_solar_R4 =  0

VSPP_quota_solar_total = 0
VSPP_quota_solar_R0 = 0
VSPP_quota_solar_R1 =  0
VSPP_quota_solar_R2 =  0
VSPP_quota_solar_R3 =   0
VSPP_quota_solar_R4 =  0

SPP_quota_biomass_total =  0
SPP_quota_biomass_R0 =  0
SPP_quota_biomass_R1 =  0
SPP_quota_biomass_R2 =  0
SPP_quota_biomass_R3 = 0
SPP_quota_biomass_R4 =  0

VSPP_quota_biomass_total = 0
VSPP_quota_biomass_R0 = 0
VSPP_quota_biomass_R1 =  0
VSPP_quota_biomass_R2 =  0
VSPP_quota_biomass_R3 =   0
VSPP_quota_biomass_R4 =  0

SPP_quota_bgec_total =  0
SPP_quota_bgec_R0 =  0
SPP_quota_bgec_R1 =  0
SPP_quota_bgec_R2 =  0
SPP_quota_bgec_R3 = 0
SPP_quota_bgec_R4 =  0

VSPP_quota_bgec_total = 0
VSPP_quota_bgec_R0 = 0
VSPP_quota_bgec_R1 =  0
VSPP_quota_bgec_R2 =  0
VSPP_quota_bgec_R3 =   0
VSPP_quota_bgec_R4 =  0

SPP_quota_msw_total =  0
SPP_quota_msw_R0 =  0
SPP_quota_msw_R1 =  0
SPP_quota_msw_R2 =  0
SPP_quota_msw_R3 = 0
SPP_quota_msw_R4 =  0

VSPP_quota_msw_total = 0
VSPP_quota_msw_R0 = 0
VSPP_quota_msw_R1 =  0
VSPP_quota_msw_R2 =  0
VSPP_quota_msw_R3 =   0
VSPP_quota_msw_R4 =  0
######## PDP ############################################################

######################## model #####################################################
m = linopy.Model()
##############################################################################
built_SPP_wind = m.add_variables(integer=True,upper=1,lower=0, coords=xr_ref.coords, name='built_SPP_wind')
cap_SPP_wind = m.add_variables(upper=90,lower=0.00, coords=xr_ref.coords, name='cap_SPP_wind')

built_VSPP_wind = m.add_variables(integer=True,upper=1,lower=0, coords=xr_ref.coords, name='built_VSPP_wind')
cap_VSPP_wind = m.add_variables(upper=10,lower=0.00, coords=xr_ref.coords, name='cap_VSPP_wind')
##############################################################################
built_SPP_solar = m.add_variables(integer=True,upper=1,lower=0, coords=xr_ref.coords, name='built_SPP_solar')
cap_SPP_solar = m.add_variables(upper=90,lower=0.00, coords=xr_ref.coords, name='cap_SPP_solar')

built_VSPP_solar = m.add_variables(integer=True,upper=1,lower=0, coords=xr_ref.coords, name='built_VSPP_solar')
cap_VSPP_solar = m.add_variables(upper=10,lower=0.00, coords=xr_ref.coords, name='cap_VSPP_solar')
##############################################################################
built_SPP_biomass = m.add_variables(integer=True,upper=1,lower=0, coords=xr_ref.coords, name='built_SPP_biomass')
cap_SPP_biomass = m.add_variables(upper=90,lower=0.00, coords=xr_ref.coords, name='cap_SPP_biomass')

built_VSPP_biomass = m.add_variables(integer=True,upper=1,lower=0, coords=xr_ref.coords, name='built_VSPP_biomass')
cap_VSPP_biomass = m.add_variables(upper=10,lower=0.00, coords=xr_ref.coords, name='cap_VSPP_biomass')
##############################################################################
built_SPP_bgec = m.add_variables(integer=True,upper=1,lower=0, coords=xr_ref.coords, name='built_SPP_bgec')
cap_SPP_bgec = m.add_variables(upper=90,lower=0.00, coords=xr_ref.coords, name='cap_SPP_bgec')

built_VSPP_bgec = m.add_variables(integer=True,upper=1,lower=0, coords=xr_ref.coords, name='built_VSPP_bgec')
cap_VSPP_bgec = m.add_variables(upper=10,lower=0.00, coords=xr_ref.coords, name='cap_VSPP_bgec')
##############################################################################

built_SPP_msw = m.add_variables(integer=True,upper=1,lower=0, coords=xr_ref.coords, name='built_SPP_msw')
cap_SPP_msw = m.add_variables(upper=90,lower=0.00, coords=xr_ref.coords, name='cap_SPP_msw')

built_VSPP_msw = m.add_variables(integer=True,upper=1,lower=0, coords=xr_ref.coords, name='built_VSPP_msw')
cap_VSPP_msw = m.add_variables(upper=10,lower=0.00, coords=xr_ref.coords, name='cap_VSPP_msw')
##############################################################################

############################################ Constraint Building Location Logic ##############################################################################
constr_built_logic =  m.add_constraints(
(
    built_SPP_wind.rolling(lat = rollingwindow_spp_wind+1,center=True,min_periods=1).sum().rolling(lon = rollingwindow_spp_wind+1,center=True,min_periods=1).sum()
    +
    built_VSPP_wind.rolling(lat = rollingwindow_vspp_wind+1,center=True,min_periods=1).sum().rolling(lon = rollingwindow_vspp_wind+1,center=True,min_periods=1).sum()
    +
    built_SPP_solar.rolling(lat = rollingwindow_spp_solar+1,center=True,min_periods=1).sum().rolling(lon = rollingwindow_spp_solar+1,center=True,min_periods=1).sum()
    +
    built_VSPP_solar.rolling(lat = rollingwindow_vspp_solar+1,center=True,min_periods=1).sum().rolling(lon = rollingwindow_vspp_solar+1,center=True,min_periods=1).sum()
    +
    built_SPP_biomass.rolling(lat = rollingwindow_spp_biomass+1,center=True,min_periods=1).sum().rolling(lon = rollingwindow_spp_biomass+1,center=True,min_periods=1).sum()
    +
    built_VSPP_biomass.rolling(lat = rollingwindow_vspp_biomass+1,center=True,min_periods=1).sum().rolling(lon = rollingwindow_vspp_biomass+1,center=True,min_periods=1).sum()
    +
    built_SPP_bgec.rolling(lat = rollingwindow_spp_bgec+1,center=True,min_periods=1).sum().rolling(lon = rollingwindow_spp_bgec+1,center=True,min_periods=1).sum()
    +
    built_VSPP_bgec.rolling(lat = rollingwindow_vspp_bgec+1,center=True,min_periods=1).sum().rolling(lon = rollingwindow_vspp_bgec+1,center=True,min_periods=1).sum()
    +
    built_SPP_msw.rolling(lat = rollingwindow_spp_msw+1,center=True,min_periods=1).sum().rolling(lon = rollingwindow_spp_msw+1,center=True,min_periods=1).sum()
    +
    built_VSPP_msw.rolling(lat = rollingwindow_vspp_msw+1,center=True,min_periods=1).sum().rolling(lon = rollingwindow_vspp_msw+1,center=True,min_periods=1).sum()
)
    <= 1
    , name='constr_built_logic'
)

############################################ Constraint Capacity ##############################################################################
constr_maxcap_SPP_wind = m.add_constraints(
    cap_SPP_wind <= (built_SPP_wind) * (mwperkm2_wind * xr_ref['AVA_Wind'].rolling(lon = rollingwindow_spp_wind, lat = rollingwindow_spp_wind, min_periods=1,center=True).sum())
    ,name = 'constr_maxcap_SPP_wind'
)
constr_mincap_SPP_wind = m.add_constraints(
    cap_SPP_wind >= (built_SPP_wind) * 16
    ,name = 'constr_mincap_SPP_wind'
)
constr_maxcap_VSPP_wind = m.add_constraints(
    cap_VSPP_wind <= (built_VSPP_wind) * (mwperkm2_wind * xr_ref['AVA_Wind'].rolling(lon = rollingwindow_vspp_wind, lat = rollingwindow_vspp_wind, min_periods=1,center=True).sum())
    ,name = 'constr_maxcap_VSPP_wind'
)
constr_mincap_VSPP_wind = m.add_constraints(
    cap_VSPP_wind >= (built_VSPP_wind) * 8
    ,name = 'constr_mincap_VSPP_wind'
)
##############################################################################
constr_maxcap_SPP_solar = m.add_constraints(
    cap_SPP_solar <= (built_SPP_solar) * (mwperkm2_solar * xr_ref['AVA_Solar'].rolling(lon = rollingwindow_spp_solar, lat = rollingwindow_spp_solar, min_periods=1,center=True).sum())
    ,name = 'constr_maxcap_SPP_solar'
)
constr_mincap_SPP_solar = m.add_constraints(
    cap_SPP_solar >= (built_SPP_solar) * 15
    ,name = 'constr_mincap_SPP_solar'
)
constr_maxcap_VSPP_solar = m.add_constraints(
    cap_VSPP_solar <= (built_VSPP_solar) * (mwperkm2_solar * xr_ref['AVA_Solar'].rolling(lon = rollingwindow_vspp_solar, lat = rollingwindow_vspp_solar, min_periods=1,center=True).sum())
    ,name = 'constr_maxcap_VSPP_solar'
)
constr_mincap_VSPP_solar = m.add_constraints(
    cap_VSPP_solar >= (built_VSPP_solar) * 1
    ,name = 'constr_mincap_VSPP_solar'
)
##############################################################################
constr_maxcap_SPP_biomass = m.add_constraints(
    cap_SPP_biomass <= (built_SPP_biomass) * (xr_ref['A_Biomass'].rolling(lon = rollingwindow_spp_biomass, lat = rollingwindow_spp_biomass, min_periods=1,center=True).sum())
    ,name = 'constr_maxcap_SPP_biomass'
)
constr_mincap_SPP_biomass = m.add_constraints(
    cap_SPP_biomass >= (built_SPP_biomass) * 15
    ,name = 'constr_mincap_SPP_biomass'
)
constr_maxcap_VSPP_biomass = m.add_constraints(
    cap_VSPP_biomass <= (built_VSPP_biomass) * (xr_ref['A_Biomass'].rolling(lon = rollingwindow_vspp_biomass, lat = rollingwindow_vspp_biomass, min_periods=1,center=True).sum())
    ,name = 'constr_maxcap_VSPP_biomass'
)
constr_mincap_VSPP_biomass = m.add_constraints(
    cap_VSPP_biomass >= (built_VSPP_biomass) * 1
    ,name = 'constr_mincap_VSPP_biomass'
)
##############################################################################
constr_maxcap_SPP_bgec = m.add_constraints(
    cap_SPP_bgec <= (built_SPP_bgec) * (xr_ref['A_BGEC'].rolling(lon = rollingwindow_spp_bgec, lat = rollingwindow_spp_bgec, min_periods=1,center=True).sum())
    ,name = 'constr_maxcap_SPP_bgec'
)
constr_mincap_SPP_bgec = m.add_constraints(
    cap_SPP_bgec >= (built_SPP_bgec) * 15
    ,name = 'constr_mincap_SPP_bgec'
)
constr_maxcap_VSPP_bgec = m.add_constraints(
    cap_VSPP_bgec<= (built_VSPP_bgec) * (xr_ref['A_BGEC'].rolling(lon = rollingwindow_vspp_bgec, lat = rollingwindow_vspp_bgec, min_periods=1,center=True).sum())
    ,name = 'constr_maxcap_VSPP_bgec'
)
constr_mincap_VSPP_bgec = m.add_constraints(
    cap_VSPP_bgec >= (built_VSPP_bgec) * 1
    ,name = 'constr_mincap_VSPP_bgec'
)
##############################################################################
constr_maxcap_SPP_msw = m.add_constraints(
    cap_SPP_msw <= (built_SPP_msw) * (xr_ref['A_MSW'].rolling(lon = rollingwindow_spp_msw, lat = rollingwindow_spp_msw, min_periods=1,center=True).sum())
    ,name = 'constr_maxcap_SPP_msw'
)
constr_mincap_SPP_msw = m.add_constraints(
    cap_SPP_msw >= (built_SPP_msw) * 15
    ,name = 'constr_mincap_SPP_msw'
)
constr_maxcap_VSPP_msw = m.add_constraints(
    cap_VSPP_msw<= (built_VSPP_msw) * (xr_ref['A_MSW'].rolling(lon = rollingwindow_vspp_msw, lat = rollingwindow_vspp_msw, min_periods=1,center=True).sum())
    ,name = 'constr_maxcap_VSPP_msw'
)
constr_mincap_VSPP_msw = m.add_constraints(
    cap_VSPP_msw >= (built_VSPP_msw) * 1
    ,name = 'constr_mincap_VSPP_msw'
)
##############################################################################

############################################ Constraint Capacity ##############################################################################

############################################################################################################################################################

constr_SPP_quota_wind = m.add_constraints((cap_SPP_wind.sum()) >= SPP_quota_wind_total , name='constr_SPP_quota_wind')

constr_SPP_quota_wind_r0 = m.add_constraints(lhs = (cap_SPP_wind).where(xr_ref['region'] == 'R0',drop=True).sum()
                                             , sign = '>=' , rhs = SPP_quota_wind_R0, name='constr_SPP_quota_wind_r0')

constr_SPP_quota_wind_r1 = m.add_constraints(lhs = (cap_SPP_wind).where(xr_ref['region'] == 'R1',drop=True).sum()
                                             , sign = '>=' , rhs = SPP_quota_wind_R1, name='constr_SPP_quota_wind_r1')

constr_SPP_quota_wind_r2 = m.add_constraints(lhs = (cap_SPP_wind).where(xr_ref['region'] == 'R2',drop=True).sum()
                                             , sign = '>=' , rhs = SPP_quota_wind_R2, name='constr_SPP_quota_wind_r2')

constr_SPP_quota_wind_r3 = m.add_constraints(lhs = (cap_SPP_wind).where(xr_ref['region'] == 'R3',drop=True).sum()
                                             , sign = '>=' , rhs = SPP_quota_wind_R3, name='constr_SPP_quota_wind_r3')

constr_SPP_quota_wind_r4 = m.add_constraints(lhs = (cap_SPP_wind).where(xr_ref['region'] == 'R4',drop=True).sum()
                                             , sign = '>=' , rhs = SPP_quota_wind_R4, name='constr_SPP_quota_wind_r4')

constr_VSPP_quota_wind = m.add_constraints((cap_VSPP_wind.sum()) >= VSPP_quota_wind_total , name='constr_VSPP_quota_wind')

constr_VSPP_quota_wind_r0 = m.add_constraints(lhs = (cap_VSPP_wind).where(xr_ref['region'] == 'R0',drop=True).sum()
                                             , sign = '>=' , rhs = VSPP_quota_wind_R0, name='constr_VSPP_quota_wind_r0')

constr_VSPP_quota_wind_r1 = m.add_constraints(lhs = (cap_VSPP_wind).where(xr_ref['region'] == 'R1',drop=True).sum()
                                             , sign = '>=' , rhs = VSPP_quota_wind_R1, name='constr_VSPP_quota_wind_r1')

constr_VSPP_quota_wind_r2 = m.add_constraints(lhs = (cap_VSPP_wind).where(xr_ref['region'] == 'R2',drop=True).sum()
                                             , sign = '>=' , rhs = VSPP_quota_wind_R2, name='constr_VSPP_quota_wind_r2')

constr_VSPP_quota_wind_r3 = m.add_constraints(lhs = (cap_VSPP_wind).where(xr_ref['region'] == 'R3',drop=True).sum()
                                             , sign = '>=' , rhs = VSPP_quota_wind_R3, name='constr_VSPP_quota_wind_r3')

constr_VSPP_quota_wind_r4 = m.add_constraints(lhs = (cap_VSPP_wind).where(xr_ref['region'] == 'R4',drop=True).sum()
                                             , sign = '>=' , rhs = VSPP_quota_wind_R4, name='constr_VSPP_quota_wind_r4')
############################################################################################################################################################

constr_SPP_quota_solar = m.add_constraints((cap_SPP_solar.sum()) >= SPP_quota_solar_total, name='constr_SPP_quota_solar')

constr_SPP_quota_solar_r0 = m.add_constraints(lhs = (cap_SPP_solar).where(xr_ref['region'] == 'R0',drop=True).sum()
                                             , sign = '>=' , rhs = SPP_quota_solar_R0, name='constr_SPP_quota_solar_r0')

constr_SPP_quota_solar_r1 = m.add_constraints(lhs = (cap_SPP_solar).where(xr_ref['region'] == 'R1',drop=True).sum()
                                             , sign = '>=' , rhs = SPP_quota_solar_R1, name='constr_SPP_quota_solar_r1')

constr_SPP_quota_solar_r2 = m.add_constraints(lhs = (cap_SPP_solar).where(xr_ref['region'] == 'R2',drop=True).sum()
                                             , sign = '>=' , rhs = SPP_quota_solar_R2, name='constr_SPP_quota_solar_r2')

constr_SPP_quota_solar_r3 = m.add_constraints(lhs = (cap_SPP_solar).where(xr_ref['region'] == 'R3',drop=True).sum()
                                             , sign = '>=' , rhs = SPP_quota_solar_R3, name='constr_SPP_quota_solar_r3')

constr_SPP_quota_solar_r4 = m.add_constraints(lhs = (cap_SPP_solar).where(xr_ref['region'] == 'R4',drop=True).sum()
                                             , sign = '>=' , rhs = SPP_quota_solar_R4, name='constr_SPP_quota_solar_r4')

constr_VSPP_quota_solar = m.add_constraints((cap_VSPP_solar.sum()) >= VSPP_quota_solar_total, name='constr_VSPP_quota_solar')

constr_VSPP_quota_solar_r0 = m.add_constraints(lhs = (cap_VSPP_solar).where(xr_ref['region'] == 'R0',drop=True).sum()
                                             , sign = '>=' , rhs = VSPP_quota_solar_R0, name='constr_VSPP_quota_solar_r0')

constr_VSPP_quota_solar_r1 = m.add_constraints(lhs = (cap_VSPP_solar).where(xr_ref['region'] == 'R1',drop=True).sum()
                                             , sign = '>=' , rhs = VSPP_quota_solar_R1, name='constr_VSPP_quota_solar_r1')

constr_VSPP_quota_solar_r2 = m.add_constraints(lhs = (cap_VSPP_solar).where(xr_ref['region'] == 'R2',drop=True).sum()
                                             , sign = '>=' , rhs = VSPP_quota_solar_R2, name='constr_VSPP_quota_solar_r2')

constr_VSPP_quota_solar_r3 = m.add_constraints(lhs = (cap_VSPP_solar).where(xr_ref['region'] == 'R3',drop=True).sum()
                                             , sign = '>=' , rhs = VSPP_quota_solar_R3, name='constr_VSPP_quota_solar_r3')

constr_VSPP_quota_solar_r4 = m.add_constraints(lhs = (cap_VSPP_solar).where(xr_ref['region'] == 'R4',drop=True).sum()
                                             , sign = '>=' , rhs = VSPP_quota_solar_R4, name='constr_VSPP_quota_solar_r4')

############################################################################################################################################################

constr_SPP_quota_biomass = m.add_constraints((cap_SPP_biomass.sum()) >= SPP_quota_biomass_total, name='constr_SPP_quota_biomass')

constr_SPP_quota_biomass_r0 = m.add_constraints(lhs = (cap_SPP_biomass).where(xr_ref['region'] == 'R0',drop=True).sum()
                                             , sign = '>=' , rhs = SPP_quota_biomass_R0, name='constr_SPP_quota_biomass_r0')

constr_SPP_quota_biomass_r1 = m.add_constraints(lhs = (cap_SPP_biomass).where(xr_ref['region'] == 'R1',drop=True).sum()
                                             , sign = '>=' , rhs = SPP_quota_biomass_R1, name='constr_SPP_quota_biomass_r1')

constr_SPP_quota_biomass_r2 = m.add_constraints(lhs = (cap_SPP_biomass).where(xr_ref['region'] == 'R2',drop=True).sum()
                                             , sign = '>=' , rhs = SPP_quota_biomass_R2, name='constr_SPP_quota_biomass_r2')

constr_SPP_quota_biomass_r3 = m.add_constraints(lhs = (cap_SPP_biomass).where(xr_ref['region'] == 'R3',drop=True).sum()
                                             , sign = '>=' , rhs = SPP_quota_biomass_R3, name='constr_SPP_quota_biomass_r3')

constr_SPP_quota_biomass_r4 = m.add_constraints(lhs = (cap_SPP_biomass).where(xr_ref['region'] == 'R4',drop=True).sum()
                                             , sign = '>=' , rhs = SPP_quota_biomass_R4, name='constr_SPP_quota_biomass_r4')

constr_VSPP_quota_biomass = m.add_constraints((cap_VSPP_biomass.sum()) >= VSPP_quota_biomass_total, name='constr_VSPP_quota_biomass')

constr_VSPP_quota_biomass_r0 = m.add_constraints(lhs = (cap_VSPP_biomass).where(xr_ref['region'] == 'R0',drop=True).sum(),
                                                 sign = '>=' , rhs = VSPP_quota_biomass_R0, name='constr_VSPP_quota_biomass_r0')

constr_VSPP_quota_biomass_r1 = m.add_constraints(lhs = (cap_VSPP_biomass).where(xr_ref['region'] == 'R1',drop=True).sum()
                                                 , sign = '>=' , rhs = VSPP_quota_biomass_R1, name='constr_VSPP_quota_biomass_r1')

constr_VSPP_quota_biomass_r2 = m.add_constraints(lhs = (cap_VSPP_biomass).where(xr_ref['region'] == 'R2',drop=True).sum()
                                                 , sign = '>=' , rhs = VSPP_quota_biomass_R2, name='constr_VSPP_quota_biomass_r2')

constr_VSPP_quota_biomass_r3 = m.add_constraints(lhs = (cap_VSPP_biomass).where(xr_ref['region'] == 'R3',drop=True).sum()
                                                 , sign = '>=' , rhs = VSPP_quota_biomass_R3, name='constr_VSPP_quota_biomass_r3')

constr_VSPP_quota_biomass_r4 = m.add_constraints(lhs = (cap_VSPP_biomass).where(xr_ref['region'] == 'R4',drop=True).sum()
                                                 , sign = '>=' , rhs = VSPP_quota_biomass_R4, name='constr_VSPP_quota_biomass_r4')

############################################################################################################################################################

constr_SPP_quota_bgec = m.add_constraints((cap_SPP_bgec.sum()) >= SPP_quota_bgec_total, name='constr_SPP_quota_bgec')

constr_SPP_quota_bgec_r0 = m.add_constraints(lhs = (cap_SPP_bgec).where(xr_ref['region'] == 'R0',drop=True).sum()
                                             , sign = '>=' , rhs = SPP_quota_bgec_R0, name='constr_SPP_quota_bgec_r0')

constr_SPP_quota_bgec_r1 = m.add_constraints(lhs = (cap_SPP_bgec).where(xr_ref['region'] == 'R1',drop=True).sum()
                                             , sign = '>=' , rhs = SPP_quota_bgec_R1, name='constr_SPP_quota_bgec_r1')

constr_SPP_quota_bgec_r2 = m.add_constraints(lhs = (cap_SPP_bgec).where(xr_ref['region'] == 'R2',drop=True).sum()
                                             , sign = '>=' , rhs = SPP_quota_bgec_R2, name='constr_SPP_quota_bgec_r2')

constr_SPP_quota_bgec_r3 = m.add_constraints(lhs = (cap_SPP_bgec).where(xr_ref['region'] == 'R3',drop=True).sum()
                                             , sign = '>=' , rhs = SPP_quota_bgec_R3, name='constr_SPP_quota_bgec_r3')

constr_SPP_quota_bgec_r4 = m.add_constraints(lhs = (cap_SPP_bgec).where(xr_ref['region'] == 'R4',drop=True).sum()
                                             , sign = '>=' , rhs = SPP_quota_bgec_R4, name='constr_SPP_quota_bgec_r4')

constr_VSPP_quota_bgec = m.add_constraints((cap_VSPP_bgec.sum()) >= VSPP_quota_bgec_total, name='constr_VSPP_quota_bgec')

constr_VSPP_quota_bgec_r0 = m.add_constraints(lhs = (cap_VSPP_bgec).where(xr_ref['region'] == 'R0',drop=True).sum()
                                            , sign = '>=' , rhs = VSPP_quota_bgec_R0, name='constr_VSPP_quota_bgec_r0')

constr_VSPP_quota_bgec_r1 = m.add_constraints(lhs = (cap_VSPP_bgec).where(xr_ref['region'] == 'R1',drop=True).sum()
                                            , sign = '>=' , rhs = VSPP_quota_bgec_R1, name='constr_VSPP_quota_bgec_r1')

constr_VSPP_quota_bgec_r2 = m.add_constraints(lhs = (cap_VSPP_bgec).where(xr_ref['region'] == 'R2',drop=True).sum()
                                            , sign = '>=' , rhs = VSPP_quota_bgec_R2, name='constr_VSPP_quota_bgec_r2')

constr_VSPP_quota_bgec_r3 = m.add_constraints(lhs = (cap_VSPP_bgec).where(xr_ref['region'] == 'R3',drop=True).sum()
                                            , sign = '>=' , rhs = VSPP_quota_bgec_R3, name='constr_VSPP_quota_bgec_r3')

constr_VSPP_quota_bgec_r4 = m.add_constraints(lhs = (cap_VSPP_bgec).where(xr_ref['region'] == 'R4',drop=True).sum()
                                            , sign = '>=' , rhs = VSPP_quota_bgec_R4, name='constr_VSPP_quota_bgec_r4')

############################################################################################################################################################

constr_SPP_quota_msw = m.add_constraints((cap_SPP_msw.sum()) >= SPP_quota_msw_total, name='constr_SPP_quota_msw')

constr_SPP_quota_msw_r0 = m.add_constraints(lhs = (cap_SPP_msw).where(xr_ref['region'] == 'R0',drop=True).sum()
                                            , sign = '>=' , rhs = SPP_quota_msw_R0, name='constr_SPP_quota_msw_r0')

constr_SPP_quota_msw_r1 = m.add_constraints(lhs = (cap_SPP_msw).where(xr_ref['region'] == 'R1',drop=True).sum()
                                            , sign = '>=' , rhs = SPP_quota_msw_R1, name='constr_SPP_quota_msw_r1')

constr_SPP_quota_msw_r2 = m.add_constraints(lhs = (cap_SPP_msw).where(xr_ref['region'] == 'R2',drop=True).sum()
                                            , sign = '>=' , rhs = SPP_quota_msw_R2, name='constr_SPP_quota_msw_r2')

constr_SPP_quota_msw_r3 = m.add_constraints(lhs = (cap_SPP_msw).where(xr_ref['region'] == 'R3',drop=True).sum()
                                            , sign = '>=' , rhs = SPP_quota_msw_R3, name='constr_SPP_quota_msw_r3')

constr_SPP_quota_msw_r4 = m.add_constraints(lhs = (cap_SPP_msw).where(xr_ref['region'] == 'R4',drop=True).sum()
                                            , sign = '>=' , rhs = SPP_quota_msw_R4, name='constr_SPP_quota_msw_r4')

constr_VSPP_quota_msw = m.add_constraints((cap_VSPP_msw.sum()) >= VSPP_quota_msw_total, name='constr_VSPP_quota_msw')

constr_VSPP_quota_msw_r0 = m.add_constraints(lhs = (cap_VSPP_msw).where(xr_ref['region'] == 'R0',drop=True).sum()
                                              , sign = '>=' , rhs = VSPP_quota_msw_R0, name='constr_VSPP_quota_msw_r0')

constr_VSPP_quota_msw_r1 = m.add_constraints(lhs = (cap_VSPP_msw).where(xr_ref['region'] == 'R1',drop=True).sum()
                                              , sign = '>=' , rhs = VSPP_quota_msw_R1, name='constr_VSPP_quota_msw_r1')

constr_VSPP_quota_msw_r2 = m.add_constraints(lhs = (cap_VSPP_msw).where(xr_ref['region'] == 'R2',drop=True).sum()
                                              , sign = '>=' , rhs = VSPP_quota_msw_R2, name='constr_VSPP_quota_msw_r2')

constr_VSPP_quota_msw_r3 = m.add_constraints(lhs = (cap_VSPP_msw).where(xr_ref['region'] == 'R3',drop=True).sum()
                                              , sign = '>=' , rhs = VSPP_quota_msw_R3, name='constr_VSPP_quota_msw_r3')

constr_VSPP_quota_msw_r4 = m.add_constraints(lhs = (cap_VSPP_msw).where(xr_ref['region'] == 'R4',drop=True).sum()
                                              , sign = '>=' , rhs = VSPP_quota_msw_R4, name='constr_VSPP_quota_msw_r4')

############################################################################################################################################################

obj = (-1000) * (
    ########################################################################
    ( 
     xr_ref['SI_Wind'].rolling(lon = rollingwindow_spp_wind, lat = rollingwindow_spp_wind, min_periods=1,center=True).sum()
     * (cap_SPP_wind) 
     / maxcapacityfor_spp
     )
    +
    ( 
     xr_ref['SI_Wind'].rolling(lon = rollingwindow_vspp_wind, lat = rollingwindow_vspp_wind, min_periods=1,center=True).sum()
     * (cap_VSPP_wind) 
     / maxcapacityfor_vspp
     )
    ########################################################################
    +
    ( 
     xr_ref['SI_Solar'].rolling(lon = rollingwindow_spp_solar, lat = rollingwindow_spp_solar, min_periods=1,center=True).sum()
     * (cap_SPP_solar) 
     / maxcapacityfor_spp
     )
    +
    ( 
     xr_ref['SI_Solar'].rolling(lon = rollingwindow_vspp_solar, lat = rollingwindow_vspp_solar, min_periods=1,center=True).sum()
     * (cap_VSPP_solar) 
     / maxcapacityfor_vspp
     )
    ########################################################################
    +
    ( 
     xr_ref['SI_Biomass'].rolling(lon = rollingwindow_spp_biomass, lat = rollingwindow_spp_biomass, min_periods=1,center=True).sum()
     * (cap_SPP_biomass) 
     / maxcapacityfor_spp
     )
    +
    ( 
     xr_ref['SI_Biomass'].rolling(lon = rollingwindow_spp_biomass, lat = rollingwindow_spp_biomass, min_periods=1,center=True).sum()
     * (cap_SPP_biomass) 
     / maxcapacityfor_vspp
     )
    ########################################################################
    +
    ( 
     xr_ref['SI_BGEC'].rolling(lon = rollingwindow_spp_bgec, lat = rollingwindow_spp_bgec, min_periods=1,center=True).sum()
     * (cap_SPP_bgec) 
     / maxcapacityfor_spp
     )
    +
    ( 
     xr_ref['SI_BGEC'].rolling(lon = rollingwindow_spp_bgec, lat = rollingwindow_spp_bgec, min_periods=1,center=True).sum()
     * (cap_SPP_bgec) 
     / maxcapacityfor_vspp
     )
    ########################################################################
    +
    ( 
     xr_ref['SI_MSW'].rolling(lon = rollingwindow_spp_msw, lat = rollingwindow_spp_msw, min_periods=1,center=True).sum()
     * (cap_SPP_msw) 
     / maxcapacityfor_spp
     )
    +
    ( 
     xr_ref['SI_MSW'].rolling(lon = rollingwindow_spp_msw, lat = rollingwindow_spp_msw, min_periods=1,center=True).sum()
     * (cap_SPP_msw) 
     / maxcapacityfor_vspp
     )
)   

m.add_objective(obj)

print("presolve = ",m)
m.solve(solver_name='highs',
        mip_abs_gap = 0.002,
        mip_rel_gap = 0.002,
        )

print('aftersolve = ',m)
solution = m.solution
solution = solution.fillna(0)
print(solution)

xr_ref['cap_SPP_wind'] = solution['cap_SPP_wind']
xr_ref['cap_VSPP_wind'] = solution['cap_VSPP_wind']
xr_ref['cap_SPP_solar'] = solution['cap_SPP_solar']
xr_ref['cap_VSPP_solar'] = solution['cap_VSPP_solar']
xr_ref['cap_SPP_biomass'] = solution['cap_SPP_biomass']
xr_ref['cap_VSPP_biomass'] = solution['cap_VSPP_biomass']
xr_ref['cap_SPP_bgec'] = solution['cap_SPP_bgec']
xr_ref['cap_VSPP_bgec'] = solution['cap_VSPP_bgec']
xr_ref['cap_SPP_msw'] = solution['cap_SPP_msw']
xr_ref['cap_VSPP_msw'] = solution['cap_VSPP_msw']

print("cap_SPP_wind = ",xr_ref['cap_SPP_wind'].sum())
print("cap_VSPP_wind = ",xr_ref['cap_VSPP_wind'].sum())
print("cap_SPP_solar = ",xr_ref['cap_SPP_solar'].sum())
print("cap_VSPP_solar = ",xr_ref['cap_VSPP_solar'].sum())
print("cap_SPP_biomass = ",xr_ref['cap_SPP_biomass'].sum())
print("cap_VSPP_biomass = ",xr_ref['cap_VSPP_biomass'].sum())
print("cap_SPP_bgec = ",xr_ref['cap_SPP_bgec'].sum())
print("cap_VSPP_bgec = ",xr_ref['cap_VSPP_bgec'].sum())
print("cap_SPP_msw = ",xr_ref['cap_SPP_msw'].sum())
print("cap_VSPP_msw = ",xr_ref['cap_VSPP_msw'].sum())

xr_ref.to_netcdf(path='Output\\xr_output_all.nc')


