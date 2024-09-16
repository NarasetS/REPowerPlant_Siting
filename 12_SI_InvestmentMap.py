import linopy
import pandas as pd
import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
from pathlib import Path
from shapely import wkt
import numpy as np

###################### Config ################################################################################################
coarsenscale = 1
lccs_resolution = 300 * coarsenscale #m
areapergrid = (lccs_resolution/1000) ** 2 ## km2
scenario_SI = 0  ## Include area where SI >= scenario_SI

mwperkm2_wind = 4.5 ## originally 9 MW/km2 but deduct by 50% of the technically available from IEA's Thailand CET
mwperkm2_solar = 15 ## originally 30 MW/km2 but deduct by 50% of the technically available from IEA's Thailand CET

mwpergrid_wind = np.round(areapergrid * mwperkm2_wind,2)
mwpergrid_solar =  np.round(areapergrid * mwperkm2_solar,2)

suitablearea_biomass = 2500 ## km2 
suitablearea_bgec = 2500 ## km2
suitablearea_msw = 2500 ## km2

rollingwindow_biomass = int(np.ceil(np.sqrt(suitablearea_biomass/areapergrid)))
rollingwindow_bgec = int(np.ceil(np.sqrt(suitablearea_bgec/areapergrid)))
rollingwindow_msw = int(np.ceil(np.sqrt(suitablearea_msw/areapergrid)))

print('areapergrid = ',areapergrid)
print('mwpergrid_wind = ',mwpergrid_wind)
print('mwpergrid_solar = ',mwpergrid_solar)
print('rollingwindow_biomass = ',rollingwindow_biomass,' * ',rollingwindow_biomass)
print('rollingwindow_bgec = ',rollingwindow_bgec,' * ',rollingwindow_bgec)
print('rollingwindow_msw = ',rollingwindow_msw,' * ',rollingwindow_msw)

xr_final_SI_raw = xr.open_dataset('Output\\xr_final_SI_all.nc')
xr_final_SI_raw = xr_final_SI_raw.drop_vars('ADM1_EN')
# xr_final_SI_raw = xr_final_SI_raw.drop_vars('SI_BGEC')
xr_final_SI_raw = xr_final_SI_raw.drop_vars('SI_BGWW')
# xr_final_SI_raw = xr_final_SI_raw.drop_vars('SI_MSW')
xr_final_SI_raw = xr_final_SI_raw.drop_vars('SI_IEW')
# xr_final_SI_raw = xr_final_SI_raw.drop_vars('A_BGEC')
xr_final_SI_raw = xr_final_SI_raw.drop_vars('A_BGWW')
# xr_final_SI_raw = xr_final_SI_raw.drop_vars('A_MSW')
xr_final_SI_raw = xr_final_SI_raw.drop_vars('A_IEW')
# xr_final_SI_raw = xr_final_SI_raw.drop_vars('AVA_BGEC')
xr_final_SI_raw = xr_final_SI_raw.drop_vars('AVA_BGWW')
# xr_final_SI_raw = xr_final_SI_raw.drop_vars('AVA_MSW')
xr_final_SI_raw = xr_final_SI_raw.drop_vars('AVA_IEW')

xr_final_SI_raw['SI_Wind'] = xr.where(xr_final_SI_raw['SI_Wind'] >= scenario_SI ,xr_final_SI_raw['SI_Wind'],0)
xr_final_SI_raw['AVA_Wind'] = xr.where(xr_final_SI_raw['SI_Wind'] == 0 , 0 , xr_final_SI_raw['AVA_Wind'])

xr_final_SI_raw['SI_Solar'] = xr.where(xr_final_SI_raw['SI_Solar'] >= scenario_SI ,xr_final_SI_raw['SI_Solar'],0)
xr_final_SI_raw['AVA_Solar'] = xr.where(xr_final_SI_raw['SI_Solar'] == 0 , 0 , xr_final_SI_raw['AVA_Solar'])

xr_final_SI_raw['SI_Biomass'] = xr.where(xr_final_SI_raw['SI_Biomass'] >= scenario_SI ,xr_final_SI_raw['SI_Biomass'],0)
xr_final_SI_raw['AVA_Biomass'] = xr.where(xr_final_SI_raw['SI_Biomass'] == 0 , 0 , xr_final_SI_raw['AVA_Biomass'])

xr_final_SI_raw['SI_BGEC'] = xr.where(xr_final_SI_raw['SI_BGEC'] >= scenario_SI ,xr_final_SI_raw['SI_BGEC'],0)
xr_final_SI_raw['AVA_BGEC'] = xr.where(xr_final_SI_raw['SI_BGEC'] == 0 , 0 , xr_final_SI_raw['AVA_BGEC'])

xr_final_SI_raw['SI_MSW'] = xr.where(xr_final_SI_raw['SI_MSW'] >= scenario_SI ,xr_final_SI_raw['SI_MSW'],0)
xr_final_SI_raw['AVA_MSW'] = xr.where(xr_final_SI_raw['SI_MSW'] == 0 , 0 , xr_final_SI_raw['AVA_MSW'])

xr_final_SI = xr_final_SI_raw.coarsen(lat = coarsenscale, lon= coarsenscale, boundary='trim').sum()
xr_final_SI['SI_Wind'] = xr_final_SI['SI_Wind'] / (coarsenscale**2)
xr_final_SI['SI_Solar'] = xr_final_SI['SI_Solar'] / (coarsenscale**2)
xr_final_SI['SI_Biomass'] = xr_final_SI['SI_Biomass'] / (coarsenscale**2)
xr_final_SI['SI_BGEC'] = xr_final_SI['SI_BGEC'] / (coarsenscale**2)
xr_final_SI['SI_MSW'] = xr_final_SI['SI_MSW'] / (coarsenscale**2)

######### Find min/max SI within cell #############
# xr_final_SI['SI_Wind_max'] = xr_final_SI_raw['SI_Wind'].coarsen(lat = coarsenscale, lon= coarsenscale, boundary='trim').max()
# xr_final_SI['SI_Wind_min'] = xr_final_SI_raw['SI_Wind'].coarsen(lat = coarsenscale, lon= coarsenscale, boundary='trim').min()
# xr_final_SI['SI_Solar_max'] = xr_final_SI_raw['SI_Solar'].coarsen(lat = coarsenscale, lon= coarsenscale, boundary='trim').max()
# xr_final_SI['SI_Solar_min'] = xr_final_SI_raw['SI_Solar'].coarsen(lat = coarsenscale, lon= coarsenscale, boundary='trim').min()

maxcap_biomass = xr_final_SI['A_Biomass'].rolling(lon = rollingwindow_biomass, lat = rollingwindow_biomass, min_periods=1,center=True).sum().where(xr_final_SI['SI_Biomass']>0).max()
maxcap_bgec = xr_final_SI['A_BGEC'].rolling(lon = rollingwindow_bgec, lat = rollingwindow_bgec, min_periods=1,center=True).sum().where(xr_final_SI['SI_BGEC']>0).max()
maxcap_msw = xr_final_SI['A_MSW'].rolling(lon = rollingwindow_msw, lat = rollingwindow_msw, min_periods=1,center=True).sum().where(xr_final_SI['SI_MSW']>0).max()
print('maxcap_biomass = ',maxcap_biomass)
print('maxcap_bgec = ',maxcap_bgec)
print('maxcap_msw = ',maxcap_msw)


xr_final_SI_raw.close()
print(xr_final_SI.data_vars)
print("AVA Wind = ",xr_final_SI['AVA_Wind'].sum())
print("AVA Solar = ",xr_final_SI['AVA_Solar'].sum())
print("A_Biomass = ",xr_final_SI['A_Biomass'].sum())
print("A_BGEC = ",xr_final_SI['A_BGEC'].sum())
print("A_MSW = ",xr_final_SI['A_MSW'].sum())
print('Max SI_Wind = ',xr_final_SI['SI_Wind'].max())
print('Max SI_Solar = ',xr_final_SI['SI_Solar'].max())
print('Max SI_Biomass = ',xr_final_SI['SI_Biomass'].max())
print('Max SI_BGEC = ',xr_final_SI['SI_BGEC'].max())
print('Max SI_MSW = ',xr_final_SI['SI_MSW'].max())
###################### Config ################################################################################################

######### Next I assign Region to xarray ################################################################################
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
#########################################################################################

###################### Summary ################################################################################################
xr_ref = xr_final_SI
print("AVA Wind = ",xr_ref['AVA_Wind'].sum())
print("AVA Solar = ",xr_ref['AVA_Solar'].sum()) 
print("A_Biomass = ",xr_ref['A_Biomass'].sum())
print("A_BGEC = ",xr_ref['A_BGEC'].sum())
print("A_MSW = ",xr_ref['A_MSW'].sum())
print('Mean SI Wind = ',xr_ref['SI_Wind'].where(xr_ref['SI_Wind']>0).mean())
print('Mean SI Solar = ',xr_ref['SI_Solar'].where(xr_ref['SI_Solar']>0).mean())
print('Mean SI_Biomass = ',xr_ref['SI_Biomass'].where(xr_ref['SI_Biomass']>0).mean())
print('Mean SI_BGEC = ',xr_ref['SI_BGEC'].where(xr_ref['SI_BGEC']>0).mean())
print('Mean SI_MSW = ',xr_ref['SI_MSW'].where(xr_ref['SI_MSW']>0).mean())
print('Max SI Wind = ',xr_ref['SI_Wind'].max())
print('Max SI Solar = ',xr_ref['SI_Solar'].max())
print('Max SI Biomass = ',xr_ref['SI_Biomass'].max())
print('Max SI BGEC = ',xr_ref['SI_BGEC'].max())
print('Max SI MSW = ',xr_ref['SI_MSW'].max())
print('coarsenscale = ',coarsenscale)
print('areapergrid = ',areapergrid)
print('mwpergrid_wind = ',mwpergrid_wind)
print('mwpergrid_solar = ',mwpergrid_solar)
###################### Summary ################################################################################################

####### PDP ############################################################
quota_wind_total =  6344.8
quota_wind_R0 = 0
quota_wind_R1 = 0
quota_wind_R2 = 6084.7
quota_wind_R3 =  260.1
quota_wind_R4 =  0

quota_solar_total =  18129
quota_solar_R0 =  23.34
quota_solar_R1 =  3255.36
quota_solar_R2 =  6072.28
quota_solar_R3 = 5602.08
quota_solar_R4 =  3175.94

# quota_biomass_total = 710.65
# quota_biomass_R0 = 0
# quota_biomass_R1 =  148.96
# quota_biomass_R2 =  260
# quota_biomass_R3 =   70.70
# quota_biomass_R4 =  230.99

# quota_bgec_total = 968.49
# quota_bgec_R0 = 0
# quota_bgec_R1 =  21.99
# quota_bgec_R2 =  420
# quota_bgec_R3 =   67.50
# quota_bgec_R4 =  262

# quota_msw_total = 843.50
# quota_msw_R0 = 195
# quota_msw_R1 =  219.52
# quota_msw_R2 =  233.65
# quota_msw_R3 =   104.33
# quota_msw_R4 =  91
####### PDP ############################################################

######################## model #####################################################
m = linopy.Model()

built_wind = m.add_variables(binary=True, coords=xr_ref.coords, name='built_wind')
cap_wind = m.add_variables(lower=0.00, coords=xr_ref.coords, name='cap_wind')

built_solar = m.add_variables(binary=True, coords=xr_ref.coords, name='built_solar')
cap_solar = m.add_variables(lower=0.00, coords=xr_ref.coords, name='cap_solar')

# built_biomass = m.add_variables(binary=True, coords=xr_ref.coords, name='built_biomass')
# cap_biomass = m.add_variables(lower=0.00, coords=xr_ref.coords, name='cap_biomass')

# built_bgec = m.add_variables(binary=True, coords=xr_ref.coords, name='built_bgec')
# cap_bgec = m.add_variables(lower=0.00, coords=xr_ref.coords, name='cap_bgec')

# built_msw = m.add_variables(binary=True, coords=xr_ref.coords, name='built_msw')
# cap_msw = m.add_variables(lower=0.00, coords=xr_ref.coords, name='cap_msw')

############################################ Constraint Building Location Logic ##############################################################################
constr_built_logic =  m.add_constraints(
    (
        built_wind
        +
        built_solar
        # +
        # built_biomass.rolling(lat = rollingwindow_biomass,min_periods=1,center=True).sum().rolling(lon = rollingwindow_biomass,min_periods=1,center=True).sum()
        # +
        # built_bgec.rolling(lat = rollingwindow_bgec,min_periods=1,center=True).sum().rolling(lon = rollingwindow_bgec,min_periods=1,center=True).sum()
        # +
        # built_msw.rolling(lat = rollingwindow_msw,min_periods=1,center=True).sum().rolling(lon = rollingwindow_msw,min_periods=1,center=True).sum()
    )
        <= 1
        , name='constr_built_logic'
)
############################################ Constraint Building Location Logic ##############################################################################

############################################ Constraint Capacity ##############################################################################
constr_maxcap_wind = m.add_constraints(
    cap_wind <= (built_wind * xr_ref['AVA_Wind'] * mwperkm2_wind)
    ,name = 'constr_maxcap_wind'
)
constr_mincap_wind = m.add_constraints(
    cap_wind >= (built_wind * 1)
    ,name = 'constr_mincap_wind'
)
constr_builtarea_wind = m.add_constraints(
    built_wind <= (xr_ref['AVA_Wind'] * 10000)
    ,name = 'constr_builtarea_wind'
)

constr_maxcap_solar = m.add_constraints(
    cap_solar <= (built_solar * xr_ref['AVA_Solar'] * mwperkm2_solar)
    ,name = 'constr_maxcap_solar'
)
constr_mincap_solar = m.add_constraints(
    cap_solar >= (built_solar * 1)
    ,name = 'constr_mincap_solar'
)
constr_builtarea_solar = m.add_constraints(
    built_solar <= (xr_ref['AVA_Solar'] * 10000)
    ,name = 'constr_builtarea_solar'
)

# constr_maxcap_biomass = m.add_constraints(
#     cap_biomass <= (built_biomass * (xr_ref['A_Biomass'].rolling(lon = rollingwindow_biomass, lat = rollingwindow_biomass, min_periods=1,center=True).sum()))
#     ,name = 'constr_maxcap_biomass'
# )
# # constr_mincap_biomass = m.add_constraints(
# #     cap_biomass >= (built_biomass * 1)
# #     ,name = 'constr_mincap_biomass'
# # )
# constr_builtarea_biomass = m.add_constraints(
#     built_biomass <= (xr_ref['AVA_Biomass'] * 10000)
#     ,name = 'constr_builtarea_biomass'
# )

# constr_maxcap_bgec = m.add_constraints(
#     cap_bgec <= (built_bgec * (xr_ref['A_BGEC'].rolling(lon = rollingwindow_bgec, lat = rollingwindow_bgec, min_periods=1,center=True).sum()))
#     ,name = 'constr_maxcap_bgec'
# )
# # constr_mincap_bgec = m.add_constraints(
# #     cap_bgec >= (built_bgec * 1)
# #     ,name = 'constr_mincap_bgec'
# # )
# constr_builtarea_bgec = m.add_constraints(
#     built_bgec <= (xr_ref['AVA_BGEC'] * 10000)
#     ,name = 'constr_builtarea_bgec'
# )

# constr_maxcap_msw = m.add_constraints(
#     cap_msw <= ((built_msw) * (xr_ref['A_MSW'].rolling(lon = rollingwindow_msw, lat = rollingwindow_msw, min_periods=1,center=True).sum()))
#     ,name = 'constr_maxcap_msw'
# )
# # constr_mincap_msw = m.add_constraints(
# #     cap_msw >= ((built_msw) * 1)
# #     ,name = 'constr_mincap_msw'
# # )
# constr_builtarea_msw = m.add_constraints(
#     built_msw <= (xr_ref['AVA_MSW'] * 10000)
#     ,name = 'constr_builtarea_msw'
# )
############################################ Constraint Capacity ##############################################################################

###########################################################################################################################################################

constr_quota_wind = m.add_constraints((cap_wind.sum()) == quota_wind_total , name='constr_quota_wind')

constr_quota_wind_r0 = m.add_constraints(lhs = (cap_wind).where((xr_ref['region'] == 'R0'),drop=True).sum()
                                             , sign = '>=' , rhs = quota_wind_R0, name='constr_quota_wind_r0')

constr_quota_wind_r1 = m.add_constraints(lhs = (cap_wind).where((xr_ref['region'] == 'R1'),drop=True).sum()
                                             , sign = '>=' , rhs = quota_wind_R1, name='constr_quota_wind_r1')

constr_quota_wind_r2 = m.add_constraints(lhs = (cap_wind).where((xr_ref['region'] == 'R2'),drop=True).sum()
                                             , sign = '>=' , rhs = quota_wind_R2, name='constr_quota_wind_r2')

constr_quota_wind_r3 = m.add_constraints(lhs = (cap_wind).where((xr_ref['region'] == 'R3'),drop=True).sum()
                                             , sign = '>=' , rhs = quota_wind_R3, name='constr_quota_wind_r3')

constr_quota_wind_r4 = m.add_constraints(lhs = (cap_wind).where((xr_ref['region'] == 'R4'),drop=True).sum()
                                             , sign = '>=' , rhs = quota_wind_R4, name='constr_quota_wind_r4')

############################################################################################################################################################

constr_quota_solar = m.add_constraints((cap_solar.sum()) == quota_solar_total, name='constr_quota_solar')

constr_quota_solar_r0 = m.add_constraints(lhs = (cap_solar).where((xr_ref['region'] == 'R0'),drop=True).sum()
                                             , sign = '>=' , rhs = quota_solar_R0, name='constr_quota_solar_r0')

constr_quota_solar_r1 = m.add_constraints(lhs = (cap_solar).where((xr_ref['region'] == 'R1'),drop=True).sum()
                                             , sign = '>=' , rhs = quota_solar_R1, name='constr_quota_solar_r1')

constr_quota_solar_r2 = m.add_constraints(lhs = (cap_solar).where((xr_ref['region'] == 'R2'),drop=True).sum()
                                             , sign = '>=' , rhs = quota_solar_R2, name='constr_quota_solar_r2')

constr_quota_solar_r3 = m.add_constraints(lhs = (cap_solar).where((xr_ref['region'] == 'R3'),drop=True).sum()
                                             , sign = '>=' , rhs = quota_solar_R3, name='constr_quota_solar_r3')

constr_quota_solar_r4 = m.add_constraints(lhs = (cap_solar).where((xr_ref['region'] == 'R4'),drop=True).sum()
                                             , sign = '>=' , rhs = quota_solar_R4, name='constr_quota_solar_r4')

###########################################################################################################################################################

# constr_quota_biomass = m.add_constraints((cap_biomass.sum()) == quota_biomass_total, name='constr_quota_biomass')

# constr_quota_biomass_r0 = m.add_constraints(lhs = (cap_biomass).where((xr_ref['region'] == 'R0'),drop=True).sum()
#                                              , sign = '>=' , rhs = quota_biomass_R0, name='constr_quota_biomass_r0')

# constr_quota_biomass_r1 = m.add_constraints(lhs = (cap_biomass).where((xr_ref['region'] == 'R1'),drop=True).sum()
#                                              , sign = '>=' , rhs = quota_biomass_R1, name='constr_quota_biomass_r1')

# constr_quota_biomass_r2 = m.add_constraints(lhs = (cap_biomass).where((xr_ref['region'] == 'R2'),drop=True).sum()
#                                              , sign = '>=' , rhs = quota_biomass_R2, name='constr_quota_biomass_r2')

# constr_quota_biomass_r3 = m.add_constraints(lhs = (cap_biomass).where((xr_ref['region'] == 'R3'),drop=True).sum()
#                                              , sign = '>=' , rhs = quota_biomass_R3, name='constr_quota_biomass_r3')

# constr_quota_biomass_r4 = m.add_constraints(lhs = (cap_biomass).where((xr_ref['region'] == 'R4'),drop=True).sum()
#                                              , sign = '>=' , rhs = quota_biomass_R4, name='constr_quota_biomass_r4')

##########################################################################################################################################################

# constr_quota_bgec = m.add_constraints((cap_bgec.sum()) == quota_bgec_total, name='constr_quota_bgec')

# constr_quota_bgec_r0 = m.add_constraints(lhs = (cap_bgec).where((xr_ref['region'] == 'R0'),drop=True).sum()
#                                              , sign = '>=' , rhs = quota_bgec_R0, name='constr_quota_bgec_r0')

# constr_quota_bgec_r1 = m.add_constraints(lhs = (cap_bgec).where((xr_ref['region'] == 'R1'),drop=True).sum()
#                                              , sign = '>=' , rhs = quota_bgec_R1, name='constr_quota_bgec_r1')

# constr_quota_bgec_r2 = m.add_constraints(lhs = (cap_bgec).where((xr_ref['region'] == 'R2'),drop=True).sum()
#                                              , sign = '>=' , rhs = quota_bgec_R2, name='constr_quota_bgec_r2')

# constr_quota_bgec_r3 = m.add_constraints(lhs = (cap_bgec).where((xr_ref['region'] == 'R3'),drop=True).sum()
#                                              , sign = '>=' , rhs = quota_bgec_R3, name='constr_quota_bgec_r3')

# constr_quota_bgec_r4 = m.add_constraints(lhs = (cap_bgec).where((xr_ref['region'] == 'R4'),drop=True).sum()
#                                              , sign = '>=' , rhs = quota_bgec_R4, name='constr_quota_bgec_r4')

##########################################################################################################################################################

# constr_quota_msw = m.add_constraints((cap_msw.sum()) == quota_msw_total, name='constr_quota_msw')

# constr_quota_msw_r0 = m.add_constraints(lhs = (cap_msw).where((xr_ref['region'] == 'R0'),drop=True).sum()
#                                              , sign = '>=' , rhs = quota_msw_R0, name='constr_quota_msw_r0')

# constr_quota_msw_r1 = m.add_constraints(lhs = (cap_msw).where((xr_ref['region'] == 'R1'),drop=True).sum()
#                                              , sign = '>=' , rhs = quota_msw_R1, name='constr_quota_msw_r1')

# constr_quota_msw_r2 = m.add_constraints(lhs = (cap_msw).where((xr_ref['region'] == 'R2'),drop=True).sum()
#                                              , sign = '>=' , rhs = quota_msw_R2, name='constr_quota_msw_r2')

# constr_quota_msw_r3 = m.add_constraints(lhs = (cap_msw).where((xr_ref['region'] == 'R3'),drop=True).sum()
#                                              , sign = '>=' , rhs = quota_msw_R3, name='constr_quota_msw_r3')

# constr_quota_msw_r4 = m.add_constraints(lhs = (cap_msw).where((xr_ref['region'] == 'R4'),drop=True).sum()
#                                              , sign = '>=' , rhs = quota_msw_R4, name='constr_quota_msw_r4')

##########################################################################################################################################################

###################### Objective function ################################################################################################
obj = (-1000) * (

    ( xr_ref['SI_Wind'] * (cap_wind  / mwpergrid_wind))
    +
    ( xr_ref['SI_Solar'] * (cap_solar / mwpergrid_solar))
#     +
#     ( 
#         xr_ref['SI_Biomass'].rolling(lon = rollingwindow_biomass, lat = rollingwindow_biomass, min_periods=1,center=True).sum() 
#         * (cap_biomass / maxcap_biomass)
#      )
#     +
#     ( 
#         xr_ref['SI_BGEC'].rolling(lon = rollingwindow_bgec, lat = rollingwindow_bgec, min_periods=1,center=True).sum() 
#         * (cap_bgec / maxcap_bgec)
#      )
#     +
#     ( 
#         xr_ref['SI_MSW'].rolling(lon = rollingwindow_msw, lat = rollingwindow_msw, min_periods=1,center=True).sum() 
#         * (cap_msw / maxcap_msw)
#      )
)

m.add_objective(obj)
###################### Objective function ################################################################################################

###################### Solver ################################################################################################
print("presolve = ",m)
m.solve(solver_name='highs',
        mip_abs_gap = 0.03,
        mip_rel_gap = 0.03,
        )

print('aftersolve = ',m)
solution = m.solution
solution = solution.fillna(0)
print(solution)

xr_ref['cap_wind'] = np.round(solution['cap_wind'] * solution['built_wind'],1)
xr_ref['cap_solar'] = np.round(solution['cap_solar'] * solution['built_solar'],1)
# xr_ref['cap_biomass'] = np.round(solution['cap_biomass'] * solution['built_biomass'],1)
# xr_ref['cap_bgec'] = np.round(solution['cap_bgec'] * solution['built_bgec'],1)
# xr_ref['cap_msw'] = np.round(solution['cap_msw'] * solution['built_msw'],1)
# xr_ref['cap_solar'] = np.round(solution['cap_solar'] * solution['built_solar'],1)


print("cap_wind = ",xr_ref['cap_wind'].sum())
print("  R0 cap_wind = ",xr_ref['cap_wind'].where(xr_ref['region'] == 'R0').sum())
print("  R1 cap_wind = ",xr_ref['cap_wind'].where(xr_ref['region'] == 'R1').sum())
print("  R2 cap_wind = ",xr_ref['cap_wind'].where(xr_ref['region'] == 'R2').sum())
print("  R3 cap_wind = ",xr_ref['cap_wind'].where(xr_ref['region'] == 'R3').sum())
print("  R4 cap_wind = ",xr_ref['cap_wind'].where(xr_ref['region'] == 'R4').sum())

print("cap_solar = ",xr_ref['cap_solar'].sum())
print("  R0 cap_solar = ",xr_ref['cap_solar'].where(xr_ref['region'] == 'R0').sum())
print("  R1 cap_solar = ",xr_ref['cap_solar'].where(xr_ref['region'] == 'R1').sum())
print("  R2 cap_solar = ",xr_ref['cap_solar'].where(xr_ref['region'] == 'R2').sum())
print("  R3 cap_solar = ",xr_ref['cap_solar'].where(xr_ref['region'] == 'R3').sum())
print("  R4 cap_solar = ",xr_ref['cap_solar'].where(xr_ref['region'] == 'R4').sum())

# print("cap_biomass = ",xr_ref['cap_biomass'].sum())
# print("  R0 cap_biomass = ",xr_ref['cap_biomass'].where(xr_ref['region'] == 'R0').sum())
# print("  R1 cap_biomass = ",xr_ref['cap_biomass'].where(xr_ref['region'] == 'R1').sum())
# print("  R2 cap_biomass = ",xr_ref['cap_biomass'].where(xr_ref['region'] == 'R2').sum())
# print("  R3 cap_biomass = ",xr_ref['cap_biomass'].where(xr_ref['region'] == 'R3').sum())
# print("  R4 cap_biomass = ",xr_ref['cap_biomass'].where(xr_ref['region'] == 'R4').sum())

# print("cap_bgec = ",xr_ref['cap_bgec'].sum())
# print("  R0 cap_bgec = ",xr_ref['cap_bgec'].where(xr_ref['region'] == 'R0').sum())
# print("  R1 cap_bgec = ",xr_ref['cap_bgec'].where(xr_ref['region'] == 'R1').sum())
# print("  R2 cap_bgec = ",xr_ref['cap_bgec'].where(xr_ref['region'] == 'R2').sum())
# print("  R3 cap_bgec = ",xr_ref['cap_bgec'].where(xr_ref['region'] == 'R3').sum())
# print("  R4 cap_bgec = ",xr_ref['cap_bgec'].where(xr_ref['region'] == 'R4').sum())

# print("cap_msw = ",xr_ref['cap_msw'].sum())
# print("  R0 cap_msw = ",xr_ref['cap_msw'].where(xr_ref['region'] == 'R0').sum())
# print("  R1 cap_msw = ",xr_ref['cap_msw'].where(xr_ref['region'] == 'R1').sum())
# print("  R2 cap_msw = ",xr_ref['cap_msw'].where(xr_ref['region'] == 'R2').sum())
# print("  R3 cap_msw = ",xr_ref['cap_msw'].where(xr_ref['region'] == 'R3').sum())
# print("  R4 cap_msw = ",xr_ref['cap_msw'].where(xr_ref['region'] == 'R4').sum())

print(xr_ref.data_vars)
xr_ref.to_netcdf(path='Output\\xr_output_all_SSI_' + str(scenario_SI) + "_CS_"+str(coarsenscale)+ '_.nc')