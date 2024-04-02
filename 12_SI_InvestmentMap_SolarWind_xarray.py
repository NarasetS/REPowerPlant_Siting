import linopy
import pandas as pd
import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
from pathlib import Path
from shapely import wkt
import numpy as np

coarsenscale = 3
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

############################################## Coarsen #########################################################################################################################
xr_final_SI = xr.open_dataset('Output\\xr_final_SI.nc')
xr_final_SI = xr_final_SI.drop_vars('ADM1_EN')

xr_final_SI = xr_final_SI.coarsen(lat = coarsenscale, lon= coarsenscale, boundary='pad').sum()
xr_final_SI['SI_Wind'] = xr_final_SI['SI_Wind'] / (coarsenscale**2)
xr_final_SI['SI_Solar'] = xr_final_SI['SI_Solar'] / (coarsenscale**2)
print("AVA Wind = ",xr_final_SI['AVA_Wind'].sum())
print("AVA Solar = ",xr_final_SI['AVA_Solar'].sum())
print('Max SI Wind = ',xr_final_SI['SI_Wind'].max())
print('Max SI Solar = ',xr_final_SI['SI_Solar'].max())
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
print("max AVA Wind = ",xr_final_SI['AVA_Wind'].max())
print("max AVA Solar = ",xr_final_SI['AVA_Solar'].max())
print('Max SI Wind = ',xr_final_SI['SI_Wind'].max())
print('Max SI Solar = ',xr_final_SI['SI_Solar'].max())
xr_ref = xr_final_SI
print('coarsenscale = ',coarsenscale)
print('areapergrid = ',areapergrid)
print('mwpergrid_wind = ',mwpergrid_wind)
print('mwpergrid_solar = ',mwpergrid_solar)
print('rollingwindow_spp_wind = ',rollingwindow_spp_wind,' * ',rollingwindow_spp_wind,' = ', mwpergrid_wind * (rollingwindow_spp_wind**2))
print('rollingwindow_vspp_wind = ',rollingwindow_vspp_wind,' * ',rollingwindow_vspp_wind,' = ', mwpergrid_wind * (rollingwindow_vspp_wind**2))
print('rollingwindow_spp_solar = ',rollingwindow_spp_solar,' * ',rollingwindow_spp_solar,' = ', mwpergrid_solar * (rollingwindow_spp_solar**2))
print('rollingwindow_vspp_solar = ',rollingwindow_vspp_solar,' * ',rollingwindow_vspp_solar,' = ', mwpergrid_solar * (rollingwindow_vspp_solar**2))
############################################## Summary #########################################################################################################################

# ######## PDP ############################################################
# SPP_quota_wind_total =  0
# SPP_quota_wind_R0 = 0
# SPP_quota_wind_R1 = 0
# SPP_quota_wind_R2 = 0
# SPP_quota_wind_R3 =  0
# SPP_quota_wind_R4 =  0

# VSPP_quota_wind_total =  0
# VSPP_quota_wind_R0 =  0
# VSPP_quota_wind_R1 =  0
# VSPP_quota_wind_R2 =  0
# VSPP_quota_wind_R3 =  0
# VSPP_quota_wind_R4 = 0

# SPP_quota_solar_total =  0
# SPP_quota_solar_R0 =  0
# SPP_quota_solar_R1 =  0
# SPP_quota_solar_R2 =  0
# SPP_quota_solar_R3 = 0
# SPP_quota_solar_R4 =  0

# VSPP_quota_solar_total = 0
# VSPP_quota_solar_R0 = 0
# VSPP_quota_solar_R1 =  0
# VSPP_quota_solar_R2 =  0
# VSPP_quota_solar_R3 =   0
# VSPP_quota_solar_R4 =  0
# ######## PDP ############################################################

######## PDP ############################################################
SPP_quota_wind_total =  6335
SPP_quota_wind_R0 = 0
SPP_quota_wind_R1 = 0
SPP_quota_wind_R2 = 6075
SPP_quota_wind_R3 =  260
SPP_quota_wind_R4 =  0

VSPP_quota_wind_total =  0
VSPP_quota_wind_R0 =  0
VSPP_quota_wind_R1 =  0
VSPP_quota_wind_R2 =  0
VSPP_quota_wind_R3 =  0
VSPP_quota_wind_R4 = 0

SPP_quota_solar_total =  18810
SPP_quota_solar_R0 =  0
SPP_quota_solar_R1 =  3150
SPP_quota_solar_R2 =  3060
SPP_quota_solar_R3 = 8100
SPP_quota_solar_R4 =  4500

VSPP_quota_solar_total = 5372
VSPP_quota_solar_R0 = 16
VSPP_quota_solar_R1 =  983
VSPP_quota_solar_R2 =  1929
VSPP_quota_solar_R3 =   725
VSPP_quota_solar_R4 =  1719
######## PDP ############################################################

print("max AVA SPP Wind = ",xr_ref['AVA_Wind'].rolling(lon = rollingwindow_spp_wind, lat = rollingwindow_spp_wind, min_periods=1,center=True).sum().max())
print("max AVA SPP Solar = ",xr_ref['AVA_Solar'].rolling(lon = rollingwindow_spp_solar, lat = rollingwindow_spp_solar, min_periods=1,center=True).sum().max())

######################## model #####################################################
m = linopy.Model()

built_SPP_wind = m.add_variables(integer=True,upper=1,lower=0, coords=xr_ref.coords, name='built_SPP_wind')
cap_SPP_wind = m.add_variables(upper=90,lower=0.00, coords=xr_ref.coords, name='cap_SPP_wind')

built_VSPP_wind = m.add_variables(integer=True,upper=1,lower=0, coords=xr_ref.coords, name='built_VSPP_wind')
cap_VSPP_wind = m.add_variables(upper=10,lower=0.00, coords=xr_ref.coords, name='cap_VSPP_wind')

built_SPP_solar = m.add_variables(integer=True,upper=1,lower=0, coords=xr_ref.coords, name='built_SPP_solar')
cap_SPP_solar = m.add_variables(upper=90,lower=0.00, coords=xr_ref.coords, name='cap_SPP_solar')

built_VSPP_solar = m.add_variables(integer=True,upper=1,lower=0, coords=xr_ref.coords, name='built_VSPP_solar')
cap_VSPP_solar = m.add_variables(upper=10,lower=0.00, coords=xr_ref.coords, name='cap_VSPP_solar')

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
)
    <= 1
    , name='constr_built_logic'
)

# constr_built_location_SPP_wind = m.add_constraints(
#     built_SPP_wind * xr_ref['AVA_Wind'] >= (areapergrid * built_SPP_wind)  , name = 'constr_built_location_SPP_wind'
# )

# constr_built_location_VSPP_wind = m.add_constraints(
#     built_VSPP_wind * xr_ref['AVA_Wind'] >= (areapergrid * built_VSPP_wind)  , name = 'constr_built_location_VSPP_wind'
# )

# constr_built_location_SPP_solar = m.add_constraints(
#     built_SPP_solar * xr_ref['AVA_Solar'] >= (areapergrid * built_SPP_solar)  , name = 'constr_built_location_SPP_solar'
# )

# constr_built_location_VSPP_solar = m.add_constraints(
#     built_VSPP_solar * xr_ref['AVA_Solar'] >= (areapergrid * built_VSPP_solar) , name = 'constr_built_location_VSPP_solar'
# )
############################################ Constraint Building Location Logic ##############################################################################

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

############################################################################################################################################################

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

############################################################################################################################################################

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

############################################## Objective Function #########################################################################################################################
obj = (-1000) * (
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
)
m.add_objective(obj)
############################################## Objective Function #########################################################################################################################

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

print("cap_SPP_wind = ",xr_ref['cap_SPP_wind'].sum())
print("cap_VSPP_wind = ",xr_ref['cap_VSPP_wind'].sum())
print("cap_SPP_solar = ",xr_ref['cap_SPP_solar'].sum())
print("cap_VSPP_solar = ",xr_ref['cap_VSPP_solar'].sum())

xr_ref.to_netcdf(path='Output\\xr_output.nc')