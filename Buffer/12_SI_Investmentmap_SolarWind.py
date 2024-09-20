import linopy
import pandas as pd
import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
from pathlib import Path
from shapely import wkt
import numpy as np

###################### Config ################################################################################################
coarsenscale = 5
lccs_resolution = 300 * coarsenscale #m
areapergrid = (lccs_resolution/1000) ** 2 ## km2
scenario_SI = 0 ## Include area where SI >= scenario_SI

mwperkm2_wind = 4.5 ## originally 9 MW/km2 but deduct by 50% of the technically available from IEA's Thailand CET
mwperkm2_solar = 15 ## originally 30 MW/km2 but deduct by 50% of the technically available from IEA's Thailand CET

mwpergrid_wind = np.round(areapergrid * mwperkm2_wind,2)
mwpergrid_solar =  np.round(areapergrid * mwperkm2_solar,2)

print('areapergrid = ',areapergrid)
print('mwpergrid_wind = ',mwpergrid_wind)
print('mwpergrid_solar = ',mwpergrid_solar)

xr_final_SI_raw = xr.open_dataset('Output\\xr_final_SI.nc')
xr_final_SI_raw = xr_final_SI_raw.drop_vars('ADM1_EN')
xr_final_SI_raw = xr_final_SI_raw.drop_vars('A_BGEC')
xr_final_SI_raw = xr_final_SI_raw.drop_vars('A_Biomass')
xr_final_SI_raw = xr_final_SI_raw.drop_vars('A_BGWW')
xr_final_SI_raw = xr_final_SI_raw.drop_vars('A_MSW')
xr_final_SI_raw = xr_final_SI_raw.drop_vars('A_IEW')
xr_final_SI_raw['SI_Wind'] = xr.where(xr_final_SI_raw['SI_Wind'] >= scenario_SI ,xr_final_SI_raw['SI_Wind'],0)
xr_final_SI_raw['SI_Solar'] = xr.where(xr_final_SI_raw['SI_Solar'] >= scenario_SI ,xr_final_SI_raw['SI_Solar'],0)
xr_final_SI_raw['AVA_Wind'] = xr.where(xr_final_SI_raw['SI_Wind'] == 0 , 0 , xr_final_SI_raw['AVA_Wind'])
xr_final_SI_raw['AVA_Solar'] = xr.where(xr_final_SI_raw['SI_Solar'] == 0 , 0 , xr_final_SI_raw['AVA_Solar'])

xr_final_SI = xr_final_SI_raw.coarsen(lat = coarsenscale, lon= coarsenscale, boundary='trim').sum()
xr_final_SI['SI_Wind'] = xr_final_SI['SI_Wind'] / (coarsenscale**2)
xr_final_SI['SI_Solar'] = xr_final_SI['SI_Solar'] / (coarsenscale**2)

######### Find min/max SI within cell #############
# xr_final_SI['SI_Wind_max'] = xr_final_SI_raw['SI_Wind'].coarsen(lat = coarsenscale, lon= coarsenscale, boundary='trim').max()
# xr_final_SI['SI_Wind_min'] = xr_final_SI_raw['SI_Wind'].coarsen(lat = coarsenscale, lon= coarsenscale, boundary='trim').min()
# xr_final_SI['SI_Solar_max'] = xr_final_SI_raw['SI_Solar'].coarsen(lat = coarsenscale, lon= coarsenscale, boundary='trim').max()
# xr_final_SI['SI_Solar_min'] = xr_final_SI_raw['SI_Solar'].coarsen(lat = coarsenscale, lon= coarsenscale, boundary='trim').min()

xr_final_SI_raw.close()
print(xr_final_SI.data_vars)
print("AVA Wind = ",xr_final_SI['AVA_Wind'].sum())
print("AVA Solar = ",xr_final_SI['AVA_Solar'].sum())
print('Max SI Wind = ',xr_final_SI['SI_Wind'].max())
print('Max SI Solar = ',xr_final_SI['SI_Solar'].max())

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
print('Mean SI Wind = ',xr_ref['SI_Wind'].where(xr_ref['SI_Wind']>0).mean())
print('Mean SI Solar = ',xr_ref['SI_Solar'].where(xr_ref['SI_Solar']>0).mean())
print('Max SI Wind = ',xr_ref['SI_Wind'].max())
print('Max SI Solar = ',xr_ref['SI_Solar'].max())
print('coarsenscale = ',coarsenscale)
print('areapergrid = ',areapergrid)
print('mwpergrid_wind = ',mwpergrid_wind)
print('mwpergrid_solar = ',mwpergrid_solar)
###################### Summary ################################################################################################

# ####### no ############################################################
# quota_wind_total =  0
# quota_wind_R0 = 0
# quota_wind_R1 = 0
# quota_wind_R2 = 0
# quota_wind_R3 =  0
# quota_wind_R4 =  0

# quota_solar_total =  0
# quota_solar_R0 =  0
# quota_solar_R1 =  0
# quota_solar_R2 =  0
# quota_solar_R3 = 0
# quota_solar_R4 =  0
# ####### no ############################################################

######## yes ############################################################
quota_wind_R0 = 0
quota_wind_R1 = 0
quota_wind_R2 = 6100 #6084.7
quota_wind_R3 =  300 #260.1
quota_wind_R4 =  0
quota_wind_total = quota_wind_R0 +  quota_wind_R1 + quota_wind_R2 + quota_wind_R3 + quota_wind_R4 #6344.8

quota_solar_R0 =  50 #23.34
quota_solar_R1 =  3300 #3255.36
quota_solar_R2 =  6100 #6072.28
quota_solar_R3 = 5650 #5602.08
quota_solar_R4 =  3200 #3175.94
quota_solar_total =  quota_solar_R0 + quota_solar_R1 + quota_solar_R2 + quota_solar_R3 + quota_solar_R4 #18129

######## yes ############################################################

######################## model #####################################################
m = linopy.Model()

built_wind = m.add_variables(binary=True, coords=xr_ref.coords, name='built_wind')
cap_wind = m.add_variables(lower=0.00, coords=xr_ref.coords, name='cap_wind')

built_solar = m.add_variables(binary=True, coords=xr_ref.coords, name='built_solar')
cap_solar = m.add_variables(lower=0.00, coords=xr_ref.coords, name='cap_solar')


############################################ Constraint Building Location Logic ##############################################################################
constr_built_logic =  m.add_constraints(
    (
        built_wind
        +
        built_solar
    )
        <= 1
        , name='constr_built_logic'
)

############################################ Constraint Building Location Logic ##############################################################################

############################################ Constraint Capacity ##############################################################################
constr_maxcap_wind = m.add_constraints(
    cap_wind <= (built_wind) * xr_ref['AVA_Wind'] * mwperkm2_wind
    ,name = 'constr_maxcap_wind'
)
constr_mincap_wind = m.add_constraints(
    cap_wind >= (built_wind) * 1
    ,name = 'constr_mincap_wind'
)

constr_builtarea_wind = m.add_constraints(
    built_wind <= xr_ref['AVA_Wind'] * 10000
    ,name = 'constr_builtarea_wind'
)

constr_maxcap_solar = m.add_constraints(
    cap_solar <= (built_solar) * xr_ref['AVA_Solar'] * mwperkm2_solar
    ,name = 'constr_maxcap_solar'
)
constr_mincap_solar = m.add_constraints(
    cap_solar >= (built_solar) * 1
    ,name = 'constr_mincap_solar'
)
constr_builtarea_solar = m.add_constraints(
    built_solar <= xr_ref['AVA_Solar'] * 10000
    ,name = 'constr_builtarea_solar'
)
############################################ Constraint Capacity ##############################################################################

###########################################################################################################################################################

constr_quota_wind = m.add_constraints((cap_wind.sum()) >= quota_wind_total , name='constr_quota_wind')

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

constr_quota_solar = m.add_constraints((cap_solar.sum()) >= quota_solar_total, name='constr_quota_solar')

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

###########################################################################################################################################################'

obj = (-1000) * (
    ( xr_ref['SI_Wind'] * (cap_wind / mwpergrid_wind) )
    +
    ( xr_ref['SI_Solar'] * (cap_solar / mwpergrid_solar) )
 )

m.add_objective(obj)

print("presolve = ",m)
m.solve(solver_name='highs',
        # mip_abs_gap = 0.002,
        # mip_rel_gap = 0.002,
        )

print('aftersolve = ',m)
solution = m.solution
solution = solution.fillna(0)
print(solution)

xr_ref['cap_wind'] = solution['cap_wind'] * solution['built_wind']
xr_ref['cap_solar'] = solution['cap_solar'] * solution['built_solar']

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


print(xr_ref.data_vars)
xr_ref.to_netcdf(path='Output\\xr_output_SSI_' + str(scenario_SI) + "_CS_"+str(coarsenscale)+ '_.nc')














