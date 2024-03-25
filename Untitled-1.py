import linopy
import pandas as pd
import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from shapely import wkt



xr_ref = xr.open_dataset('Output\\xr_final_SI_exc_r.nc')

mwperkm2_wind = 4.5 ## originally 9 MW/km2 but deduct by 50% of the technically available from IEA's Thailand CET
mwperkm2_solar = 15 ## originally 30 MW/km2 but deduct by 50% of the technically available from IEA's Thailand CET

lccs_resoultion = (300/1000) ** 2 ## meters
print(lccs_resoultion)
## max capacity of vspp 10 MW
maxcapacityinpolygon = 10

maxnoofgrid_wind = np.ceil((maxcapacityinpolygon/(mwperkm2_wind * lccs_resoultion)))
maxnoofgrid_solar = np.ceil((maxcapacityinpolygon/(mwperkm2_solar * lccs_resoultion)))

rolling_latlon_wind = np.ceil(np.sqrt(maxnoofgrid_wind)).astype('int16')
rolling_latlon_solar = np.ceil(np.sqrt(maxnoofgrid_solar)).astype('int16')

print("maxnoofgrid_wind = ",maxnoofgrid_wind,"sqrt = ",np.ceil(np.sqrt(maxnoofgrid_wind)))
print("maxnoofgrid_solar = ",maxnoofgrid_solar,"sqrt = ",np.ceil(np.sqrt(maxnoofgrid_solar)))

######################## model #####################################################
m = linopy.Model()

built_solar = m.add_variables(integer=True,upper=1,lower=0, coords=xr_ref.coords, name='built_solar')

constr_built_logic = m.add_constraints(built_solar <=  1, name='built_logic')
constr_rolling_lat = m.add_constraints(built_solar.rolling(lat=rolling_latlon_solar,min_periods=1,center=True).sum() <=1 ,name='constr_rolling_lat')
constr_rolling_lon = m.add_constraints(built_solar.rolling(lon=rolling_latlon_solar,min_periods=1,center=True).sum() <=1 ,name='constr_rolling_lon')

obj = (-1000) * (
    (built_solar*xr_ref['SI_Solar']).sum()
)

m.solve(solver_name='highs',
        mip_abs_gap = 0.0018,
        mip_rel_gap = 0.0018,
        )

print('aftersolve = ',m)
solution = m.solution
solution = solution.fillna(0)

print(solution)