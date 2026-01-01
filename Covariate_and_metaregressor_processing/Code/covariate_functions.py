import pandas as pd
import numpy as np


def find_nearest_pt(latitude, longitude, ds):
  '''Identifies the closest grid cell in a raster to a given lat-lon coordinate'''
  d_lat = ds.lat - latitude
  d_lon = ds.lon - longitude
  r2 = d_lat**2 + d_lon**2
  i_j_loc = np.where(r2 == np.min(r2))
  nearest_point = ds.isel(lat=i_j_loc[0], lon=i_j_loc[1])
  return nearest_point

def find_plot_gridcells(ds, plots):
    '''Takes in a dataframe of plots and find the nearest raster grid cell to the plot centroid'''
    df_dict =  {}
    for i in plots.index:
        site = plots.loc[i, 'unique_nm']
        nrst_pt = find_nearest_pt(plots.loc[i, 'Latitude'], plots.loc[i, 'Longitude'], ds)
        df = nrst_pt.to_dataframe()
        df['Site'] = plots.loc[i, 'Site']
        df['UnitID'] = plots.loc[i, 'UnitID']
        df['PlotID'] = plots.loc[i, 'PlotID']
        df['unique_nm'] = plots.loc[i, 'unique_nm']
        df_dict[site] = df
    out = pd.concat(df_dict.values())
    out.reset_index(inplace=True, drop=False)
    return out