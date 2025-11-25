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

def cvar_mean_by_plot(cvar, df):
    '''Loops through sites; for each site, finds the mean value of the climate variable for each 
    measurement unit (plot) over the study timeframe'''
    site_list = []
    plot_list = []
    val_list = []
    unit_nms = []
    
    for site in df['Site'].unique():
        ## Open tree-level information
        pft_df = pd.read_csv('Outputs/Processed_data_for_statespace/{}/{}_pft_df.csv'.format(site,site))
        
        # Subset by initial and final year of the study
        yri = pft_df['start_year'].min() + 1 # Include the next climate year but not the previous
        yrf = pft_df['end_year'].max()
        data_subset = df.loc[(df['Site']==site) & (df['year']>=yri) & (df['year']<=yrf)]

        ## Loop through units
        for u in data_subset['UnitID'].unique():
            ## Loop through subplots
          for p in data_subset.loc[data_subset['UnitID']==u, 'PlotID'].unique():
              site_list = site_list + [site]
              unit_nms = unit_nms + [u]
              plot_list = plot_list + data_subset.loc[(data_subset['UnitID']==u) & (data_subset['PlotID']==p), 'unique_nm'].unique().tolist()
              # Mean value of climate variable
              val_list = val_list + [data_subset.loc[(data_subset['UnitID']==u) & (data_subset['PlotID']==p), cvar].mean()]
        
    mean_plot_cvar = pd.DataFrame({'Site':site_list, 'UnitID':unit_nms, 'unique_nm':plot_list, cvar:val_list})
    return mean_plot_cvar