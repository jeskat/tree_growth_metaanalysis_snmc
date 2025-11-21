import pandas as pd
import numpy as np

# Functions
def melt_inventory_df(df,vars_to_melt,id_vars,var_name):
  return df.loc[:,id_vars+vars_to_melt].melt(id_vars=id_vars,var_name=var_name)

def tree_thru_time(tree_id, df, treeid_col, date_col, dbh_col, status_col,dead_ind,date_incr):
    '''This function takes in a dataframe `df` with repeat forest inventory data.
    It returns a dataframe for a single tree ('tree_id') with diameter measurements arranged 
    chronologically, and growth increment calculated for each measurement period.'''

    # Subset tree and sort records chronologically
    one_tree = df.loc[df[treeid_col]==tree_id].sort_values(by=date_col)
    
    time_incr = []
    dbh_incr=[]
    start_yr = []

    # Loop through measurements. For each measurement after the first one, calculate growth as
    # the difference between DBH measurements, divided by the length of time between measurements.
    for i in range(len(one_tree.index)):
        # The first measurement gets null values for growth.
        if i == 0:
            time_incr.append(np.nan)
            dbh_incr.append(np.nan)
            start_yr.append(np.nan)
        else:
            ## Calculate the number of years between measurements
            if date_incr == 'day': # Calculate fraction of year based on number of days
                time_incr.append((one_tree[date_col].iloc[i]-one_tree[date_col].iloc[i-1]).days/365)
                start_yr.append(one_tree[date_col].iloc[i-1].year)
            elif date_incr=='year':
                time_incr.append((one_tree[date_col].iloc[i]-one_tree[date_col].iloc[i-1]))
                start_yr.append(one_tree[date_col].iloc[i-1])
            # If the record is a repeat measurement on a dead tree, the DBH increment recevies a null value.
            if (one_tree[status_col].iloc[i]==dead_ind) & (one_tree[status_col].iloc[i-1] == dead_ind):
                dbh_incr.append(np.nan)
            else: # Otherwise, calculate the change in DBH between measurements
                dbh_incr.append(one_tree[dbh_col].iloc[i]-one_tree[dbh_col].iloc[i-1]) 

    # Add new columns to tree timeseries
    one_tree['time_increment_yr']=time_incr
    one_tree['dbh_delta_cm'] = dbh_incr
    one_tree['dbh_increment_cm/yr'] = np.divide(one_tree['dbh_delta_cm'],one_tree['time_increment_yr'])
    one_tree['start_year'] = start_yr
    if date_incr == 'day':
        one_tree['end_year'] = one_tree[date_col].dt.year
    elif date_incr == 'year':
        one_tree['end_year'] = one_tree[date_col]
    
    return one_tree


def mismatched_species(tree_df, species_col):
    '''This function takes in a dataframe with timeseries inventory records for a 
    single tree (`tree_df`) and returns a value of 0 if all of the species records 
    are identical. Otherwise, it returns a value of 1.'''
    flag = 0
    if all(elem==tree_df[species_col].iloc[0] for elem in tree_df[species_col]):
        pass
    else:
        flag = 1
    return flag


def resurrected_trees(tree_df, status_col, live_ind, dead_ind):
    '''This function takes in a dataframe with timeseries inventory records for a 
    single tree (`tree_df`) and returns a value of 1 if the status of the tree is 
    indicated to be alive at a time following a record of the same tree being dead.
    Otherwise, it returns a value of 0.'''
    flag = 0
    if len(tree_df) > 1:
        for i in range(1,len(tree_df)):
            if (tree_df[status_col].iloc[i] == live_ind) & (tree_df[status_col].iloc[i-1]==dead_ind):
                flag = 1
                break
            else:
                continue
    else:
        pass
        
    return flag

def id_outliers(df, species, species_col, calculation_col='dbh_increment_cm/yr'):
    '''For a given species `species`, subsets forest inventory data and calculates upper and 
    lower bounds for annual growth increment as +/- 3 standard deviations from the species mean.
    Returns a dataframe with records associated with growth outliers, as well as the bounds.'''
    sub_df = df.loc[df[species_col]==species] # subset df to include only the chosen species
    # calculate mean and stdev annual dbh increment for species
    mean_ddbh = sub_df[calculation_col].mean()
    std_ddbh = sub_df[calculation_col].std()
    # define outliers: annual dbh increment as measurements  outside of mean +/- 3*sd
    ulim = mean_ddbh + 3*std_ddbh
    llim = mean_ddbh - 3*std_ddbh
    # Return a dataframe containing the outlier records, along with the upper and lower bound for growth
    return sub_df.loc[(sub_df[calculation_col]>ulim) | (sub_df[calculation_col]<llim)], llim, ulim
    

def manage_outliers(master_df, outlier_df, outlier_lookup, 
                    species_col, treeid_col, dbh_col, 
                    year_col, date_col, status_col, dead_ind, 
                    date_incr, calculation_col='dbh_increment_cm/yr'):
    '''
    Identifies DBH records associated with growth outliers using the following process:
    
    1) If a tree has one or more outliers in its growth timeseries, flag the record associated with the first outlier in the tree’s growth timeseries.
    2) Recalculate the mean annual DBH increment for the remaining records
    3) Check the adjusted timeseries for outliers
    4) Repeat (1-3) until all the outliers have been flagged
    5) Replace original growth timeseries with adjusted one
    6) Log which records were flagged
    
    If the first growth increment is the outlier, decide whether first or second DBH record needs to be dropped.
    If mean growth in first AND third interval exceed that of the mean growth between interval 1 and 3, the second record is the outlier. 
    Otherwise, it’s the first.
    '''
    ddbh_noOutliers = master_df.copy() # this is the dataframe that will be modified and ultimately returned
    ddbh_noOutliers['Outlier'] = 0

    # Loop through tree with one or more outlier ddbh increments recorded
    for i in outlier_df[treeid_col].unique(): 
        tree = ddbh_noOutliers.loc[ddbh_noOutliers[treeid_col]==i].copy() # locate all (ordered) records for the tree with outlier(s) 
        tree_outliers = outlier_df.loc[outlier_df[treeid_col]==i].copy() # find the outliers for that tree
        
        # Identify the upper and lower bounds for record inclusion based on the species of the tree
        species = tree.iloc[0][species_col]
        llim = outlier_lookup.loc[outlier_lookup['Species']==species,'Lower_limit'].values[0]
        ulim = outlier_lookup.loc[outlier_lookup['Species']==species, 'Upper_limit'].values[0]
        
        while len(tree_outliers)>0: # repeat for every outlier in a tree's timeseries
            if tree_outliers.index[0] == tree.index[1]: # Is the first growth interval the outlier?
                if len(tree) > 2: # If we drop an outlier, will we still have at least one growth interval?
                # If mean growth between the first and second AND the second and third growth intervals BOTH exceed mean growth between the first and third intervals, the second record is the outlier. 
                # Otherwise, the first record is the outlier.
                    # Calculate mean growth between the first and third measurement year
                    ddbh_0to2 = abs((tree.iloc[2][dbh_col]-tree.iloc[0][dbh_col])/(tree.iloc[2][year_col]-tree.iloc[0][year_col]))

                    # If mean growth between the first and second AND the second and third growth intervals BOTH exceed mean growth 
                    # between the first and third intervals, then the second record is the outlier. 
                    if (abs(tree.iloc[1][calculation_col])>ddbh_0to2) & (abs(tree.iloc[2][calculation_col])>ddbh_0to2):
                        tree.drop(tree_outliers.index[0], inplace = True) # drop the outlier (second record)
                        # ddbh_noOutliers.loc[tree_outliers.index[0],calculation_col] = np.nan
                        ddbh_noOutliers.loc[tree_outliers.index[0],'Outlier'] = 1
                        tree = tree_thru_time(tree_id=i,
                                              df=tree,
                                              date_col=date_col,
                                              dbh_col=dbh_col,
                                              status_col=status_col,
                                              dead_ind=dead_ind,
                                              treeid_col=treeid_col,
                                              date_incr=date_incr) # recalculate dbh increment with remaining trees (don't use tree_thru_time_seq!)
                        tree_outliers = tree.loc[(tree[calculation_col]>ulim) | (tree[calculation_col]<llim)]  # get an updated list of outliers for tree
                    else:
                        # ddbh_noOutliers.loc[tree.index[0],calculation_col] = np.nan
                        ddbh_noOutliers.loc[tree.index[0],'Outlier'] = 1
                        tree.drop(tree.index[0], inplace = True) # drop the outlier (first record)
                        tree = tree_thru_time(tree_id=i,
                                              df=tree,
                                              date_col=date_col,
                                              dbh_col=dbh_col,
                                              status_col=status_col,
                                              dead_ind=dead_ind,
                                              treeid_col=treeid_col,
                                              date_incr=date_incr) # recalculate dbh increment with remaining records
                        tree_outliers = tree.loc[(tree[calculation_col]>ulim) | (tree[calculation_col]<llim)] # get an updated list of outliers for tree
                else: # if there are only two growth records
                  # check whether one of the growth records is zero. If so, that record is flagged as the outlier.
                    if 0 in tree[dbh_col].tolist():
                        out = tree.loc[tree[dbh_col]==0]
                        # ddbh_noOutliers.loc[out.index,calculation_col] = np.nan
                        ddbh_noOutliers.loc[out.index,'Outlier'] = 1
                        tree.drop(out.index, inplace=True)
                        tree = tree_thru_time(tree_id=i,
                                              df=tree,
                                              date_col=date_col,
                                              dbh_col=dbh_col,
                                              status_col=status_col,
                                              dead_ind=dead_ind,
                                              treeid_col=treeid_col,
                                              date_incr=date_incr) # recalculate dbh increment with remaining trees
                        tree_outliers = tree.loc[(tree[calculation_col]>ulim) | (tree[calculation_col]<llim)]  # get an updated list of outliers for tree
                    # otherwise, flag both records as outliers. 
                    else:
                        tree['Outlier'] = 1 
                        tree_outliers = pd.DataFrame()   # have to do this to keep the loop moving     
            else: # for trees where the first growth interval is not the outlier
                tree.drop(tree_outliers.index[0], inplace = True) # drop the first record that is an outlier
                # ddbh_noOutliers.loc[tree_outliers.index[0],calculation_col] = np.nan
                ddbh_noOutliers.loc[tree_outliers.index[0],'Outlier'] = 1
                tree = tree_thru_time(tree_id=i,
                                      df=tree,
                                      date_col=date_col,
                                      dbh_col=dbh_col,
                                      status_col=status_col,
                                      dead_ind=dead_ind,
                                      treeid_col=treeid_col,
                                      date_incr=date_incr) # recalculate dbh increment with remaining trees
                tree_outliers = tree.loc[(tree[calculation_col]>ulim) | (tree[calculation_col]<llim)]  # get an updated list of outliers for tree
        ddbh_noOutliers.loc[tree.index, :]= tree # add the updated timeseries for the tree back to the master dataframe
    return ddbh_noOutliers 


def single_records(df, treeid_col):
    '''Identify the trees for which there is only a single DBH record (i.e., no way to calculate ddbh)'''
    treeID_single_measurement = []
    # Loop through all unique tree IDs
    for t in df[treeid_col].unique():
        if len(df.loc[df[treeid_col]==t])<2:
            treeID_single_measurement.append(t)
    return treeID_single_measurement

def make_statespace_df(no_outlier_df, 
                       year_col, 
                       unit_col, 
                       treeid_col, 
                       species_col, 
                       status_col, 
                       dbh_col, 
                       pft='all'): # Default is to include all PFTs
    '''Reformats inventory dataframe `no_outlier_df` into two new dataframes.
    The first (tree_level_df) is the dataframe of tree attributes that is needed as input to the statespace models.
    The second (pft_df) is a master dataframe for this site, which is useful for additional analysis and visualization.
    '''
    
    if pft!='all':
        pft_df = no_outlier_df.loc[no_outlier_df['PFT']==pft,:].copy() # subset the growth df for the pft of choice
    else:
        pft_df = no_outlier_df.copy()

    # Determine the first and last year of measurements across the whole experiment
    first_yr = pft_df[year_col].min()
    last_yr = pft_df[year_col].max()

    # Assign the year number within the experiment timeseries to each measurement
    pft_df['timeseries_year'] = pft_df[year_col] - first_yr + 1 # Add the one so that the first year is Year 1, not Year 0

    # Assign a count number to the observation within each tree's timeseries
    pft_df['counter'] = np.nan

    tree_id = []
    st_yr = []
    end_yr = []
    plot = []
    n_obs = []
    pft_lst = []

    # Loop through tree IDs
    for tree in pft_df.sort_values(treeid_col)[treeid_col].unique():
        tree_id.append(tree) # get the TreeIDs in the correct order
        tree_df = pft_df.loc[pft_df[treeid_col]==tree] # df for the tree
        st_yr.append(tree_df['timeseries_year'].min()) # first year of measurements for this tree
        end_yr.append(tree_df['timeseries_year'].max()) # final year of measurements for this tree
        pft_lst.append(tree_df['PFT'].iloc[0]) # tree pft
        plot.append(tree_df[unit_col].iloc[0]) # treatment unit associated with tree
        n_obs.append(len(tree_df)) # number of times tree was measured

        # add a count associated with each tree measurement (1:n indexing instead of Python indexing)
        pft_df.loc[pft_df[treeid_col]==tree,'counter'] = range(1,len(tree_df)+1) 

    # Construct the tree-level attributes table
    d = {'TreeID':tree_id,'start_year':st_yr, 'end_year':end_yr, 'PlotID':plot, 'n_obs':n_obs, 'PFT':pft_lst}
    tree_level_df = pd.DataFrame(d)

    # Standardize some column names across sites
    pft_df.rename(columns={unit_col: "UnitID", treeid_col: "TreeID", 
                         species_col: "Species", year_col: "Year", 
                         status_col: "Status", dbh_col: "DBH"}, inplace=True)
    return tree_level_df, pft_df
