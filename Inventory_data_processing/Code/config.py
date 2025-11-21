## Minimum DBH threshold to include in analysis
glbl_dbh_min = 15

## Species to PFT mapping
# Keys are the species we are focusing on in this analysis. Values are their PFT shorthand. 
pft_dict ={'ABCO':'fir',
           'ABMA':'fir',
           'ACMA':'other_hardwood',
           'ALRH':'other_hardwood',
           'ARME':'other_hardwood',
           'CADE':'cedar',
           'CHCH':'other_hardwood',
           'CONU':'other_hardwood',
           'NODE':'other_hardwood',
           'OTHER':'unknown',
           'PIAL':'white_pine',
           'PIBA':'white_pine',
           'PICO':'other_conifer',
           'PIJE':'yellow_pine',
           'PILA':'white_pine',
           'PIMO':'white_pine',
           'PIPO':'yellow_pine',
           'PIPO/PIJE':'yellow_pine',
           'PISP':'other_conifer',
           'PREM':'other_hardwood',
           'PSME':'other_conifer',
           'QUCH':'evergreen_oak',
           'QUKE':'deciduous_oak',
           'QUWI':'evergreen_oak',
           'QUSP': 'other_hardwood', 
           'SASC':'other_hardwood',
           'SASP':'other_hardwood',
           'SEGI':'sequoia','UNKN':'unknown'}

## Unit ID to treatment mapping
treat_dict = {
    # Teakettle
    'uc1':'Thin','uc2':'Thin','uc3':'Thin', 
    'un1':'None', 'un2':'None', 'un3':'None', 
    'us1':'Thin', 'us2':'Thin', 'us3':'Thin', 
    'bc1':'Burn+Thin','bc2':'Burn+Thin', 'bc3':'Burn+Thin', 
    'bn1':'Burn', 'bn2':'Burn', 'bn3':'Burn', 
    'bs1':'Burn+Thin', 'bs2':'Burn+Thin', 'bs3':'Burn+Thin',
    # Blodgett
    'control':'None','burn':'Burn','mechburn':'Burn+Thin','mech':'Thin', 
    # W. Lake Tahoe
    'Control':'None', 'Treated':'Thin', 
    # Tharp's Creek
    'LOTHAR':'Burn', 'UPTHAR':'Burn', 'LOLOG':'None', 'UPLOG':'None', 
    # Sequoia
    'FFS2BURN':'Burn', 'FFS5BURN':'Burn', 'FFS6BURN':'Burn', 'FFS7CONTROL':'None',
    # STEF
    1: 'Burn+Thin',
    2: 'Burn+Thin',
    3: 'Burn',
    4: 'Thin',
    5: 'Thin',
    6: 'None',
    7: 'None',
    8: 'Thin',
    9: 'Thin',
    10: 'Thin',
    11: 'Thin',
    12: 'None',
    13: 'Burn+Thin',
    14: 'Burn',
    15: 'Burn+Thin',
    16: 'Burn',
    17: 'Burn+Thin',
    18: 'Burn+Thin',
    19: 'Burn',
    20: 'None',
    21: 'Thin',
    22: 'Burn+Thin',
    23: 'Burn+Thin',
    24: 'Thin'
}