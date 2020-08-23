import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle


def fracWithHosts(transient_dict):
    count = 0
    for name, host in transient_dict.items():
        # only do matching if there's a found host
        if isinstance(host, list):
            if len(host) > 0 and np.any(np.array(host == host)):
                count += 1
        elif isinstance(host, np.ndarray):
            if len(host) > 0 and np.any(np.array(host == host)):
                count += 1
        else:
            if host == host:
                count+=1
            else:
                print(host)
                #count += 1
    return count/len(transient_dict.keys())


def build_ML_df(dic, hostDF, transientDF):
    hostDF = hostDF.reset_index(drop=True)
    hostDF["TransientClass"] = ""
    hostDF["Transient Name"] = ""
    for name, host in dic.items():
        # only do matching if there's a found host
        chosenHost = ""
        if (host == host):
            if isinstance(host, np.ndarray):
                if host:
                    chosenHost = host[0]
            else:
                chosenHost = host
        if chosenHost:
            #find host in df
            idx = hostDF['objID'] == chosenHost
            idx_transient = transientDF['Name'] == str(name)
            if hostDF.loc[idx, "TransientClass"].values[0] != "":
                print("Found a double!")
                hostDF = hostDF.append([hostDF[idx]], ignore_index=True)
                idx = hostDF.index[-1]
            hostDF.loc[idx, "TransientClass"] = transientDF.loc[idx_transient, 'Obj. Type'].to_string(index=False).strip()
            hostDF.loc[idx, "Transient Name"] = transientDF.loc[idx_transient, 'Name'].to_string(index=False).strip()
            transCoord = SkyCoord(transientDF.loc[idx_transient, 'RA'], transientDF.loc[idx_transient, 'DEC'], unit=(u.hourangle, u.deg))
            if len(transCoord) > 1:
                transCoord = transCoord[0]
            hostDF.loc[idx, "TransientRA"] = transCoord.ra.deg
            hostDF.loc[idx, "TransientDEC"] = transCoord.dec.deg
            hostDF.loc[idx, "TransientDiscoveryDate"] = transientDF.loc[idx_transient, 'Discovery Date (UT)'].to_string(index=False).strip()
            hostDF.loc[idx, "TransientDiscoveryMag"] = transientDF.loc[idx_transient, 'Discovery Mag'].to_string(index=False).strip()
            hostDF.loc[idx, "TransientRedshift"] = transientDF.loc[idx_transient, 'Redshift'].to_string(index=False).strip()
    #hostDF = hostDF[hostDF["TransientClass"] != ""]
    #hostDF = hostDF.drop_duplicates(subset=hostDF.columns.difference(['distance']))
    hostDF = hostDF.reset_index(drop=True)
    return hostDF
