path = '/home/alexgagliano/Documents/Research/Transient_ML/scripts'
import os
os.chdir(path)
from PS1QueryFunctions import *
from TNSQueryFunctions import *
from NEDQueryFunctions import *
from starSeparation import *
from sourceCleaning import *
from hostMatching import *
from stellarLocus import *
from DLR import *
from datetime import datetime
from astropy import units as u
from astropy.coordinates import SkyCoord
import pandas as pd
import numpy as np
import pickle
from collections import Counter
path = '/home/alexgagliano/Documents/Research/Transient_ML/'
os.chdir(path)
fullData = pd.read_csv("/home/alexgagliano/Documents/Research/Transient_ML_Box/tables/FinalAssociationTable_0406_Cleaned_FullNED_Fixed.csv")

now = datetime.now()
dateStr = "%i%.02i%.02i" % (now.year,now.month,now.day)

dateStr = "20200402"
rad = 30 #arcsec
#fn_SN = "SNe_TNS_%s.csv" % dateStr
fn_SN  ='./tables/allSNe_2020.tar.gz'
fn_Host = "SNe_TNS_%s_PS1Hosts_%iarcsec.csv" % (dateStr, rad)
fn_Dict = fn_Host.replace(".csv", "") + ".p"

#begin doing the heavy lifting to associate transients with hosts
#get_transients(path, fn_SN)
#get_hosts(path, fn_SN, fn_Host, rad)
host_DF = pd.read_csv(path+'../Transient_ML_Box/tables/'+fn_Host)

cuts = ["n", "quality", "coords", "primary", "best", "duplicate"]
#find dict - we'll be using this later
transient_dict = pickle.load(open(path+"dictionaries/"+fn_Dict, "rb"))

transient_dict = {k.replace(' ', ''): v for k, v in transient_dict.items()}
# check how many supernovae we have potential hosts for - should be nearly all of em
fracWithHosts(transient_dict)

host_DF = makeCuts(host_DF, cuts, transient_dict)
transient_dict = clean_dict(transient_dict, host_DF, [])
host_DF = getColors(host_DF)
host_DF = removePS1Duplicates(host_DF)

host_DF_noNED = getNEDInfo(host_DF_noNED)
host_DF = pd.concat([host_DF_wNED, host_DF_noNED])
#host_DF.to_csv("checkpoint_NEDInfo.tar.gz")
host_DF = pd.read_csv(path+"../Transient_ML_Box/tables/checkpoint_NEDInfo_wColors_0318_FinalSearch.tar.gz")

#~25k found - we'll take it!
host_DF = getColors(host_DF)
host_DF = calc_7DCD(host_DF)
host_gals_DF, stars_DF = separateStars(host_DF, plot=1)

transient_dict = clean_dict(transient_dict, host_gals_DF, [])

host_gals_DF= host_gals_DF[host_gals_DF['NED_type'] != '!*']
host_gals_DF= host_gals_DF[host_gals_DF['NED_type'] != '!V*']
host_gals_DF= host_gals_DF[host_gals_DF['NED_type'] != 'V*']
host_gals_DF= host_gals_DF[host_gals_DF['NED_type'] != 'PN']

from collections import Counter
Counter(host_gals_DF['NED_type'])

host_gals_DF.to_csv("gals_checkpoint_StarSeparation_0402_OldLocus.tar.gz",index=False)
stars_DF.to_csv("stars_checkpoint_StarSeparation_0402_OldLocus.tar.gz",index=False)

plotLocus(host_gals_DF, color=1, save=1, type="Gals", timestamp=dateStr)
plotLocus(stars_DF, color=1, save=1, type="Stars", timestamp=dateStr)

host_dict_nospace = {k.replace(' ', ''): v for k, v in transient_dict.items()}
fn = "RedoneSupernovae_CataloguedByName.txt"
transients = pd.read_csv(path+fn_SN)
clean_dict(host_dict_nospace, host_gals_DF, [])
fracWithHosts(host_dict_nospace)

with open("./dictionaries/" + "gals_checkpoint_preDLR.p", 'wb') as fp:
    pickle.dump(host_dict_nospace, fp, protocol=pickle.HIGHEST_PROTOCOL)

host_DF, host_dict_nospace = chooseByDLR(host_gals_DF, transients, fn, host_dict_nospace, host_dict_nospace, todo="r")

host_DF = build_ML_df(host_dict_nospace, host_DF, transients)

with open("./dictionaries/" + "Final_Dictionary_0402_OldLocus.p", 'wb') as fp:
    pickle.dump(host_dict_nospace, fp, protocol=pickle.HIGHEST_PROTOCOL)

host_DF.drop_duplicates(subset=['Transient Name'], inplace=True)
host_DF = host_DF[host_DF['Transient Name'] != ""]
host_DF.reset_index(inplace=True, drop=True)

host_DF['Transient Name'] = [x.replace(" ", "") for x in host_DF['Transient Name']]

host_DF['TransientDiscoveryYear'] = [str(x).split("/")[0] for x in host_DF['TransientDiscoveryDate']]
host_DF['TransientDiscoveryYear'] = [str(x).split("-")[0] for x in host_DF['TransientDiscoveryYear']]
host_DF.to_csv("FinalAssociationTable_0402_OldLocus.csv", index=False)

hSpecPath = "/home/alexgagliano/Documents/Research/Transient_ML_Box/hostSpectra"
tSpecPath = "/home/alexgagliano/Documents/Research/Transient_ML/tables/SNspectra"
#getAllPostageStamps(host_DF, 120, psPath)
host_DF = pd.read_csv("/home/alexgagliano/Documents/Research/Transient_ML/tables/FinalAssociationTable_20200225_Cleaned_FullNED_ApMag.csv")
getNEDSpectra(host_DF, hSpecPath)
getTNSSpectra(host_DF, tSpecPath)
clean_spectra('/home/alexgagliano/Documents/Research/Transient_ML/tables/SNspectra')
