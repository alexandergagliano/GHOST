import numpy as np
from astropy.table import Table
import requests
from PIL import Image
from io import BytesIO
import pylab
import os
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
import re
from PS1QueryFunctions import *
import glob
def remove_prefix(text, prefix):
	return text[text.startswith(prefix) and len(prefix):]

# Runs through pipeline to find host for new event,
# and returns all data for transient host.
# input - TransientCoord, the coordinates of
#         a new transient being queried
# output - A data frame of the host
#          associated with the transient
#          being queried, with PS1 information.
#          If the supernova is currently in the
#          database, return the database object.
#          If not, complete the pipeline for
#          host association and return the database
#          object.
#          If no host is found, return None
def getHostFromTransientCoords(transientCoord):
	fullTable = fullData()
	c2 = SkyCoord(fullTable['TransientRA']*u.deg, fullTable['TransientDEC']*u.deg, frame='icrs')
	sep = np.array(transientCoord.separation(c2).arcsec)
	host = None
	if np.nanmin(sep) <= 1:
		host_idx = np.where(sep == np.nanmin(sep))[0][0]
		host = fullTable.iloc[[host_idx]]
	else:
		print("Sorry, that supernova was not found in our database! The closest supernova is %.2f arcsec away.\n" %np.nanmin(sep))
	return host

# Returns all data for transient host,
# based on supernova name
# input - TransientCoord, the coordinates of
#         a new transient being queried
# output - A data frame of the host
#          associated with the transient
#          being queried, with PS1 information.
#          If the supernova is currently in the
#          database, return the database object.
#          If no host is found, return None
def getHostFromTransientName(SNName):
	fullTable = fullData()
	SNName = re.sub(r"\s+", "", str(SNName))
	host = None
	possibleNames = [SNName, SNName.upper(), SNName.lower(), "SN"+SNName]
	for name in possibleNames:
		if len(fullTable[fullTable['Transient Name'] == name])>0:
			host = fullTable[fullTable['Transient Name'] == name]
			break
	else:
		print("Sorry, that supernova was not found in our database!\n")
	return host

# Returns all data for transient and host,
# based on host name
# input - hostName
# output - A data frame of the host
#          associated with the transient
#          being queried, with PS1 information.
#          If the supernova is currently in the
#          database, return the database object.
#          If no host is found, return None
def getHostFromHostName(hostName):
	fullTable = fullData()
	hostName = re.sub(r"\s+", "", str(hostName))
	allHosts = np.array([re.sub(r"\s+", "", str(x)) for x in fullTable['NED_name']])
	possibleNames = [hostName, hostName.upper(), hostName.lower()]
	for name in possibleNames:
		if name in allHosts:
			host = fullTable.iloc[[np.where(name == allHosts)[0][0]]]
			return host
	print("Sorry, that host was not found in our database!\n")
	return None

def getHostFromHostCoords(hostCoords):
	fullTable = fullData()
	c2 = SkyCoord(fullTable['raMean']*u.deg, fullTable['decMean']*u.deg, frame='icrs')
	sep = np.array(hostCoords.separation(c2).arcsec)
	host = None
	if np.nanmin(sep) <= 1:
		host_idx = np.where(sep == np.nanmin(sep))[0][0]
		host = fullTable.iloc[[host_idx]]
	else:
		print("Sorry, that host was not found in our database! The closest host is %.2f arcsec away.\n" %np.nanmin(sep))
	return host

# Returns basic statistics for transient,
# based on a query of the coordinates of
# its host
# inputs - the position of the host
# outputs - A printout of name,
#           discovery year and discovery mag
#           for all transients associated
#           with a host at that coordinates

def getTransientStatsFromHostCoords(hostCoord):
	host = getHostFromHostCoords(hostCoord)
	i = 0
	if len(host) > 0:
		print("Found info for host %s.\n"%host['NED_name'].values[0])
		print("Associated supernovae: ")
		for SN in np.array(host['Transient Name'].values):
			SN_frame = host.loc[host['Transient Name'] == SN]
			print("%i. %s"%(i+1, SN_frame['Transient Name'].values[0]))
			print("RA, DEC (J2000): %f, %f"%(SN_frame['TransientRA'].values[0], SN_frame['TransientDEC'].values[0]))
			print("Redshift: %f"%SN_frame['TransientRedshift'].values[0])
			print("Discovery Date: %s"%SN_frame['TransientDiscoveryDate'].values[0].split(" ")[0])
			print("Discovery Mag: %.2f"%SN_frame['TransientDiscoveryMag'].values[0])
			i+= 1
	return

# Returns basic statistics for transient,
# based on a query of the coordinates of
# its host
# inputs - the position of the host
# outputs - A printout of name,
#           discovery year and discovery mag
#           for all transients associated
#           with a host at that coordinates
def getTransientStatsFromHostName(hostName):
	host = getHostFromHostName(hostName)
	hostCoord = SkyCoord(np.unique(host['raMean'])*u.deg, np.unique(host['decMean'])*u.deg, frame='icrs')
	getTransientStatsFromHostCoords(hostCoord)
	return

# Returns basic statistics for the most likely
# host of a previously identified transient, using
# the pipeline outlined above.
# inputs - the coordinates of the transient to search
def getHostStatsFromTransientCoords(transientCoords):
	fullTable = fullData()
	c2 = SkyCoord(fullTable['TransientRA']*u.deg, fullTable['TransientDEC']*u.deg, frame='icrs')
	sep = np.array(transientCoords.separation(c2).arcsec)
	host = None
	if np.nanmin(sep) <= 1:
		host_idx = np.where(sep == np.nanmin(sep))[0][0]
		host = fullTable.iloc[[host_idx]]
	getHostStatsFromTransientName(host['Transient Name'].values[0])

# Returns basic statistics for the most likely
# host of a previously identified transient, using
# the pipeline outlined above.
# inputs - the coordinates of the transient to search
def getHostStatsFromTransientName(SNName):
	host = getHostFromTransientName(SNName)
	if len(host) > 0:
		if np.unique(host['NED_name']) != "":
			print("Found host %s.\n"%host['NED_name'].values[0])
		else:
			print("Found host %s"%host['objName'])
		print("RA, DEC (J2000): %f, %f"%(host['raMean'], host['decMean']))
		if host['NED_redshift'].values[0] != '':
			print("Redshift: %f"%host['NED_redshift'].values[0])
		print("PS1 rMag: %.2f"%host['rApMag'].values[0])
		print("g-r    r-i    i-z")
		print("%.2f   %.2f   %.2f"%(host['g-r'].values[0], host['r-i'].values[0], host['i-z'].values[0]))
	#    print("g-r: %.2f."%host['g-r'].values[0])
	#    print("r-i: %.2f."%host['r-i'].values[0])
	#    print("i-z: %.2f."%host['i-z'].values[0])
		#print("There are %i supernovae associated with this object.\n"%len(host))
		print("Associated supernovae: ")
		for SN in np.array(host['Transient Name']):
			print(SN)
	return

# Returns a postage stamp of the most likely
# host in one of the PS1 bands - g,r,i,z,y - as a
# fits file with radius rad, and plots the image.
# inputs
# outputs
def getHostImage(transientName='', band="grizy", rad=60, save=0):
	if transientName == '':
		print("Error! Please enter a supernova!\n")
		return
	fullTable = fullData()
	host = getHostFromTransientName(transientName)
	host.reset_index(drop=True, inplace=True)
	tempSize = int(4*float(rad))
	fn_save = host['objID'].values[0]
	if np.unique(host['NED_name']) != "":
		fn_save = host['NED_name'].values[0]
		print("Showing postage stamp for %s"%np.unique(host['NED_name'])[0])
	ra = np.unique(host['raMean'])[0]
	dec = np.unique(host['decMean'])[0]
	tempSize = int(4*rad)
	img = getcolorim(ra, dec, output_size=tempSize, size=tempSize, filters=band, format="png")
	plt.figure(figsize=(10,10))
	plt.imshow(img)
	if save:
		img.save("%s.png" % fn_save)
	return

# description
# inputs
# outputs
def getTransientSpectra(path, SNname):#
	SNname = remove_prefix(SNname, 'SN')
	files = glob.glob(path+"*%s*"%SNname)
	specFiles = []
	if len(files) > 0:
		for file in files:
			if remove_prefix(file, path).startswith("osc_"):
				spectra = pd.read_csv(file, header=None, names=['Wave(A)', 'I', 'Ierr'])
			else:
				spectra = pd.read_csv(file, delim_whitespace=True, header=None, names=['x', 'y', 'z'])
			#spectra = spectra.astype('float')
			specFiles.append(spectra)
	else:
		print("Sorry! No spectra found.")
	return specFiles

# description
# inputs
# outputs
def getHostSpectra(SNname, path):#
	SNname = remove_prefix(SNname, 'SN')
	files = glob.glob(path+"*%s_hostSpectra.csv*"%SNname)
	specFiles = []
	if len(files) > 0:
		for file in files:
			spectra = pd.read_csv(file)
			#spectra = spectra.astype('float')
			specFiles.append(spectra)
	else:
		print("Sorry! No spectra found.")
	return specFiles

# A cone search for all transient-host pairs
# within a certain radius, returned as a pandas dataframe
# inputs location, in astropy coordinates, and radius in arcsec
# outputs the data frame of all nearby pairs
def coneSearchPairs(coord, radius):
	fullTable = fullData()
	c2 = SkyCoord(fullTable['TransientRA']*u.deg, fullTable['TransientDEC']*u.deg, frame='icrs')
	sep = np.array(coord.separation(c2).arcsec)
	hosts = None
	if np.nanmin(sep) < radius:
		host_idx = np.where(sep < radius)[0]
		hosts = fullTable.iloc[host_idx]
		#hosts = pd.concat(hosts)
	return hosts

# Returns the full table of data
# inputs
# outputs
def fullData():
	fullTable = pd.read_csv("../database/GHOST.tar.gz")
	return fullTable
