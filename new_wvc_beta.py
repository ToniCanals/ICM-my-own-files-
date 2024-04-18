from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pylab as pylab


import math as m

import numpy.ma as ma
#from geopy.distance import geodesic
#import cartopy.crs as ccrs
#import pygrib

from matplotlib.pyplot import figure,show
from matplotlib import cm

#from scipy.interpolate import griddata

#import xarray as xr
#import metpy as mp
#from metpy.calc import divergence,lat_lon_grid_deltas

import timeit

import netCDF4 as nc

start = timeit.default_timer()


#fpath = "C:\\Users\\T.C\\Desktop\\Pràctiques en Empresa ICM\\Codis\\IRIS\\antics nc anar moguent pel oscat_proves i set_spherical.py\\oscat_20120106_135109_ocsat2_12113_o_500_ovw_l2.nc"
fpath = 'oscat_20120106_135109_ocsat2_12113_o_500_ovw_l2_COPY.nc' # Si funciona
#fpath = 'oscat_20121201_011232_ocsat2_16891_o_250_ovw_l2_COPY.nc' #NO fnciona ?????
#fpath = 'nwp_20130205_06_03_125.nc'

fnetcdf = Dataset(fpath)

print(fnetcdf.variables.keys())
#print(fnetcdf.dimensions)

#-------------------------------------------------------------------------------

phiN = fnetcdf.variables['wind_dir'].__array__()


#u = fnetcdf.variables['northward_wind'].__array__()
#v = fnetcdf.variables['eastward_wind'].__array__()

lon = fnetcdf.variables['lon'].__array__()
lat = fnetcdf.variables['lat'].__array__()

U10 = fnetcdf.variables['wind_speed'].__array__()

v = U10*np.cos(np.deg2rad(phiN))
u = U10*np.sin(np.deg2rad(phiN))

#------------------------------------------------------------------------------
"""
Leer varables lon lt de netcdf
lamar la funcion call_wvc..
"""

#---------------------------------------------------------------------

slice_column_of_first_column = np.expand_dims(lon[:,0],axis =1)
#slice_column_of_first_column = slice_row_of_first_column[:, np.newaxis]
# b = np.expand_dims(a, axis=1) # Equivalent to x[:,np.newaxis]
#print(slice_column_of_first_column,slice_column_of_first_column.shape)
#print(slice_column_of_first_column,slice_column_of_first_column.shape)
i_half = int(lat.shape[1]/2)
print(i_half)
print(lat[:,i_half+1:].shape)
#---------------------------------------------------------------------
def WVC_Orientation(lat1,lat2,lon1,lon2):        # Cada un hauria de ser una columna per fer els càlculs??
    # En graus sexagesimals
	DLam = lon2 - lon1                           # lon2ndim != lon1, dim lon2 = Nx1 i lon1 = NxM/2, volem que iteras al llarg de la les columnes la resta?
#Ignoram per ara
#---------------------------------------------------------------

	DLam_zero = np.equal(lon2,lon1)              # lon2.shape ha de ser igual a lon1.shape, però amb lon_left,lon_right_margin és fals

	Alfa1[lat2>lat1 & DLam_zero] =  90
	Alfa2[lat2>lat1 & DLam_zero] =  90
	Alfa1[lat2<lat1 & DLam_zero] = -90
	Alfa2[lat2<lat1 & DLam_zero] = -90           # A nes final? Ho fa per cada element?? S'hauria de canviar, nosaltres volem que ho faci quan es compleix
	                                             # la condició per un element no per tota la matriu , seria estúpid

#--------------------------------------------------------------

	SinPhi1 = np.sin(np.deg2rad(lat1)) # DEg o rad dubte, pensar de no emplear numpy ja que serien escalars
	SinPhi2 = np.sin(np.deg2rad(lat2)) # Potser posar sin(Deg2Rad*lat2)
	CosPhi1 = np.cos(np.deg2rad(lat1))
	CosPhi2 = np.cos(np.deg2rad(lat2))


	SinDLam = np.sin(np.deg2rad(DLam))
	CosDLam = np.cos(np.deg2rad(DLam))

	CosC = SinPhi1*SinPhi2 + CosPhi1*CosPhi2*CosDLam
	SinC_inv = 1.0/(np.sqrt(1.0 - CosC**2))

	Aux = CosPhi2*SinDLam*SinC_inv

	Aux[Aux >= 1] = 1
	Aux[Aux <= -1] = -1

	Alfa1=Rad2Deg*np.arccos(Aux)


	Aux = CosPhi1*SinDLam*SinC_inv

	Aux[Aux >= 1] = 1
	Aux[Aux <= -1] = -1

	Alfa2=Rad2Deg*np.arccos(Aux)


	return Alfa1,Alfa2

def WVC_Orientation__dubte(lat1,lat2,lon1,lon2): # Cada un hauria de ser una columna per fer els càlculs
    # En graus sexagesimals
	DLam = lon2 - lon1                           # lon2ndim != lon1, dim lon2 = Nx1 i lon1 = NxM/2, volem que iteras al llarg de la les columnes la resta?

		#Ignoram per ara
    #---------------------------------------------------------------

	DLam_zero = np.equal(lon2,lon1)              # lon2.shape ha de ser igual a lon1.shape, però amb lon_left,lon_right_margin és fals
	m=int(lat1.shape[0])

	Alfa1 = np.empty((m,1)) #!!!!
	Alfa2 = np.empty((m,1)) #!!!!

	Alfa1[(lat2>lat1) & DLam_zero] =  90.0
	Alfa2[(lat2>lat1) & DLam_zero] =  90.0
	Alfa1[(lat2<lat1) & DLam_zero] = -90.0
	Alfa2[(lat2<lat1) & DLam_zero] = -90.0          # A nes final? Ho fa per cada element?? S'hauria de canviar, nosaltres volem que ho faci quan es compleix
	                                                # la condició per un element no per tota la matriu , seria estúpid
	#--------------------------------------------------------------


	SinPhi1 = np.sin(np.deg2rad(lat1)) # DEg o rad dubte, pensar de no emplear numpy ja que serien escalars
	SinPhi2 = np.sin(np.deg2rad(lat2)) # Potser posar sin(Deg2Rad*lat2)
	CosPhi1 = np.cos(np.deg2rad(lat1))
	CosPhi2 = np.cos(np.deg2rad(lat2))


	SinDLam = np.sin(np.deg2rad(DLam))
	CosDLam = np.cos(np.deg2rad(DLam))

	CosC = SinPhi1*SinPhi2 + CosPhi1*CosPhi2*CosDLam
	SinC_inv = 1.0/(np.sqrt(1.0 - CosC**2))

	Aux = CosPhi2*SinDLam*SinC_inv

	Aux[Aux >= 1] = 1
	Aux[Aux <= -1] = -1

	Alfa1=np.rad2deg(np.arccos(Aux))

	Aux = CosPhi1*SinDLam*SinC_inv

	Aux[Aux >= 1] = 1
	Aux[Aux <= -1] = -1

	Alfa2 = np.rad2deg(np.arccos(Aux)) # deg2rad rad2deg


# 	#Ignoram per ara
# #---------------------------------------------------------------

# 	DLam_zero = np.equal(lon2,lon1)              # lon2.shape ha de ser igual a lon1.shape, però amb lon_left,lon_right_margin és fals

# 	Alfa1[(lat2>lat1) & DLam_zero] =  90.0
# 	Alfa2[(lat2>lat1) & DLam_zero] =  90.0
# 	Alfa1[(lat2<lat1) & DLam_zero] = -90.0
# 	Alfa2[(lat2<lat1) & DLam_zero] = -90.0          # A nes final? Ho fa per cada element?? S'hauria de canviar, nosaltres volem que ho faci quan es compleix
# 	                                                # la condició per un element no per tota la matriu , seria estúpid
# 	#--------------------------------------------------------------
	#Passam a deg, ra2deg
	#Alfa1_deg = np.rad2deg(Alfa1)
	#Alfa2_deg = np.rad2deg(Alfa2)


	return Alfa1,Alfa2 #Alfa1,Alfa2




"""

def call_WVC_Orientation(lons, lats):

	lat_left_margin  = np.expand_dims(lats[:,0],axis=1)  # mirar como expandir dimensiones en numpy, tiene que ser una columa nx1 # --> Ya tene las dimensiones correctas de Nx1 ,no? [a,b,c,d,e,f, ... ]
	lat_right_margin = np.expand_dims(lats[:,-1],axis=1) # mirar como expandir dimensiones en numpy, tiene que ser una columa nx1
	lon_left_margin  = np.expand_dims(lons[:,0],axis=1)  # mirar como expandir dimensiones en numpy, tiene que ser una columa nx1
	lon_right_margin = np.expand_dims(lons[:,-1],axis=1) # mirar como expandir dimensiones en numpy, tiene que ser una columa nx1
	i_half    = int(lats.shape[1]/2)                     # Utilitzar numpy.int????
	lat_left  = lats[:,0:i_half]
	lat_right = lats[:,i_half+1:]
	lon_left  = lons[:,0:i_half]                         # En teoria tenen lon 76 pertant separam en 38 i 38 (0,37) i (38,75) 
	lon_right = lons[:,-1,i_half+1:]


	alfa1_fh, alfa2_fh = WVC_Orientation_calc(lat_left,lat_right_margin,lon_left,lon_right_margin)
	alfa1_sh, alfa2_sh = WVC_Orientation_calc(lat_left_margin,lat_right, lon_left_margin,lon_right)
	#concatenar 2 mitades con np.concatenate, alfa1_fh con alfa1_sh y viceversa

	alfa1 = np.concatenate((alfa1_fh, alfa2_fh), axis=1)
	alfa2 = np.concatenate((alfa2_fh,alfa2_sh), axis=1)

	#N = np.nanmax(np.shape(lat))

	# alfa1_fh, alfa2_fh = WVC_Orientation_calc1(lat_left,lat_right_margin,lon_left,lon_right_margin)
	# alfa1_sh, alfa2_sh = WVC_Orientation_calc2(lat_left_margin,lat_right, lon_left_margin,lon_right)


#--------------------------------------------------------------------------------------------------------------
	Alpha = np.full((N,1), np.nan)
    AlphaSet = np.full((N, 1), False)
#---------------------------------------------------------------------------------------------------------------
    SinPhi1 = np.sin(np.deg2rad(lat1)) # DEg o rad dubte, pensar de no emplear numpy ja que serien escalars
	SinPhi2 = np.sin(np.deg2rad(lat2)) # Potser posar sin(Deg2Rad*lat2)
	CosPhi1 = np.cos(np.deg2rad(lat1))
	CosPhi2 = np.cos(np.deg2rad(lat2))


	SinDLam = np.sin(np.deg2rad(DLam))
	CosDLam = np.cos(np.deg2rad(DLam))

	CosC = SinPhi1*SinPhi2 + CosPhi1*CosPhi2*CosDLam
	SinC_inv = 1.0/(np.sqrt(1.0 - CosC**2))

	Aux = CosPhi2*SinDLam*SinC_inv

	Aux[Aux >= 1] = 1
	Aux[Aux <= -1] = -1

	Alfa1=Rad2Deg*np.arccos(Aux)


	Aux = CosPhi1*SinDLam*SinC_inv

	Aux[Aux >= 1] = 1
	Aux[Aux <= -1] = -1

	Alfa2 = Rad2Deg*np.arccos(Aux)




"""

def sdo2uv(S,phio): # U10,phiN

	phiu = 90 - phio
	u = S*np.cos(np.deg2rad(phiu))
	v = S*np.sin(np.deg2rad(phiu))

	return u,v

def oce2ca(S,phio,Alfad): #U10,phiN

	u,v = sdo2uv(S,phio)

	CosA = np.cos(np.deg2rad(Alfad))
	SinA = np.sin(np.deg2rad(Alfad))

	uc = CosA*u + SinA*v 
	ua = -SinA*u + CosA*v

	return uc,ua




def WVC_Orientation_obtain_alfa(lat,lon):

	half = int(lat.shape[1]/2)

	# Alfa1_final_left= np.array([])
	# Alfa2_final_left=np.array([])
	# Alfa1_final_right=np.array([])
	# Alfa2_final_right=np.array([])


	# Alfa1_final_left= np.empty(N,1)
	# Alfa2_final_left=np.empty(N,1)
	# Alfa1_final_right=np.empty(N,1)
	# Alfa2_final_right=np.empty(N,1)


	# lat_left  = lats[:,0:i_half]
	# lat_right = lats[:,i_half+1:]
	# lon_left  = lons[:,0:i_half]                         # En teoria tenen lon 76 pertant separam en 38 i 38 (0,37) i (38,75) 
	# lon_right = lons[:,-1,i_half+1:]

	# alfa1_fh, alfa2_fh = WVC_Orientation_calc1(lat_left,lat_right_margin,lon_left,lon_right_margin)   #First
	# alfa1_sh, alfa2_sh = WVC_Orientation_calc2(lat_left_margin,lat_right, lon_left_margin,lon_right)

	Alfa1_anterior_l,Alfa2_anterior_l = WVC_Orientation__dubte(lat[:,0],lat[:,-1],lon[:,0],lon[:,-1]) # First iteration, is it right??
	Alfa1_anterior_l = np.expand_dims(Alfa1_anterior_l,axis =1)
	#print(Alfa1_anterior_l)
	# FER UN PRINT PER VEURE SI ESTÀN BÉ AQUESTS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Alfa2_anterior_l = np.expand_dims(Alfa2_anterior_l,axis =1)

	print(np.shape(Alfa1_anterior_l))

	Alfa2_anterior_r,Alfa1_anterior_r = WVC_Orientation__dubte(lat[:,0],lat[:,half],lon[:,0],lon[:,half]) # First iteration, is it right??
	Alfa1_anterior_r = np.expand_dims(Alfa1_anterior_r,axis =1)
	Alfa2_anterior_r = np.expand_dims(Alfa2_anterior_r,axis =1)


	#print(Alfa1_anterior_r)
	#print('Imprimim Alfa1_anterior_r de la fila {} de python, {} - 1 matlab:'.format(half,half))
	#for x in np.nditer(Alfa2_anterior_r):
	#	print(x)


	for num in range(half-1): #-1 pq rang va fins a N-1 i -1 una altra vegada per la primera iteració

		Alfa1,Alfa2 = WVC_Orientation__dubte(np.expand_dims(lat[:,num + 1],axis = 1),np.expand_dims(lat[:,-1],axis = 1),np.expand_dims(lon[:,num + 1],axis =1),np.expand_dims(lon[:,-1],axis =1)) #Sumam + 1 pq no comenci a zero
		# Ja tenim la primera iteració
		# Volem treure una matriu de la primera opració  Nx1 guardarla i fer apend, conservar la concatenada de la iteracio anterior i fer concatenate de la iteració actual

		print(np.shape(Alfa1))

		Alfa1_final_left = np.hstack((Alfa1_anterior_l,Alfa1))       # hem de guardar la iteració anterior i fer concatenate amb la actual
		Alfa2_final_left = np.hstack((Alfa2_anterior_l,Alfa2))
		Alfa1_anterior_l = Alfa1_final_left
		Alfa2_anterior_l = Alfa2_final_left
		print('Dimensions de Alfa1_final_left : {}'.format(np.shape(Alfa1_final_left)))

		Alfa2,Alfa1 = WVC_Orientation__dubte(np.expand_dims(lat[:,0],axis = 1),np.expand_dims(lat[:,half + num + 1],axis =1),np.expand_dims(lon[:,0],axis=1),np.expand_dims(lon[:,half + num + 1],axis=1)) 
		# Suma +1 pq ja hem calculat l'anterior
		Alfa1_final_right = np.hstack((Alfa1_anterior_r,Alfa1))
		Alfa2_final_right = np.hstack((Alfa2_anterior_r,Alfa2))
		Alfa1_anterior_r = Alfa1_final_right
		Alfa2_anterior_r = Alfa2_final_right


	Alfa1_final = np.hstack((Alfa1_final_left,Alfa1_final_right))
	Alfa2_final = np.hstack((Alfa2_final_left,Alfa2_final_right))

	#Problema de machaque de la matriu alfa1 i 2 inicials de la primera iteració respecte a la segona

	Alpha_left = Alfa1_final_left
	Alpha = np.hstack((Alpha_left,Alfa1_final_right))


	# Alpha_left = Alfa1_final[:,0:half]                    #0  a 37
	# print(Alpha_left,np.shape(Alpha_left))
	# Alpha_right = Alfa2_final[:,0:half]
	# Alpha = np.hstack((Alpha_left,Alpha_right))   #38 a 76
	# print(Alpha_left,np.shape(Alpha))

	# for num in range(half):
	# 	Alpha_left = np.expand_dims(Alfa1_final[:,num],axis=1)
	# 	print(Alpha_left,np.shape(Alpha_left))
	# Alpha = np.hstack((Alpha_left,Alfa2[:,half-]))

		#Amb alfa1_final i alfa2_final com treure Alpha?

	print('Dimensions de Alfa1_final : {}'.format(np.shape(Alfa1_final)))


	return Alpha # Alfa1_final,Alfa2_final,

def div_numpy_2D_vector_field(Fx,Fy,gridkm): #Easy one
	#divergence computes the partial derivatives in its definition by using finite differences. For interior data points, the partial derivatives are calculated using central difference.
	#For data points along the edges, the partial derivatives are calculated using single-sided (forward) difference.
	#For example, consider a 2-D vector field F that is represented by the matrices Fx and Fy at locations X and Y with size m-by-n.
	#The locations are 2-D grids created by [X,Y] = meshgrid(x,y), where x is a vector of length n and y is a vector of length m.
	# Divergence then computes the partial derivatives ∂Fx / ∂x and ∂Fy / ∂y as


	#PYTHON INDEXES GO FROM 0 TO NP.SHAPE(ARRAY[j])-1

	#[x,y] = np.meshgrid(x,y, indexing='ij')

	#gridkm ha de ser #25km en m

	m = int(Fx.shape[0])
	n = int(Fx.shape[1])

	print('m: {} i n: {}'.format(m,n))

	#Must have equal dimensions both Fx and Fy


	dFx = np.empty((m,n))
	dFy = np.empty((m,n))	


	#for data points at the left and right edges.

	dFx[:,0] = (Fx[:,1] - Fx[:,0])/(gridkm) #and
	dFx[:,-1] = (Fx[:,n-1] - Fx[:,n-2])/(gridkm) 

	#for interior data points.

	for i in range(1,n-3):
	    dFx[:,i] = (Fx[:,i] - Fx[:,i-1])/(2*gridkm) #and

	for j in range(1,m-3):
		dFy[j,:] = (Fy[j+1,:] - Fy[j-1,:])/(2*gridkm)
	


	#for data points at the top and bottom edges.

	dFy[0,:] = (Fy[1,:] - Fy[0,:])/(gridkm)#/(y[1,:] - y[0,:]) #and
	dFy[m-1,:] = (Fy[m-1,:] - Fy[m-2,:])/(gridkm)#/(y[m-1,:] - y[m-2:])   # last element is m-1 because it starts at 0 andd has len m


	#The numerical divergence of the vector field is equal to div = dFx + dFy.

	div = dFx + dFy
	#rot = dFx - dFy

	return div

def div_numpy_2D_findiff_4_elements(Fx,Fy,gridkm):

	# This code makes use of a convolution style algorithm to accelerate the calculation, with a 2x2 convolution kernel or filter.
	# Bottom and right edges are calculated with single differences 

	# n-1 column needs a look (correction required)

	m = int(Fx.shape[0])
	n = int(Fx.shape[1])

	print('m: {} i n: {}'.format(m,n))

	#Must have equal dimensions both Fx and Fy

	#gridkm ha de ser #25km en m

	dc = gridkm*1000 # in meters
	da = gridkm*1000

	dFx = np.empty((m,n))
	dFy = np.empty((m,n))

	roll_1 = np.roll(Fx,(-1,-1),axis=(0,1)) # u[i,j] <== u[i+1,j+1]
	roll_2 = np.roll(Fx,(-1,0),axis=(0,1))  # u[i,j] <== u[i+1,j]
	roll_3 = np.roll(Fx, (0,-1), axis=(0,1)) # u[i,j] <== u[i,j+1]

	roll_4 = np.roll(Fy,(-1,-1),axis=(0,1)) # u[i,j] <== u[i+1,j+1]
	roll_5 = np.roll(Fy,(-1,0),axis=(0,1))  # u[i,j] <== u[i+1,j]
	roll_6 = np.roll(Fy, (0,-1), axis=(0,1)) # u[i,j] <== u[i,j+1]

	dFx_dx = 0.5*(roll_1 - roll_2 + roll_3 - Fx)/dc

	dFy_dy = 0.5*(roll_4 - roll_6 + roll_5 - Fy)/da

	#Edges (Furthest right column and furthest down colum)


	dFx_dx[m-1,:] = (Fx[m-1,:] - Fx[m-2,:])/dc  # Comprovar si està girat 
	dFx_dx[:,n-1] = (Fx[:,n-1] - Fx[:,n-2])/dc

	dFy_dy[m-1,:] = (Fy[m-1,:] - Fy[m-2,:])/da
	dFy_dy[:,n-1] = (Fy[:,n-1] - Fy[:,n-2])/da


	div = dFx_dx + dFy_dy


	return div


def div_straightouttamatlab(uc,ua,gridkm):

	na = int(ua.shape[0])
	nc = int(ua.shape[1])

	

	dc = gridkm*1000 # in meters
	da = gridkm*1000

	duc_dc = np.empty((na,nc))
	duc_da = np.empty((na,nc))
	dua_dc = np.empty((na,nc))
	dua_da = np.empty((na,nc))


	for j0 in range(0,na-1):
		j1 = j0 + 1;
		for i0 in range(0,nc-1):
			i1 = i0 + 1
			duc_dc[j0,i0] = 0.5*((uc[j1,i1] - uc[j1,i0]) + (uc[j0,i1] - uc[j0,i0]))/dc
			duc_da[j0,i0] = 0.5*((uc[j1,i0] - uc[j0,i0]) + (uc[j1,i1] - uc[j0,i1]))/da
			dua_dc[j0,i0] = 0.5*((ua[j1,i1] - ua[j1,i0]) + (ua[j0,i1] - ua[j0,i0]))/dc
			dua_da[j0,i0] = 0.5*((ua[j1,i1] - ua[j0,i1]) + (ua[j1,i0] - ua[j0,i0]))/da

	duc_dc[na-1,:] = (uc[na-1,:] - uc[na-2,:])/dc
	duc_dc[:,nc-1] = (uc[:,nc-1] - uc[:,nc-2])/da

	dua_da[na-1,:] = (ua[na-1,:] - ua[na-2,:])/dc
	dua_da[:,nc-1] = (ua[:,nc-1] - ua[:,nc-2])/da 


	Div  = duc_dc + dua_da
	Curl = dua_dc - duc_da

	return Div


def add_variable_2nc_file_out(nc_file_path, variable, metadata):
    """
    Adds variable into an existing NetCDF file
    :param nc_file_path: Filepath to the NetCDF file to be updated
    :param variable: Variable 2d numpy array to be added
    :param metadata: Metadata dictionary that contains variable description that has the following fields
    variable_metadata = {
            "short_name": "t",
            "long_name": "Surface temperature at 10 m",
            "units": "K",
            "type": "float"}
    :return: Confirmation that the variable has been added to print out to the log output
    """
    ncfile = Dataset(nc_file_path, mode='r+')
    lon_nc = ncfile.dimensions['lon'].name
    lat_nc = ncfile.dimensions['lat'].name
    time_nc = ncfile.dimensions['time'].name
    var = ncfile.createVariable(metadata['short_name'], metadata['type'], (time_nc, lat_nc, lon_nc))
    var.standard_name = metadata['short_name']
    var.long_name = metadata['long_name']
    var.units = metadata['units']
    if "scale_factor" in metadata:
        var.scale_factor = metadata['scale_factor']
    ncfile.variables[metadata['short_name']][:] = variable
    ncfile.close()
    return "Added variable " + metadata['short_name'] + " to " +  nc_file_path





#Problema ara de com feim quin elements de alfa1,alfa2 acabin en quins elements de Alpha
#
#print(lat[2,3],lon[2,3],lat[2,-1],lon[2,-1])


#Alfa1,Alfa2= WVC_Orientation__dubte(lat[:,3],lat[:,-1],lon[:,3],lon[:,-1])

Alpha = WVC_Orientation_obtain_alfa(lat,lon)
#print(Alpha) 


uc,ua = oce2ca(U10,phiN,Alpha)


# print(uc,ua)

#div = div_numpy_2D_vector_field(ua,uc,25000)
# print('Això és div : {}'.format(div))

# div_to_numpy = div.flatten()
# div_restraint = div[abs(div) <= 1E-3]



#div1 = div_straightouttamatlab(uc,ua,25)
div = div_numpy_2D_findiff_4_elements(uc,ua,25)
#---------------------------------------------------------------------------------
# nc_file_path = fpath
# variable = div
# metadata = {
#             "short_name": "div",
#             "long_name": "Wind divergence at 10 m",
#             "units": "s-1",
#             "type": "float"}
# add_variable_2nc_file_out(nc_file_path, variable, metadata)
#--------------------------------------------------------------------------------
#div = np.array(div)
#div =div.data

#print('La diferència entre els algoritmes és: {}'.format(np.sum(div1[:-1,:-1] - div[:-1,:-1])))

# #Let's set alfa1 then to  Alpha to test it out

#n=464

# #n=2
# print('Dimensions de alfa : {}'.format(np.shape(Alpha))) #Dimensions problem we have[a,b,c,d,e,f,g,...] and we want [[a,b,c,d,e,f,g,...]]
# #print('Elements de la {} columna :'.format(n),Alfa1[n,:])
# print('Elements de la fila {}:'.format(n),Alpha[n,:])

# # print('Elements de la quarta columna :',Alfa1[:,3])
# # #print('Elements de la quarta columna :',Alpha[:,3])




# print('Imprimim Alpha de la fila {} de python, {} - 1 matlab:'.format(n,n))
# for x in np.nditer(Alpha[n,:]):
# 	print(x)

# print('Imprimim uc de la fila {} de python, {} - 1 matlab:'.format(n,n))
# for x in np.nditer(ua[n,:]):
# 	print(x)

# print('Imprimim div de la fila {} de python, {} - 1 matlab:'.format(n,n))
# for x in np.nditer(div[n,:]):
# 	print(x)





# print('Dimensions de phiN : {}'.format(np.shape(phiN)))
# print('Dimensions de U10 : {}'.format(np.shape(U10)))
# print('Dimensions de Alfa1 : {}'.format(np.shape(Alfa1)))
# print('Dimensions de Alfa2: {}'.format(np.shape(Alfa2)))




# uc,ua = oce2ca(U10,phiN,Alpha)

# div = div_numpy_2D_vector_field(ua,uc,lon,lat)

# print(div)

# print(uc,ua)


# #-------------------------------------------

# xr_uc = xr.DataArray(uc)
# xr_ua = xr.DataArray(ua)



# dx,dy = mp.calc.lat_lon_grid_deltas(lon,lat)



# div = mp.calc.divergence(xr_uc,xr_ua,dx= dx,dy= dy)

#fig, ax = plt.subplots()

# #x_values = range(len(months))
# div_numpy = div.to_numpy()

#div_numpy = div.flatten()
#div_numpy1 = div1.flatten()

# #----------------------------------------------------------------------------------------------
# n,x,_ = plt.hist(div_numpy, bins = np.linspace(-2E-3,2E-3,60), histtype=u'step',log = True )
# plt.clf()
# bin_centers = 0.5*(x[1:]+x[:-1])
# maxim = np.amax(n)
# #print(maxim)
# #y= range(len(n/maxim))
# #x1=range(len(bin_centers))
# plt.plot(bin_centers,n/maxim) ## using bin_centers rather than edges   ,n
# plt.yscale('log')


# plt.ylim(5E-7,1)


# plt.title('Tropical Atlantic WDIV (Optimised)')  
# plt.xlabel('WDIV s-1')  
# plt.ylabel('Probability') 
# plt.ticklabel_format(style='sci', axis='x',scilimits=(0,0))

# fig.tight_layout()


# plt.show()


# n,x,_ = plt.hist(div_numpy1, bins = np.linspace(-2E-3,2E-3,60), histtype=u'step',log = True )
# plt.clf()
# bin_centers = 0.5*(x[1:]+x[:-1])
# maxim = np.amax(n)
# #print(maxim)
# #y= range(len(n/maxim))
# #x1=range(len(bin_centers))
# plt.plot(bin_centers,n/maxim) ## using bin_centers rather than edges   ,n
# plt.yscale('log')


# plt.ylim(5E-7,1)


# plt.title('Tropical Atlantic WDIV')  
# plt.xlabel('WDIV s-1')  
# plt.ylabel('Probability') 
# plt.ticklabel_format(style='sci', axis='x',scilimits=(0,0))

# fig.tight_layout()


# plt.show()








# stop = timeit.default_timer()
# time = stop - start

# print('Time: {} seconds '.format(time))

#-------------------------------------------------------
"""
fig, ax = plt.subplots()



fig = plt.figure(figsize = (15, 20))
# scipy interp. cubic
f = interp2d(lon, lat, div, kind='cubic')

# xnew = np.arange(323, 350, .1)
# ynew = np.arange(0, 50, .1)
# data1 = f(xnew,ynew)
# Xn, Yn = np.meshgrid(xnew, ynew)
p=plt.pcolormesh(lon, lat, f, cmap='RdBu')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
fig.colorbar(p, ax=ax,format='%.0e')
plt.title(r'Divergence map', fontsize=25)
plt.show()
"""
#------------------------------------------------------

"""
fn = "C:\\Users\\T.C\\Desktop\\Pràctiques en Empresa ICM\\Codis\\IRIS\\nous ncs\\new_beta.nc"
ds = nc.Dataset(fn,'w', format = 'NETCDF4')

lat = ds.createDimension('lat')
lon = ds.createDimension('lon')

#div = ds.createDimension('div',10)

lats = ds.createVariable('lat','f4',('lat',))
lons = ds.createVariable('lon','f4',('lon',))

div = ds.createVariable('div','f4',())


lats = ds.createVariable('lat', 'f4', ('lat',))
lons = ds.createVariable('lon', 'f4', ('lon',))
value = ds.createVariable('value', 'f4', ('lat', 'lon',))
value.units = 'Unknown'

# lats[:] = lat
# lons[:] = lon

lats = lat
lons = lon

print('var size before adding data', value.shape)

value[0, :, :] = div
value = div

# print('var size after adding first data', value.shape)
# xval = np.linspace(0.5, 5.0, 10)
# yval = np.linspace(0.5, 5.0, 10)
# value[1, :, :] = np.array(xval.reshape(-1, 1) + yval)

# print('var size after adding second data', value.shape)

ds.close()
"""
# #----------------------------------------------------------------------------------------------
"""
# plot for smoothing kind of like an interpolation
fig = plt.figure(figsize = (15, 20))
ax = plt.axes()
p=plt.contourf(lon,lat,div,1000,cmap =cm.jet,vmin = -1E-4,vmax =1E-4)
p.set_clim(-0.0001, 0.0001)
plt.colorbar(p,boundaries=np.linspace(-0.0001, 0.0001, 10))
#fig.colorbar(p,ax=ax,format='%.0e') # =cbar
#cbar.set_ticks(np.arange(-1E-4,1E-4,1E-5))
#plt.imshow(zi)#, extent=(0,1,0,1), origin='lower')
#plt.contourf(lon_int,lat_int,zi,cmap = 'jet', vmin = -1E-4,vmax = 1E-4)
#p = ax.pcolor(lon,lat,zi,cmap = 'jet', vmin = -1E-4,vmax = 1E-4)
#plt.plot(lon,lat,zi)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title(r'Divergence map', fontsize=25)
plt.show()
"""
#----------------------------------------------------------------------------------------------



#-------------------------------------------------------
#np.set_printoptions(threshold=100)
#print(div)


#n=465
#print('Imprimim div de la fila {} de python, {} - 1 matlab:'.format(n,n))
#for x in np.nditer(div[n,:]):
#	print(x)
#fig, ax = plt.subplots()

"""
fig = plt.figure(figsize = (15, 20))
ax = plt.axes()
plt.xlabel('Longitude')
plt.ylabel('Latitude')

#lon[lon  180] = lon[lon > 180] - 360
#lon[]
p = ax.pcolor(lon[280:520],lat[280:520],div[280:520],cmap = 'jet', vmin = -4E-4,vmax = 4E-4)
#p = ax.pcolormesh(lon,lat,div,cmap = 'jet', vmin = -1E-4,vmax = 1E-4,shading='gouraud')
#p = ax.pcolor(lon,lat,div,cmap = 'jet', vmin = -1E-4,vmax = 1E-4)#,shading='gouraud')
#p = ax.plot(div_hist )
fig.colorbar(p, ax=ax,format='%.0e')
plt.title(r'Divergence map', fontsize=25)
plt.show()

"""