from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np


import math as m



from matplotlib.pyplot import figure,show
#from matplotlib import cm


import timeit

import netCDF4 as nc
"""
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
Both div_numpy_2D_findiff_4_elements(Fx,Fy,gridkm) and div_straightoutofmatlab(uc,ua,gridkm) do not have yet the half-grid and the borders are not well calculated

In order to be called one must do as in this script file: 

1) Have all the same libraries
2) Call the functions as done below:
fpath ="....nc"
fnetcdf = Dataset(fpath)

phiN = fnetcdf.variables['wind_dir'].__array__()

#To be uncommented or commented wether we use 500 or 250 type files
#u = fnetcdf.variables['northward_wind'].__array__()
#v = fnetcdf.variables['eastward_wind'].__array__()
lon = fnetcdf.variables['lon'].__array__()
lat = fnetcdf.variables['lat'].__array__()

U10 = fnetcdf.variables['wind_speed'].__array__()

#To be commented when using 250 type files
v = U10*np.cos(np.deg2rad(phiN))
u = U10*np.sin(np.deg2rad(phiN))

Alpha = WVC_Orientation_obtain_alfa(lat,lon)   
uc,ua = oce2ca(U10,phiN,Alpha)

div1 = div_straightoutofmatlab(uc,ua,25)
div = div_numpy_2D_findiff_4_elements(uc,ua,25)  Call only one of these two*


Note*: For performance reasons it is more convenient to use the div_numpy_2D_findiff_4_elements(uc,ua,25)
Author: Antoni Canals Gayà
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""

start = timeit.default_timer()


#fpath = "C:\\Users\\T.C\\Desktop\\Pràctiques en Empresa ICM\\Codis\\IRIS\\antics nc anar moguent pel oscat_proves i set_spherical.py\\oscat_20120106_135109_ocsat2_12113_o_500_ovw_l2.nc"
#fpath = 'oscat_20120106_135109_ocsat2_12113_o_500_ovw_l2_COPY.nc'  #500 example
fpath = 'oscat_20121201_011232_ocsat2_16891_o_250_ovw_l2_COPY.nc' #250 example
#fpath = 'nwp_20130205_06_03_125.nc'

fnetcdf = Dataset(fpath)

#print(fnetcdf.variables.keys())
#print(fnetcdf.dimensions)

#-------------------------------------------------------------------------------

phiN = fnetcdf.variables['wind_dir'].__array__()

#To be uncommented or commented wether we use 500 or 250 type files
u = fnetcdf.variables['northward_wind'].__array__()
v = fnetcdf.variables['eastward_wind'].__array__()

lon = fnetcdf.variables['lon'].__array__()
lat = fnetcdf.variables['lat'].__array__()

U10 = fnetcdf.variables['wind_speed'].__array__()

#To be commented when using 250 type files

#v = U10*np.cos(np.deg2rad(phiN))
#u = U10*np.sin(np.deg2rad(phiN))



def WVC_Orientation_optimized(lat1,lat2,lon1,lon2): # Each argument must be a Nx1 column, check how WVC_Orientation_obtain_alfa(lat,lon) calls theis function

    # In degrees
	DLam = lon2 - lon1
    #---------------------------------------------------------------

	DLam_zero = np.equal(lon2,lon1)              # If both longitudes are the same we require same matrix shapes
	m=lat1.shape[0]

	Alfa1 = np.empty((m,1)) #!!!!
	Alfa2 = np.empty((m,1)) #!!!!

	Alfa1[(lat2>lat1) & DLam_zero] =  90.0
	Alfa2[(lat2>lat1) & DLam_zero] =  90.0
	Alfa1[(lat2<lat1) & DLam_zero] = -90.0
	Alfa2[(lat2<lat1) & DLam_zero] = -90.0          

	#--------------------------------------------------------------


	SinPhi1 = np.sin(np.deg2rad(lat1)) 
	SinPhi2 = np.sin(np.deg2rad(lat2)) 
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


	return Alfa1,Alfa2 #Alfa1,Alfa2




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


	Alfa1_anterior_l,Alfa2_anterior_l = WVC_Orientation_optimized(lat[:,0],lat[:,-1],lon[:,0],lon[:,-1]) # First iteration
	Alfa1_anterior_l = np.expand_dims(Alfa1_anterior_l,axis =1)

	Alfa2_anterior_l = np.expand_dims(Alfa2_anterior_l,axis =1)

	#print(np.shape(Alfa1_anterior_l)) Check proper shape

	Alfa2_anterior_r,Alfa1_anterior_r = WVC_Orientation_optimized(lat[:,0],lat[:,half],lon[:,0],lon[:,half]) # First iteration
	Alfa1_anterior_r = np.expand_dims(Alfa1_anterior_r,axis =1)
	Alfa2_anterior_r = np.expand_dims(Alfa2_anterior_r,axis =1)


	for num in range(half-1): # -1 because the range goes up to N-1 and -1 again for the first iteration #-1 pq rang va fins a N-1 i -1 una altra vegada per la primera iteració

		Alfa1,Alfa2 = WVC_Orientation_optimized(np.expand_dims(lat[:,num + 1],axis = 1),np.expand_dims(lat[:,-1],axis = 1),np.expand_dims(lon[:,num + 1],axis =1),np.expand_dims(lon[:,-1],axis =1)) #Sumam + 1 pq no comenci a zero
		# We have the first iteration\Ja tenim la primera iteració 
		# We want to obtain a matrix of the first operation Nx1, save it concatenate, save the concatenated matrix of the previous iteration and concatenate the current one 
		#\Volem treure una matriu de la primera opració  Nx1 guardarla i fer apend, conservar la concatenada de la iteracio anterior i fer concatenate de la iteració actual

		#print('Dimensions of Alfa1: {}\n'.format(np.shape(Alfa1)))#CAN BE COMMENTED OR REMOVED TO INCREASE SPEED. Dimensions of the arrays that we concatenate.(Columns)

		Alfa1_final_left = np.hstack((Alfa1_anterior_l,Alfa1))      #We have to save the current iteration output and concatenate with the current iteration # Hem de guardar la iteració anterior i fer concatenate amb la actual
		Alfa2_final_left = np.hstack((Alfa2_anterior_l,Alfa2))
		Alfa1_anterior_l = Alfa1_final_left
		Alfa2_anterior_l = Alfa2_final_left
		#print('Dimensions of Alfa1_final_left : {}'.format(np.shape(Alfa1_final_left))) #We are checking for porper dimensions . CAN BE COMMENTED OR REMOVED TO ACCELERATE SPEED. Dimensions of the concatenated arrrays

		Alfa2,Alfa1 = WVC_Orientation_optimized(np.expand_dims(lat[:,0],axis = 1),np.expand_dims(lat[:,half + num + 1],axis =1),np.expand_dims(lon[:,0],axis=1),np.expand_dims(lon[:,half + num + 1],axis=1)) 
		# Suma +1 pq ja hem calculat l'anterior
		Alfa1_final_right = np.hstack((Alfa1_anterior_r,Alfa1))
		Alfa2_final_right = np.hstack((Alfa2_anterior_r,Alfa2))
		Alfa1_anterior_r = Alfa1_final_right
		Alfa2_anterior_r = Alfa2_final_right


	Alfa1_final = np.hstack((Alfa1_final_left,Alfa1_final_right))
	Alfa2_final = np.hstack((Alfa2_final_left,Alfa2_final_right))


	Alpha_left = Alfa1_final_left
	Alpha = np.hstack((Alpha_left,Alfa1_final_right))



	#print('Dimensions of Alfa1_final : {}'.format(np.shape(Alfa1_final)))


	return Alpha # Alfa1_final,Alfa2_final


def div_numpy_2D_findiff_4_elements(Fx,Fy,gridkm):

	# This code makes use of a convolution style algorithm to accelerate the calculation, with a 2x2 convolution kernel or filter.
	# Bottom and right edges are calculated with single differences 

	# n-1 column needs a look (correction required)

	m = int(Fx.shape[0])
	n = int(Fx.shape[1])

	#print('m: {} i n: {}'.format(m,n))

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


	dFx_dx[m-1,:] = (Fx[m-1,:] - Fx[m-2,:])/dc  #Should be changed
	dFx_dx[:,n-1] = (Fx[:,n-1] - Fx[:,n-2])/dc

	dFy_dy[m-1,:] = (Fy[m-1,:] - Fy[m-2,:])/da
	dFy_dy[:,n-1] = (Fy[:,n-1] - Fy[:,n-2])/da


	div = dFx_dx + dFy_dy


	return div

"""
def div_straightoutofmatlab(uc,ua,gridkm):

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

	duc_dc[na-1,:] = (uc[na-1,:] - uc[na-2,:])/dc  #Should be changed
	duc_dc[:,nc-1] = (uc[:,nc-1] - uc[:,nc-2])/da

	dua_da[na-1,:] = (ua[na-1,:] - ua[na-2,:])/dc
	dua_da[:,nc-1] = (ua[:,nc-1] - ua[:,nc-2])/da 


	Div  = duc_dc + dua_da
	Curl = dua_dc - duc_da

	return Div

"""






Alpha = WVC_Orientation_obtain_alfa(lat,lon)

uc,ua = oce2ca(U10,phiN,Alpha)






#div1 = div_straightoutofmatlab(uc,ua,25)
div = div_numpy_2D_findiff_4_elements(uc,ua,25)



#Plotting examples
#-------------------------------------------



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



stop = timeit.default_timer()
time = stop - start

print('Time: {} seconds '.format(time))

"""
# #----------------------------------------------------------------------------------------------

# plot for smoothing kind of like an interpolation
fig = plt.figure(figsize = (15, 20))
ax = plt.axes()
p=plt.contourf(lon,lat,div,1000,cmap =cm.jet,vmin = -1E-4,vmax =1E-4)
p.set_clim(-0.0001, 0.0001)
#plt.colorbar(p,boundaries=np.linspace(-0.0001, 0.0001, 10))
fig.colorbar(p,ax=ax,format='%.0e') # =cbar
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
#---------------------------------------------------------------------------------------------


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
