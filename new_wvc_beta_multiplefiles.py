

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pylab as pylab


import math as m

import numpy.ma as ma


from matplotlib.pyplot import figure,show
import os

import timeit
 
start = timeit.default_timer()

#------------------------------------------------------------------------------
"""
Leer varables lon lt de netcdf
lamar la funcion call_wvc..
"""

#---------------------------------------------------------------------
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


#   #Ignoram per ara
# #---------------------------------------------------------------

#   DLam_zero = np.equal(lon2,lon1)              # lon2.shape ha de ser igual a lon1.shape, però amb lon_left,lon_right_margin és fals

#   Alfa1[(lat2>lat1) & DLam_zero] =  90.0
#   Alfa2[(lat2>lat1) & DLam_zero] =  90.0
#   Alfa1[(lat2<lat1) & DLam_zero] = -90.0
#   Alfa2[(lat2<lat1) & DLam_zero] = -90.0          # A nes final? Ho fa per cada element?? S'hauria de canviar, nosaltres volem que ho faci quan es compleix
#                                                   # la condició per un element no per tota la matriu , seria estúpid
#   #--------------------------------------------------------------
    #Passam a deg, ra2deg
    #Alfa1_deg = np.rad2deg(Alfa1)
    #Alfa2_deg = np.rad2deg(Alfa2)


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
    #   print(x)


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
    #   Alpha_left = np.expand_dims(Alfa1_final[:,num],axis=1)
    #   print(Alpha_left,np.shape(Alpha_left))
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

    dFx_dx[m-1,:] = (Fx[m-1,:] - Fx[m-2,:])/dc
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







pathfiles_directory = "C:\\Users\\T.C\\Desktop\\Pràctiques en Empresa ICM\\Pràctiques redacció\\Practiques redacció\\Codi_practiques\\Carpeta unica\\arxius_oscat"

hist = []
div_hist = np.array(hist)
#div_hist = np.empty()

for filename in os.scandir(pathfiles_directory):
    if filename.is_file():
        if '.nc' in filename.path:

            fnetcdf = Dataset(filename.path)

            phiN = fnetcdf.variables['wind_dir'].__array__()
            U10 = fnetcdf.variables['wind_speed'].__array__()
            lon = fnetcdf.variables['lon'].__array__()
            lat = fnetcdf.variables['lat'].__array__()
            u = fnetcdf.variables['northward_wind'].__array__()
            v = fnetcdf.variables['eastward_wind'].__array__()

            Alpha = WVC_Orientation_obtain_alfa(lat,lon)
            uc,ua = oce2ca(U10,phiN,Alpha)
            div = div_numpy_2D_findiff_4_elements(uc,ua,25)

            div_final = np.concatenate(div)

            stop = timeit.default_timer()
            time = stop - start
            print('Time: {} seconds '.format(time))


            print('Esteim al fitxer {}'.format(filename.path))















div_numpy = div_final.flatten()

n,x,_ = plt.hist(div_numpy, bins = np.linspace(-2E-3,2E-3,120), histtype=u'step',log = True )
plt.clf()
bin_centers = 0.5*(x[1:]+x[:-1])
maxim = np.amax(n)
#print(maxim)
#y= range(len(n/maxim))
#x1=range(len(bin_centers))
plt.plot(bin_centers,n/maxim) ## using bin_centers rather than edges   ,n
plt.yscale('log')


plt.ylim(5E-7,1)


plt.title('Tropical Atlantic WDIV')  
plt.xlabel('WDIV s-1')  
plt.ylabel('Probability') 
plt.ticklabel_format(style='sci', axis='x',scilimits=(0,0))

plt.tight_layout()


plt.show()
































