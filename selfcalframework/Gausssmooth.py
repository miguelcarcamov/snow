import numpy as np
from astropy.io import fits
from astropy.convolution import interpolate_replace_nans
from astropy.convolution import convolve as astropy_convolve
from scipy.signal import convolve as scipy_convolve
from scipy.signal import choose_conv_method

def Gauss_filter(inarray, stdev_x, stdev_y, PA, Plot=False, Spike=False,Boundary='extend',ConvoMethod='auto',UseAstropy=False):
	''' 
        inarray: input data array
        stdev_x: float
	1sigma std dev. along BMAJ in pixels 
        stdev_y: float
	1sigma std dev. along BMIN in pixels
        PA: East of North in degrees  '''
	
	IsCube=False
	LargerThanACube=False
	datashape=inarray.shape
	if (len(datashape)>3): 
		print("larger than a cube, looping over  assumming first 3 are the relevant ones")
		(nz0,ny0,nx0) = (datashape[-3],datashape[-2],datashape[-1])
		print( "looping over 3rd dimmension with nz0=",nz0)
		LargerThanACube=True
		inarray0=inarray.copy()
		inarray=inarray0[0,:]
		IsCube=True
	elif (len(datashape)==3):
		IsCube=True
		(nz0,ny0,nx0) = datashape
		print( "this is a cube, looping over 3rd dimmension with nz0=",nz0)
	elif (len(datashape)==2):
		(ny0,nx0) = datashape
	else:
		sys.exit("smaller than an image")


	nx = min(ny0,min( nx0, int(7.*stdev_x)))
	nx = nx + ((nx+1) % 2)
	ny = nx

	#x=np.arange(1,nx+1)
	#y=np.arange(1,ny+1)
	x=np.arange(0,nx)
	y=np.arange(0,ny)
	X, Y = np.meshgrid(x, y)
	X0 = (float(nx)-1)/2.
	Y0 = (float(ny)-1)/2.

	print("Kernel center: x0",X0,"Y0",Y0," Kernel dimmensions:",nx,ny)

	#----------
	theta =   np.pi * (90.-PA) / 180.  # 
	A = 1
	a = np.cos(theta)**2/(2.*stdev_x**2) + np.sin(theta)**2/(2.*stdev_y**2)
	b = np.sin(2*theta)/(4.*stdev_x**2) - np.sin(2.*theta)/(4.*stdev_y**2)
	c = np.sin(theta)**2/(2.*stdev_x**2) + np.cos(theta)**2/(2. *stdev_y**2)

	Z=A*np.exp(-(a*(X-X0)**2-2.*b*(X-X0)*(Y-Y0)+c*(Y-Y0)**2))

	Z /= np.sum(Z)


	if Plot==True:
		import matplotlib.pyplot as plt
		import matplotlib.cm as cm
		plt.imshow(Z, cmap=cm.jet,origin='lower')
		plt.show()

	print( "data dims",inarray.shape)
	print( "kernel  dims",Z.shape)




	if (IsCube):
		smootharray=inarray.copy()
		for k in range(nz0):
			print( "k",k)
			im=np.double(inarray[k,:,:])
			if UseAstropy:
                                im_smooth = astropy_convolve(im, Z)
			else:
                                ##im = interpolate_replace_nans(im, Z)
                                bestmethod=choose_conv_method(im, Z)
                                #if (bestmethod != ConvoMethod):
                                #        print("WARNING: not using the best ConvoMethod, which is ",bestmethod)
                                #        #im_smooth = scipy_convolve(im, Z,method=ConvoMethod,mode='same')
                                ##print("im ",im.shape,im.dtype)
                                ##print("Z",Z.shape,Z.dtype)
                                if ( (bestmethod != ConvoMethod) and (ConvoMethod != 'auto')):
                                        print("WARNING not using default  ConvoMethod, which is ",bestmethod)
                                        im_smooth = scipy_convolve(im, Z, mode='same', method=ConvoMethod)
                                else:
                                        im_smooth = scipy_convolve(im, Z, mode='same')
                                print("im_smooth.shape",im_smooth.shape,"im.shape",im.shape)

			smootharray[k,:,:]=im_smooth
			if Plot==True:
				plt.imshow(im_smooth, cmap=cm.magma, origin='lower')
				plt.show()
	elif (len(datashape)==2):
                if UseAstropy:
                        smootharray = astropy_convolve(inarray, Z)
                else:
                        # im = interpolate_replace_nans(im, Z)
                        im=np.double(inarray)
                        bestmethod=choose_conv_method(im, Z)
                        if ( (bestmethod != ConvoMethod) and (ConvoMethod != 'auto')):
                                print("WARNING not using default  ConvoMethod, which is ",bestmethod)
                                smootharray = scipy_convolve(im, Z, mode='same', method=ConvoMethod)
                        else:
                                smootharray = scipy_convolve(im, Z, mode='same')
                if Plot==True:
                        plt.imshow(smootharray, cmap=cm.magma,origin='lower')
                        plt.show()
                if (Spike):
                        im_spike=np.zeros(inarray.shape)
                        im_spike[int(nx0/2),int(ny0/2)]=1.
                        print("spike at",int(nx0/2),int(ny0/2))
                        smooth_spike = convolve(im_spike, Z,boundary='extend')
                        hdu = fits.PrimaryHDU()
                        hdu.data=im_spike
                        hdu.writeto('spike.fits',overwrite=True)
                        hdu.data=smooth_spike
                        hdu.writeto('spike_smooth.fits',overwrite=True)


				 
	if LargerThanACube:
		inarray0[0,:]=smootharray
		return inarray0
	else:
		return smootharray

