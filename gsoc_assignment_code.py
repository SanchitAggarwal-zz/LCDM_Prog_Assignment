__author__ = 'saggarwal'
__date__ = '21 March 2014'

# Programming Assignment
from astropy.io import fits
from matplotlib import colors, cm, pyplot as plt
import numpy

# Reading Sample Data
image_data = fits.open('M86.fits')[0].data
image_size =  image_data.shape
image_mean = numpy.mean(image_data)
image_std  = numpy.std(image_data)
# Removing the Noise
Threshold  = image_mean + image_std
Moments = numpy.zeros(shape=(3,3))
M00 = 0
M01 = 0
M10 = 0
M11 = 0
M20 = 0
M02 = 0
for i in range(1, image_size[0]):
	for j in range(1, image_size[1]):
		if image_data[i][j] > Threshold:
			M00 = M00 + image_data[i][j]
			M10 = M10 + i * image_data[i][j]
			M01 = M01 + j * image_data[i][j]
			M11 = M11 + i*j*image_data[i][j]
			M02 = M02 + j*j*image_data[i][j]
			M20 = M20 + i*i*image_data[i][j]

Centroid_X = M10/M00
Centroid_Y = M01/M00

print repr(Centroid_X)+','+repr(Centroid_Y)

mu00 = M00
mu01 = 0
mu10 = 0
mu11 = M11 - Centroid_X*M01
mu20 = M20 - Centroid_X*M10
mu02 = M02 - Centroid_Y*M01
scale = 5    
theta = 0.5 * numpy.arctan2((2*mu11),(mu20 - mu02))
covarianceMatrix = numpy.array([[mu20/mu00,mu11/mu00],[mu11/mu00,mu02/mu00]])
#covarianceMatrix = numpy.array(mu20/mu00,mu11/mu00,mu11/mu00,mu02/mu00)
norm = colors.LogNorm(image_data.mean() + 0.5 * image_data.std(), image_data.max(), clip='True')
imgplot = plt.matshow(image_data, cmap=cm.gray, norm=norm, origin="lower")
eigvals = numpy.linalg.eigvals(covarianceMatrix)
axixLength = eigvals[1]/scale
axixLength2 = eigvals[0]/scale
Point_Ax = [Centroid_Y-axixLength*numpy.cos(theta), Centroid_Y+axixLength*numpy.cos(theta)]
Point_Bx = [Centroid_X-axixLength*numpy.sin(theta), Centroid_X+axixLength*numpy.sin(theta)]
Point_Ay = [Centroid_Y-axixLength2*numpy.cos(theta+numpy.pi/2), Centroid_Y+axixLength2*numpy.cos(theta+numpy.pi/2)]
Point_By = [Centroid_X-axixLength2*numpy.sin(theta+numpy.pi/2), Centroid_X+axixLength2*numpy.sin(theta+numpy.pi/2)]
plt.plot(Point_Ax,Point_Bx , 'g-')
plt.plot(Point_Ay,Point_By , 'r-')
circle=plt.Circle([Centroid_Y,Centroid_X,-1],5,color='g',fill=False);
fig = plt.gcf()
fig.gca().add_artist(circle)
plt.show()
fig.savefig('output.png')
