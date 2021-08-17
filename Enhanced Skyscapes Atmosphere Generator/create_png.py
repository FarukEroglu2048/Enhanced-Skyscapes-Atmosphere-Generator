import numpy
import PIL.Image

array = numpy.fromfile("transmittance.bin", dtype=numpy.uint8)
array = numpy.reshape(array, (1024, 1024, 3))

image = PIL.Image.fromarray(array)
image.save("transmittance.png")