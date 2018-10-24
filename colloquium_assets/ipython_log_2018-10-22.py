########################################################
# Started Logging At: 2018-10-22 13:09:15
########################################################

########################################################
# # Started Logging At: 2018-10-22 13:09:15
########################################################
import pyavm
pyavm.AVM('Taurus_NIR.jpg')
rslt = pyavm.AVM('Taurus_NIR.jpg')
rslt.Image
rslt.to_wcs()
get_ipython().system('wget https://cdn.eso.org/images/large/eso1423c.jpg')
pyavm.AVM('eso1423c.jpg')
pyavm.AVM('eso1423c.jpg').to_wcs()
get_ipython().system('wget https://www.eso.org/public/archives/images/original/eso1423c.tif')
pyavm.AVM('eso1423c.tif').to_wcs()
get_ipython().magic('rm eso1423c.tif')
get_ipython().system('open eso1423c.jpg')
get_ipython().magic('ls ')
get_ipython().system('open ESA_Herschel_Taurus_1280w.jpg')
pyavm.AVM('ESA_Herschel_Taurus_1280w.jpg').to_wcs()
