"""
    PlotAnimation.py 0.1 - creates Animation with Titles for matplotlib

    Copyright (C) 2008 Paul Hilscher

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""



from pylab import *
from numpy import *
import string
import glob

from PIL import Image


prefixImages1 = sys.argv[1]
prefixImages2 = sys.argv[2]

# Get all pictures and sort to get right order
images1 =  glob.glob(prefixImages1+'*')
images2 =  glob.glob(prefixImages2+'*')

images1.sort()
images2.sort()


if len(images1) != len(images2):
  print len(images1)
  print len(images2)
  print "Unequal number of Files"

max_len = max(len(images1), len(images2))
for i in range(max_len):
  # Use the last picture if images have unequal size
  if i < len(images1):
    image1 = images1[i];
  else:
    image1 = images1[-1];

  if i < len(images2):
    image2 = images2[i];
  else:
    image2 = images2[-1];

  im_1frame = Image.open(image1)
  im_2frame = Image.open(image2)

  im_1frame = im_1frame.resize((800,800))
  im_2frame = im_2frame.resize((800,800))

  # Create new, white picture
  image = Image.new("RGB", (1750,900), (255,255,255))

  image.paste(im_1frame, (50, 50, 850,850))
  image.paste(im_2frame, (900, 50, 1700,850))


  # If we provide an Helios HDF5 file, print timestep to title
  # TODO
  #if(sys.argv[3] != None):





  image_save_str = "Image_" + string.zfill(i, int(log(max_len))) + ".png"
  image.save(image_save_str)
  print "Saved (", i, "/", max_len, ') : ', image1, '   ', image2, '  ---->   ' , image_save_str



