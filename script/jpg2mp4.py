

import cv2
import numpy as np
import glob

# https://blog.csdn.net/m0_51004308/article/details/125439217

img_array = []
filenamelist = []
for filename in glob.glob('*.jpg'):
   img = cv2.imread(filename)
   height, width, layers = img.shape
   size = (width,height)
   img_array.append(img)
   filenamelist.append(filename)


print(filenamelist)
fourcc = int(cv2.VideoWriter_fourcc(*'avc1'))
out = cv2.VideoWriter('video.mp4',fourcc, 15, size)

for i in range(len(img_array)):
   out.write(img_array[i])
out.release()