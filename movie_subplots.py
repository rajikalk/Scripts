import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

maximage = 906
dir1 = ""
dir2 = "Density_xy_plane/"
image_name = "movie_frame_"
format = ".png"

for it in range(906):
    img1 = mpimg.imread(dir1 + image_name + ("%06d" % it) + format)
    img2 = mpimg.imread(dir2 + image_name + ("%06d" % it) + format)
    fig = plt.figure()
    xy = fig.add_subplot(1,2,1)
    implot = plt.imshow(img1)
    plt.axis('off')
    xz = fig.add_subplot(1,2,2)
    implot = plt.imshow(img2)
    plt.axis('off')
    plt.savefig("frame_" + ("%06d" % it) + ".png", bbox_inches='tight')