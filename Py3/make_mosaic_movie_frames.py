import sys
from PIL import Image
import glob

left_frames = sorted(glob.glob(sys.argv[1]+'*.jpg'))
right_frames = sorted(glob.glob(sys.argv[2]+'*.jpg'))
save_dir = sys.argv[3]

if len(left_frames) < len(right_frames):
    range_number = len(left_frames)
else:
    range_number = len(right_frames)

for image_ind in range(range_number):
    file_name = save_dir + "movie_frame_" + ("%06d" % image_ind)
    images = [Image.open(x) for x in [left_frames[image_ind], right_frames[image_ind]]]
    
    widths, heights = zip(*(i.size for i in images))

    if images[0].size[1] > images[1].size[1]:
        ratio = images[1].size[1]/images[0].size[1]
        images[0] = images[0].resize((int(images[0].size[0]*ratio), int(images[0].size[1]*ratio)))
    else:
        ratio = images[0].size[1]/images[1].size[1]
        images[1] = images[1].resize((int(images[1].size[0]*ratio), int(images[1].size[1]*ratio)))
    
    widths, heights = zip(*(i.size for i in images))

    total_width = sum(widths)
    max_height = max(heights)

    new_im = Image.new('RGB', (total_width, max_height))

    x_offset = 0
    for im in images:
        new_im.paste(im, (x_offset,0))
        x_offset += im.size[0]
        
    new_im.save(file_name+'.jpg')
    print('Created frame:',file_name)
