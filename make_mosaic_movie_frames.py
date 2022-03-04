import argparse
from PIL import Image
import glob
import os

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-left", "--left_files", help="left files", type=str)
    parser.add_argument("-right", "--right_files", help="right files", type=str)
    parser.add_argument("-save_dir", "--save_directory", help="save_directory", type=str)
    args = parser.parse_args()
    return args
    
args = parse_inputs()

left_files = args.left_files
right_files = args.right_files
left_frames = sorted(glob.glob(left_files))
right_frames = sorted(glob.glob(right_files))
save_dir = args.save_directory
if os.path.exists(save_dir) == False:
    os.makedirs(save_dir)

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

print("Finished making mosaic frames")
