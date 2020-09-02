#!/usr/bin/env python
import glob
import sys
import os
from PIL import Image
import PIL.ImageDraw as ImageDraw
import PIL.ImageFont as ImageFont

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-logos", "--logo_directory", default="./", type=str)
    parser.add_argument("-EU_logo_scale", "--EU_logo_scale", default=0.5, type=float)
    parser.add_argument("-KU_logo_scale", "--KU_logo_scale", default=0.5, type=float)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

#=======MAIN=======
args = parse_inputs()

input_dir = sys.argv[1]
save_dir = sys.argv[2]
if os.path.exists(save_dir) == False:
    os.makedirs(save_dir)
    
#load logos:
#EU logo:
logo_dir = args.logo_directory
EU_flag = Image.open(logo_dir+"EU_flag.png")
EU_flag_width = int(EU_flag.size[0]*args.EU_logo_scale)
EU_flag_height = int(EU_flag.size[1]*args.EU_logo_scale)
EU_flag.thumbnail((EU_flag_width, EU_flag_height), Image.ANTIALIAS)
EU_text = Image.open(logo_dir+"EU_text.png")
EU_text_width = int(EU_text.size[0]*args.EU_logo_scale)
EU_text_height = int(EU_text.size[1]*args.EU_logo_scale)
EU_text.thumbnail((EU_text_width, EU_text_height), Image.ANTIALIAS)
font = ImageFont.truetype("Arial Bold.ttf",27)

#KU_png
KU_logo = Image.open(logo_dir+"KU_icon.png")
KU_logo_width = int(KU_logo.size[0]*args.KU_logo_scale)
KU_logo_height = int(KU_logo.size[1]*args.KU_logo_scale)
KU_logo.thumbnail((KU_logo_width, KU_logo_height), Image.ANTIALIAS)
KU_text = Image.open(logo_dir+"KU_text.png")
KU_text_width = int(KU_text.size[0]*args.KU_logo_scale*1.1)
KU_text_height = int(KU_text.size[1]*args.KU_logo_scale*1.1)
KU_text.thumbnail((KU_text_width, KU_text_height), Image.ANTIALIAS)

movie_frames = sorted(glob.glob(input_dir+"movie_frame*.jpg"))
for frame in movie_frames:
    #Open image wiht PIL:
    try:
        frame_pil = Image.open(frame)
        frame_pil.paste(EU_flag, (10,(frame_pil.size[1] - EU_flag.size[1]-10)), EU_flag)
        frame_pil.paste(EU_text, (10 + EU_flag.size[0] + 10,(frame_pil.size[1] - EU_flag.size[1]+43)), EU_text)
        draw = ImageDraw.Draw(frame_pil)
        draw.text((10 + EU_flag.size[0] + 9,(frame_pil.size[1] - EU_flag.size[1]+68)), "INTERACTIONS Fellowship", font=font, fill='rgb(0, 51, 153)')
        frame_pil.paste(KU_logo, (frame_pil.size[0]-KU_logo.size[0]-5,(frame_pil.size[1] - KU_logo.size[1])), KU_logo)
        frame_pil.paste(KU_text, (frame_pil.size[0]-KU_logo.size[0]-15-KU_text.size[0],(frame_pil.size[1] - KU_text.size[1]-20)), KU_text)
        file_name = save_dir + frame.split('/')[-1]
        frame_pil.save(file_name)
        print("Added logos to frame:", frame)
    except:
        print("Couldn't open frame:", frame)
