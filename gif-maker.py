#!/usr/bin/env python3
import glob
from PIL import Image
import argparse
import textwrap
import sys
from collections import defaultdict


parser = argparse.ArgumentParser(
    prog="gif-maker.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    ************************************************************************



    Developed by Arkadiy Garber
    Please send comments and inquiries to ark@midauthobio.edu
    ************************************************************************
    '''))

parser.add_argument('-d', type=str, help="directory with images", default="NA")

parser.add_argument('-o', type=str, help="output gif", default="NA")

parser.add_argument('-x', type=str, help="filename extension", default="png")

parser.add_argument('-s', type=int, help="frame speed", default=100)

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_known_args()[0]


def lastItem(ls):
    x = ''
    for i in ls:
        if i != "":
            x = i
    return x

# key=lambda i: ((i.split(".")[0]), float(lastItem(i.split(".")[0].split("/"))) )

Dict = defaultdict(lambda: 'EMPTY')
for i in glob.glob(args.d + "/*" + args.x):
    POS = lastItem(i.split("/"))
    POS = POS.split(".")[0]
    Dict[int(POS)] = i

sortedList = []
for i in sorted(Dict.keys()):
    sortedList.append(Dict[i])


def make_gif(frame_folder):
    frames = [Image.open(image) for image in sortedList]
    # frames.reverse()
    frame_one = frames[0]
    frame_one.save(args.o, format="GIF", append_images=frames,
                   save_all=True, duration=args.s, loop=0)


if __name__ == "__main__":
    make_gif(args.d)
