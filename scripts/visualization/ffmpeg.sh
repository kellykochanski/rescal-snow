#!/bin/bash

# Use ffmpeg to make a video from png files

#ffmpeg -framerate 10 -i DUN0000-%02d.png -s:v 1280x720  \
#	-pix_fmt yuv420p rescal_demo.mp4
#ffmpeg -r 1 -i DUN0000-%2d.png -pix_fmt yuv420p -r 10 rescal_demo.mp4

ffmpeg -r 20 -f image2 -i DUN%06d.png -f mp4 -q:v 0 -vcodec mpeg4 -r 20 rescal_demo.mp4
