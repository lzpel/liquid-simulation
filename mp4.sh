#!/usr/bin/env bash
set -eux

cd `dirname $0`/cmake-build

./ffmpeg -y -r 5 -i output%05d.prof.png -vcodec libx264 -pix_fmt yuv420p -r 30 out.mp4