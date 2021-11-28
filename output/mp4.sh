set -eux

cd `dirname $0`
cd ../cmake-build-debug

./ffmpeg -y -r 5 -i output_%04d.prof.png -vcodec libx264 -pix_fmt yuv420p -r 30 out.mp4