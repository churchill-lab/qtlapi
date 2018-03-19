set -ex

# SET THE FOLLOWING VARIABLES
# docker hub username
USERNAME=mattjvincent
# image name
IMAGE=qtlapi

docker build -t $USERNAME/$IMAGE:latest .



