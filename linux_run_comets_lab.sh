#!/bin/bash

set -e

IMAGE_TAR="comets-lab.tar"
IMAGE_NAME="comets-lab"
CONTAINER_NAME="comets-lab-container"
PORT=8888
NOTEBOOK_DIR="$(pwd)/notebooks"

echo "Checking Docker installation..."
if ! command -v docker &> /dev/null; then
    echo "Docker is not installed."
    exit 1
fi

echo "Checking Docker daemon..."
if ! docker info &> /dev/null; then
    echo "Docker daemon not running. Attempting to start..."
    if command -v systemctl &> /dev/null; then
        sudo systemctl start docker
    elif command -v service &> /dev/null; then
        sudo service docker start
    elif command -v dockerd &> /dev/null; then
        sudo dockerd > /tmp/dockerd.log 2>&1 &
    else
        echo "Could not determine how to start Docker."
        exit 1
    fi

    echo "Waiting for Docker to start..."
    while ! docker info &> /dev/null; do
        sleep 2
    done
fi

# Create notebooks directory if missing
[ ! -d "$NOTEBOOK_DIR" ] && mkdir -p "$NOTEBOOK_DIR"

echo "Loading Docker image..."
docker load -i "$IMAGE_TAR" > /dev/null

echo "Removing old container (if exists)..."
docker rm -f "$CONTAINER_NAME" > /dev/null 2>&1 || true

echo "Starting JupyterLab container..."
docker run -d \
    --name "$CONTAINER_NAME" \
    -p $PORT:8888 \
    -v "$NOTEBOOK_DIR":/workspace/notebooks \
    -w /workspace/notebooks \
    "$IMAGE_NAME" \
    jupyter lab \
        --ip=0.0.0.0 \
        --port=8888 \
        --no-browser \
        --allow-root \
        --ServerApp.token='' \
        --ServerApp.password=''

sleep 5

# Open browser automatically if Linux GUI available
if command -v xdg-open &> /dev/null; then
    echo "Opening browser..."
    xdg-open http://localhost:$PORT
else
    echo ""
    echo "Could not detect browser opener. Open manually:"
    echo "http://localhost:$PORT"
fi

echo ""
echo "JupyterLab is running at: http://localhost:$PORT"
echo "To stop it later, run: docker stop $CONTAINER_NAME"

