#!/bin/bash
#
# Start the COMETS JupyterLab container as a Jupyter *server* that you connect
# to from VS Code (instead of using JupyterLab in the browser).
#
# In VS Code: open a notebook -> click the kernel picker (top-right) ->
#   "Select Another Kernel..." -> "Existing Jupyter Server..." ->
#   enter the URL printed at the end of this script (http://localhost:8888) ->
#   then pick the "Python 3" kernel. That kernel runs inside the container,
#   with COMETS / cometspy available.
#
# macOS only.

set -e

IMAGE_TAR="comets-lab.tar"
IMAGE_NAME="comets-lab"
CONTAINER_NAME="comets-lab-container"
PORT=8888

# Package dir (this script's location) and the project root one level up.
PACKAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$PACKAGE_DIR/.." && pwd)"

echo "Checking Docker installation..."
if ! command -v docker &> /dev/null; then
    echo "Docker is not installed."
    exit 1
fi

echo "Checking Docker daemon..."
if ! docker info &> /dev/null; then
    echo "Docker daemon not running. Starting Docker Desktop..."
    open --background -a Docker
    echo "Waiting for Docker to start..."
    while ! docker info &> /dev/null; do
        sleep 2
    done
fi

# Load the image only if it isn't already present (saves time on restarts).
if ! docker image inspect "$IMAGE_NAME" &> /dev/null; then
    echo "Loading Docker image (first run only)..."
    docker load -i "$IMAGE_TAR" > /dev/null
else
    echo "Docker image '$IMAGE_NAME' already loaded."
fi

echo "Removing old container (if exists)..."
docker rm -f "$CONTAINER_NAME" > /dev/null 2>&1 || true

echo "Starting Jupyter server container..."
# The project folder is mounted at the SAME absolute path inside the container,
# so absolute and relative file paths (SBML models, notebooks, COMETS logs)
# resolve identically whether code runs on the host or in the container.
docker run -d \
    --name "$CONTAINER_NAME" \
    -p $PORT:8888 \
    -v "$PROJECT_DIR":"$PROJECT_DIR" \
    -w "$PROJECT_DIR" \
    "$IMAGE_NAME" \
    jupyter lab \
        --ip=0.0.0.0 \
        --port=8888 \
        --no-browser \
        --allow-root \
        --ServerApp.token='' \
        --ServerApp.password='' \
        --ServerApp.allow_origin='*' \
        --ServerApp.allow_remote_access=True \
        --ServerApp.disable_check_xsrf=True

echo "Waiting for the Jupyter server to come up..."
for i in $(seq 1 30); do
    if curl -sf "http://localhost:$PORT/api" > /dev/null 2>&1; then
        break
    fi
    sleep 1
done

echo ""
echo "============================================================"
echo " COMETS Jupyter server is running."
echo ""
echo " Connect from VS Code:"
echo "   1. Open a notebook (e.g. in $PROJECT_DIR)."
echo "   2. Kernel picker (top-right) -> 'Select Another Kernel...'"
echo "   3. 'Existing Jupyter Server...'  ->  enter:"
echo ""
echo "        http://localhost:$PORT"
echo ""
echo "      (leave the token empty when prompted)"
echo "   4. Choose the 'Python 3' kernel that appears."
echo ""
echo " Project folder is mounted at: $PROJECT_DIR"
echo " To stop the server later:  docker stop $CONTAINER_NAME"
echo "============================================================"
