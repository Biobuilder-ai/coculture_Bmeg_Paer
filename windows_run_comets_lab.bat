@echo off
setlocal ENABLEDELAYEDEXPANSION

set IMAGE_TAR=comets-lab.tar
set IMAGE_NAME=comets-lab
set CONTAINER_NAME=comets-lab-container
set PORT=8888

echo Checking Docker installation...
where docker >nul 2>nul
if %errorlevel% neq 0 (
    echo Docker is not installed or not in PATH.
    pause
    exit /b 1
)

echo Checking Docker daemon...
docker info >nul 2>nul
if %errorlevel% neq 0 (
    echo Starting Docker Desktop...
    start "" "C:\Program Files\Docker\Docker\Docker Desktop.exe"

    echo Waiting for Docker to start...
    :waitloop
    timeout /t 2 >nul
    docker info >nul 2>nul
    if %errorlevel% neq 0 goto waitloop
)

echo Docker is running.

echo Loading Docker image...
docker load -i %IMAGE_TAR% >nul

echo Removing old container (if exists)...
docker rm -f %CONTAINER_NAME% >nul 2>nul

echo Starting JupyterLab container...

docker run -d ^
    --name %CONTAINER_NAME% ^
    -p %PORT%:8888 ^
    -v %cd%\notebooks:/workspace/notebooks ^
    -w /workspace/notebooks ^
    %IMAGE_NAME% ^
    jupyter lab ^
        --ip=0.0.0.0 ^
        --port=8888 ^
        --no-browser ^
        --allow-root ^
        --ServerApp.token='' ^
        --ServerApp.password=''

echo Waiting for Jupyter to become available...
timeout /t 5 >nul

echo Opening browser...
start http://localhost:%PORT%

echo.
echo JupyterLab is running at:
echo http://localhost:%PORT%
echo.
echo To stop it later, run:
echo docker stop %CONTAINER_NAME%
echo.
pause

