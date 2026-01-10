@echo off
REM Build the Docker container for RFdiffusion
REM Usage: build.cmd [--no-cache]

echo Checking Docker daemon...
wsl bash -c "docker ps > /dev/null 2>&1"
if errorlevel 1 (
    echo Docker daemon is not running. Starting it...
    echo NOTE: You may need to enter your WSL password
    wsl sudo service docker start
    timeout /t 3 /nobreak > nul
)

echo.
echo Building RFdiffusion Docker container...
echo.

if "%1"=="--no-cache" (
    wsl docker compose -f docker-compose.local.yml build --no-cache
) else (
    wsl docker compose -f docker-compose.local.yml build
)

echo.
echo Build complete. Run 'start.cmd' to start the container.
