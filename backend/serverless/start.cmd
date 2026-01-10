@echo off
REM Start the RFdiffusion Docker container
REM Usage: start.cmd

echo Checking Docker daemon...
wsl bash -c "docker ps > /dev/null 2>&1"
if errorlevel 1 (
    echo Docker daemon is not running. Starting it...
    echo NOTE: You may need to enter your WSL password
    wsl sudo service docker start
    timeout /t 3 /nobreak > nul
)

echo.
echo Starting RFdiffusion container...
echo.

wsl docker compose -f docker-compose.local.yml up -d

echo.
echo Container started. API available at http://localhost:8000
echo.
echo Test with:
echo   curl http://localhost:8000/runsync -X POST -H "Content-Type: application/json" -d "{\"input\": {\"task\": \"health\"}}"
echo.
echo View logs with: logs.cmd
