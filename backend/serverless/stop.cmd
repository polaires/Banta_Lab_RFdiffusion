@echo off
REM Stop the RFdiffusion Docker container
REM Usage: stop.cmd

echo Stopping RFdiffusion container...
wsl docker compose -f docker-compose.local.yml down
echo Container stopped.
