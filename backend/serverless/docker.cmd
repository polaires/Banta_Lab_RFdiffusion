@echo off
REM Docker helper script - runs Docker commands via WSL
REM Usage: docker.cmd [docker command]
REM Example: docker.cmd compose up -d
REM          docker.cmd build .

wsl docker %*
