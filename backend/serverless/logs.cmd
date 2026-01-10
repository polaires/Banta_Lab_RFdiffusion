@echo off
REM View logs from the RFdiffusion Docker container
REM Usage: logs.cmd [-f]
REM   -f : Follow log output (like tail -f)

if "%1"=="-f" (
    wsl docker compose -f docker-compose.local.yml logs -f rfdiffusion
) else (
    wsl docker compose -f docker-compose.local.yml logs --tail=100 rfdiffusion
)
