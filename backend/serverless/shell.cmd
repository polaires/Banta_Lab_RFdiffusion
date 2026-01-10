@echo off
REM Open a shell in the running RFdiffusion container
REM Usage: shell.cmd

echo Opening shell in RFdiffusion container...
wsl docker compose -f docker-compose.local.yml exec rfdiffusion /bin/bash
