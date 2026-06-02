@echo off
cd /d E:\sunlab_JYLee\_github\structural-time-integrators
git add -A
set /p msg="Commit message (e.g. update files): "
git commit -m "%msg%"
git push
echo.
echo Done. Press any key to close.
pause >nul