@echo off
echo Installing dependencies...
"C:\Users\albey\AppData\Local\Python\pythoncore-3.14-64\python.exe" -m pip install -r requirements.txt --quiet
echo.
echo Starting CPB Analyzer at http://localhost:8888
echo Press Ctrl+C to stop.
echo.
"C:\Users\albey\AppData\Local\Python\pythoncore-3.14-64\python.exe" app.py
pause
