@ECHO OFF
cd classifier
python classifier.py > results.txt
cd..
move .\classifier\results.txt .\
PAUSE
