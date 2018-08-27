./demo.exe 1 1 -v &> /tmp/11.txt
diff 11.txt /tmp/11.txt 
./demo.exe 1 2 -v &> /tmp/12.txt
diff 12.txt /tmp/12.txt 
./demo.exe 2 1 -v &> /tmp/21.txt
diff 21.txt /tmp/21.txt 
./demo.exe 2 2 -v &> /tmp/22.txt
diff 22.txt /tmp/22.txt 
