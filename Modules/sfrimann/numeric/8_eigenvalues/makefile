all: jacobi

jacobi: main.py jacobi.py
	python main.py > output.txt
	cat output.txt

clean: 
	rm -f output.txt
	rm -f *.pyc
