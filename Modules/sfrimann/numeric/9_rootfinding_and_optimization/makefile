all: root_optimization

root_optimization: main.py optimization.py root.py
	python main.py > output.txt
	cat output.txt

clean:
	rm -f output.txt
	rm -f *.pyc
