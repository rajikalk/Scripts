all: ODE

ODE: main.py ODE.py
	python main.py > output.txt
	cat output.txt
	eog damped_oscillator.png

clean:
	rm -f *.png
	rm -f *.txt
	rm -f *.pyc
