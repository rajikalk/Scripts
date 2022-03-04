all: fft

fft: main.py fft.py
	python main.py
	eog FourierTransform.png
	eog CompositePower.png
	eog Signal.png

clean:
	rm -f *.png
	rm -f *.pyc
