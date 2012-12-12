f2p = f2py

f90s   = ./wrapping/mac.f90
opts   = --opt='-fopenmp' -lgomp
pylib  = ./source/macMethod.py     \
         ./source/sliderWindow.py


main:
	$(f2p) -c -m mac $(f90s) $(opts)
	mv mac.so ./includes/
	python -m compileall -l ./source 
	mv ./source/*.pyc ./includes

clean:
	rm -rf ./includes/*.pyc
	rm -rf ./includes/*.so
	rm -rf *.pyc