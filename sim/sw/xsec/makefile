
default: build

clean:
	rm -f xsec-dbg

xsec-dbg: xsec_calc.cxx
	g++ -o xsec-dbg xsec_calc.cxx

build: xsec-dbg

run: xsec-dbg
	mkdir $(label)
	./xsec-dbg --energy 5 \
		--ap-mass 0.1 \
		--2d $(label)/dsdxtheta.csv \
		--1d $(label)/dsdx.csv \
		--total $(label)/xsec.csv	

