
default: build

clean:
	rm -f xsec-dbg

xsec-dbg: xsec_calc.cxx
	g++ -o xsec-dbg xsec_calc.cxx

build: xsec-dbg

run: xsec-dbg
	mkdir $(label)
	./xsec-dbg --energy 100 \
		--ap-mass 0.1 \
		--2d $(label)/dsdxtheta.csv \
		--1d $(label)/dsdx.csv \
		--chi $(label)/chi.csv \
		--total $(label)/total.csv	

