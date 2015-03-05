all:
	 cd ./Scripts; tar xvzf weblogo-3.2.tar.gz ; cd ../Programs; make;
clean:
	cd ./Programs; make clean
