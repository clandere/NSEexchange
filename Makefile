rebuild:
	R CMD build nseExchangeSource
	R CMD INSTALL NSEexchange_0.1.tar.gz --library=~/cubfits/preston/

clean:
	rm -rf NSEexchange
	rm -rf NSEexchange_0.1.tar.gz
