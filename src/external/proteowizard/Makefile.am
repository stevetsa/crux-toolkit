# Automake template for external code
EXTRA_DIST=pwiz-src-without-tv-3_0_3533.tar.bz2 create_links.sh
PWIZ_LIB_PATH=install

all:
	tar -jxvf pwiz-src-without-tv-3_0_3533.tar.bz2;\
	./quickbuild.sh --without-mz5 --prefix=$(PWIZ_LIB_PATH) libraries pwiz/data/msdata/ pwiz/data/identdata/ 1> out 2> er;\
	./create_links.sh --prefix=$(PWIZ_LIB_PATH)







