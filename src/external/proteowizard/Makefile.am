# Automake template for external code
#PWIZ_TARBALL=pwiz-src-without-v-3_0_3949.tar.bz2
PWIZ_TARBALL=pwiz-src-without-tv-3_0_3533.tar.bz2


EXTRA_DIST=${PWIZ_TARBALL} create_links.sh
PWIZ_LIB_PATH=install

all:
	tar -jxvf ${PWIZ_TARBALL};\
	./quickbuild.sh --without-mz5 --prefix=$(PWIZ_LIB_PATH) libraries pwiz/data/msdata/ pwiz/data/identdata/ 1> out 2> er;\
	./create_links.sh --prefix=$(PWIZ_LIB_PATH)







