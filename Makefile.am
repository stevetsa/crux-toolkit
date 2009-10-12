SUBDIRS = src doc

if SEARCH
  SEARCH_ENABLED = true
else
  SEARCH_ENABLED =
endif

distdir = $(PACKAGE)_$(VERSION)
distclean-local: clean-local

clean-local: clean-ac
clean-ac:
	rm -rf autom4te.cache

# makes a binary only distribution, using zip.
# make sure to 'make' before you make a 'bdist'
ar=`arch`
name=`uname`
bdist:
	rm -fr bdist
	mkdir bdist
	mkdir bdist/doc
	mkdir bdist/bin
	cp -R \
		doc/user/data/demo.ms2 \
		doc/user/data/small-yeast.fasta \
		doc/user/*.html \
		doc/user/*.png \
		doc/user/*.gif \
		doc/user/default.params \
		doc/user/*.css \
		bdist/doc/
	cp 	src/c/crux \
		src/c/crux-get-ms2-spectrum \
		src/c/crux-generate-peptides \
		src/c/crux-predict-peptide-ions \
		bin/spectacle.pl \
		bdist/bin/
	zip -r $(distdir)-$(ar)-$(name).zip bdist
	rm -fr bdist/

checksum:
	md5sum $(distdir)-$(ar)-$(name).zip > $(distdir)-$(ar)-$(name).zip.md5sum.txt
	md5sum $(distdir).tar.gz > $(distdir).tar.gz.md5sum.txt
	
ignore-files:
	svn update .svnignore
	svn propset svn:ignore -F .svnignore .

# This is needed to keep 'make distcheck' from complaining about
# the file version.h, which doesn't get cleaned up (and doesn't
# need to be  cleaned up) in 'make distclean'
distcleancheck_listfiles = $(top_srcdir)/src/c/version.h

# Copy the appropriate copyright statement, COPYING for 
# limited license distribution, COPYRIGHT for open source distribution.
dist-hook:
	if test -z "$(SEARCH_ENABLED)"; \
    then cp COPYRIGHT $(distdir)/COPYING; fi
