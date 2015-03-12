R=/home/jweile/bin/Rscript
NOW := $(shell date +"%y-%m-%d")
STABLE=v0.1

#default target: unit test and build
all: build_stable

# #run all unit tests
# test: setup_test
# 	cd testbench;\
# 	for test in *_test.R; do \
# 		$(R) $$test >$${test}.log;\
# 	done

# #setup a test directory
# setup_test:
# 	mkdir -p testbench/lib/
# 	mkdir -p testbench/res/
# 	mkdir -p testbench/html/
# 	cp src/main/R/* testbench/lib/
# 	cp src/main/resources/* testbench/res/
# 	cp src/main/html/* testbench/html
# 	cp src/test/R/* testbench/
# 	cp src/test/bash/* testbench/
# 	cp src/test/resources/* testbench/res/

# #delete test data
# clean:
# 	rm -r testbench/

#load the latest version of the code
latest:
	hg update

#load the last stable version of the code
stable:
	hg update $(STABLE)

#build a zip file with the final product
_build: 
	#create build directory
	mkdir -p kiloseq/lib/
	# mkdir -p kiloseq/html/
	mkdir -p kiloseq/res/
	#copy code into dir structure
	cp src/main/R/* kiloseq/lib/
	cp src/main/resources/* kiloseq/res/
	# cp src/main/html/* kiloseq/html/
	cp src/main/bash/* kiloseq/
	mv kiloseq/demuxStats.sh kiloseq/lib/
	#interpolate binary locations
	src/make/interpolate.sh bins.cfg kiloseq/
	#setup permissions
	chmod u+x kiloseq/lib/*.R
	chmod u+x kiloseq/*.sh
	#zip build
	zip kiloseq_$(NOW).zip -r kiloseq/ 
	#delete build directory
	rm -r kiloseq/

build_latest: latest _build

build_stable: stable _build

#install the built contents as well as required dependencies
install: 
	$(eval R=`grep Rscript bins.cfg|cut -d, -f2`)
	$(R) src/make/install_dependencies.R
	unzip kiloseq_$(NOW).zip -d $${HOME}
	# mkdir -p $${HOME}/www/html/
	# if ! [ -h $${HOME}/www/html/kiloseq ]; then \
	# 	ln -s $${HOME}/kiloseq/ $${HOME}/www/html/;\
	# fi
