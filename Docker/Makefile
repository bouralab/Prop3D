# Definitions
build_output = runtime/molmimic
runtime_fullpath = $(realpath runtime)
build_tool = runtime-container.DONE
#git_commit ?= $(shell git log --pretty=oneline -n 1 -- ../molmimic | cut -f1 -d " ")
name = edraizen/molmimic
tag = 0.0.1

# Steps
build:
	#Uncomment if custom toil
	rm -rf ../toil
	cp -r ../../toil ../toil

	docker build --no-cache -t ${name}:${tag} -f ./Dockerfile ../
	-docker rmi -f ${name}:latest
	docker tag ${name}:${tag} ${name}:latest
	touch ${build_tool}

push: build
	# Requires ~/.dockercfg
	sudo docker push ${name}:${tag}
	sudo docker push ${name}:latest

# test: build
# 	python test.py

clean:
	-rm ${build_tool}
	-rm ${build_output}
