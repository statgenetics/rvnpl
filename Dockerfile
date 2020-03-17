FROM gaow/seqlink

MAINTAINER Linhai Zhao (Linhai.Zhao@bcm.edu)

# RV-NPL
ARG DUMMY=unknown
RUN cd /tmp && wget https://github.com/statgenetics/rvnpl/archive/master.tar.gz && tar xzvf master.tar.gz && cd rvnpl-master && python setup.py install && rm -rf /tmp/*

# To build,
# docker build --build-arg DUMMY=`date +%s` -t statisticalgenetics/rvnpl .
# docker push statisticalgenetics/rvnpl

# To use,
# alias rvnpl='docker run --rm --security-opt label:disable -t '\
#	'-P -h RV-NPL -w $PWD -v /tmp:/tmp -v $PWD:$PWD '\
#	'-u $UID:${GROUPS[0]} -e HOME=/seqlink -e USER=$USER statisticalgenetics/rvnpl rvnpl'
