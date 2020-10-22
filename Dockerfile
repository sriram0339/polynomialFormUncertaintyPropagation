FROM ubuntu:20.10 AS stage1

COPY ./src /usr/polyforms/src
COPY CMakeLists.txt /usr/polyforms/CMakeLists.txt
COPY ./cmake-includes /usr/polyforms/cmake-includes
COPY ./test /usr/polyforms/test
COPY ./Dockerfile /usr/polyforms/Dockerfile
RUN apt-get update && apt-get install -y cmake libmpfr-dev libmpfi-dev flex bison libglpk-dev gcc g++ libboost-dev

FROM stage1 AS buildStage
COPY ./README.txt /usr/polyforms/README.txt
WORKDIR /usr/polyforms/
RUN cmake ./CMakeLists.txt && make
ENV LD_LIBRARY_PATH=/usr/local/lib
ENV DISPLAY=host.docker.internal:0
