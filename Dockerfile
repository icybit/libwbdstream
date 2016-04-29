FROM ubuntu:14.04

ENV PKG_BUILD_DIR /tmp
ENV CZMQ_VERSION 3.0.2
ENV SODIUM_VERSION 1.0.8
ENV ZMQ_VERSION 4.1.4

RUN apt-get update && apt-get install -y \
	autoconf \
	automake \
	build-essential \
	curl \
	gdb \
	libtool \
	pkg-config \
	uuid-dev \
&& rm -rf /var/lib/apt/lists/*

WORKDIR $PKG_BUILD_DIR

RUN curl -sL https://download.libsodium.org/libsodium/releases/libsodium-$SODIUM_VERSION.tar.gz | tar xz \
	&& cd libsodium-$SODIUM_VERSION \
	&& ./autogen.sh && ./configure && make && make install \
	&& cd $PKG_BUILD_DIR \
	&& rm -rf libsodium-$SODIUM_VERSION

RUN curl -sL http://download.zeromq.org/zeromq-$ZMQ_VERSION.tar.gz | tar xz \
	&& cd zeromq-$ZMQ_VERSION \
	&& ./configure && make && make install && ldconfig \
	&& cd $PKG_BUILD_DIR \
	&& rm -rf zeromq-$ZMQ_VERSION

RUN curl -sL https://github.com/zeromq/czmq/archive/v$CZMQ_VERSION.tar.gz | tar xz \
	&& cd czmq-$CZMQ_VERSION \
	&& ./autogen.sh && ./configure && make && make install && ldconfig \
	&& cd $PKG_BUILD_DIR \
	&& rm -rf czmq-$CZMQ_VERSION

#ENV GDB_VERSION 7.7.1
#RUN curl -sL http://ftp.gnu.org/gnu/gdb/gdb-$GDB_VERSION.tar.gz | tar xz \
#	&& cd gdb-$GDB_VERSION \
#	&& ./configure && make && make install \
#	&& cd $PKG_BUILD_DIR \
#	&& rm -rf gdb-$GDB_VERSION

RUN apt-get purge -y \
	autoconf \
	automake \
	build-essential \
	curl \
	libtool \
	pkg-config \
	uuid-dev

RUN groupadd -r d-cell && useradd -r -m -g d-cell d-cell
WORKDIR /home/d-cell/d-cell
COPY . .
RUN chown -R d-cell:d-cell /home/d-cell && chmod g+s /home/d-cell
USER d-cell
RUN make

EXPOSE 5556

CMD ["build/d-cell"]
