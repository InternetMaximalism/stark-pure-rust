FROM rust

RUN rustup toolchain install nightly --allow-downgrade
RUN rustup default nightly

WORKDIR /root
