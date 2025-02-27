run:
    make -C build
    ./build/bridge
setup:
    mkdir -p build
    cmake -B build
