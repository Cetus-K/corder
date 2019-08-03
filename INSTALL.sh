
export PREFIX=$(pwd)

# install dependencies
cd $PREFIX/deps/voro++; make; make install;
cd $PREFIX/deps/qhull; make;
cd $PREFIX/deps/includes; make;
cd $PREFIX/deps/linalg; make;

# install corder
cd $PREFIX/src; make;
cd $PREFIX;
mkdir ./bin; mv ./src/corder ./bin
