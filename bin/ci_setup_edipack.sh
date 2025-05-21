#Building EDIpack
#Errors
set -e

cd edipack
mkdir build
cd build

echo "cmake .."
cmake ..

echo "make"
make -j

echo "make install"
make install

echo "source ~/opt/edipack/gnu/*/bin/edipack_config_user.sh" >> ~/.edipack_config_user
echo -e "\e[32m EDIpack installed and sourced \e[0m"

