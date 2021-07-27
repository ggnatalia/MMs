MMS_DIR=$PREFIX/MMs
#$PREFIX is the MMs directory if there is nothing to build. Otherwise, it is the common conda env structure
mkdir $MMS_DIR
cp -r $SRC_DIR/* $MMS_DIR
cd $PREFIX/bin
ln -s $MMS_DIR/makemocks.py .
ln -s $MMS_DIR/utils/* .
# Download NanoSim
cd $PREFIX
git clone https://github.com/bcgsc/NanoSim.git
cd NanoSim/pre-trained_models
# Extract metagenome model
tar -xzvf metagenome_ERR3152364_Even.tar.gz
cd ../..
#rm $PREFIX/lib/python3.7/site-packages/iss
#ln -s $MMS_DIR/extlibs/iss/ $PREFIX/lib/python3.7/site-packages/
ln -s $MMS_DIR/bin/* $PREFIX/bin/
cd $MMS_DIR/libs
#echo pwd
python3 setup.py build_ext --inplace
cd ..
