MMS_DIR=$PREFIX/MMs
#$PREFIX is the MMs directory if there is nothing to build. Otherwise, it is the common conda env structure
mkdir $MMS_DIR
cp -r $SRC_DIR/* $MMS_DIR
cd $PREFIX/bin
ln -s $MMS_DIR/makemocks.py .
ln -s $MMS_DIR/utils/* .
#rm $PREFIX/lib/python3.7/site-packages/iss
#ln -s $MMS_DIR/extlibs/iss/ $PREFIX/lib/python3.7/site-packages/
ln -s $MMS_DIR/bin/* $PREFIX/bin/
cd $MMS_DIR/libs
#echo pwd
python setup.py build_ext --inplace
cd ..
