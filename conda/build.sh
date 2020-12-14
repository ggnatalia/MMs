MMS_DIR=$PREFIX/MMs
echo $PREFIX
echo $MMS_DIR
mkdir $MMS_DIR
cp -r $SRC_DIR/* $MMS_DIR
echo $SRC_DIR
cd $PREFIX/bin
ln -s $MMS_DIR/makemocks.py .
ln -s $MMS_DIR/utils/* .
ln -s $MMS_DIR/extlibs/iss/ ./../lib/python3.7/site-packages/
ln -s $MMS_DIR/bin/iss .

