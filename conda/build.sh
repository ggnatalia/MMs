MMS_DIR=$PREFIX/MMs
mkdir $MMS_DIR
cp -r $SRC_DIR/* $MMS_DIR
cd $PREFIX/bin
ln -s $MMS_DIR/makemocks.py .
ln -s $MMS_DIR/utils/* .
