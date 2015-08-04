cd build/bin

./particle2dtest

./unicycletest

./body3dtest

./qrotortest ../../bin/qrotor.cfg

./hrotortest

./helicartest

./body2dtracktest2 ../../bin/body2dtrack2.cfg

./cartest

./rccartest

##./sddpbulletrccartest

./airmopttest

# @mk that below got broken at some point: need to investigate
#./fixedchaintest ../../bin/fixedchain1.cfg

./airbottest ../../bin/airbot.cfg

# ./imuekftest

./uuvtest ../../bin/jhurov.cfg

#./cecartest ../../bin/cecar.cfg

./body3dstab

cd ../..
build/bin/body3drhc bin/body3drhc.cfg

build/bin/body3ddemstab bin/body3ddemstab.cfg


build/bin/dynvisest xxx bin/dynvisest.cfg
