# !bin/bash
#
# dist_bin.sh
# assemble .ZIP files for Frame3dd binary distribution release

export VERSION=20091022
echo $VERSION

# clean out prior distribution files
echo "cleaning out prior distribution files ... "
rm -rf dist/Frame3DD

# make directory structure for distribution .ZIP files in the ./dist directory
echo "creating directory structure for new distribution ... "
mkdir dist/Frame3DD
mkdir dist/Frame3DD/doc
mkdir dist/Frame3DD/doc/img
mkdir dist/Frame3DD/examples
mkdir dist/Frame3DD/matlab

# copy documentation
echo "copying documentation ... "
cp --preserve=mode,timestamps doc/*.html            dist/Frame3DD/doc/.
cp --preserve=mode,timestamps doc/template.3dd      dist/Frame3DD/doc/.
cp --preserve=mode,timestamps doc/img/*.png         dist/Frame3DD/doc/img/.
cp --preserve=mode,timestamps doc/img/*.jpg         dist/Frame3DD/doc/img/.
cp --preserve=mode,timestamps doc/img/*.gif         dist/Frame3DD/doc/img/.
cp --preserve=mode,timestamps doc/LinuxTerminal.pdf dist/Frame3DD/doc/.	

# copy ChangeLog, license, and README files
cp --preserve=mode,timestamps ChangeLog.txt         dist/Frame3DD/.
cp --preserve=mode,timestamps LICENSE.txt           dist/Frame3DD/.
cp --preserve=mode,timestamps README.txt            dist/Frame3DD/.
cp --preserve=mode,timestamps README-win32.txt      dist/Frame3DD/.

# copy example files
echo "copying example files ... "
cp --preserve=mode,timestamps examples/ex*.3dd      dist/Frame3DD/examples/.
cp --preserve=mode,timestamps examples/ex*.out      dist/Frame3DD/examples/.
cp --preserve=mode,timestamps examples/saveplot     dist/Frame3DD/examples/.
cp --preserve=mode,timestamps examples/saveplot_w32 dist/Frame3DD/examples/.

# copy matlab files
echo "copying matlab files ... "
cp --preserve=mode,timestamps matlab/*              dist/Frame3DD/matlab/.

# make an empty temp directory
echo "making an empty Frame3DD/temp directory ... "
mkdir dist/Frame3DD/temp/

# assemble the .zip file
echo "assembling .zip file ... "
cd dist                                           # change to trunk/dist
#zip -r Frame3DD_$(date +%Y%m%d).zip Frame3DD/*
zip -r Frame3DD_$(echo $VERSION).zip Frame3DD/*

# make copies of .zip files for each operating system
#cp Frame3DD_$(date +%Y%m%d).zip Frame3DD_$(date +%Y%m%d)_linux.zip
#cp Frame3DD_$(date +%Y%m%d).zip Frame3DD_$(date +%Y%m%d)_osx.zip
#cp Frame3DD_$(date +%Y%m%d).zip Frame3DD_$(date +%Y%m%d)_win32.zip

cp Frame3DD_$(echo $VERSION).zip Frame3DD_$(echo $VERSION)_linux.zip
cp Frame3DD_$(echo $VERSION).zip Frame3DD_$(echo $VERSION)_osx.zip
cp Frame3DD_$(echo $VERSION).zip Frame3DD_$(echo $VERSION)_win32.zip

# add Linux executable
echo "adding Linux executable ... "
cp --preserve=mode,timestamps ../build/frame3dd        Frame3DD/.
# zip Frame3DD_$(date +%Y%m%d)_linux.zip Frame3DD/frame3dd
zip Frame3DD_$(echo $VERSION)_linux.zip Frame3DD/frame3dd
rm Frame3DD/frame3dd

# add OS X executable
echo "adding OS X executable ... "
cp --preserve=mode,timestamps ../build/frame3ddosx     Frame3DD/frame3dd
#zip Frame3DD_$(date +%Y%m%d)_osx.zip Frame3DD/frame3dd
zip Frame3DD_$(echo $VERSION)_osx.zip Frame3DD/frame3dd
rm Frame3DD/frame3dd

# add Windows executable
echo "adding Windows executable ... "
cp --preserve=mode,timestamps ../build/frame3dd.exe    Frame3DD/.
#zip Frame3DD_$(date +%Y%m%d)_win32.zip Frame3DD/frame3dd.exe
zip Frame3DD_$(echo $VERSION)_win32.zip Frame3DD/frame3dd.exe
rm Frame3DD/frame3dd.exe

#rm Frame3DD_$(date +%Y%m%d).zip
rm Frame3DD_$(echo $VERSION).zip
rm -rf Frame3DD

echo "Frame3DD .zip archive complete"

# Uploading distribution files to sourceforge.net ...
#
# rsync -uav *.zip hpgavin,frame3dd@frs.sourceforge.net:/home/frs/project/f/fr/frame3dd/frame3dd/0.20091022/.
#
# rsync -uav *.bz2  hpgavin,frame3dd@frs.sourceforge.net:/home/frs/project/f/fr/frame3dd/frame3dd/0.20091022/.

# ----------------------------------------------------------- dist_bin.sh
# Henri P. Gavin  2009.10.20
# updated: 2009.10.22, 2009.10.25, 2009.10.27

