#/bin/bash

#tmp folder - $dir/tmp
#input folder - $dir/input

#input parameters
inputfile='test.txt'

dir=$PWD

if [ -f $dir/sample/$inputfile ] 
then
	echo 'File found';
else
	echo 'File Not Found - Exiting program';
	exit;
fi

#clean tmp folder
rm -rf $dir/tmp
mkdir -p tmp

#test
cd $dir/tmp
for i in {1..100}
do
	m4 -DFILE_NUMBER=1.${i} $dir/sample/test.txt > ${i}_${inputfile}
done
cd ..

echo 'Script Done'
