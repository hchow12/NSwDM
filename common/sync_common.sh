#! /bin/bash

file_to_sync=$(ls)

echo "syncing common files..."

for file in $file_to_sync
do
	echo $file
	if [ $file = "sync_common.sh" ]
	then
		break
	fi
	cp $file ../binary/
	cp $file ../gw_ft/
done

echo "Done syncing common files"
