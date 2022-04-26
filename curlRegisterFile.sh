#!/bin/bash

if [ $# -ne 3 ]; then
  echo "Please input userID, API key and WURCS file path!" 1>&2
  exit 1
fi

while read line
do
  echo ""
  curl -X POST --header 'Content-Type: application/json' --header 'Accept: application/json' --user $1:$2 -d '{ "sequence": "'$line'" }' 'https://api.glytoucan.org/glycan/register'
done < $3
