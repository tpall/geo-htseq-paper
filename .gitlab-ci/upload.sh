#!/bin/bash
# Upload files to Zenodo.
# Adapted from jhpoelen/zenodo-upload
#
# usage: ./zenodo_upload.sh [deposition id] [filename]
#

set -xe

DEPOSITION=$1
FILEPATH=$2
FILENAME=$(echo $FILEPATH | sed 's+.*/++g')

BUCKET=$(curl -H "Accept: application/json" -H "Authorization: Bearer $ACCESS_TOKEN" "https://www.zenodo.org/api/deposit/depositions/$DEPOSITION" | jq --raw-output .links.bucket)

curl --progress-bar -o /dev/null --upload-file $FILEPATH $BUCKET/$FILENAME?access_token=$ACCESS_TOKEN
