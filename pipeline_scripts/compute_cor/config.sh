#!/bin/bash
istr="/path/to"
ostr="$PWD"
cat compute_correlation_KHV.sh > pipeline_compute_correlation_KHV.sh
eval "sed -i -e 's#"$istr"#"$ostr"#g' pipeline_compute_correlation_KHV.sh"

chmod -R 777 *.R

echo "Done!"
