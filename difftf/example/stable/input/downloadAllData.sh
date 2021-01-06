# This script can be used to download all required data for the example analysis.
# Simply execute with "sh downloadAllData.sh"

filename="exampleData.tar.gz"
curl -0 https://www.embl.de/download/zaugg/diffTF/$filename -o $filename  && tar xvzf $filename  && rm $filename

filename="mm10.fa.tar.gz"
curl -0 https://www.embl.de/download/zaugg/diffTF/referenceGenome/$filename -o $filename  && mkdir referenceGenome && tar xvzf $filename -C referenceGenome  && rm $filename

filename="TFBS_mm10_PWMScan_HOCOMOCOv10.tar.gz"
curl -0 https://www.embl.de/download/zaugg/diffTF/TFBS/$filename -o $filename  && tar xvzf $filename  && rm $filename
