DEMULT_DIR="/Users/ramon.gallegosimon/Projects/Nextera_Dada2/data"

FILEPATH=$( find "${DEMULT_DIR}" -name *".fastq"  )
#FILEPATH2=$( ls "${DEMULT_DIR}/"*".fastq" )
FILEPATH3="${DEMULT_DIR}"/*.fastq
#

echo "${FILEPATH[@]}"
#cho "${FILEPATH2[@]}"
echo "${FILEPATH3}"[@]
exit
#This is the good way of doing it, skipping the ls
for file in "${DEMULT_DIR}"/*.fastq; do

#file in "${FILEPATH2}"; do
  echo "${file}"
  cat "${file}"  |\
  #Now find the sequence line of each read (starting from the
  #second line, every fourlines),
  sed -n '2~4p' |
  #And now extract the 6 characters of the barcode
  cut -c 4-9 |
  #Now sort | uniq  to get the counts
  sort | uniq -c | sort -r | head -n 5

done
