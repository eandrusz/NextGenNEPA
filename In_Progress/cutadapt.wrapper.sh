#!/usr/bin/env bash
echo "This is the working directory"
pwd
# Usage bash cutadapt.wrapper.sh {folder with fastqs} {metadata_file} {output_folder}
echo "This was the command executed"
echo "bash cutadapt.wrapper.sh ${1} ${2} ${3} ${4}"
MAIN_DIR="$(dirname "$0")"
echo "${MAIN_DIR}"
echo "That is the directory with the scripts"
SCRIPT_DIR="${MAIN_DIR}"/core
for file in "${SCRIPT_DIR}"/* ; do
	source "${file}"
done

# Capture input folder
INPUT_DIR=${1}
#primers_file=${2}
# Make Output Folder and logfile
OUTPUT_DIR=${3}/noprimers
mkdir "${OUTPUT_DIR}"
LOGFILE="${OUTPUT_DIR}"/log.cutadapt.txt

# Capture Minimum sequence length
MIN_LENGTH=${4}

# Capture metadata_file
SEQUENCING_METADATA=${2}
OUTPUT_METADATA="${OUTPUT_DIR}"/output.metadata.csv
## Thanks to Jimmy O'Donnell - fix metadata lines if needed
if [[ $( file "${SEQUENCING_METADATA}" ) == *"CRLF"* ]]; then

  echo "The file has CRLF endings. Let me fix that for you..."

  BASE="${SEQUENCING_METADATA%.*}"

  EXT="${SEQUENCING_METADATA##*.}"

  NEWLINES_FIXED="${BASE}"_fix."${EXT}"

  tr -d '\r' < "${SEQUENCING_METADATA}" > "${NEWLINES_FIXED}"

  echo "the old file was: ${SEQUENCING_METADATA}"

  echo "The new file is here:"

  echo "${NEWLINES_FIXED}"

else

  echo "The file passes test for CRLF. Everybody dance!"
  echo

fi

if [[ -s "${NEWLINES_FIXED}" ]]; then
	SEQUENCING_METADATA="${NEWLINES_FIXED}"
fi

## READ METADATA ###########

METADATA_DIM=($( awk -F, 'END{print NR, NF}' "${SEQUENCING_METADATA}" ))
echo "Metadata has" "${METADATA_DIM[0]}" "rows and" "${METADATA_DIM[1]}" "columns including header."
N_SAMPLES=$((METADATA_DIM[0] - 1))
echo "Expecting" "${N_SAMPLES}" "samples total."
echo
# Filnames
COLNUM_FILE1=$( get_colnum "file1" "${SEQUENCING_METADATA}")
COLNUM_FILE2=$( get_colnum "file2" "${SEQUENCING_METADATA}")

# Sample names
COLNUM_SAMPLE=$( get_colnum "Sample_name" "${SEQUENCING_METADATA}")

# Primers
COLNUM_PRIMER1=$( get_colnum "PrimerF" "${SEQUENCING_METADATA}")
COLNUM_PRIMER2=$( get_colnum "PrimerR" "${SEQUENCING_METADATA}")

# Locus ID

COLNUM_LOCUS=$( get_colnum "Locus" "${SEQUENCING_METADATA}")

all_columns=(COLNUM_FILE1 COLNUM_FILE2 COLNUM_LOCUS \
COLNUM_SAMPLE COLNUM_PRIMER1 COLNUM_PRIMER2)

echo "Checking that all columns in metadata are there"

for column in "${all_columns[@]}" ; do

 if [ "${!column}" -gt 0 ]; then
	 echo "looking good, ${column}"
 else
  echo "Something went wrong with column name ${column}"
	echo "exiting script"
	exit
fi
done
echo "All columns passed test"

### GET FILENAMES and PRIMERS and Sample SAMPLE_NAMES and LocusInfo

FILE1=($(awk -F',' -v COLNUM=$COLNUM_FILE1 \
  'NR>1 {  print $COLNUM }' $SEQUENCING_METADATA \
   ))

FILE2=($(awk -F',' -v COLNUM=$COLNUM_FILE2 \
  'NR>1 {print $COLNUM}' $SEQUENCING_METADATA \
    ))

PRIMER1=($(awk -F',' -v COLNUM=$COLNUM_PRIMER1 \
    'NR > 1 { print $COLNUM }' $SEQUENCING_METADATA \
      ))

PRIMER2=($(awk -F',' -v COLNUM=$COLNUM_PRIMER2 \
    'NR > 1 { print $COLNUM }' $SEQUENCING_METADATA \
     ))

SAMPLE_NAMES=($(awk -F',' -v COLNUM=$COLNUM_SAMPLE \
      'NR>1 { print $COLNUM }' "${SEQUENCING_METADATA}" \
       ))

LOCUS=($(awk -F',' -v COLNUM=$COLNUM_LOCUS \
      'NR>1 { print $COLNUM }' "${SEQUENCING_METADATA}"\
        ))


# Start metadata file
echo "Sample_name,Locus,file1,file2" > "${OUTPUT_METADATA}"


# For loop to make cutadapt cut all pcr primers

for (( i=0; i < "${#FILE1[@]}"; i++ )); do    # for each fwd file

echo "${FILE1[i]}"

short_file_F=$(basename "${FILE1[i]%.*}")
short_file_R=$(basename "${FILE2[i]%.*}")

OUTPUT_FILE_F="${OUTPUT_DIR}"/Locus_"${LOCUS[i]}"_"${short_file_F}".fastq
OUTPUT_FILE_R="${OUTPUT_DIR}"/Locus_"${LOCUS[i]}"_"${short_file_R}".fastq

cutadapt -g "${PRIMER1[i]}" \
         --discard-untrimmed \
         -G "${PRIMER2[i]}" \
         -o "${OUTPUT_FILE_F}" \
         -p "${OUTPUT_FILE_R}" \
         -m "${MIN_LENGTH}" \
         "${INPUT_DIR}"/"${FILE1[i]}" "${INPUT_DIR}"/"${FILE2[i]}"  2>> "${LOGFILE}"

         #-o "${OUTPUT_DIR}"/Locus_"${LOCUS}"_"${short_file_F}".fastq \

         #-p "${OUTPUT_DIR}"/Locus_"${LOCUS}"_"${short_file_R}".fastq \


printf "%s,%s,%s,%s" \
"${SAMPLE_NAMES[i]}" "${LOCUS[i]}" \
Locus_"${LOCUS[i]}"_"${short_file_F}".fastq \
Locus_"${LOCUS[i]}"_"${short_file_R}".fastq >> "${OUTPUT_METADATA}"
#Locus_"${LOCUS[i]}"_"${short_file_F[i]}".fastq Locus_"${LOCUS[i]}"_"${short_file_R[i]}".fastq \

printf "\n" >> "${OUTPUT_METADATA}"

 done
