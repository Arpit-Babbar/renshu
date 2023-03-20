# Run with bash, not sh

FILENAME=$1  # Include .pdf extension
RAND=$RANDOM # Needed to avoid overwriting errors
NEW_FILE=${RAND}_${FILENAME}
echo $RAND
echo $FILENAME
cp $FILENAME $NEW_FILE
gs -o /tmp/$NEW_FILE -dNoOutputFonts -sDEVICE=pdfwrite $FILENAME
soffice --infilter=impress_pdf_import --convert-to pptx /tmp/$NEW_FILE
rm /tmp/$NEW_FILE
PPTX_NAME_WITH_RAND="${NEW_FILE%.*}"
# Trim extension https://stackoverflow.com/questions/12152626/how-can-i-remove-the-extension-of-a-filename-in-a-shell-script
PPTX_NAME_WITHOUT_RAND="${FILENAME%.*}"
mv $PPTX_NAME_WITH_RAND.pptx $PPTX_NAME_WITHOUT_RAND.pptx
echo "Final file is ${PPTX_NAME_WITHOUT_RAND}.pptx"
