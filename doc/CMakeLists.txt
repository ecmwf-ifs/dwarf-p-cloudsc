set(USERGUIDESRC ${CMAKE_CURRENT_SOURCE_DIR})
set(USERGUIDE ${CMAKE_CURRENT_BINARY_DIR})



# HTML GENERATION =============================================================
file(MAKE_DIRECTORY ${USERGUIDE}/html)

add_custom_target(dwarf2-user-guide-html
   export TEXINPUTS=${CMAKE_SOURCE_DIR}//:${USERGUIDESRC}//: &&
   ${HTLATEX} ${USERGUIDESRC}/user-guide.tex
   "${USERGUIDESRC}/styling.cfg,html,3,next,NoFonts"
   WORKING_DIRECTORY ${USERGUIDE}/html )


# If tex4ht successful, create img dir and copy images across
# .png and .jpg images --------------------------------------------------------
FILE(GLOB_RECURSE imgfiles
    "introduction/imgs/*.png" "introduction/imgs/*.jpg"
    "getting-started/*/imgs/*.png" "getting-started/*/imgs/*jpg"
    "core-functionalities/*/imgs/*.png" "core-functionalities/*/imgs/*jpg"
    "applications/*/imgs/*.png" "applications/*/imgs/*jpg")

add_custom_command(
   TARGET dwarf2-user-guide-html
   POST_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory ${USERGUIDE}/html/imgs)

FOREACH(img ${imgfiles})
   ADD_CUSTOM_COMMAND(TARGET dwarf2-user-guide-html
       POST_BUILD COMMAND
       ${CMAKE_COMMAND} -E copy ${img} ${USERGUIDE}/html/imgs)
ENDFOREACH()
# -----------------------------------------------------------------------------

# .pdf images -----------------------------------------------------------------
FILE(GLOB_RECURSE pdffiles
    "introduction/imgs/*.pdf"
    "getting-started/*/imgs/*.pdf"
    "core-functionalities/*/imgs/*.pdf"
    "applications/*/imgs/*.pdf")

FIND_PROGRAM(CONVERT convert)

FOREACH(pdf ${pdffiles})
   GET_FILENAME_COMPONENT(BASENAME ${pdf} NAME_WE)
   ADD_CUSTOM_COMMAND(
       TARGET dwarf2-user-guide-html
       POST_BUILD COMMAND
       ${CONVERT} ${pdf} ${USERGUIDE}/html/imgs/${BASENAME}.png)
ENDFOREACH()

# =============================================================================

FIND_PROGRAM(BIBTEX bibtex)
MARK_AS_ADVANCED(BIBTEX)

# PDF GENERATION ==============================================================
add_custom_target(dwarf2-user-guide-pdf
   export TEXINPUTS=${CMAKE_SOURCE_DIR}//: &&
   ${PDFLATEX} --output-directory ${USERGUIDE} ${USERGUIDESRC}/user-guide.tex
   COMMAND TEXMFOUTPUT=${USERGUIDE} ${BIBTEX} ${USERGUIDE}/user-guide.aux
   COMMAND export TEXINPUTS=${CMAKE_SOURCE_DIR}//: &&
   ${PDFLATEX} --output-directory ${USERGUIDE} ${USERGUIDESRC}/user-guide.tex
   ${PDFLATEX} --output-directory ${USERGUIDE} ${USERGUIDESRC}/user-guide.tex
   WORKING_DIRECTORY ${USERGUIDESRC} )
# =============================================================================

add_custom_target(dwarf2_doc_pdf)
add_dependencies(dwarf2_doc_pdf dwarf2-user-guide-pdf)

add_custom_target(dwarf2_doc_html)
add_dependencies(dwarf2_doc_html dwarf2-user-guide-html)

add_custom_target(dwarf2_doc)
add_dependencies(dwarf2_doc dwarf2_doc_pdf)

if( PROJECT_NAME STREQUAL CMAKE_PROJECT_NAME )
  add_custom_target(doc)
  add_dependencies(doc dwarf2_doc)
endif()

