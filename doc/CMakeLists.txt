# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

set(DOCUMENTATIONSRC ${CMAKE_CURRENT_SOURCE_DIR})
set(DOCUMENTATION ${CMAKE_CURRENT_BINARY_DIR})

set( VERSION_STR ${${PROJECT_NAME}_VERSION_STR} )
configure_file( version.tex.in version.tex )


# HTML GENERATION =============================================================
file(MAKE_DIRECTORY ${DOCUMENTATION}/html)

add_custom_target(dwarf-P-cloudMicrophysics-IFSScheme-documentation-html
   export TEXINPUTS=${DOCUMENTATIONSRC}//:${CMAKE_CURRENT_BINARY_DIR}//: &&
   ${HTLATEX} ${DOCUMENTATIONSRC}/documentation.tex
   "${DOCUMENTATIONSRC}/styling.cfg,html,3,next,NoFonts"
   WORKING_DIRECTORY ${DOCUMENTATION}/html )


# If tex4ht successful, create img dir and copy images across
# .png and .jpg images --------------------------------------------------------
FILE(GLOB_RECURSE imgfiles
    "introduction/imgs/*.png" "introduction/imgs/*.jpg"
    "getting-started/*/imgs/*.png" "getting-started/*/imgs/*jpg"
    "core-functionalities/*/imgs/*.png" "core-functionalities/*/imgs/*jpg"
    "applications/*/imgs/*.png" "applications/*/imgs/*jpg")

add_custom_command(
   TARGET dwarf-P-cloudMicrophysics-IFSScheme-documentation-html
   POST_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory ${DOCUMENTATION}/html/imgs)

FOREACH(img ${imgfiles})
   ADD_CUSTOM_COMMAND(TARGET dwarf-P-cloudMicrophysics-IFSScheme-documentation-html
       POST_BUILD COMMAND
       ${CMAKE_COMMAND} -E copy ${img} ${DOCUMENTATION}/html/imgs)
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
       TARGET dwarf-P-cloudMicrophysics-IFSScheme-documentation-html
       POST_BUILD COMMAND
       ${CONVERT} ${pdf} ${DOCUMENTATION}/html/imgs/${BASENAME}.png)
ENDFOREACH()

# =============================================================================

FIND_PROGRAM(BIBTEX bibtex)
MARK_AS_ADVANCED(BIBTEX)

# PDF GENERATION ==============================================================
add_custom_target(dwarf-P-cloudMicrophysics-IFSScheme-documentation-pdf
   export TEXINPUTS=${DOCUMENTATIONSRC}//:${CMAKE_CURRENT_BINARY_DIR}//: &&
   ${PDFLATEX} --output-directory ${DOCUMENTATION} ${DOCUMENTATIONSRC}/documentation.tex
   COMMAND TEXMFOUTPUT=${DOCUMENTATION} ${BIBTEX} ${DOCUMENTATION}/documentation.aux
   COMMAND export TEXINPUTS=${CMAKE_SOURCE_DIR}//: &&
   ${PDFLATEX} --output-directory ${DOCUMENTATION} ${DOCUMENTATIONSRC}/documentation.tex
   ${PDFLATEX} --output-directory ${DOCUMENTATION} ${DOCUMENTATIONSRC}/documentation.tex
   WORKING_DIRECTORY ${DOCUMENTATIONSRC} )
# =============================================================================

add_custom_target(dwarf-P-cloudMicrophysics-IFSScheme_doc_pdf)
add_dependencies(dwarf-P-cloudMicrophysics-IFSScheme_doc_pdf dwarf-P-cloudMicrophysics-IFSScheme-documentation-pdf)

add_custom_target(dwarf-P-cloudMicrophysics-IFSScheme_doc_html)
add_dependencies(dwarf-P-cloudMicrophysics-IFSScheme_doc_html dwarf-P-cloudMicrophysics-IFSScheme-documentation-html)

add_custom_target(dwarf-P-cloudMicrophysics-IFSScheme_doc)
add_dependencies(dwarf-P-cloudMicrophysics-IFSScheme_doc dwarf-P-cloudMicrophysics-IFSScheme_doc_pdf)

if( PROJECT_NAME STREQUAL CMAKE_PROJECT_NAME )
  add_custom_target(doc)
  add_dependencies(doc dwarf-P-cloudMicrophysics-IFSScheme_doc)
endif()

