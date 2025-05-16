# LaTeX documentation
find_package(LATEX COMPONENTS PDFLATEX BIBTEX)
if(LATEX_FOUND)
include(${CMAKE_SOURCE_DIR}/CMakeModules/UseLATEX.cmake)
set(LATEX_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR})
add_latex_document(Userguide.tex
    IMAGE_DIRS images
    BIBFILES TBsd.bib
)
add_latex_document(Formalism.tex
    IMAGE_DIRS images
    BIBFILES TBsd.bib
)
else()
message("LaTeX needs to be installed to generate the LaTeX documentation")
endif()