include (../../shared.pri)

HEADERS       += PhucFilter.h \
    PhucMesh.h \
    Point.h \
    Vertice.h \
    types.h \
    Edge.h \
    Facet.h \
    utilities.h \
    Mesh.h \
    Mesher.h \
    Vertice.h \
    utilities.h \
    types.h \


SOURCES       += PhucFilter.cpp \
    PhucMesh.cpp \
    Point.cpp \
    Vertice.cpp \
    Edge.cpp \
    Facet.cpp \
    Mesh.cpp \
    Mesher.cpp \

		
TARGET        = PhucFilter

RESOURCES += \
    PhucFilter.qrc \
