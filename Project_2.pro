TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp

LIBS += -llapack -lblas -larmadillo -lcurl
QMAKE_CXXFLAGS += -std=c++0x -lrt

