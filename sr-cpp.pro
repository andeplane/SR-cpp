TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++11

LIBS += -L/usr/local/Cellar/armadillo/5.200.1/lib -larmadillo
INCLUDEPATH += /usr/local/Cellar/armadillo/5.200.1/include
SOURCES += main.cpp \
    sr.cpp \
    lorentzgenerator.cpp


include(deployment.pri)
qtcAddDeployment()

HEADERS += \
    sr.h \
    lorentzgenerator.h

