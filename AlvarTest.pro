#-------------------------------------------------
#
# Project created by QtCreator 2015-10-15T14:54:30
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = AlvarTest
TEMPLATE = app

QT      += multimedia multimediawidgets

include(c:/Projects/Workspace/ProjectsCommon/Druai.pri)
include(c:/Projects/Workspace/ProjectsCommon/Boost.pri)
include(c:/Projects/Workspace/ProjectsCommon/OpenCV2411.pri)

INCLUDEPATH += C:\Projects\Libraries\OpenCV2.4.8\sources\include\opencv

SOURCES += \
   AlvarTestApp.cpp \
   AlvarTestMainWindow.cpp \
   Camera.cpp \
   Line.cpp \
   Marker.cpp \
   Pose.cpp \
   Rotation.cpp \
   ConnectedComponents.cpp \
   Draw.cpp \
   Bitset.cpp \
   AlvarVersion.cpp


HEADERS  += \
   AlvarTestMainWindow.h \
   Camera.h \
   Line.h \
   Marker.h \
   Pose.h \
   Rotation.h \
   Util.h \
   Alvar.h \
   AlvarVersion.h \
   MarkerDetector.h \
   ConnectedComponents.h \
   Draw.h \
   Bitset.h \
   ../Druai/CameraCalibration.h

CONFIG += c++11 warn_on
