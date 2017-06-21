# Alvar
Update to the original ALVAR Augmented Reality Library by VTT Technical Research Centre of Finland by Druai Consulting.

Converted the Alvar library to use the C++ interface to the OpenCV library, Version 2.4.11.  Cleaned up some
issues, worked on const correctness, changed some argument passing for efficiency, updated to some C++ 11 features.

So far only the detection of single markers has been done since I'm mostly interested in the library to use in
robot localization and not its augmented reality features.  If and when I get around to updating the rest of 
the library is unknown at this point.

Currently, the updated modules are compiling as standalone with a test application and haven't been reincorporated
into a library build.

Tim Craig 03/15/16
