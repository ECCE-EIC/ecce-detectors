Imported from https://gitlab.com/eic/detector/pid/-/tree/master/dRICH 

# dRICH parameterisation

The main file is `dualRICH.h` that contains the definition and implementation of the `dualRICH_Aerogel` and `dualRICH_C2F6` classes.
Notice that the aerogel and the C2F6 radiators are implemented in two separate classes, therefore the full dRICH PID parameterisation is split over the two radiators.

Other files are used to define generic interfaces, such a `genericRICH.h` and `genericDetector.h` and are an integral part of this software.

Both `dualRICH_Aerogel` and `dualRICH_C2F6` are by default operated in `RICH mode`. The extension of the operation in `threshold mode` can be achieved by requesting it explicitely with `setThresholdMode(true)`.

A test example can be found in `dRICH.C`.

For details and information about usage, please contact Roberto Preghenella (preghenella@bo.infn.it).
